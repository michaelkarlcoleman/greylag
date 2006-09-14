

//	$Id$


#include "cxtpy.h"


#include <cerrno>
#include <cstdlib>
#include <cmath>
#include <stdexcept>

#include <iostream>


// speed up reading, at some expense to portability
// fast version assumes no threads and single-precision good enough
// (mass and intensity values currently under 7 digits)
// XXX: even faster: mmap + strtol only??

//#define FAST

#ifndef FAST
#define strtodX std::strtod
#define fgetsX std::fgets
#define ferrorX std::ferror
#else
#define strtodX std::strtof
#define fgetsX fgets_unlocked
#define ferrorX ferror_unlocked
#endif


parameters parameters::the;


char *
peak::__repr__() const {
  static char temp[128];
  sprintf(temp, "<peak mz=%.4f intensity=%.4f>", mz, intensity);
  return &temp[0];
}


int spectrum::next_id = 1;
int spectrum::next_physical_id = 1;

std::vector<spectrum> spectrum::searchable_spectra;
// parent mass -> searchable_spectra index
std::multimap<double,
	      std::vector<spectrum>::size_type> spectrum::spectrum_mass_index;


char *
spectrum::__repr__() const {
  static char temp[1024];
  sprintf(temp,
	  "<spectrum #%d (phys #%d) '%s' %.4f/%+d"
	  " [%d peaks, maxI=%f, sumI=%f]>",
	  id, physical_id, name.c_str(), mass, charge, peaks.size(),
	  max_peak_intensity, sum_peak_intensity);
  return &temp[0];
}



// This is an error exit for read_spectra.
static void
read_error(FILE *f, const char *message="") {
  if (ferrorX(f))
    message = "I/O error while reading ms2 file";
  throw std::ios_base::failure(message);
}


// Read spectra from a file in ms2 format.  Multiply charge spectra (e.g.,
// +2/+3) are split into separate spectra having the same physical id.  The
// peak list is initially sorted by mz.  Throws an exception on invalid
// input.

std::vector<spectrum>
spectrum::read_spectra(FILE *f, const int file_id) {
  std::vector<spectrum> spectra;

  const int bufsiz = 1024;
  char buf[bufsiz];
  char *endp;

  char *result = fgetsX(buf, bufsiz, f);
  while (true) {
    std::vector<std::string> names;
    std::vector<double> masses;
    std::vector<int> charges;
    std::vector<peak> peaks;
    assert(peaks.empty());

    // read headers
    while (true) {
      if (ferrorX(f))
	read_error(f);
      if (not result or buf[0] != ':')
	break;
      names.push_back(std::string(buf+1, buf+std::strlen(buf)-1));
      if (not fgetsX(buf, bufsiz, f))
	read_error(f, "bad ms2 format: mass/charge line expected");
      errno = 0;
      masses.push_back(std::strtod(buf, &endp)); // need double accuracy here
      if (errno or endp == buf)
	read_error(f, "bad ms2 format: bad mass");
      charges.push_back(std::strtol(endp, &endp, 10));
      if (errno or endp == buf)
	read_error(f, "bad ms2 format: bad charge");
      result = fgetsX(buf, bufsiz, f);
    }
    if (not result)
      break;
    // read peaks
    while (true) {
      peak p;
      errno = 0;
      p.mz = strtodX(buf, &endp);
      if (errno or endp == buf)
	read_error(f, "bad ms2 format: bad peak mz");
      p.intensity = strtodX(endp, &endp);
      if (errno or endp == buf)
	read_error(f, "bad ms2 format: bad peak intensity");
      peaks.push_back(p);

      result = fgetsX(buf, bufsiz, f);
      if (ferrorX(f))
	read_error(f);
      if (not result or buf[0] == ':')
	break;
    }

    // add spectra to vector
    if (names.empty() and not peaks.empty())
      read_error(f, "bad ms2 format: missing header lines?");

    std::sort(peaks.begin(), peaks.end(), peak::less_mz);

    for (std::vector<double>::size_type i=0; i<names.size(); i++) {
      spectrum sp(masses[i], charges[i]);
      sp.peaks = peaks;
      sp.name = names[i];
      sp.file_id = file_id;
      sp.physical_id = next_physical_id;
      sp.calculate_intensity_statistics();
      spectra.push_back(sp);
    }
    spectrum::next_physical_id++;
  }
  return spectra;
}


// Sets max/sum_peak_intensity, according to peaks and normalization_factor.
void
spectrum::calculate_intensity_statistics() {
  sum_peak_intensity = 0;
  max_peak_intensity = -1;

  for (std::vector<peak>::const_iterator it=peaks.begin(); it != peaks.end();
       it++) {
    sum_peak_intensity += it->intensity;
    max_peak_intensity = std::max<double>(it->intensity, max_peak_intensity);
  }
  max_peak_intensity *= normalization_factor;
  sum_peak_intensity *= normalization_factor;
}


// For each *non-overlapping* group of peaks with mz within an interval of
// (less than) length 'limit', keep only the most intense peak.  (The
// algorithm works from the right, which matter here.)

// FIX: Is this really the best approach?  The obvious alternative would be to
// eliminate peaks so that no interval of length 'limit' has more than one
// peak (eliminating weak peaks to make this so).

static void
remove_close_peaks(std::vector<peak> &peaks, double limit) {
  for (std::vector<peak>::size_type i=0; i+1<peaks.size(); i++) {
    double ub_mz = peaks[i].mz + limit;
    while (i+1 < peaks.size() and peaks[i+1].mz < ub_mz)
      peaks.erase(peaks.begin()
		  + (peaks[i+1].intensity > peaks[i].intensity
		     ? i : i+1));
  }
}

// static void
// remove_close_peaks(std::vector<peak> &peaks, double limit) {
//   for (std::vector<peak>::size_type i=0; i+1<peaks.size();)
//     if (peaks[i+1].mz - peaks[i].mz < limit)
//       if (peaks[i+1].intensity > peaks[i].intensity)
// 	peaks.erase(peaks.begin() + i);
//       else
// 	peaks.erase(peaks.begin() + i+1);
//     else
//       i++;
// }



// Examine this spectrum to see if it should be kept.  If so, return true and
// sort the peaks by mz, normalize them and limit their number.
bool
spectrum::filter_and_normalize(double minimum_fragment_mz,
			       double dynamic_range, int minimum_peaks,
			       int total_peaks) {
  // "spectrum, maximum parent charge", noise suppression check,
  // "spectrum, minimum parent m+h" not implemented

  if (minimum_fragment_mz < 0 or dynamic_range <= 0 or minimum_peaks < 0
      or total_peaks < 0)
    throw std::invalid_argument("invalid argument value");

  // remove_isotopes
  // >>> removes multiple entries within 0.95 Da of each other, retaining the
  // highest value. this is necessary because of the behavior of some peak
  // finding routines in commercial software

  // peaks already sorted by mz upon read
  remove_close_peaks(peaks, 0.95);

  // remove parent
  // >>> set up m/z regions to ignore: those immediately below the m/z of the
  // parent ion which will contain uninformative neutral loss ions, and those
  // immediately above the parent ion m/z, which will contain the parent ion
  // and its isotope pattern

  const double parent_mz = 1.00727 + (mass - 1.00727) / charge;
  for (std::vector<peak>::size_type i=0; i<peaks.size();)
    if (parent_mz >= peaks[i].mz
	? std::abs(parent_mz - peaks[i].mz) < (50.0 / charge)
	: std::abs(parent_mz - peaks[i].mz) < (5.0 / charge))
      peaks.erase(peaks.begin() + i);
    else
      i++;
      

  // remove_low_masses
  std::vector<peak>::iterator pi;
  for (pi=peaks.begin(); pi != peaks.end() and pi->mz <= minimum_fragment_mz;
       pi++);
  peaks.erase(peaks.begin(), pi);


  // normalize spectrum
  // >>> use the dynamic range parameter to set the maximum intensity value
  // for the spectrum.  then remove all peaks with a normalized intensity < 1
  if (peaks.empty())
    return false;
  std::vector<peak>::iterator most_intense_peak
    = std::max_element(peaks.begin(), peaks.end(), peak::less_intensity);
  normalization_factor
    = std::max<double>(1.0, most_intense_peak->intensity) / dynamic_range;
  for (std::vector<peak>::iterator p=peaks.begin(); p != peaks.end(); p++)
    p->intensity /= normalization_factor;
  for (std::vector<peak>::size_type i=0; i<peaks.size();)
    if (peaks[i].intensity < 1.0)
      peaks.erase(peaks.begin() + i);
    else
      i++;

  if (peaks.size() < (unsigned int) minimum_peaks)
    return false;
    
  //     # check is_noise (NYI)
  //     # * is_noise attempts to determine if the spectrum is simply
  //     # noise. if the spectrum
  //     # * does not have any peaks within a window near the parent ion mass,
  //     # it is considered
  //     # * noise.
  //     #limit = sp.mass / sp.charge
  //     #if sp.charge < 3:
  //     #    limit = sp.mass - 600.0
  //     #if not [ True for mz, i in sp.peaks if mz > limit ]:
  //     #    return "looks like noise"

  // clean_isotopes removes peaks that are probably C13 isotopes
  remove_close_peaks(peaks, 1.5);

  // keep only the most intense peaks
  if (peaks.size() > (unsigned int) total_peaks) {
    std::sort(peaks.begin(), peaks.end(), peak::less_intensity);
    peaks.erase(peaks.begin(), peaks.end() - total_peaks);
    std::sort(peaks.begin(), peaks.end(), peak::less_mz);
  }

  // update statistics, now that we've removed some peaks
  calculate_intensity_statistics();
  return true;
}


// Store the list of spectra that search_peptide will search against, and also
// build spectrum_mass_index.
void
spectrum::set_searchable_spectra(const std::vector<spectrum> &spectra) {
  searchable_spectra = spectra;
  for (std::vector<spectrum>::size_type i=0; i<spectra.size(); i++)
    spectrum_mass_index.insert(std::make_pair(spectra[i].mass, i));
}


int 
spectrum::search_peptide(int idno, int offset, int begin,
			 const std::string &peptide_seq, bool is_N,
			 bool is_C, int missed_cleavage_count,
			 score_stats stats) {
  return 0;
}



// Return the similarity score between this spectrum and that, and also a
// count of common peaks in *peak_count.

// xtandem quantizes mz values into bins of width fragment_mass_error first,
// and adds a bin on either side.  This makes the effective fragment mass
// error for each peak an arbitrary value between 1x and 2x the parameter,
// based on quantization of that peak.  If quirks_mode is on, we'll try to
// emulate this behavior.
double
spectrum::score_similarity(const spectrum &x, const spectrum &y,
			   int *peak_count) {
  const double frag_err = parameters::the.fragment_mass_error;
  const bool quirks_mode = parameters::the.quirks_mode;

  double score = 0;
  *peak_count = 0;
  std::vector<peak>::const_iterator x_it = x.peaks.begin();
  std::vector<peak>::const_iterator y_it = y.peaks.begin();
  
  while (x_it != x.peaks.end() and y_it != y.peaks.end()) {
    double x_mz = x_it->mz, y_mz = y_it->mz;
    if (quirks_mode) {
      const float ffe = frag_err;
      x_mz = std::floor(static_cast<float>(x_mz) / ffe) * ffe;
      y_mz = std::floor(static_cast<float>(y_mz) / ffe) * ffe;
    }
    // in quirks mode, must be within one frag_err-wide bin
    if (std::abs(x_mz - y_mz) < (quirks_mode ? frag_err*1.5 : frag_err)) {
      *peak_count += 1;
      score += (x_it->intensity * y_it->intensity);
      //std::cerr << "hit " << x_mz << " vs " << y_mz << ", i = " << x_it->intensity * y_it->intensity << std::endl;
#if 0
      // assume peaks aren't "close" together
      x_it++;
      y_it++;
      continue;
#endif
    }
    if (x_it->mz < y_it->mz)
      x_it++;
    else
      y_it++;
  }

  return score;
}


// returns the monoisotopic mass
double
get_peptide_mass(const std::string &peptide_seq, bool is_N, bool is_C) {
  const parameters &CP = parameters::the;
  double NC_mods = 0;
  if (is_N)
    NC_mods += CP.modification_mass['['];
  if (is_C)
    NC_mods += CP.modification_mass[']'];

  double residue_sum = 0;
  for (std::string::const_iterator it=peptide_seq.begin();
       it != peptide_seq.end(); it++)
    residue_sum += (CP.monoisotopic_residue_mass[*it]
		    + CP.modification_mass[*it]);
  return (NC_mods + residue_sum
	  + CP.cleavage_N_terminal_mass_change
	  + CP.cleavage_C_terminal_mass_change
	  + CP.proton_mass);
}


// returns log-scaled hyperscore
double
scale_hyperscore(double hyper_score) {
  return 4 * std::log10(hyper_score);
}

