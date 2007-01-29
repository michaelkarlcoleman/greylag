// C++ module for greylag

//     Copyright (C) 2006-2007, Stowers Institute for Medical Research
//
//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License along
//     with this program; if not, write to the Free Software Foundation, Inc.,
//     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


#include "cgreylag.hpp"


#include <cerrno>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>

#include <iostream>


#ifdef __GNU__
#define NOTHROW __attribute__ ((nothrow))
// use faster, non-threadsafe versions, since we don't use threads
#define fgetsX fgets_unlocked
#define fputsX fputs_unlocked
#define fprintfX __builtin_fprintf_unlocked
#define ferrorX ferror_unlocked
#else
#define NOTHROW
#define fgetsX std::fgets
#define fputsX std::fputs
#define fprintfX std::fprintf
#define ferrorX std::ferror
#endif


parameters parameters::the;


char *
peak::__repr__() const {
  static char temp[128];
  sprintf(temp, "<peak mz=%.4f intensity=%.4f>", mz, intensity);
  return &temp[0];
}

const int FIRST_SPECTRUM_ID = 1;
int spectrum::next_id = FIRST_SPECTRUM_ID;
int spectrum::next_physical_id = FIRST_SPECTRUM_ID;

std::vector<spectrum> spectrum::searchable_spectra;
// parent mass -> searchable_spectra index
std::multimap<double,
	      std::vector<spectrum>::size_type> spectrum::spectrum_mass_index;

typedef
  std::multimap<double, std::vector<spectrum>::size_type>::const_iterator
  spmi_c_it;			// ugh


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



// This is an error exit for read_spectra and filter_ms2_by_mass.
static inline void
io_error(FILE *f, const char *message="") {
  if (ferrorX(f))
    message = "I/O error while reading (or writing) spectrum file";
  throw std::ios_base::failure(message);
}


// Return sorted list of parent masses present in the given ms2 files.
std::vector<double>
spectrum::read_ms2_spectrum_masses(std::vector<int> fds) {
  std::vector<double> masses;

  for (std::vector<int>::const_iterator fdit=fds.begin(); fdit != fds.end();
       fdit++) {
    FILE *f = fdopen(*fdit, "r");

    const int bufsiz = 1024;
    char buf[bufsiz];
    char *endp;
    char *result = fgetsX(buf, bufsiz, f);
    while (true) {
      // read headers
      while (true) {
	if (ferrorX(f))
	  io_error(f);
	if (not result or buf[0] != ':')
	  break;
	if (not fgetsX(buf, bufsiz, f))
	  io_error(f, "bad ms2 format: mass/charge line expected");
	errno = 0;
	masses.push_back(std::strtod(buf, &endp)); // need double accuracy here
	if (errno or endp == buf)
	  io_error(f, "bad ms2 format: bad mass");
	result = fgetsX(buf, bufsiz, f);
      }
      if (not result)
	break;
      // read peaks
      std::vector<peak> peaks;
      while (true) {
	result = fgetsX(buf, bufsiz, f);
	if (ferrorX(f))
	  io_error(f);
	if (not result or buf[0] == ':')
	  break;
      }
    }
  }
  std::sort(masses.begin(), masses.end());
  return masses;
}


// Read spectra from file in ms2 format, tagging them with file_id.  Spectra
// with charge zero are omitted from the result.  (All of them are read,
// though, to keep spectrum ids in sync.)  If file_id == -1, the ms2 file is
// an annotated file produced by split_ms2_by_mass_band.

// Multiply charge spectra (e.g., +2/+3) are split into separate spectra
// having the same physical id.  The peak list is initially sorted by mz.
// Throws an exception on invalid input.

std::vector<spectrum>
spectrum::read_spectra_from_ms2(FILE *f, const int file_id) {
  std::vector<spectrum> spectra;

  const int bufsiz = 1024;
  char buf[bufsiz];
  char *endp;

  char *result = fgetsX(buf, bufsiz, f);
  while (true) {
    std::vector<std::string> names;
    std::vector<double> masses;
    std::vector<int> charges;
    // for annotated ms2
    std::vector<int> file_ids;
    std::vector<int> physical_ids;
    std::vector<int> ids;

    // read headers
    while (true) {
      if (ferrorX(f))
	io_error(f);
      if (not result or buf[0] != ':')
	break;
      if (file_id != -1)
	names.push_back(std::string(buf+1, buf+std::strlen(buf)-1));
      else {
	char *anno_hash = std::strrchr(buf, '#'); // const, but C++ doesn't
						  // like it declared thus--yuck
	if (not anno_hash)
	  io_error(f, "bad ms2+ format: '#' not found");
	names.push_back(std::string(buf+1, anno_hash));
	errno = 0;
	file_ids.push_back(std::strtol(anno_hash+1, &endp, 10));
	physical_ids.push_back(std::strtol(endp+1, &endp, 10));
	ids.push_back(std::strtol(endp+1, &endp, 10));
	if (errno)		// FIX: improve error check
	  io_error(f, "bad ms2+ format: bad values");
      }
      if (not fgetsX(buf, bufsiz, f))
	io_error(f, "bad ms2 format: mass/charge line expected");
      errno = 0;
      masses.push_back(std::strtod(buf, &endp)); // need double accuracy here
      if (errno or endp == buf)
	io_error(f, "bad ms2 format: bad mass");
      const char *endp0 = endp;
      charges.push_back(std::strtol(endp0, &endp, 10));
      if (errno or endp == endp0)
	io_error(f, "bad ms2 format: bad charge");
      result = fgetsX(buf, bufsiz, f);
    }
    if (not result)
      break;
    // read peaks
    std::vector<peak> peaks;
    while (true) {
      peak p;
      errno = 0;
      p.mz = strtod(buf, &endp);
      if (errno or endp == buf)
	io_error(f, "bad ms2 format: bad peak mz");
      const char *endp0 = endp;
      p.intensity = strtod(endp0, &endp);
      if (errno or endp == endp0)
	io_error(f, "bad ms2 format: bad peak intensity");
      peaks.push_back(p);

      result = fgetsX(buf, bufsiz, f);
      if (ferrorX(f))
	io_error(f);
      if (not result or buf[0] == ':')
	break;
    }

    // add spectra to vector
    if (names.empty() and not peaks.empty())
      io_error(f, "bad ms2 format: missing header lines?");

    std::sort(peaks.begin(), peaks.end(), peak::less_mz);

    for (std::vector<double>::size_type i=0; i<names.size(); i++) {
      spectrum sp(masses[i], charges[i]);
      sp.peaks = peaks;
      sp.name = names[i];
      if (file_id != -1) {
	sp.file_id = file_id;
	sp.physical_id = next_physical_id;
      } else {
	sp.file_id = file_ids[i];
	sp.physical_id = physical_ids[i];
	sp.id = ids[i];
	spectrum::next_id = std::max(spectrum::next_id, sp.id+1);
	spectrum::next_physical_id = std::max(spectrum::next_physical_id,
					      sp.physical_id+1);
      }
      sp.calculate_intensity_statistics();
      spectra.push_back(sp);
    }
    if (file_id != -1)
      spectrum::next_physical_id += 1;
  }
  return spectra;
}


// Copy spectra from an ms2 file to a set of output ms2 files, one for each
// band in the set of mass bands described by their upper bounds.  Extra
// information is written in the first header line of each spectra so that
// id's may be recovered.  So, for example, ':0002.0002.1' might become
// ':0002.0002.1 # 0 45 78', where 0, 45, 78 are the spectrum's file_id,
// physical_id, and id, respectively.  Multiply charged spectra will be
// split into separate spectra (having the same physical_id).

// FIX: cleaner to just pass in a vector of filenames, rather than fds?
void
spectrum::split_ms2_by_mass_band(FILE *inf, const std::vector<int> &outfds,
				 const int file_id,
				 const std::vector<double>
				   &mass_band_upper_bounds) {
  std::map<double, FILE *> mass_file_index;

  if (outfds.size() < 1
      or outfds.size() != mass_band_upper_bounds.size())
    throw std::invalid_argument("invalid argument value");

  for (std::vector<int>::size_type i=0; i<outfds.size(); i++) {
    FILE *f = fdopen(outfds[i], "w");
    if (not f)
      throw std::ios_base::failure("fdopen failed: fd limit exceeded?");
    mass_file_index.insert(std::make_pair(mass_band_upper_bounds[i], f));
  }

  int o_next_id = FIRST_SPECTRUM_ID;
  int o_next_physical_id = FIRST_SPECTRUM_ID;

  const int bufsiz = 1024;
  char buf0[bufsiz], buf1[bufsiz];
  char *endp;
  char *result = fgetsX(buf0, bufsiz, inf);
  while (true) {
    double mass;
    std::set<FILE *> sp_files;

    // copy headers
    while (true) {
      if (ferrorX(inf))
	io_error(inf);
      if (not result or buf0[0] != ':')
	break;
      if (not fgetsX(buf1, bufsiz, inf))
	io_error(inf, "bad ms2 format: mass/charge line expected");
      errno = 0;
      mass = std::strtod(buf1, &endp); // need double accuracy here
      if (errno or endp == buf1)
	io_error(inf, "bad ms2 format: bad mass");
      errno = 0;
      std::map<double, FILE *>::const_iterator bit
	= mass_file_index.lower_bound(mass);
      if (bit == mass_file_index.end())
	throw std::invalid_argument("internal error: mass out of range");
      FILE *bandf = bit->second;
      sp_files.insert(bandf);
      assert(buf0[strlen(buf0)-1] == '\n');
      buf0[strlen(buf0)-1] = 0;	// chop newline
      fprintfX(bandf, "%s # %d %d %d\n", buf0, file_id, o_next_physical_id,
	       o_next_id++);
      fputsX(buf1, bandf);
      if (errno)
	io_error(bandf, "error writing ms2 file");
      result = fgetsX(buf0, bufsiz, inf);
    }
    if (not result)
      break;
    if (sp_files.empty())
      io_error(inf, "bad ms2 format: missing header lines?");
    // copy peaks
    while (true) {
      peak p;
      errno = 0;
      for (std::set<FILE *>::const_iterator fit = sp_files.begin();
	   fit != sp_files.end(); fit++) {
	fputsX(buf0, *fit);
	if (errno)
	  io_error(*fit, "error writing ms2 file");
      }
      result = fgetsX(buf0, bufsiz, inf);
      if (ferrorX(inf))
	io_error(inf);
      if (not result or buf0[0] == ':')
	break;
    }
    o_next_physical_id += 1;
  }

  // flush all output files (is this necessary?)
  for (std::map<double, FILE *>::const_iterator it = mass_file_index.begin();
       it != mass_file_index.end(); it++)
    std::fflush(it->second);
}


// For each *non-overlapping* group of peaks with mz within an interval of
// (less than) length 'limit', keep only the most intense peak.  (The
// algorithm works from the right, which matter here.)

// FIX: Is this really the best approach?  The obvious alternative would be to
// eliminate peaks so that no interval of length 'limit' has more than one
// peak (eliminating weak peaks to make this so).

static inline void
remove_close_peaks(std::vector<peak> &peaks, double limit) {
  for (std::vector<peak>::size_type i=0; i+1<peaks.size(); i++) {
    double ub_mz = peaks[i].mz + limit;
    while (i+1 < peaks.size() and peaks[i+1].mz < ub_mz)
      peaks.erase(peaks.begin()
		  + (peaks[i+1].intensity > peaks[i].intensity
		     ? i : i+1));
  }
}


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


// Multiply peak intensities by residue-dependent factors to generate a more
// realistic spectrum.  If is_XYZ, walk through the peptide sequence in
// reverse.
static inline void
synthesize_ladder_intensities(std::vector<peak> &mass_ladder,
			      const std::string &peptide_seq,
			      const bool is_XYZ) NOTHROW {
  const parameters &CP = parameters::the;

  assert(mass_ladder.size() == peptide_seq.size() - 1);
  assert(mass_ladder.size() >= 2);

  // An intensity boost is given for the second peptide bond/peak (*) if the
  // N-terminal residue is 'P':
  //      0 * 2 3 (peaks)
  //     A P C D E
  const int proline_bonus_position = is_XYZ ? mass_ladder.size()-2 : 1;
  mass_ladder[proline_bonus_position].intensity = 3.0;
  if (peptide_seq[1] == 'P')
    mass_ladder[proline_bonus_position].intensity = 10.0;

  if (CP.spectrum_synthesis)
    for (unsigned int i=0; i<mass_ladder.size(); i++) {
      // FIX: there must be a better way
      const char peptide_index = is_XYZ ? mass_ladder.size()-i-1 : i;
      switch (peptide_seq[peptide_index]) {
      case 'D':
	mass_ladder[i].intensity *= 5.0;
	break;
      case 'V':
      case 'E':
      case 'I':
      case 'L':
	mass_ladder[i].intensity *= 3.0;
	break;
      case 'N':
      case 'Q':
	mass_ladder[i].intensity *= 2.0;
	break;
      }
      switch (peptide_seq[peptide_index+1]) {
      case 'P':
	mass_ladder[i].intensity *= 5.0;
	break;
      }
    }
}


// Calculate peaks for a synthesized mass (not mz) ladder.
static inline void
get_synthetic_Y_mass_ladder(std::vector<peak> &mass_ladder,
			    const std::string &peptide_seq,
			    const std::vector<double> &mass_list,
			    const double fragment_C_fixed_mass) NOTHROW {
  assert(peptide_seq.size() == mass_list.size());
  double m = fragment_C_fixed_mass;

  const int ladder_size = mass_ladder.size();
  for (int i=ladder_size-1; i>=0; i--) {
    m += mass_list[i+1];
    mass_ladder[ladder_size-1-i] = peak(m, 1.0);
  }

  synthesize_ladder_intensities(mass_ladder, peptide_seq, true);
}


// Calculate peaks for a synthesized mass (not mz) ladder.
static inline void
get_synthetic_B_mass_ladder(std::vector<peak> &mass_ladder,
			    const std::string &peptide_seq,
			    const std::vector<double> &mass_list,
			    const double fragment_N_fixed_mass) NOTHROW {
  double m = fragment_N_fixed_mass;

  const int ladder_size = mass_ladder.size();
  for (int i=0; i<=ladder_size-1; i++) {
    m += mass_list[i];
    mass_ladder[i] = peak(m, 1.0);
  }

  synthesize_ladder_intensities(mass_ladder, peptide_seq, false);
}


static inline void
synthetic_spectra(spectrum synth_sp[/* max_fragment_charge+1 */][ION_MAX],
		  const std::string &peptide_seq,
		  const std::vector<double> &mass_list,
		  const double fragment_N_fixed_mass,
		  const double fragment_C_fixed_mass,
		  const double max_fragment_charge) NOTHROW {
  std::vector<peak> mass_ladder(mass_list.size()-1);

  for (ion ion_type=ION_MIN; ion_type<ION_MAX; ion_type++) {
    switch (ion_type) {
    case ION_Y:
      get_synthetic_Y_mass_ladder(mass_ladder, peptide_seq, mass_list,
				  fragment_C_fixed_mass);
      break;
    case ION_B:
      get_synthetic_B_mass_ladder(mass_ladder, peptide_seq, mass_list,
				  fragment_N_fixed_mass);
      break;
    default:
      assert(false);
    }

    //for (unsigned int i=0; i<mass_ladder.size(); i++)
    //  std::cerr << "ladder step " << peptide_seq[i] << " " << mass_ladder[i].mz << std::endl;

    for (int charge=1; charge<=max_fragment_charge; charge++) {
      spectrum &sp = synth_sp[charge][ion_type];
      //sp.mass = fragment_peptide_mass;
      sp.charge = charge;
      sp.peaks.resize(mass_ladder.size());
      for (std::vector<peak>::size_type i=0; i<mass_ladder.size(); i++) {
	sp.peaks[i].intensity = mass_ladder[i].intensity;
	sp.peaks[i].mz = peak::get_mz(mass_ladder[i].mz, charge);
      }
    }
  }
}


// Return the similarity score between this spectrum and that, and also a
// count of common peaks in *peak_count.

// xtandem quantizes mz values into bins of width fragment_mass_error first,
// and adds a bin on either side.  This makes the effective fragment mass
// error for each peak an arbitrary value between 1x and 2x the parameter,
// based on quantization of that peak.  If quirks_mode is on, we'll try to
// emulate this behavior.
static inline double
score_similarity_(const spectrum &x, const spectrum &y,
		  int *peak_count) NOTHROW {
  const double frag_err = parameters::the.fragment_mass_error;
  const bool quirks_mode = parameters::the.quirks_mode;
  const double error_limit = quirks_mode ? frag_err*1.5 : frag_err;

  double score = 0;
  *peak_count = 0;
  std::vector<peak>::const_iterator x_it = x.peaks.begin();
  std::vector<peak>::const_iterator y_it = y.peaks.begin();

  while (x_it != x.peaks.end() and y_it != y.peaks.end()) {
    double x_mz = x_it->mz, y_mz = y_it->mz;
    if (quirks_mode) {
      const float ffe = frag_err;
      x_mz = int(static_cast<float>(x_mz) / ffe) * ffe;
      y_mz = int(static_cast<float>(y_mz) / ffe) * ffe;
    }
    // in quirks mode, must be within one frag_err-wide bin
    const double delta = y_mz - x_mz;
    if (std::abs(delta) < error_limit) {
      *peak_count += 1;
      score += (x_it->intensity * y_it->intensity);
      //std::cerr << "hit " << x_mz << " vs " << y_mz << ", i = " << x_it->intensity * y_it->intensity << std::endl;
    }
    if (delta > 0)
      x_it++;
    else
      y_it++;
  }

  return score;
}


double
spectrum::score_similarity(const spectrum &x, const spectrum &y, int
			   *peak_count) {
  return score_similarity_(x, y, peak_count);
}


// IDEA: could we combine all these spectra and just do the correlation once?

// Score the specified spectrum against a group of synthetic fragment
// spectra.
static inline void
score_spectrum(double &hyper_score, double &convolution_score,
	       std::vector<double> &ion_scores, std::vector<int> &ion_peaks,
	       spectrum synth_sp[/* max_fragment_charge+1 */][ION_MAX],
	       const int spectrum_id, const int max_fragment_charge) NOTHROW {
  const parameters &CP = parameters::the;
  const spectrum sp = spectrum::searchable_spectra[spectrum_id];
  const int spectrum_charge = sp.charge;
  hyper_score = 1.0;
  convolution_score = 0;

  for (ion ion_type=ION_MIN; ion_type<ION_MAX; ion_type++) {
    int i_peaks = 0;
    double i_scores = 0.0;
    const int charge_limit = std::max<int>(1, spectrum_charge-1);
    for (int charge=1; charge<=charge_limit; charge++) {
      assert(charge <= max_fragment_charge);

      // convolution score is just sum over charges/ions
      // hyperscore is product of p! over charges/ions (where p is corr peak
      // count) times the convolution score (clipped to FLT_MAX)
      // > blurred!
      int common_peak_count;	// FIX!
      double conv_score = score_similarity_(synth_sp[charge][ion_type],
					    sp, &common_peak_count);
      i_peaks += common_peak_count;
      i_scores += conv_score;
      hyper_score *= CP.factorial[common_peak_count];
      convolution_score += conv_score;
    }

    ion_peaks[ion_type] = i_peaks;
    ion_scores[ion_type] = i_scores;
  }
  hyper_score *= convolution_score;
}


// FIX: does this actually help inlining?
// returns log-scaled hyperscore
static inline double
scale_hyperscore_(double hyper_score) NOTHROW {
  assert(hyper_score >= 0);	// otherwise check for EDOM, ERANGE
  if (hyper_score == 0)
    return -DBL_MAX;		// FIX: this shouldn't be reached?
  return 4 * std::log10(hyper_score);
}
double
scale_hyperscore(double hyper_score) { return scale_hyperscore_(hyper_score); }


struct mass_trace_list {
  mass_trace_item item;
  const mass_trace_list *next;

  mass_trace_list(const mass_trace_list *p=0) : next(p) { }
};


// Search for matches of this particular peptide modification variation
// against the spectra.  Updates score_stats and returns the number of
// candidate spectra found.
static inline void
evaluate_peptide(const search_context &context, match &m,
		 const mass_trace_list *mtlp,
		 const std::vector<double> &mass_list,
		 const spmi_c_it &candidate_spectra_info_begin,
		 const spmi_c_it &candidate_spectra_info_end,
		 score_stats &stats) NOTHROW {
  const parameters &CP = parameters::the;
  assert(m.peptide_sequence.size() == mass_list.size());

  int max_candidate_charge = 0;
  for (spmi_c_it it=candidate_spectra_info_begin;
       it != candidate_spectra_info_end; it++)
    max_candidate_charge
      = std::max<int>(max_candidate_charge,
		      spectrum::searchable_spectra[it->second].charge);
  assert(max_candidate_charge >= 1);

  int max_fragment_charge = max_candidate_charge;
  if (not CP.check_all_fragment_charges)
    max_fragment_charge = std::max<int>(1, max_candidate_charge-1);
  assert(max_fragment_charge <= spectrum::max_supported_charge);
  // e.g. +2 -> 'B' -> spectrum
  // FIX: should this be a std::vector??
  static spectrum synth_sp[spectrum::max_supported_charge+1][ION_MAX];
  synthetic_spectra(synth_sp, m.peptide_sequence, mass_list,
		    context.fragment_N_fixed_mass,
		    context.fragment_C_fixed_mass, max_fragment_charge);

  for (spmi_c_it candidate_it = candidate_spectra_info_begin;
       candidate_it != candidate_spectra_info_end; candidate_it++) {
    m.spectrum_index = candidate_it->second;
    stats.candidate_spectrum_count += 1;
    if (CP.estimate_only)
      continue;

    //std::cerr << "sp " << m.spectrum_index << std::endl;
    score_spectrum(m.hyper_score, m.convolution_score, m.ion_scores,
		   m.ion_peaks, synth_sp, m.spectrum_index,
		   max_fragment_charge);
    //std::cerr << "score " << m.convolution_score << std::endl;
    if (m.convolution_score <= 2.0)
      continue;

    ////std::cerr << "histod: " << m.spectrum_index << " " << m.peptide_sequence
    //<< " " << m.peptide_mass << " " <<
    //spectrum::searchable_spectra[m.spectrum_index].mass << std::endl;

    // update spectrum histograms
    std::vector<int> &hh = stats.hyperscore_histogram[m.spectrum_index];
    int binned_scaled_hyper_score = int(0.5 + scale_hyperscore_(m.hyper_score));
    assert(binned_scaled_hyper_score >= 0);
    // FIX
    if (static_cast<int>(hh.size())-1 < binned_scaled_hyper_score) {
      assert(hh.size() < INT_MAX/2);
      hh.resize(std::max<unsigned int>(binned_scaled_hyper_score+1,
				       2*hh.size()));
    }
    hh[binned_scaled_hyper_score] += 1;
    // incr m_tPeptideScoredCount
    // nyi: permute stuff
    // FIX
    const int sp_ion_count = m.ion_peaks[ION_B] + m.ion_peaks[ION_Y];
    if (sp_ion_count < CP.minimum_ion_count)
      continue;
    const bool has_b_and_y = m.ion_peaks[ION_B] and m.ion_peaks[ION_Y];
    if (not has_b_and_y)
      continue;

    // check that parent masses are within error range (isotope ni) already
    // done above (why does xtandem do it here?  something to do with isotope
    // jitter?)  if check fails, only eligible for 2nd-best record

    // Remember all of the highest-hyper-scoring matches against each
    // spectrum.  These might be in multiple domains (pos, length, mods) in
    // multiple proteins.

    // To avoid the problems that xtandem has with inexact float-point
    // calculations, we consider hyperscores within the "epsilon ratio" to be
    // "equal".  The best_match vector will need to be re-screened against
    // this ratio, because it may finally contain entries which don't meet our
    // criterion.  (We don't take the time to get it exactly right here,
    // because we will often end up discarding these later anyway.)

    const double current_ratio = (m.hyper_score
				  / stats.best_score[m.spectrum_index]);

    if (CP.quirks_mode and (m.hyper_score < stats.best_score[m.spectrum_index])
	or (not CP.quirks_mode
	    and (current_ratio < CP.hyper_score_epsilon_ratio))) {
      // This isn't a best score, so see if it's a second-best score.
      if (m.hyper_score > stats.second_best_score[m.spectrum_index])
	stats.second_best_score[m.spectrum_index] = m.hyper_score;
      continue;
    }
    if (CP.quirks_mode and (m.hyper_score > stats.best_score[m.spectrum_index])
	or (not CP.quirks_mode
	    and (current_ratio > (1/CP.hyper_score_epsilon_ratio)))) {
      // This score is significantly better, so forget previous matches.
      stats.improved_candidates += 1;
      //if (stats.best_match[m.spectrum_index].empty()) why doesn't this work?
      if (stats.best_score[m.spectrum_index] == 100.0)
	stats.spectra_with_candidates += 1;
      stats.best_score[m.spectrum_index] = m.hyper_score;
      stats.best_match[m.spectrum_index].clear();
    }
    // Now remember this match because it's at least as good as those
    // previously seen.
    m.mass_trace.clear();
    for (const mass_trace_list *p=mtlp; p; p=p->next)
      m.mass_trace.push_back(p->item);
    stats.best_match[m.spectrum_index].push_back(m);
  }
}


// IDEA FIX: Use a cost parameter to define the iterative search front.


// Choose one possible residue modification position.  Once they're all
// chosen, then evaluate.
static inline void
choose_residue_mod(const search_context &context, match &m,
		   const mass_trace_list *mtlp, std::vector<double> &mass_list,
		   const spmi_c_it &candidate_spectra_info_begin,
		   const spmi_c_it &candidate_spectra_info_end,
		   score_stats &stats,
		   std::vector<int> &db_remaining,
		   const unsigned remaining_positions_to_choose,
		   const unsigned next_position_to_consider) NOTHROW {
  const parameters &CP = parameters::the;
  assert(remaining_positions_to_choose
	 <= m.peptide_sequence.size() - next_position_to_consider);

  if (CP.maximum_modification_combinations_searched
      and (stats.combinations_searched
	   >= CP.maximum_modification_combinations_searched))
    return;

  if (remaining_positions_to_choose == 0) {
    stats.combinations_searched += 1;
    evaluate_peptide(context, m, mtlp, mass_list, candidate_spectra_info_begin,
		     candidate_spectra_info_end, stats);
  } else {
    mass_trace_list mtl(mtlp);

    // consider all of the positions where we could next add a mod
    for (unsigned int i=next_position_to_consider;
	 i <= m.peptide_sequence.size()-remaining_positions_to_choose; i++) {
      mtl.item.position = i;
      const double save_pos_mass=mass_list[i];
      const char pos_res = m.peptide_sequence[i];

      // consider the possibilities for this position
      for (std::vector<int>::const_iterator
	     it=context.delta_bag_lookup[pos_res].begin();
	   it != context.delta_bag_lookup[pos_res].end(); it++) {
	const int db_index = *it;
	if (db_remaining[db_index] < 1)
	  continue;
	db_remaining[db_index] -= 1;
	mass_list[i] = save_pos_mass + context.delta_bag_delta[db_index];
	mtl.item.delta = context.delta_bag_delta[db_index];
	choose_residue_mod(context, m, &mtl, mass_list,
			   candidate_spectra_info_begin,
			   candidate_spectra_info_end, stats, db_remaining,
			   remaining_positions_to_choose-1, i+1);
	db_remaining[db_index] += 1;
      }
      mass_list[i] = save_pos_mass;
    }
  }
}


// p_begin (p_end) to begin_index (end_index), updating p_mass in the process.
// The sign determines whether mass increases or decreases as p_begin (p_end)
// is moved forward--so it should be -1 for the begin case and +1 for the end
// case.
// FIX: is this worth its complexity?  kill or else add reset
static inline void
update_p_mass(double &p_mass, int &p_begin, int begin_index, int sign,
	      const std::string &run_sequence,
	      const std::vector<double> &fixed_residue_mass) {
  assert(sign == +1 or sign == -1);
  assert(begin_index >= 0);
  const bool moving_forward = begin_index >= p_begin;
  int p0=p_begin, p1=begin_index;
  if (not moving_forward) {
    std::swap(p0, p1);
    sign = -sign;
  }
  for (int i=p0; i<p1; i++)
    p_mass += sign * fixed_residue_mass[run_sequence[i]];
  p_begin = begin_index;
}



// Search a sequence run for matches according to the context against the
// spectra.  Updates score_stats and the number of candidate spectra found.

// FIX: examine carefully for signed/unsigned problems

void inline
search_run(const search_context &context, const sequence_run &sequence_run,
	   score_stats &stats) {
  const parameters &CP = parameters::the;
  const int min_peptide_length = std::max<int>(CP.minimum_peptide_length,
					       context.mod_count);
  const std::vector<double> &fixed_parent_mass \
    = CP.parent_mass_regime[context.mass_regime_index].fixed_residue_mass;

  const std::string &run_sequence = sequence_run.sequence;
  const std::vector<int> &cleavage_points = sequence_run.cleavage_points;

  // This match will be passed inward and used to record information that we
  // need to remember about a match when we finally see one.  At that point, a
  // copy of this match will be saved.
  match m;			// FIX: move inward?
  m.sequence_index = sequence_run.sequence_index;
  m.sequence_offset = sequence_run.sequence_offset;

  assert(context.delta_bag_delta.size()
	 == context.delta_bag_count.size());
  // counts remaining as mod positions are chosen
  std::vector<int> db_remaining = context.delta_bag_count;

  // the rightmost 'end' seen when all spectra masses were too high
  // (0 means none yet encountered.)
  unsigned int next_end = 0;

  // p_mass is the parent mass of the peptide run_sequence[p_begin:p_end]
  double p_mass = context.parent_fixed_mass;
  int p_begin=0, p_end=0;

  // FIX: optimize the non-specific cleavage case?
  for (unsigned int begin=0; begin<cleavage_points.size()-1; begin++) {
    const int begin_index = m.peptide_begin = cleavage_points[begin];
    if (not context.pca_residues.empty())
      if (context.pca_residues.find(run_sequence[begin_index])
	  == std::string::npos)
	continue;
    update_p_mass(p_mass, p_begin, begin_index, -1, run_sequence,
		  fixed_parent_mass);
    unsigned int end = std::max<unsigned int>(begin + 1, next_end);
    for (; end<cleavage_points.size(); end++) {
      m.missed_cleavage_count = end - begin - 1;
      if (m.missed_cleavage_count > context.maximum_missed_cleavage_sites)
	break;

      const int end_index = cleavage_points[end];
      assert(end_index - begin_index > 0);
      const int peptide_size = end_index - begin_index;
      assert(peptide_size > 0);
      if (peptide_size < min_peptide_length)
	continue;
      update_p_mass(p_mass, p_end, end_index, +1, run_sequence,
		    fixed_parent_mass);
      //std::cerr << "peptide: " << std::string(run_sequence, begin_index, peptide_size) << " p_mass: " << p_mass << std::endl;


      double sp_mass_lb = p_mass + CP.parent_monoisotopic_mass_error_minus;
      double sp_mass_ub = p_mass + CP.parent_monoisotopic_mass_error_plus;

      const spmi_c_it candidate_spectra_info_begin
	= spectrum::spectrum_mass_index.lower_bound(sp_mass_lb);
      if (candidate_spectra_info_begin == spectrum::spectrum_mass_index.end()) {
	// spectrum masses all too low to match peptide (peptide too long)
	//std::cerr << "peptide: sp low" << std::endl;
	break;
      }

      const spmi_c_it candidate_spectra_info_end
	= spectrum::spectrum_mass_index.upper_bound(sp_mass_ub);
      if (candidate_spectra_info_end == spectrum::spectrum_mass_index.begin()) {
	// spectrum masses all too high to match peptide
	// (peptide is too short, so increase end, if possible, else return)
	//std::cerr << "peptide: sp high" << std::endl;
	if (end == cleavage_points.size() - 1)
	  return;
	next_end = end;
	continue;
      }
      if (candidate_spectra_info_begin == candidate_spectra_info_end) {
	// no spectrum with close-enough parent mass
	//std::cerr << "peptide: sp none nearby" << std::endl;
	continue;
      }

      m.parent_peptide_mass = p_mass;
      //m.fragment_peptide_mass = ???;  unneeded???

      m.peptide_sequence.assign(run_sequence, begin_index, peptide_size);

      std::vector<double> mass_list(peptide_size);
      for (std::vector<double>::size_type i=0; i<m.peptide_sequence.size();
	   i++)
	mass_list[i] = (CP.fragment_mass_regime[context.mass_regime_index]
			.fixed_residue_mass[m.peptide_sequence[i]]);

      choose_residue_mod(context, m, NULL, mass_list,
			 candidate_spectra_info_begin,
			 candidate_spectra_info_end, stats, db_remaining,
			 context.mod_count, 0);
    }
  }
}


// Search all sequence runs for matches according to the context against the
// spectra.  Updates score_stats and the number of candidate spectra found.
void
spectrum::search_runs(const search_context &context, score_stats &stats) {
  const parameters &CP = parameters::the;

  const int no_runs = context.sequence_runs.size();
  for (int i=0; i<no_runs; i++) {
    search_run(context, context.sequence_runs[i], stats);

    if (CP.show_progress)
      std::cerr << i+1 << " of " << no_runs << " sequences, "
		<< stats.candidate_spectrum_count << " cand for "
		<< stats.spectra_with_candidates << " sp, "
		<< stats.improved_candidates << "++\r" << std::flush;
  }

  if (CP.show_progress)
    std::cerr << std::endl;
}
