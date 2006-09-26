

//	$Id$


#include "cxtpy.hpp"


#include <cerrno>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include <cmath>
#include <stdexcept>

#include <iostream>


//#define FAST_LOWER_READ_ACCURACY

#ifdef __GNU__
#define FAST_GNU_ONLY
#endif

#ifdef FAST_GNU_ONLY
// use non-threadsafe versions, since we don't use threads
#define fgetsX fgets_unlocked
#define ferrorX ferror_unlocked
#else
#define fgetsX std::fgets
#define ferrorX std::ferror
#endif

#ifdef FAST_LOWER_READ_ACCURACY
// probably faster, at some cost in accuracy
#define strtodX std::strtof
#else
#define strtodX std::strtod
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


// LATER
//bool
//spectrum::generate_mod_patterns(potential_mod_pattern, peptide_seq, begin) {
//}


// Set peaks to the synthesized mass (not mz) ladder.
void
spectrum::set_synthetic_Y_peaks(const std::string &peptide_seq,
				const std::vector<double> *potential_mod_pattern,
				const unsigned mass_regime) {
  const parameters &CP = parameters::the;
  const bool is_C = false;	// FIX
  double m = (CP.fragment_mass_regime[mass_regime].water_mass
	      + (CP.cleavage_C_terminal_mass_change
		 - CP.fragment_mass_regime[mass_regime].hydroxyl_mass));
  if (is_C) {
    m += CP.fragment_mass_regime[mass_regime].modification_mass[']'];
    //if False:                       # ] is diff modded
    //        #m += CP.potential_modification_mass[ord(']')]...
  }

  const int ladder_size = peptide_seq.size() - 1;
  peaks.resize(ladder_size);
  for (int i=ladder_size-1; i>=0; i--) {
    m += CP.fragment_mass_regime[mass_regime].residue_mass[peptide_seq[i+1]];
    // m += the mod delta, too
    peaks[ladder_size-1-i] = peak(m, 1.0);
  }
  assert(peptide_seq.size() >= 3);
  peaks[ladder_size-2].intensity = 3.0;
  if (peptide_seq[1] == 'P')
    peaks[ladder_size-2].intensity = 10.0;
}

//  0 * 2 3 (peaks)
// A P C D E

// Set peaks to the synthesized mass (not mz) ladder.
void
spectrum::set_synthetic_B_peaks(const std::string &peptide_seq,
				const std::vector<double> *potential_mod_pattern,
				const unsigned mass_regime) {
  const parameters &CP = parameters::the;
  const bool is_N = false;	// FIX
  double m = 0.0 + (CP.cleavage_N_terminal_mass_change - CP.hydrogen_mass);
  if (is_N) {
    m += CP.fragment_mass_regime[mass_regime].modification_mass['['];
    //if False:                       # [ is diff modded
    //        #m += CP.potential_modification_mass[ord('[')]...
  }

  const int ladder_size = peptide_seq.size() - 1;
  peaks.resize(ladder_size);
  for (int i=0; i<=ladder_size-1; i++) {
    m += CP.fragment_mass_regime[mass_regime].residue_mass[peptide_seq[i]];
    // m += the mod delta, too
    peaks[i] = peak(m, 1.0);
  }
  assert(peptide_seq.size() >= 2);
  peaks[1].intensity = 3.0;
  if (peptide_seq[1] == 'P')
    peaks[1].intensity = 10.0;
}


// Multiply peak intensities by residue-dependent factors to generate a more
// realistic spectrum.
void
spectrum::synthesize_peak_intensities(const std::string &peptide_seq) {
  assert(peaks.size() + 1 == peptide_seq.size());

  for (unsigned int i=0; i<peaks.size(); i++) {
    switch (peptide_seq[i]) {
    case 'D':
      peaks[i].intensity *= 5.0;
      break;
    case 'V':
    case 'E':
    case 'I':
    case 'L':
      peaks[i].intensity *= 3.0;
      break;
    case 'N':
    case 'Q':
      peaks[i].intensity *= 2.0;
      break;
    }
    switch (peptide_seq[i+1]) {
    case 'P':
      peaks[i].intensity *= 5.0;
      break;
    }
  }
}


// Construct a synthetic spectrum from these parameters.
// 'potential_mod_pattern' may be null if there are no mods.
spectrum::spectrum(ion ion_type, int charge, const std::string &peptide_seq,
		   const std::vector<double> *peptide_mod_pattern,
		   double peptide_mod_mass, bool synthesize_intensities) {
  init(peptide_mod_mass, charge);

  switch (ion_type) {
  case ION_Y:
    set_synthetic_Y_peaks(peptide_seq, peptide_mod_pattern);
    break;
  case ION_B:
    set_synthetic_B_peaks(peptide_seq, peptide_mod_pattern);
    break;
  default:
    assert(false);
  }

  if (synthesize_intensities)
    synthesize_peak_intensities(peptide_seq);

  //std::cerr << "synspect" << std::endl;
  //for (std::vector<peak>::iterator it=peaks.begin(); it != peaks.end(); it++)
  //  std::cerr << "synpeak: " << it->mz << " " << charge << std::endl;

  // FIX: move outward?  B1 and B2 are the same until this point?
  // FIX: move this into set_synth* functions?
  for (std::vector<peak>::iterator it=peaks.begin(); it != peaks.end(); it++)
    it->mz = peak::get_mz(it->mz, charge);
}


static void
synthetic_spectra(spectrum synth_sp[/* max_fragment_charge+1 */][ION_MAX],
		  const std::string &peptide_seq,
		  const std::vector<double> &peptide_mod_pattern,
		  double peptide_mod_mass, double max_fragment_charge) {
  const parameters &CP = parameters::the;

  for (ion ion_type=ION_MIN; ion_type<ION_MAX; ion_type++)
    for (int charge=1; charge<=max_fragment_charge; charge++)
      synth_sp[charge][ion_type] = spectrum(ion_type, charge, peptide_seq,
					    &peptide_mod_pattern,
					    peptide_mod_mass,
					    CP.spectrum_synthesis);
}


// IDEA: could we combine all these spectra and just do the correlation once?

// Score the specified spectrum against a group of synthetic fragment
// spectra.
static void
score_spectrum(double &hyper_score, double &convolution_score,
	       std::vector<double> &ion_scores, std::vector<int> &ion_peaks,
	       spectrum synth_sp[][ION_MAX], int spectrum_id,
	       int max_fragment_charge) {
  const parameters &CP = parameters::the;
  const spectrum sp = spectrum::searchable_spectra[spectrum_id];
  const int spectrum_charge = sp.charge;
  hyper_score = 1.0;
  convolution_score = 0;

  // FIX: verify that we only access the ion_peaks/ion_scores slots we store
  // here?
  for (ion ion_type=ION_MIN; ion_type<ION_MAX; ion_type++) {
    int i_peaks = 0;
    double i_scores = 0.0;
    for (int charge=1; charge<=max_fragment_charge; charge++) {
      // FIX: move this test outward???
      if (spectrum_charge == 1 and charge > spectrum_charge
	  or spectrum_charge > 1 and charge > spectrum_charge - 1)
	continue;

      // convolution score is just sum over charges/ions
      // hyperscore is product of p! over charges/ions (where p is corr peak
      // count) times the convolution score (clipped to FLT_MAX)
      // > blurred!
      int common_peak_count;	// FIX!
      double conv_score = spectrum::score_similarity(synth_sp[charge][ion_type],
						     sp, &common_peak_count);
      i_peaks += common_peak_count;
      i_scores += conv_score;
      hyper_score *= CP.factorial[common_peak_count];
      convolution_score += conv_score;
    }

    ion_peaks[ion_type] = i_peaks;
    ion_scores[ion_type] = i_scores;

    // want to know:
    // - total peaks (over charges) for each ion_type
    // - total convolution score (over charges) for each ion_type
    // - maybe grand total convolution score
    // - hyper score: product of all factorial(peak count) times
    //                grand total convolution score
    //   [xtandem clips hyper score to FLTMAX]

  }
  hyper_score *= convolution_score;
}


#if 0

struct generator;


struct generator {
  // fx = &generator::generate_static_aa_mod;

  typedef void generator::generator_fx(score_stats &stats, const match &m,
				       std::vector<double> position_mass,
				       const std::vector<double> &terminal_mass,
				       const std::vector<generator>::iterator g_it);

  generator() : fx(0) { }

  generator_fx generate_static_aa_mod;

  generator_fx generator::*fx;
  std::vector<double> args;
};


void
generator::generate_static_aa_mod(score_stats &stats, const match &m,
				  std::vector<double> position_mass,
				  const std::vector<double> &terminal_mass,
				  const std::vector<generator>::iterator g_it) {
  for (unsigned int i=0; i<m.peptide_sequence.size(); i++)
    position_mass[i] += args[m.peptide_sequence[i]];

  ((*g_it).*(g_it->fx))(stats, m, position_mass, terminal_mass, g_it+1);
}

#endif		    



// Search for matches of this peptide (including all specified modification)
// against the spectra.  Updates score_stats and returns the number of
// candidate spectra found for this peptide.
int 
spectrum::search_peptide(int idno, int offset, int begin,
			 const std::string &peptide_seq,
			 int missed_cleavage_count, score_stats &stats) {
  const parameters &CP = parameters::the;
  int peptide_candidate_spectrum_count = 0;
  double peptide_mass = get_parent_peptide_mass(peptide_seq);

  // pyrroline?
  std::vector<double> potential_mod_pattern(peptide_seq.size());
  //while (generate_mod_patterns(potential_mod_pattern, peptide_seq, begin))
  bool once = true;		// placeholder
  while (once) {
    once = false;
    double peptide_mod_mass = peptide_mass; // + get_peptide_mod_mass(peptide_seq,
                                            //              potential_mod_pattern,
                                            //                         is_N, is_C)
    double sp_mass_lb = (peptide_mod_mass
			 + CP.parent_monoisotopic_mass_error_minus);
    double sp_mass_ub = (peptide_mod_mass
			 + CP.parent_monoisotopic_mass_error_plus);
    std::multimap<double,
      std::vector<spectrum>::size_type>::const_iterator
      candidate_spectra_info_begin = spectrum_mass_index.lower_bound(sp_mass_lb);
    std::multimap<double,
      std::vector<spectrum>::size_type>::const_iterator
      candidate_spectra_info_end = spectrum_mass_index.upper_bound(sp_mass_ub);
    if (candidate_spectra_info_begin == candidate_spectra_info_end)
      continue;

    int max_candidate_charge = 0;
    for (std::multimap<double, std::vector<spectrum>::size_type>::const_iterator
	   it = candidate_spectra_info_begin;
	 it != candidate_spectra_info_end; it++)
      max_candidate_charge = std::max<int>(max_candidate_charge,
					   searchable_spectra[it->second].charge);

    const int max_fragment_charge = std::max<int>(1, max_candidate_charge-1);
    // e.g. +2 -> 'B' -> spectrum
    // FIX: should this be a std::vector??
    spectrum synth_sp[max_fragment_charge+1][ION_MAX];
    synthetic_spectra(synth_sp, peptide_seq, potential_mod_pattern,
		      peptide_mod_mass, max_fragment_charge);
    for (std::multimap<double, std::vector<spectrum>::size_type>::const_iterator
	   candidate_it = candidate_spectra_info_begin;
	 candidate_it != candidate_spectra_info_end; candidate_it++) {
      const int spectrum_id = candidate_it->second;
      peptide_candidate_spectrum_count++;
      double hyper_score, convolution_score;
      std::vector<double> ion_scores(ION_MAX, 0.0);
      std::vector<int> ion_peaks(ION_MAX, 0);
      score_spectrum(hyper_score, convolution_score, ion_scores, ion_peaks,
		     synth_sp, spectrum_id, max_fragment_charge);
      if (convolution_score <= 2.0)
	continue;

      // update spectrum histograms
      std::vector<int> &hh = stats.hyperscore_histogram[spectrum_id];
      int binned_scaled_hyper_score = int(0.5 + scale_hyperscore(hyper_score));
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
      const int sp_ion_count = ion_peaks[ION_B] + ion_peaks[ION_Y];
      if (sp_ion_count < CP.minimum_ion_count)
	continue;
      const bool has_b_and_y = ion_peaks[ION_B] and ion_peaks[ION_Y];
      if (not has_b_and_y)
	continue;

      // check that parent masses are within error range (isotope ni)
      // already done above (why does xtandem do it here?)
      // if check fails, only eligible for 2nd-best record

      // Remember all of the highest-hyper-scoring matches against each
      // spectrum.  These might be in multiple domains (pos, length, mods) in
      // multiple proteins.

      // To avoid the problems that xtandem has with inexact float-point
      // calculations, we consider hyperscores within the "epsilon ratio" to
      // be "equal".  The best_match vector will need to be re-screened
      // against this ratio, because it may finally contain entries which
      // don't meet our criterion.  (We don't take the time to get it exactly
      // right here, because we may end up discarding it all anyway.)

      const double current_ratio = hyper_score / stats.best_score[spectrum_id];

      if (CP.quirks_mode and (hyper_score < stats.best_score[spectrum_id])
	  or (not CP.quirks_mode
	      and (current_ratio < CP.hyper_score_epsilon_ratio))) {
	// This isn't a best score, so see if it's a second-best score.
	if (hyper_score > stats.second_best_score[spectrum_id])
	  stats.second_best_score[spectrum_id] = hyper_score;
	continue;
      }
      if (CP.quirks_mode and (hyper_score > stats.best_score[spectrum_id])
	  or (not CP.quirks_mode
	      and (current_ratio > (1/CP.hyper_score_epsilon_ratio)))) {
	// This score is significantly better, so forget previous matches.
	stats.best_score[spectrum_id] = hyper_score;
	stats.best_match[spectrum_id].resize(1);
      }	else {
	// This score is similar to best previously seen, so add it to the
	// list.
	stats.best_match[spectrum_id].push_back(match());
      }
      match &m = stats.best_match[spectrum_id].back();
      // remember this best match
      // FIX: better to create one of these at the top!(?)
      m.hyper_score = hyper_score;
      m.convolution_score = convolution_score;
      m.ion_scores = ion_scores;
      m.ion_peaks = ion_peaks;
      m.missed_cleavage_count = missed_cleavage_count;
      m.spectrum_index = spectrum_id;
      m.peptide_begin = begin;
      m.peptide_sequence = peptide_seq;
      m.potential_mod_pattern = potential_mod_pattern;
      m.peptide_mod_mass = peptide_mod_mass;
      m.sequence_index = idno;
      m.sequence_offset = offset;
    }
  }  
  return peptide_candidate_spectrum_count;
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


double
get_parent_peptide_mass(const std::string &peptide_seq,
			const unsigned mass_regime) {
  const parameters &CP = parameters::the;
  const bool is_N=false, is_C=false; // FIX
  double NC_mods = 0;
  if (is_N)
    NC_mods += CP.parent_mass_regime[mass_regime].modification_mass['['];
  if (is_C)
    NC_mods += CP.parent_mass_regime[mass_regime].modification_mass[']'];

  double residue_sum = 0;
  for (std::string::const_iterator it=peptide_seq.begin();
       it != peptide_seq.end(); it++)
    residue_sum += (CP.parent_mass_regime[mass_regime].residue_mass[*it]
		    + CP.parent_mass_regime[mass_regime].modification_mass[*it]);
  return (NC_mods + residue_sum
	  + CP.cleavage_N_terminal_mass_change
	  + CP.cleavage_C_terminal_mass_change
	  + CP.proton_mass);
}


// returns log-scaled hyperscore
double
scale_hyperscore(double hyper_score) {
  assert(hyper_score >= 0);	// otherwise check for EDOM, ERANGE
  if (hyper_score == 0)
    return -DBL_MAX;		// FIX: this shouldn't be reached?
  return 4 * std::log10(hyper_score);
}

