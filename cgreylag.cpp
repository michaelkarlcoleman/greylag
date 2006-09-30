

//	$Id$


#include "cgreylag.hpp"


#include <cerrno>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include <cmath>
#include <numeric>
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


// Multiply peak intensities by residue-dependent factors to generate a more
// realistic spectrum.  If is_XYZ, walk through the peptide sequence in
// reverse.
static inline void
synthesize_ladder_intensities(std::vector<peak> &mass_ladder,
			      const std::string &peptide_seq,
			      const bool is_XYZ) {
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
			    const double C_terminal_mass,
			    const unsigned mass_regime=0) {
  const parameters &CP = parameters::the;
  assert(peptide_seq.size() == mass_list.size());
  double m = (CP.fragment_mass_regime[mass_regime].water_mass
	      + (CP.cleavage_C_terminal_mass_change
		 - CP.fragment_mass_regime[mass_regime].hydroxyl_mass)
	      + C_terminal_mass);

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
			    const double N_terminal_mass,
			    const unsigned mass_regime=0) {
  const parameters &CP = parameters::the;
  double m = (0.0 + (CP.cleavage_N_terminal_mass_change - CP.hydrogen_mass)
	      + N_terminal_mass);

  const int ladder_size = mass_ladder.size();
  for (int i=0; i<=ladder_size-1; i++) {
    m += mass_list[i];
    mass_ladder[i] = peak(m, 1.0);
  }

  synthesize_ladder_intensities(mass_ladder, peptide_seq, false);
}


// Construct a synthetic spectrum from these parameters.
spectrum::spectrum(double peptide_mod_mass, int charge,
		   const std::vector<peak> &mass_ladder) {
  init(peptide_mod_mass, charge);

  peaks = mass_ladder;
  for (std::vector<peak>::size_type i=0; i<mass_ladder.size(); i++)
    peaks[i].mz = peak::get_mz(peaks[i].mz, charge);

  //std::cerr << "synth charge " << charge << std::endl;
  //for (unsigned int i=0; i<mass_ladder.size(); i++)
  //  std::cerr << "peak mz " << peaks[i].mz << " i " << peaks[i].intensity << std::endl;
}


static inline void
synthetic_spectra(spectrum synth_sp[/* max_fragment_charge+1 */][ION_MAX],
		  const std::string &peptide_seq, const double peptide_mass,
		  const std::vector<double> &mass_list,
		  const double N_terminal_mass, const double C_terminal_mass,
		  const double max_fragment_charge) {
  for (ion ion_type=ION_MIN; ion_type<ION_MAX; ion_type++) {
    std::vector<peak> mass_ladder(mass_list.size()-1);

    switch (ion_type) {
    case ION_Y:
      get_synthetic_Y_mass_ladder(mass_ladder, peptide_seq, mass_list,
				  C_terminal_mass);
      break;
    case ION_B:
      get_synthetic_B_mass_ladder(mass_ladder, peptide_seq, mass_list,
				  N_terminal_mass);
      break;
    default:
      assert(false);
    }

    //for (unsigned int i=0; i<mass_ladder.size(); i++)
    //  std::cerr << "ladder step " << peptide_seq[i] << " " << mass_ladder[i].mz << std::endl;

    for (int charge=1; charge<=max_fragment_charge; charge++)
      synth_sp[charge][ion_type] = spectrum(peptide_mass, charge,
					    mass_ladder);
  }
}


// IDEA: could we combine all these spectra and just do the correlation once?

// Score the specified spectrum against a group of synthetic fragment
// spectra.
static inline void
score_spectrum(double &hyper_score, double &convolution_score,
	       std::vector<double> &ion_scores, std::vector<int> &ion_peaks,
	       spectrum synth_sp[/* max_fragment_charge+1 */][ION_MAX],
	       const int spectrum_id, const int max_fragment_charge) {
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
      double conv_score = spectrum::score_similarity(synth_sp[charge][ion_type],
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


// FIX: currently only one mass_list is implemented, meaning that parent and
// fragment masses must be the same (e.g., mono/mono).  A workaround is to
// just increase the parent error plus a bit (the xtandem default seems to
// already do this).

// Search for matches of this particular peptide modification variation
// against the spectra.  Updates score_stats and returns the number of
// candidate spectra found.
static inline void
evaluate_peptide_mod_variation(match &m,
			       const std::vector<double> &mass_list,
			       const double N_terminal_mass,
			       const double C_terminal_mass,
			       score_stats &stats) {
  const parameters &CP = parameters::the;
  assert(m.peptide_sequence.size() == mass_list.size());

  m.peptide_mass = get_peptide_mass(mass_list, N_terminal_mass,
				    C_terminal_mass);
  double sp_mass_lb = m.peptide_mass + CP.parent_monoisotopic_mass_error_minus;
  double sp_mass_ub = m.peptide_mass + CP.parent_monoisotopic_mass_error_plus;
  std::multimap<double, std::vector<spectrum>::size_type>::const_iterator
    candidate_spectra_info_begin
    = spectrum::spectrum_mass_index.lower_bound(sp_mass_lb);
  // FIX: is there a better way?
  std::multimap<double, std::vector<spectrum>::size_type>::const_iterator
    candidate_spectra_info_end
    = spectrum::spectrum_mass_index.upper_bound(sp_mass_ub);
  if (candidate_spectra_info_begin == candidate_spectra_info_end)
    return;

  int max_candidate_charge = 0;
  for (std::multimap<double, std::vector<spectrum>::size_type>::const_iterator
	 it = candidate_spectra_info_begin;
       it != candidate_spectra_info_end; it++)
    max_candidate_charge
      = std::max<int>(max_candidate_charge,
		      spectrum::searchable_spectra[it->second].charge);

  const int max_fragment_charge = std::max<int>(1, max_candidate_charge-1);
  assert(max_fragment_charge <= spectrum::max_supported_charge);
  // e.g. +2 -> 'B' -> spectrum
  // FIX: should this be a std::vector??
  static spectrum synth_sp[spectrum::max_supported_charge+1][ION_MAX];
  synthetic_spectra(synth_sp, m.peptide_sequence, m.peptide_mass, mass_list,
		    N_terminal_mass, C_terminal_mass, max_fragment_charge);

  for (std::multimap<double, std::vector<spectrum>::size_type>::const_iterator
	 candidate_it = candidate_spectra_info_begin;
       candidate_it != candidate_spectra_info_end; candidate_it++) {
    m.spectrum_index = candidate_it->second;
    stats.candidate_spectrum_count++;
    score_spectrum(m.hyper_score, m.convolution_score, m.ion_scores,
		   m.ion_peaks, synth_sp, m.spectrum_index,
		   max_fragment_charge);
    if (m.convolution_score <= 2.0)
      continue;

    //std::cerr << "histod: " << m.spectrum_index << " " << m.peptide_sequence
    //<< " " << m.peptide_mass << " " <<
    //spectrum::searchable_spectra[m.spectrum_index].mass << std::endl;

    // update spectrum histograms
    std::vector<int> &hh = stats.hyperscore_histogram[m.spectrum_index];
    int binned_scaled_hyper_score = int(0.5 + scale_hyperscore(m.hyper_score));
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
      stats.best_score[m.spectrum_index] = m.hyper_score;
      stats.best_match[m.spectrum_index].clear();
    }
    // Now remember this match because it's at least as good as those
    // previously seen.
    stats.best_match[m.spectrum_index].push_back(m);
  }
}


// Two separate functions for N- and C-terminus potential mods.  The former
// should include the PCA possibility.

// Choose a possible modification to account for PCA (pyrrolidone carboxyl
// acid) circularization of the peptide N-terminal.  This is excluded if a
// static N-terminal mod was specified.  (FIX: Does a PCA mod exclude other
// potential mods?)
static inline void
choose_potential_mod(match &m, std::vector<double> &mass_list,
		     const double N_terminal_mass, const double C_terminal_mass,
		     score_stats &stats) {
  const parameters &CP = parameters::the;
  evaluate_peptide_mod_variation(m, mass_list, N_terminal_mass,
				 C_terminal_mass, stats);

  const int mass_regime=0;	// FIX

}


// Choose a possible modification to account for PCA (pyrrolidone carboxyl
// acid) circularization of the peptide N-terminal.  This is excluded if a
// static N-terminal mod was specified.  (FIX: Does a PCA mod exclude other
// potential mods?)
static inline void
choose_PCA_mod(match &m, std::vector<double> &mass_list,
	       const double N_terminal_mass, const double C_terminal_mass,
	       score_stats &stats) {
  const parameters &CP = parameters::the;
  evaluate_peptide_mod_variation(m, mass_list, N_terminal_mass,
				 C_terminal_mass, stats);

  const int mass_regime=0;	// FIX

  if (CP.fragment_mass_regime[mass_regime].modification_mass['['] == 0) {
    for (int bug=0; bug<(CP.quirks_mode? 2 : 1); bug++)	// FIX
      switch (m.peptide_sequence[0]) {
      case 'E':
	evaluate_peptide_mod_variation(m, mass_list,
				       N_terminal_mass
				       -CP.fragment_mass_regime[mass_regime].water_mass,
				       C_terminal_mass, stats);
	break;
      case 'C': {
	// FIX: need better test for C+57? (symbolic?)
	const double C_mod
	  = CP.fragment_mass_regime[mass_regime].modification_mass['C'];
	if (CP.quirks_mode ? int(C_mod) != 57 : std::abs(C_mod - 57) > 0.5)
	  break;		// skip unless C+57
      }
      case 'Q':
	evaluate_peptide_mod_variation(m, mass_list,
				       N_terminal_mass
				       -CP.fragment_mass_regime[mass_regime].ammonia_mass,
				       C_terminal_mass, stats);
	break;
      }
  }
}


// Choose among the requested mass regimes (e.g., isotopes), and also choose
// the static mods (including N- and C-terminal mods), which are determined by
// the mass regime choice.
static inline void
choose_mass_regime(match &m, std::vector<double> &mass_list,
		   const double N_terminal_mass, const double C_terminal_mass,
		   score_stats &stats) {
  const parameters &CP = parameters::the;

  const int mass_regime=0;	// FIX
  for (std::vector<double>::size_type i=0; i<m.peptide_sequence.size(); i++) {
    const char res = m.peptide_sequence[i];
    mass_list[i] = (CP.fragment_mass_regime[mass_regime].residue_mass[res]
		    + CP.fragment_mass_regime[mass_regime].modification_mass[res]);
  }
  choose_PCA_mod(m, mass_list,
		 N_terminal_mass
		 + CP.fragment_mass_regime[mass_regime].modification_mass['['],
		 C_terminal_mass
		 + CP.fragment_mass_regime[mass_regime].modification_mass[']'],
		 stats);
}


// Search for matches of all modification variations of this peptide against
// the spectra.  Updates score_stats and returns the number of candidate
// spectra found.
void
spectrum::search_peptide_all_mods(int idno, int offset, int begin,
				  const std::string &peptide_seq,
				  int missed_cleavage_count,
				  score_stats &stats) {
  std::vector<double> mass_list(peptide_seq.size());
  double N_terminal_mass = 0.0;
  double C_terminal_mass = 0.0;

  // This match will be passed inward and used to record information that we
  // need to remember about a match when we finally see one.  At that point, a
  // copy of this match will be saved.
  match m;
  m.sequence_index = idno;
  m.sequence_offset = offset;
  m.peptide_begin = begin;
  m.peptide_sequence = peptide_seq;
  m.missed_cleavage_count = missed_cleavage_count;

  choose_mass_regime(m, mass_list, N_terminal_mass, C_terminal_mass, stats);
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


// FIX: obsolete?
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


double
get_peptide_mass(const std::vector<double> &mass_list,
		 const double N_terminal_mass, const double C_terminal_mass) {
  const parameters &CP = parameters::the;

  double m = (N_terminal_mass + C_terminal_mass
	      + CP.cleavage_N_terminal_mass_change
	      + CP.cleavage_C_terminal_mass_change
	      + CP.proton_mass);
  return std::accumulate(mass_list.begin(), mass_list.end(), m);
}


// returns log-scaled hyperscore
double
scale_hyperscore(double hyper_score) {
  assert(hyper_score >= 0);	// otherwise check for EDOM, ERANGE
  if (hyper_score == 0)
    return -DBL_MAX;		// FIX: this shouldn't be reached?
  return 4 * std::log10(hyper_score);
}

