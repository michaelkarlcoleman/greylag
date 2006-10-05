// cgreylag.hpp

//	$Id$

#ifndef CGREYLAG_H
#define CGREYLAG_H

#include <cassert>
#include <cstdio>
#include <list>
#include <map>
#include <vector>
#include <stdexcept>
#include <string>
#include <utility>


class score_stats;


class mass_regime_parameters {
public:
  double hydroxyl_mass;
  double water_mass;
  double ammonia_mass;

  // these are indexed by residue char (plus '[' and ']' for N and C termini)
  std::vector<double> residue_mass;
  std::vector<double> modification_mass;
  std::vector< std::vector<double> > potential_modification_mass;
  std::vector< std::vector<double> > potential_modification_mass_refine;

  mass_regime_parameters(unsigned size=0) : hydroxyl_mass(0), water_mass(0),
					    ammonia_mass(0) {
    residue_mass.resize(size);
    modification_mass.resize(size);
    potential_modification_mass.resize(size);
    potential_modification_mass_refine.resize(size);
  }
};

// FIX: want these class members to all be static, but had SWIG trouble
class parameters {
public:
  bool quirks_mode;		// try to produce results closer to xtandem's

  // hyperscores within this ratio are considered "equal"
  double hyper_score_epsilon_ratio;

  std::vector<double> factorial; // factorial[n] == n!

  // deuterium not yet implemented
  double proton_mass;
  double hydrogen_mass;
  std::vector<mass_regime_parameters> parent_mass_regime;
  std::vector<mass_regime_parameters> fragment_mass_regime;

  // from XTP
  double cleavage_N_terminal_mass_change;
  double cleavage_C_terminal_mass_change;

  double parent_monoisotopic_mass_error_plus;
  double parent_monoisotopic_mass_error_minus;
  double fragment_mass_error;
  int minimum_ion_count;
  bool spectrum_synthesis;

  static parameters the;
};


enum ion { ION_MIN, ION_Y=ION_MIN, ION_B, ION_MAX }; // NB: [ION_MIN, ION_MAX)

inline void operator++(ion &i) { i = ion(i + 1); }
inline void operator++(ion &i, int) { i = ion(i + 1); }


class peak {
 public:
  double mz;
  double intensity;

  peak(double mz=0, double intensity=0) : mz(mz), intensity(intensity) { }

  char *__repr__() const;

  static bool less_mz(peak x, peak y) { return x.mz < y.mz; }
  static bool less_intensity(peak x, peak y) {
    return x.intensity < y.intensity;
  }

  static double get_mz(double mass, int charge) {
    const parameters &CP = parameters::the;
    return mass/charge + CP.proton_mass;
  }
};


// A spectrum searched on multiple charges will be turned into separate
// spectra having differing id's, but the same physical id.

class spectrum {
public:
  double mass;
  int charge;
  static const int max_supported_charge = 10;
  std::vector<peak> peaks;
  double max_peak_intensity;
  double sum_peak_intensity;
  double normalization_factor;
  std::string name;
  int file_id;			// file # of spectra's source
  int id;			// unique for all spectra
  int physical_id;

private:
  void init(double mass_=0, int charge_=0) {
    mass = mass_;
    charge = charge_;
    if (charge > max_supported_charge)
      throw std::invalid_argument("attempt to create a spectrum with greater"
				  " than supported charge");
    set_id();
    max_peak_intensity = -1;
    sum_peak_intensity = -1;
    normalization_factor = 1;
    file_id = -1;
    physical_id = -1;
  }

public:
  // Construct an empty spectrum.
  spectrum(double mass=0, int charge=0) { init(mass, charge); }

  // Construct a (synthetic) spectrum from these parameters.
  spectrum(double peptide_mod_mass, int charge,
	   const std::vector<peak> &mass_ladder) {
    init(peptide_mod_mass, charge);
    peaks = mass_ladder;
    for (std::vector<peak>::size_type i=0; i<mass_ladder.size(); i++)
      peaks[i].mz = peak::get_mz(peaks[i].mz, charge);
  }

  char *__repr__() const;

  // for tinkering
  void set_peaks_from_matrix(const std::vector< std::vector<double> > &m) {
    peaks.resize(m.size());
    std::vector<peak>::iterator p_it = peaks.begin();
    for (std::vector< std::vector<double> >::const_iterator it = m.begin();
	 it != m.end(); it++, p_it++) {
      if (it->size() != 2)
	throw std::invalid_argument("invalid matrix (must be size N x 2)");
      p_it->mz = (*it)[0];
      p_it->intensity = (*it)[1];
    }
  }

  // Read spectra from file in ms2 format, tagging them with file_id.  Only
  // spectra with parent masses in the range [min_mass, max_mass) are
  // included.  (All of them are read, though, so the ids for any spectrum
  // will be the same regardless of the range used.)
  static std::vector<spectrum> read_spectra(FILE *f, const int file_id,
					    const double min_mass,
					    const double max_mass);

  // Sets max/sum_peak_intensity, according to peaks and
  // normalization_factor.
  void calculate_intensity_statistics() {
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


  // Examine this spectrum to see if it should be kept.  If so, return true
  // and sort the peaks by mz, normalize them and limit their number.
  bool filter_and_normalize(double minimum_fragment_mz, double dynamic_range,
			    int minimum_peaks, int total_peaks);

  // Store the list of spectra that search_peptide will search against, and
  // also build spectrum_mass_index.
  static void set_searchable_spectra(const std::vector<spectrum> &spectra);

  static void search_peptide_all_mods(int idno, int offset, int begin,
				      const std::string &peptide_seq,
				      int missed_cleavage_count,
				      score_stats &stats);

  // exposed for tinkering
  static double score_similarity(const spectrum &x, const spectrum &y,
				 int *peak_count);

  // conceptually these are 'protected:', but we're taking it easy here
public:
  void set_id() {
    id = next_id++;
    assert(next_id > 0);
  }

  static int next_id;
  static int next_physical_id;

  static std::vector<spectrum> searchable_spectra;
  // parent mass -> searchable_spectra index
  static std::multimap<double,
		       std::vector<spectrum>::size_type> spectrum_mass_index;
};


enum position { POSITION_NTERM=-3, POSITION_CTERM=-2, POSITION_WHOLE=-1 };

struct mass_trace_item {
  int position;			// enum position or 0..N
  double delta;			// 0 means delta is implied (i.e., for the
				// mass regime)
  const char *description;
};


// This is everything we want to remember about a match, so that we can report
// it later.
// FIX: Are any of these fields unneeded?
class match {
public:
  double hyper_score;
  double convolution_score;
  // FIX: what would be better for these next two?
  std::vector<double> ion_scores;
  std::vector<int> ion_peaks;
  int missed_cleavage_count;
  int spectrum_index;
  int peptide_begin;
  std::string peptide_sequence;
  double peptide_mass;
  int sequence_index;
  int sequence_offset;
  std::vector<mass_trace_item> mass_trace;

  match() : hyper_score(-1), convolution_score(-1), missed_cleavage_count(-1),
	    spectrum_index(-1), peptide_begin(-1), peptide_mass(-1),
	    sequence_index(-1), sequence_offset(-1) {
    ion_scores.resize(ION_MAX);
    ion_peaks.resize(ION_MAX);
  }

};


// This holds information about the best peptide/spectrum matches found, and
// some statistics.  Essentially, it holds the results of the search process.
class score_stats {
public:
  explicit score_stats(int n) {
    hyperscore_histogram.resize(n);
    best_score.resize(n, 100.0);
    second_best_score.resize(n, 100.0);
    best_match.resize(n);
    candidate_spectrum_count = 0;
  }
  
  // spectrum index -> (scaled, binned hyperscore -> count)
  std::vector< std::vector<int> > hyperscore_histogram;
  // spectrum index -> best_hyperscore (and 2nd best, respectively)
  std::vector<double> best_score;
  std::vector<double> second_best_score;
  // spectrum index -> [ <match info>, ... ]
  std::vector< std::vector<match> > best_match;
  unsigned long long candidate_spectrum_count; // may be > 2^32

  // This is a temp variable to implement xtandem's search limit (quirks mode
  // only)
  int combinations_searched;
};


// returns log-scaled hyperscore
double scale_hyperscore(double hyper_score);

#endif
