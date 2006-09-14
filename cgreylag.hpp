
//	$Id$


#ifndef CXTPY_H
#define CXTPY_H


#include <cassert>
#include <cstdio>
#include <list>
#include <map>
#include <vector>
#include <stdexcept>
#include <string>
#include <utility>


class score_stats;


// FIX: want these class members to all be static, but had SWIG trouble
class parameters {
public:
  bool quirks_mode;		// try to produce results closer to xtandem's

  std::vector<double> factorial; // factorial[n] == n!

  // these are indexed by atom char
  // (slightly bogus, but all the atoms we use have single-char names)
  std::vector<double> average_atomic_mass;
  std::vector<double> monoisotopic_atomic_mass;

  double proton_mass;
  // monoisotopic, for now (?)
  double hydrogen;
  double hydroxide;
  double water_mass;

  // these are indexed by residue char (plus '[' and ']' for N and C termini)
  std::vector<double> average_residue_mass;
  std::vector<double> monoisotopic_residue_mass;
  std::vector<double> modification_mass;
  std::vector< std::vector<double> > potential_modification_mass;
  std::vector< std::vector<double> > potential_modification_mass_refine;

  // from XTP
  double cleavage_N_terminal_mass_change;
  double cleavage_C_terminal_mass_change;
  double parent_monoisotopic_mass_error_plus;
  double parent_monoisotopic_mass_error_minus;
  double fragment_mass_error;
  int minimum_ion_count;

  static parameters the;
};


class peak {
 public:
  double mz;
  double intensity;

  peak(double mz=0, double intensity=0)
    : mz(mz), intensity(intensity) { }

  char *__repr__() const;

  static bool less_mz(peak x, peak y) { return x.mz < y.mz; }
  static bool less_intensity(peak x, peak y) {
    return x.intensity < y.intensity;
  }
};


// A spectrum searched on multiple charges will be turned into separate
// spectra having differing id's, but the same physical id.

class spectrum {
public:
  double mass;
  int charge;
  std::vector<peak> peaks;
  double max_peak_intensity;
  double sum_peak_intensity;
  double normalization_factor;
  std::string name;
  int file_id;			// file # of spectra's source
  int id;			// unique for all spectra
  int physical_id;

  spectrum(double mass=0, int charge=0) : mass(mass), charge(charge) {
    id = next_id++;
    assert(next_id > 0);
    max_peak_intensity = -1;
    sum_peak_intensity = -1;
    normalization_factor = 1;
    file_id = -1;
    physical_id = -1;
  }

  char *__repr__() const;

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

  // Read spectra from file in ms2 format, tagging them with file_id.
  static std::vector<spectrum> read_spectra(FILE *f, const int file_id);

  // Sets max/sum_peak_intensity, according to peaks and
  // normalization_factor.
  void calculate_intensity_statistics();

  // Examine this spectrum to see if it should be kept.  If so, return true
  // and sort the peaks by mz, normalize them and limit their number.
  bool filter_and_normalize(double minimum_fragment_mz, double dynamic_range,
			    int minimum_peaks, int total_peaks);

  // Store the list of spectra that search_peptide will search against, and
  // also build spectrum_mass_index.
  static void set_searchable_spectra(const std::vector<spectrum> &spectra);

  static int search_peptide(int idno, int offset, int begin,
			    const std::string &peptide_seq, bool is_N,
			    bool is_C, int missed_cleavage_count,
			    score_stats stats);


  // Return the similarity score between this spectrum and that, and also a
  // count of common peaks in *peak_count.
  static double score_similarity(const spectrum &x, const spectrum &y,
				 int *peak_count);


public:				// conceptually protected:
  static int next_id;
  static int next_physical_id;

  static std::vector<spectrum> searchable_spectra;
  // parent mass -> searchable_spectra index
  static std::multimap<double,
		       std::vector<spectrum>::size_type> spectrum_mass_index;
};


// This is everything we want to remember about a match, so that we can report
// it later.
class match {
public:
  double convolution_score;
  // FIX: Will these next two be too slow?  We only need to hold two values
  // here...
  std::map<char, double> ion_scores;
  std::map<char, int> ion_peaks;
  int missed_cleavage_count;
  int spectrum_index;
  int peptide_begin;
  std::string peptide_sequence;
  std::vector<double> potential_mod_pattern; // sufficient?
  double peptide_mod_mass;	// unneeded?
  int sequence_index;
  int sequence_offset;

  match() : convolution_score(-1), missed_cleavage_count(-1),
	    spectrum_index(-1), peptide_begin(-1), peptide_mod_mass(-1),
	    sequence_index(-1), sequence_offset(-1) {
  }

};


// This holds information about the best peptide/spectrum matches found, and
// some statistics.  Essentially, it holds the results of the search process.
class score_stats {
public:
  // spectrum index -> (scaled, binned hyperscore -> count)
  std::vector< std::vector<int> > hyperscore_histogram;
  // spectrum index -> best_hyperscore (and 2nd best, respectively)
  std::vector<double> best_score, second_best_score;
  // spectrum index -> [ <match info>, ... ]
  std::vector< std::vector<match> > best_match;
};


// returns the monoisotopic mass
double get_peptide_mass(const std::string &peptide_seq, bool is_N, bool is_C);

// returns log-scaled hyperscore
double scale_hyperscore(double hyper_score);

#endif
