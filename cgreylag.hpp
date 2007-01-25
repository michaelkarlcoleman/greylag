// cgreylag.hpp

// $Id$

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


#ifndef CGREYLAG_H
#define CGREYLAG_H

#include <cassert>
#include <cstdio>
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

  // residue char (incl '['/']') -> mass (per regime) with any fixed mod
  std::vector<double> fixed_residue_mass;
};

// FIX: want these class members to all be static, but had SWIG trouble.  The
// member 'the' is a handle to the singleton instance, for now.
class parameters {
public:
  bool quirks_mode;		// try to produce results closer to xtandem's
  bool estimate_only;		// just estimate work required
  bool show_progress;		// running progress on stderr

  // hyperscores within this ratio are considered "equal"
  double hyper_score_epsilon_ratio;

  std::vector<double> factorial; // factorial[n] == n!

  // deuterium not (yet) implemented
  double proton_mass;
  double hydrogen_mass;
  std::vector<mass_regime_parameters> parent_mass_regime;
  std::vector<mass_regime_parameters> fragment_mass_regime;

  // from XTP
  double parent_monoisotopic_mass_error_plus;
  double parent_monoisotopic_mass_error_minus;
  double fragment_mass_error;
  int minimum_ion_count;
  int minimum_peptide_length;

  bool spectrum_synthesis;
  bool check_all_fragment_charges;

  long maximum_modification_combinations_searched;
  unsigned int maximum_simultaneous_modifications_searched;

  static parameters the;
};


class sequence_run {
public:
  int sequence_index;		// index of this sequence in the database
  int sequence_offset;		// offset of this run's start in sequence
  std::string sequence;
  std::vector<int> cleavage_points;

  explicit sequence_run() { }
  explicit sequence_run(const int sequence_index, const int sequence_offset,
			const std::string sequence,
			const std::vector<int> cleavage_points)
    : sequence_index(sequence_index), sequence_offset(sequence_offset),
      sequence(sequence), cleavage_points(cleavage_points) {
  }
};


// information about the search that's passed in to--but not changed by--the
// C++ search code
class search_context {
public:
  int mod_count;
  int mass_regime_index;
  std::string pca_residues;
  double pca_delta;

  double N_delta;		// 0 if none
  double C_delta;		// 0 if none
  
  // maps residue char to list of indices for _count/_delta
  std::vector< std::vector<int> > delta_bag_lookup;
  std::vector<double> delta_bag_delta;
  std::vector<int> delta_bag_count;

  double parent_fixed_mass;
  double fragment_N_fixed_mass;
  double fragment_C_fixed_mass;

  // general search parameters
  int maximum_missed_cleavage_sites;

  std::vector<sequence_run> sequence_runs;
};


enum ion { ION_MIN, ION_Y=ION_MIN, ION_B, ION_MAX }; // NB: [ION_MIN, ION_MAX)

inline void operator++(ion &i) { i = ion(i + 1); }
inline void operator++(ion &i, int) { i = ion(i + 1); }


class peak {
public:
  double mz;
  double intensity;

  explicit peak(double mz=0, double intensity=0)
    : mz(mz), intensity(intensity) {
  }

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

  // Construct an empty spectrum.
  explicit spectrum(double mass=0, int charge=0) : mass(mass), charge(charge) {
    assert(charge >= 0);
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



  // Read spectra from file in ms2 format, tagging them with file_id.  Spectra
  // with charge zero are omitted from the result.  (All of them are read,
  // though, so the ids for any spectrum will be the same regardless of the
  // range used.)  If file_id == -1, the ms2 file is an annotated file
  // produced by split_ms2_by_mass_band.
  static std::vector<spectrum> read_spectra_from_ms2(FILE *f,
						     const int file_id);

  // Return sorted list of parent masses present in the given ms2 files.
  static std::vector<double> read_ms2_spectrum_masses(std::vector<int> fds);

  // Copy spectra from an ms2 file to a set of output ms2 files, one for each
  // band in the set of mass bands described by their upper bounds.  Extra
  // information is written in the first header line of each spectra so that
  // id's may be recovered.  So, for example, ':0002.0002.1' might become
  // ':0002.0002.1 # 0 45 78', where 0, 45, 78 are the spectrum's file_id,
  // physical_id, and id, respectively.  Multiply charged spectra will be
  // split into separate spectra (having the same physical_id).
  static void split_ms2_by_mass_band(FILE *inf, const std::vector<int> &outfds,
				     const int file_id, 
				     const std::vector<double>
				       &mass_band_upper_bounds);


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

  // Search sequence runs for matches according to the context against the
  // spectra.  Updates score_stats and the number of candidate spectra found.
  static void search_runs(const search_context &context, score_stats &stats);

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


struct mass_trace_item {
  int position;
  double delta;
  int id;			// FIX
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
  double parent_peptide_mass;
  //double fragment_peptide_mass;
  int sequence_index;
  int sequence_offset;
  std::vector<mass_trace_item> mass_trace;

  match() : hyper_score(-1), convolution_score(-1), missed_cleavage_count(-1),
	    spectrum_index(-1), peptide_begin(-1), parent_peptide_mass(-1),
	    sequence_index(-1), sequence_offset(-1) {
    ion_scores.resize(ION_MAX);
    ion_peaks.resize(ION_MAX);
  }

};


// This holds information about the best peptide/spectrum matches found, and
// some statistics.  Essentially, it holds the results of the search process.
class score_stats {
public:
  explicit score_stats(int n) : candidate_spectrum_count(0),
				spectra_with_candidates(0),
				improved_candidates(0),
				combinations_searched(0) {
    hyperscore_histogram.resize(n);
    best_score.resize(n, 100.0);
    second_best_score.resize(n, 100.0);
    best_match.resize(n);
  }
  
  // spectrum index -> (scaled, binned hyperscore -> count)
  std::vector< std::vector<int> > hyperscore_histogram;
  // spectrum index -> best_hyperscore (and 2nd best, respectively)
  std::vector<double> best_score;
  std::vector<double> second_best_score;
  // spectrum index -> [ <match info>, ... ]
  std::vector< std::vector<match> > best_match;

  // statistics for reporting
  unsigned long long candidate_spectrum_count; // may be > 2^32
  unsigned long spectra_with_candidates;
  unsigned long long improved_candidates;
  // to implement something like xtandem's search limit
  int combinations_searched;
};


// returns log-scaled hyperscore
double scale_hyperscore(double hyper_score);

#endif
