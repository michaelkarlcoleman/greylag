// cgreylag.hpp

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
#include <cmath>
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
  bool estimate_only;		// just estimate work required
  bool show_progress;		// running progress on stderr

  std::vector<double> ln_factorial; // factorial[n] == ln(n!)

  // deuterium not (yet) implemented
  double proton_mass;
  double hydrogen_mass;
  std::vector<mass_regime_parameters> parent_mass_regime;
  std::vector<mass_regime_parameters> fragment_mass_regime;

  // from GLP
  double parent_mass_tolerance_1; // for +1
  double parent_mass_tolerance_max; // for +N (typically +3)

  unsigned int minimum_peptide_length;

  double fragment_mass_tolerance;
  int intensity_class_count;

  static parameters the;
};


class sequence_run {
public:
  int sequence_index;		// index of this sequence in the database FIX:nn
  int sequence_offset;		// offset of this run's start in sequence
  std::string sequence;
  std::vector<int> cleavage_points;
  std::string name;		// e.g., initial defline word

  sequence_run() { }
  sequence_run(const int sequence_index, const int sequence_offset,
	       const std::string sequence, std::vector<int> cleavage_points,
	       const std::string name)
    : sequence_index(sequence_index), sequence_offset(sequence_offset),
      sequence(sequence), cleavage_points(cleavage_points),
      name(name) {
  }
};


// information about the search that's passed in to--but not changed by--the
// C++ search code
class search_context {
public:
  unsigned int mod_count;
  int mass_regime_index;
  int conjunct_index;
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
  bool nonspecific_cleavage;

  std::vector<sequence_run> sequence_runs;
};


class peak {
public:
  double mz;
  double intensity;
  int intensity_class;		// 0..N; -1 == unclassified

  explicit peak(double mz=0, double intensity=0, int intensity_class=-1)
    : mz(mz), intensity(intensity), intensity_class(intensity_class) {
  }

  char *__repr__() const;

  static bool less_mz(peak x, peak y) { return x.mz < y.mz; }

  static bool greater_intensity(peak x, peak y) {
    return x.intensity > y.intensity;
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
  double mass;			// [M + H+], aka "neutral mass" (protonated)
  int charge;
  static const int MAX_SUPPORTED_CHARGE = 15;
  static const int MAX_THEORETICAL_FRAGMENTS = 256;
  std::vector<peak> peaks;	// always ordered by increasing m/z!

  std::vector<int> intensity_class_counts; // number of peaks in each class

  // This is the mz range of peaks, but with a buffer on each end of size
  // fragment_mass_tolerance.
  double min_peak_mz;
  double max_peak_mz;

  // These are caches of the data in peaks above.  They speed up search.
  double *peak_mz_cache;
  int *peak_intensity_class_cache;

  double total_ion_current;

  // This is the total number of bins, 2*fragment_tolerance mz wide in peak
  // range.  So, for example, if the fragment tolerance is 1 and the min and
  // max mz are 200 and 1800, total peak bins would be 800.
  int total_peak_bins;
  int empty_peak_bins;	      // number of bins that aren't occupied by a peak

  std::string name;
  int file_id;			// file # of spectra's source
  int id;			// unique for all spectra
  int physical_id;		// FIX: unneeded?

  unsigned long comparisons;	// times scored against theoretical spectra

  // Construct an empty spectrum.
  explicit spectrum(double mass=0,
		    int charge=0) : mass(mass), charge(charge), min_peak_mz(0),
				    max_peak_mz(0), peak_mz_cache(0),
				    peak_intensity_class_cache(0),
				    total_ion_current(0),
				    total_peak_bins(0), empty_peak_bins(0),
				    file_id(-1), physical_id(-1),
				    comparisons(0) {
    assert(charge >= 0);
    if (charge > MAX_SUPPORTED_CHARGE)
      throw std::invalid_argument("attempt to create a spectrum with greater"
				  " than supported charge");
    set_id();
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


  // Read spectra from file in ms2 format, tagging them with file_id.  Before
  // reading, seek to absolute position offset_begin.  If offset_end != -1,
  // any spectra that begin after position offset_end in the file are not
  // read.
  static std::vector<spectrum>
  read_spectra_from_ms2(FILE *f, const int file_id, const long offset_begin,
			const long offset_end);


  // Filter peaks to limit their number according to TIC_cutoff_proportion,
  // and to remove those too large to be fragment products.  (Peaks are
  // assumed to be ordered by increasing mz.)
  void filter_peaks(double TIC_cutoff_proportion,
		    double parent_mass_tolerance_max);

  // Classify peaks and update class_counts.
  void classify(int intensity_class_count, double intensity_class_ratio,
		double fragment_mass_tolerance);

  // Update the *_cache fields from the peaks field.
  void update_peak_cache();

  void clear_peak_cache() {
    if (peak_mz_cache) {
      //delete[] peak_mz_cache;
      peak_mz_cache = 0;
    }
    if (peak_intensity_class_cache) {
      //delete[] peak_intensity_class_cache;
      peak_intensity_class_cache = 0;
    }
  }

  // Strictly speaking, we should release the cache on destruction.  This
  // would mean implementing bug-prone copy constructor and assignment
  // constructor functions.  Instead, skip it.  This leaks memory if
  // update_peak_cache() is called multiple times, but we don't do that.
  //
  // ~spectrum() { clear_peak_cache(); }

  // Store the list of spectra that search_peptide will search against, and
  // also build spectrum_mass_index.
  static void set_searchable_spectra(const std::vector<spectrum> &spectra);

  // Search sequence runs for matches according to the context against the
  // spectra.  Updates score_stats and the number of candidate spectra found.
  static void search_runs(const search_context &context, score_stats &stats);

  // conceptually these are 'protected:', but we're taking it easy here
public:
  void set_id() { id = next_id++; assert(next_id > 0); }

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

  // grrrr
  bool operator==(const mass_trace_item &x) const {
    return this->position == x.position and this->delta == x.delta;
  }
};


// This is everything we want to remember about a match, so that we can report
// it later.
class match {
public:
  double score;
  int spectrum_index;
  std::string peptide_sequence;
  char N_peptide_flank, C_peptide_flank; // adjacent residues; '-' if none
  double predicted_parent_mass;

  int mass_regime_index;
  int conjunct_index;
  double pca_delta;
  std::vector<mass_trace_item> mass_trace;

  // these are vectors of results for storing peptides found in multiple
  // sequence locations
  std::vector<int> peptide_begin; // absolute position within locus
  std::vector<std::string> sequence_name;

  match() : score(0), spectrum_index(-1), N_peptide_flank('-'),
	    C_peptide_flank('-'), predicted_parent_mass(0),
	    mass_regime_index(-1), conjunct_index(-1), pca_delta(0) { }
};


// This holds information about the best peptide/spectrum matches found, and
// some statistics.  Essentially, it holds the results of the search process.
class score_stats {
public:
  score_stats(int spectrum_count, int best_result_count)
    : candidate_spectrum_count(0), combinations_searched(0) {
    best_matches.resize(spectrum_count);
    assert(best_result_count >= 1);
    for (int i=0; i<spectrum_count; i++)
      best_matches[i].resize(best_result_count);
  }

  // spectrum index -> list of match (only top N matches per spectrum) The
  // best match lists are of fixed length (typically 5), and ordered
  // best-first.  Trailing entries with score 0 are null matches (to simplify
  // the code).
  std::vector< std::vector<match> > best_matches;

  // Statistics:
  // How many times was a spectrum scored against a peptide?
  unsigned long long candidate_spectrum_count; // may be > 2^32
  // How many different mod combinations were searched? (e.g., A*ABC and AA*BC
  // are two distinct combinations)
  int combinations_searched;	// FIX: could be > 2^31?
};


#endif
