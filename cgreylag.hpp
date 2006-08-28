
//	$Id$


#ifndef SPECTRUM_H
#define SPECTRUM_H


#include <cassert>
#include <cstdio>
#include <list>
#include <vector>
#include <string>


class parameters {
public:
  // these are indexed by residue char (plus '[' and ']' for N and C termini)
  std::vector<double> average_residue_mass;
  std::vector<double> monoisotopic_residue_mass;
  std::vector<double> average_atomic_mass;
  std::vector<double> monoisotopic_atomic_mass;
  std::vector<double> modification_mass;
  std::vector< std::vector<double> > potential_modification_mass;
  std::vector< std::vector<double> > potential_modification_mass_refine;

  double cleave_N_terminal_mass_change;
  double cleave_C_terminal_mass_change;

  double proton_mass;
  double water_mass;		// monoisotopic, for now (?)

  static parameters the;
};


class peak {
 public:
  double mz;
  double intensity;

  peak() : mz(0), intensity(0) { }

  char *__repr__() const;

  static bool less_mz(peak x, peak y) { return x.mz < y.mz; }
  static bool less_intensity(peak x, peak y) {
    return x.intensity < y.intensity;
  }
};


class spectrum {
public:
  double mass;
  int charge;
  double secondary_mass;
  int secondary_charge;		// charge+1 if present, else 0
  std::vector<peak> peaks;
  int id;
  std::string name;

  spectrum() : mass(0), charge(0), secondary_mass(0), secondary_charge(0) {
    id = next_id++;
    assert(next_id > 0);
  }

  char *__repr__() const;

  // this doesn't quite work (error over 0.1)
  double predicted_secondary_mass() const {
    const double proton = 1.008;
    return proton + (mass-proton) / charge * (secondary_charge);
  }

  // Read a spectrum, in ms2 format, returning true on success.
  bool read(FILE *f);

  // Examine this spectrum to see if it should be kept.  If so, return true
  // and sort the peaks by mz, normalize them and limit their number.
  bool filter_and_normalize(double minimum_fragment_mz, double dynamic_range,
			    int minimum_peaks, int total_peaks);

protected:
  static int next_id;		// start at 0

};



#endif
