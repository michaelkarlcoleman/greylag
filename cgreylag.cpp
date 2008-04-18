// C++ module for greylag

//     greylag, a collection of programs for MS/MS protein analysis
//     Copyright (C) 2006-2008  Stowers Institute for Medical Research
//
//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//     Contact: Mike Coleman
//              Stowers Institute for Medical Research
//              1000 East 50th Street
//              Kansas City, Missouri  64110
//              USA


#include "cgreylag.hpp"


#include <cerrno>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>



#ifdef __GNU__
// Use faster, non-threadsafe versions, since we don't use threads.  This is
// maybe 10% faster?
#define fgetsX fgets_unlocked
#define fputsX fputs_unlocked
#define fprintfX __builtin_fprintf_unlocked
#define ferrorX ferror_unlocked
#define getcX getc_unlocked
#else
#define fgetsX std::fgets
#define fputsX std::fputs
#define fprintfX std::fprintf
#define ferrorX std::ferror
#define getcX std::getc
#endif


parameters parameters::the;
const parameters &CP = parameters::the;


char *
peak::__repr__() const {
  static char temp[128];
  sprintf(temp, "<peak mz=%.4f intensity=%.4f intensity_class=%d>",
	  mz, intensity, intensity_class);
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
	  "<spectrum #%d (phys #%d) '%s' %.4f/%+d [%zd peaks]>",
	  id, physical_id, name.c_str(), mass, charge, peaks.size());
  return &temp[0];
}


// Return ln of n_C_k.
static inline double
ln_combination(unsigned int n, unsigned int k) {
  // Occasionally happens due to problems in scoring function. (FIX)
  if (n < k)
    return 0;

  assert(0 <= k and k <= n and n < CP.ln_factorial.size());
  return CP.ln_factorial[n] - CP.ln_factorial[k] - CP.ln_factorial[n-k];
}


// FIX: Currently we're reading spectrum files in Python, which is fast enough?
// // This is an error exit for read_spectra*.
// static inline void
// io_error(FILE *f, const char *message="") {
//   if (ferrorX(f))
//     message = "I/O error while reading (or writing) spectrum file";
//   throw std::ios_base::failure(message);
// }


// // true iff s consists entirely of whitespace
// static inline bool
// at_eol(const char *s) {
//   return s[strspn(s, " \t\r\n\f\v")] == 0;
// }


// // Read spectra from file in ms2 format, tagging them with file_id.

// // Multiply charged spectra (e.g., +2/+3) are split into separate spectra
// // having the same physical id.  Note that depending on
// // offset_begin/offset_end, we may end up with the first charge and not the
// // second, or vice versa.

// // The peak list is initially sorted by mz.  Throws an exception on invalid
// // input.  Error checking is stringent in this function.

// // This is pretty hideous.  Is there a simpler way to read, reasonably
// // efficiently, catching any input error?

// std::vector<spectrum>
// spectrum::read_spectra_from_ms2(FILE *f, const int file_id) {
//   std::vector<spectrum> spectra;

//   const int bufsiz = 1024;
//   char buf[bufsiz];
//   char *endp;

//   char *result = fgetsX(buf, bufsiz, f);
//   while (true) {
//     std::vector<std::string> names;
//     std::vector<double> masses;
//     std::vector<int> charges;

//     // read headers
//     while (true) {
//       if (ferrorX(f))
// 	io_error(f);
//       if (not result or buf[0] != ':')
// 	break;
//       names.push_back(std::string(buf+1, buf+std::strlen(buf)-1));
//       if (not fgetsX(buf, bufsiz, f))
// 	io_error(f, "bad ms2 format: mass/charge line expected");
//       errno = 0;
//       double mass = std::strtod(buf, &endp); // need double accuracy here
//       masses.push_back(mass);
//       if (errno or endp == buf or mass <= 0)
// 	io_error(f, "bad ms2 format: bad mass");
//       const char *endp0 = endp;
//       int charge = std::strtol(endp0, &endp, 10);
//       charges.push_back(charge);
//       if (errno or endp == endp0 or charge <= 0)
// 	io_error(f, "bad ms2 format: bad charge");
//       if (not at_eol(endp))
// 	io_error(f, "bad ms2 format: junk at end of mass/charge line");
//       result = fgetsX(buf, bufsiz, f);
//     }
//     if (not result) {
//       if (not names.empty())
// 	io_error(f, "bad ms2 format: spectrum has no peaks (file truncated?)");
//       break;
//     }
//     // read peaks
//     std::vector<peak> peaks;
//     while (true) {
//       peak p;
//       errno = 0;
//       p.mz = strtod(buf, &endp);
//       if (errno or endp == buf or p.mz <= 0)
// 	io_error(f, "bad ms2 format: bad peak mz");
//       const char *endp0 = endp;
//       p.intensity = strtod(endp0, &endp);
//       if (errno or endp == endp0 or p.intensity <= 0)
// 	io_error(f, "bad ms2 format: bad peak intensity");
//       if (not at_eol(endp))
// 	io_error(f, "bad ms2 format: junk at end of peak line");
//       peaks.push_back(p);

//       result = fgetsX(buf, bufsiz, f);
//       if (ferrorX(f))
// 	io_error(f);
//       if (not result or buf[0] == ':')
// 	break;
//     }

//     // add spectra to vector
//     if (names.empty() and not peaks.empty())
//       io_error(f, "bad ms2 format: missing header lines?");

//     std::sort(peaks.begin(), peaks.end(), peak::less_mz);

//     for (std::vector<double>::size_type i=0; i<names.size(); i++) {
//       spectrum sp(masses[i], charges[i]);
//       sp.peaks = peaks;
//       sp.name = names[i];
//       sp.file_id = file_id;
//       sp.physical_id = next_physical_id;
//       spectra.push_back(sp);
//     }
//     if (file_id != -1)
//       spectrum::next_physical_id += 1;
//   }
//   return spectra;
// }


// Filter peaks to limit their number according to TIC_cutoff_proportion, and
// to remove those too large to be fragment products.  (Peaks are assumed to
// be ordered by increasing mz.)
void
spectrum::filter_peaks(double TIC_cutoff_proportion,
		       double parent_mass_tolerance_max) {
  if (not (0 <= TIC_cutoff_proportion and TIC_cutoff_proportion <= 1)
      or parent_mass_tolerance_max < 0)
    throw std::invalid_argument("invalid argument value");

  // remove peaks too large to be fragment products
  // (i.e., peak m/z > parent [M+H+] + parent_mass_tolerance * charge_limit)
  const double peak_mz_limit = mass + parent_mass_tolerance_max;

  std::vector<peak>::iterator pi;
  for (pi=peaks.begin(); pi != peaks.end() and pi->mz <= peak_mz_limit; pi++);
  peaks.erase(pi, peaks.end());

  // FIX: note mzLowerBound..mzUpperBound here

  // now filter by TIC
  std::sort(peaks.begin(), peaks.end(), peak::greater_intensity);

  total_ion_current = 0;
  for (pi=peaks.begin(); pi != peaks.end(); pi++)
    total_ion_current += pi->intensity;

  const double i_limit = total_ion_current * TIC_cutoff_proportion;
  double accumulated_intensity = 0;
  for (pi=peaks.begin();
       pi != peaks.end() and accumulated_intensity <= i_limit; pi++)
    accumulated_intensity += pi->intensity;

  peaks.erase(pi, peaks.end());

  // restore order
  std::sort(peaks.begin(), peaks.end(), peak::less_mz);
}


// Classify peaks and update class_counts and peak bin counts.
// FIX: hoist some of this up into Python?  use (1,2,4) as parameter
void
spectrum::classify(int intensity_class_count, double intensity_class_ratio,
		   double fragment_mass_tolerance) {
  if (intensity_class_count < 1 or intensity_class_ratio <= 1.0)
    throw std::invalid_argument("invalid argument value");

  intensity_class_counts.clear();
  intensity_class_counts.resize(intensity_class_count);

  // e.g., 7 = 1 + 2 + 4
  // FIX: does this fail for non-integer ratios?
  int min_count = int((pow(intensity_class_ratio, intensity_class_count) - 1)
		      / (intensity_class_ratio - 1));

  std::sort(peaks.begin(), peaks.end(), peak::greater_intensity);
  std::vector<peak>::iterator pi=peaks.begin();
  for (int i_class=0; i_class<intensity_class_count; i_class++) {
    int peaks_this_class = int((pow(intensity_class_ratio, i_class) * peaks.size()
				/ min_count) + 0.5);
    for (int cp=0; cp < peaks_this_class and pi != peaks.end(); cp++, pi++) {
      pi->intensity_class = i_class;
      intensity_class_counts[i_class]++;
    }
  }
  std::sort(peaks.begin(), peaks.end(), peak::less_mz);

  if (not peaks.empty()) {
    min_peak_mz = peaks.front().mz - CP.fragment_mass_tolerance;
    max_peak_mz = peaks.back().mz + CP.fragment_mass_tolerance;
    total_peak_bins = (int((max_peak_mz - min_peak_mz)
			   / (2 * fragment_mass_tolerance) + 0.5));
    // This is the MyriMatch estimate.  Wouldn't it be better to simply count
    // how many are empty?
    empty_peak_bins = std::max<int>(total_peak_bins - peaks.size(), 0);
  }
}


// Update the *_cache fields from the peaks field.
void
spectrum::update_peak_cache() {
  clear_peak_cache();

  const unsigned int peak_count = peaks.size();
  peak_mz_cache = new double[peak_count+1];
  peak_intensity_class_cache = new int[peak_count+1];

  for (unsigned int i=0; i<peak_count; i++) {
    peak_mz_cache[i] = peaks[i].mz;
    peak_intensity_class_cache[i] = peaks[i].intensity_class;
  }

  // negative value is terminator
  peak_mz_cache[peak_count] = -1;
  peak_intensity_class_cache[peak_count] = -1;
}


// Store the list of spectra that search_peptide will search against, and also
// build spectrum_mass_index.
void
spectrum::set_searchable_spectra(const std::vector<spectrum> &spectra) {
  for (std::vector<spectrum>::size_type i=0; i<spectra.size(); i++)
    spectra[i].clear_peak_cache();
  searchable_spectra = spectra;

  spectrum_mass_index.clear();
  for (std::vector<spectrum>::size_type i=0; i<spectra.size(); i++) {
    spectrum_mass_index.insert(std::make_pair(spectra[i].mass, i));
    searchable_spectra[i].update_peak_cache();
  }
}


// Calculate peaks for a synthesized mass (not mz) ladder.
static inline void
get_synthetic_Y_mass_ladder(double *mass_ladder, const unsigned int ladder_size,
			    const double *mass_list,
			    const double fragment_C_fixed_mass) {
  double m = fragment_C_fixed_mass;

  for (unsigned int i=ladder_size; i>0; i--) {
    m += mass_list[i];
    mass_ladder[ladder_size-i] = m;
  }
}


// Calculate peaks for a synthesized mass (not mz) ladder.
static inline void
get_synthetic_B_mass_ladder(double *mass_ladder, const unsigned int ladder_size,
			    const double *mass_list,
			    const double fragment_N_fixed_mass) {
  double m = fragment_N_fixed_mass;

  for (unsigned int i=0; i<ladder_size; i++) {
    m += mass_list[i];
    mass_ladder[i] = m;
  }
}


// Generate synthetic spectra for a set of charges.  Only the charge and peaks
// of these spectra are initialized.
static inline void
synthetic_spectra(double synth_sp_mz[/* max_fragment_charge+1 */]
		                    [spectrum::MAX_THEORETICAL_FRAGMENTS],
		  const double *mass_list, const unsigned int mass_list_size,
		  const double fragment_N_fixed_mass,
		  const double fragment_C_fixed_mass,
		  const int max_fragment_charge) {
  assert(mass_list_size >= 1);
  const unsigned int ladder_size = mass_list_size - 1;
  double Y_mass_ladder[ladder_size];
  double B_mass_ladder[ladder_size];

  for (int charge=1; charge<=max_fragment_charge; charge++) {
    get_synthetic_Y_mass_ladder(Y_mass_ladder, ladder_size, mass_list,
				fragment_C_fixed_mass);
    get_synthetic_B_mass_ladder(B_mass_ladder, ladder_size, mass_list,
				fragment_N_fixed_mass);

    unsigned int y_i=0, b_i=0;
    double *sp_mz_p = synth_sp_mz[charge];

    // merge the two (already sorted) mass lists, converting them to mz values
    // in the process
    // FIX: try algorithm::merge?
    while (y_i < ladder_size and b_i < ladder_size)
      if (Y_mass_ladder[y_i] < B_mass_ladder[b_i])
	*(sp_mz_p++) = peak::get_mz(Y_mass_ladder[y_i++], charge);
      else
	*(sp_mz_p++) = peak::get_mz(B_mass_ladder[b_i++], charge);
    while (y_i < ladder_size)
      *(sp_mz_p++) = peak::get_mz(Y_mass_ladder[y_i++], charge);
    while (b_i < ladder_size)
      *(sp_mz_p++) = peak::get_mz(B_mass_ladder[b_i++], charge);
    *sp_mz_p = -1;		// invalid mz as terminator
  }
}


// Return the MyriMatch-style score of a (real) spectrum (x) versus a
// theoretical spectrum (y) generated from a peptide candidate.
//
// FIX: This code seems complicated.  Is the MyriMatch idea of using maps (or
// better yet, sets) for peak lists simpler?  As fast or faster?
//
// This is the innermost loop, so it's worthwhile to optimize this some.
// FIX const declaration
static inline double
score_spectrum(const spectrum &x,
	       const double synth_mzs[spectrum::MAX_THEORETICAL_FRAGMENTS]) {
  // FIX this
  unsigned int y_peak_count=0;
  for (const double *p=synth_mzs; *p >= 0; p++, y_peak_count++);

  double peak_best_delta[y_peak_count];
  int peak_best_class[y_peak_count];
  for (unsigned int i=0; i<y_peak_count; i++) {
    peak_best_delta[i] = CP.fragment_mass_tolerance;
    peak_best_class[i] = -1;
  }

  // Iterate over all closest pairs of peaks, collecting info on the class of
  // the best (meaning closest mz) real peak matching each theoretical peak
  // (for real/theoretical peak pairs that are no farther apart than
  // fragment_mass_tolerance).

  const unsigned int x_peak_count=x.peaks.size();
  assert(x_peak_count > 0 and y_peak_count > 0);

  unsigned int x_index=0, y_index=0;
  while (true) {
    const double delta = synth_mzs[y_index] - x.peak_mz_cache[x_index];
    if (std::abs(delta) <= peak_best_delta[y_index]) {
      peak_best_class[y_index] = x.peak_intensity_class_cache[x_index];
      peak_best_delta[y_index] = std::abs(delta);
    }
    if (delta > 0) {
      if (++x_index >= x_peak_count)
	break;
    } else {
      if (++y_index >= y_peak_count)
	break;
    }
  }

  assert(CP.intensity_class_count < INT_MAX);
  const unsigned int intensity_class_count = CP.intensity_class_count;
  unsigned int peak_hit_histogram[intensity_class_count];
  for (unsigned int i=0; i<intensity_class_count; i++)
    peak_hit_histogram[i] = 0;
  for (unsigned int i=0; i<y_peak_count; i++) {
    const int bc = peak_best_class[i];
    if (bc >= 0) {
      assert(bc < static_cast<int>(intensity_class_count));
      peak_hit_histogram[bc] += 1;
    }
  }

  // How many theoretical peaks overlap the real peak range?
  int valid_theoretical_peaks = 0;
  for (unsigned int i=0; i<y_peak_count; i++)
    if (x.min_peak_mz <= synth_mzs[i] and synth_mzs[i] <= x.max_peak_mz)
      valid_theoretical_peaks++;

  unsigned int total_peak_hits = 0;
  for (unsigned int i=0; i<intensity_class_count; i++)
    total_peak_hits += peak_hit_histogram[i];

  // How many valid theoretical peaks were misses?
  const int peak_misses = valid_theoretical_peaks - total_peak_hits;
  assert(peak_misses >= 0);

  double score = ln_combination(x.empty_peak_bins, peak_misses);
  for (unsigned int i=0; i<intensity_class_count; i++)
    score += ln_combination(x.intensity_class_counts[i],
			    peak_hit_histogram[i]);
  score -= ln_combination(x.total_peak_bins, valid_theoretical_peaks);

  return score;
}


// Return the score of the specified spectrum against a group of synthetic
// fragment spectra.

// For +1 and +2 precursors, this does what MyriMatch does.  For +3 (and
// above), this currently takes the sum of +1b/+1y or +2b/+2y, which may or
// may not make much sense.  (Wouldn't it be better to take the max of +1b/+2y
// or +2b/+1y?)
// FIX: implement MM smart +3 model, and/or come up with something better.
// FIX const decl
static inline double
score_spectrum_over_charges(const double
			    synth_sp_mz[/* max_fragment_charge+1 */]
			    [spectrum::MAX_THEORETICAL_FRAGMENTS],
			    const spectrum &sp,
			    const int max_fragment_charge) {
  double best_score = score_spectrum(sp, synth_sp_mz[1]);

  const int charge_limit = std::max<int>(1, sp.charge-1);
  for (int charge=2; charge<=charge_limit; charge++) {
    assert(charge <= max_fragment_charge);
    best_score += score_spectrum(sp, synth_sp_mz[charge]);
  }
  return best_score;
}


struct mass_trace_list {
  mass_trace_item item;
  const mass_trace_list *next;

  mass_trace_list(const mass_trace_list *p=0) : next(p) { }
};


// fuzzy comparison might not be necessary, but it's safer
static inline bool
score_equal(double s1, double s2) { return std::abs(s1-s2) < 1e-6; }


// Search for matches of this particular peptide modification variation
// against the spectra.  Updates score_stats and returns the number of
// candidate spectra found.
static inline void
evaluate_peptide(const search_context &context, match &m,
		 const mass_trace_list *mtlp, const double *mass_list,
		 const std::vector<std::vector<spectrum>::size_type>
		     &candidate_spectra,
		 score_stats &stats) {
  int max_candidate_charge = 0;
  typedef std::vector<std::vector<spectrum>::size_type>::const_iterator c_it;
  for (c_it it=candidate_spectra.begin(); it != candidate_spectra.end(); it++)
    max_candidate_charge
      = std::max<int>(max_candidate_charge,
		      spectrum::searchable_spectra[*it].charge);
  assert(max_candidate_charge >= 1);

  const int max_fragment_charge = std::max<int>(1, max_candidate_charge-1);
  assert(max_fragment_charge <= spectrum::MAX_SUPPORTED_CHARGE);

  const unsigned int peptide_size = m.peptide_sequence.size();
  // FIX: "2" assumes just two ion series (e.g., B and Y)
  assert(peptide_size <= 2*(spectrum::MAX_THEORETICAL_FRAGMENTS-1));

  // This array is preallocated to avoid the overhead of allocation and
  // deallocation within the inner loop.
  // FIX
  // charge -> mz array
  static double synth_sp_mz[spectrum::MAX_SUPPORTED_CHARGE+1]
                           [spectrum::MAX_THEORETICAL_FRAGMENTS];

  synthetic_spectra(synth_sp_mz, mass_list, peptide_size,
		    context.fragment_N_fixed_mass,
		    context.fragment_C_fixed_mass, max_fragment_charge);

  for (c_it it=candidate_spectra.begin(); it != candidate_spectra.end(); it++) {
    m.spectrum_index = *it;
    spectrum &sp = spectrum::searchable_spectra[m.spectrum_index];

    sp.comparisons += 1;
    stats.candidate_spectrum_count += 1;
    if (CP.estimate_only)
      continue;

    //std::cerr << "sp " << m.spectrum_index << std::endl;
    m.score = score_spectrum_over_charges(synth_sp_mz, sp, max_fragment_charge);
    //std::cerr << "score " << m.score << std::endl;

    std::vector<match> &sp_best_matches = stats.best_matches[m.spectrum_index];
    // Is this score good enough to be in the best score list?  (Better scores
    // are more negative.)
    if (m.score > sp_best_matches.back().score)
      continue;

    // We're saving this match, so remember the mass trace info, too.
    m.mass_trace.clear();
    for (const mass_trace_list *p=mtlp; p; p=p->next)
      m.mass_trace.push_back(p->item);

    // If this is duplicate, just append to the existing match and return
    for (std::vector<match>::reverse_iterator rit
	   = sp_best_matches.rbegin(); rit != sp_best_matches.rend(); rit++)
      if (score_equal(rit->score, m.score)
	  and rit->peptide_sequence == m.peptide_sequence
	  and rit->mass_regime_index == m.mass_regime_index
	  and rit->conjunct_index == m.conjunct_index
	  and rit->mass_trace == m.mass_trace) {
	rit->peptide_begin.push_back(m.peptide_begin[0]);
	rit->sequence_name.push_back(m.sequence_name[0]);
	return;
      }

    // Otherwise, insert this match in the correct position
    assert(sp_best_matches.size() >= 1);
    std::vector<match>::reverse_iterator rit = sp_best_matches.rbegin();
    for (;rit+1 != sp_best_matches.rend() and (rit+1)->score > m.score; rit++)
      *rit = *(rit+1);
    *rit = m;
  }
}


// IDEA FIX: Use a cost parameter to define the iterative search front.


// Choose one possible residue modification position.  Once they're all
// chosen, then evaluate.
static inline void
choose_residue_mod(const search_context &context, match &m,
		   const mass_trace_list *mtlp, double *mass_list,
		   const std::vector<std::vector<spectrum>::size_type>
		       &candidate_spectra,
		   score_stats &stats,
		   std::vector<int> &db_remaining,
		   const unsigned int remaining_positions_to_choose,
		   const unsigned int next_position_to_consider) {
  assert(remaining_positions_to_choose
	 <= m.peptide_sequence.size() - next_position_to_consider);

  if (remaining_positions_to_choose == 0) {
    stats.evaluation_count += 1;
    evaluate_peptide(context, m, mtlp, mass_list, candidate_spectra, stats);
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
	mtl.item.conjunct_item_index = db_index;
	choose_residue_mod(context, m, &mtl, mass_list, candidate_spectra,
			   stats, db_remaining,
			   remaining_positions_to_choose-1, i+1);
	db_remaining[db_index] += 1;
      }
      mass_list[i] = save_pos_mass;
    }
  }
}


static inline int
search_peptide(const search_context &context, const double &p_mass, match &m,
	       const std::string &run_sequence, const int &begin_index,
	       const unsigned int &peptide_size, score_stats &stats,
	       std::vector<int> &db_remaining) {
  // We'd like to consider only the spectra whose masses are within the
  // tolerance range of theoretical mass of the peptide (p_mass).
  // Unfortunately, the tolerance range of each spectrum depends on its
  // charge, since the tolerance is specified in m/z units.  Furthermore,
  // when deciding how to adjust begin/end to generate the next peptide,
  // we need to be conservative, as a high-charge spectrum may come into
  // range, even though only +1 spectra are currently in-range.  So, we
  // first screen using the maximum tolerance range, then check the actual
  // spectra below.

  const double sp_mass_lb = p_mass - CP.parent_mass_tolerance_max;
  const double sp_mass_ub = p_mass + CP.parent_mass_tolerance_max;

  const spmi_c_it candidate_spectra_info_begin
    = spectrum::spectrum_mass_index.lower_bound(sp_mass_lb);
  if (candidate_spectra_info_begin == spectrum::spectrum_mass_index.end()) {
    // spectrum masses all too low to match peptide (peptide too long)
    return -1;
  }

  const spmi_c_it candidate_spectra_info_end
    = spectrum::spectrum_mass_index.upper_bound(sp_mass_ub);
  if (candidate_spectra_info_end == spectrum::spectrum_mass_index.begin()) {
    // spectrum masses all too high to match peptide (peptide too short)
    return 1;
  }
  if (candidate_spectra_info_begin == candidate_spectra_info_end) {
    // no spectrum with close-enough parent mass
    return 0;
  }

  // Now generate the list of candidate spectra that are actually
  // in-range, considering the charge of each spectrum.
  std::vector<std::vector<spectrum>::size_type> candidate_spectra_0;
  for (spmi_c_it it=candidate_spectra_info_begin;
       it != candidate_spectra_info_end; it++) {
    const spectrum &csp = spectrum::searchable_spectra[it->second];
    if (std::abs(csp.mass - p_mass) <= (csp.charge
					* CP.parent_mass_tolerance_1))
      candidate_spectra_0.push_back(it->second);
  }
  if (candidate_spectra_0.empty())
    return 0;

  m.predicted_parent_mass = p_mass;
  m.peptide_sequence.assign(run_sequence, begin_index, peptide_size);

  double mass_list[peptide_size];
  for (unsigned int i=0; i<peptide_size; i++)
    mass_list[i] = (CP.fragment_mass_regime[context.mass_regime_index]
		    .fixed_residue_mass[m.peptide_sequence[i]]);

  choose_residue_mod(context, m, NULL, mass_list, candidate_spectra_0,
		     stats, db_remaining, context.mod_count, 0);
  return 0;
}


// Search a sequence run for matches according to the context against the
// spectra.  Updates score_stats and the number of candidate spectra found.

// FIX: examine carefully for signed/unsigned problems

static inline void
search_run(const search_context &context, const sequence_run &sequence_run,
	   score_stats &stats) {
  const unsigned int min_peptide_length
    = std::max<unsigned int>(CP.minimum_peptide_length, context.mod_count);
  const std::vector<double> &fixed_parent_mass
    = CP.parent_mass_regime[context.mass_regime_index].fixed_residue_mass;

  const std::string &run_sequence = sequence_run.sequence;

  // This match will be passed inward and used to record information that we
  // need to remember about a match when we finally see one.  At that point, a
  // copy of this match will be saved.
  // FIX: can this be done better?
  match m;
  m.sequence_name.push_back(sequence_run.name);
  m.peptide_begin.push_back(0);
  int &peptide_begin = m.peptide_begin[0];
  m.mass_regime_index = context.mass_regime_index;
  m.conjunct_index = context.conjunct_index;
  m.pca_delta = context.pca_delta;

  assert(context.delta_bag_delta.size() == context.delta_bag_count.size());
  // counts remaining as mod positions are chosen
  std::vector<int> db_remaining = context.delta_bag_count;

  // the rightmost 'end' seen when all spectra masses were too high
  // (0 means none yet encountered.)
  unsigned int next_end = 0;

  // p_mass is the parent mass of the peptide run_sequence[p_begin:p_end]
  double p_mass = context.parent_fixed_mass;
  int p_begin=0, p_end=0;

  const std::vector<int> *cleavage_points = &sequence_run.cleavage_points;
  std::vector<int>::size_type cleavage_points_size = cleavage_points->size();
  // "fake" cleavage point vector used for nonspecific cleavage (v[i]==i)
  static std::vector<int> nonspecific_cleavage_points;
  if (context.nonspecific_cleavage) {
    cleavage_points_size = run_sequence.size() + 1;
    for (unsigned int i=nonspecific_cleavage_points.size();
	 i<cleavage_points_size; i++)
      nonspecific_cleavage_points.push_back(i);
    cleavage_points = &nonspecific_cleavage_points;
  }
  assert(cleavage_points_size >= 2);	// endpoints must be present

  for (unsigned int begin=0; begin<cleavage_points_size-1; begin++) {
    const int begin_index = (*cleavage_points)[begin];
    peptide_begin = begin_index + sequence_run.sequence_offset;

    // If pca_residues specified, require one of them to be the peptide
    // N-terminal residue
    switch (context.pca_residues.size()) {
    case 0:
      break;
    case 2:
      if (context.pca_residues[1] == run_sequence[begin_index])
	break;
    case 1:
      if (context.pca_residues[0] == run_sequence[begin_index])
	break;
      continue;
    default:
      assert(false);
    }

    assert(begin_index >= p_begin and begin_index >= 0);
    for (int i=p_begin; i<begin_index; i++)
      p_mass -= fixed_parent_mass[run_sequence[i]];
    p_begin = begin_index;

    unsigned int end = std::max<unsigned int>(begin + 1, next_end);
    for (; end<cleavage_points_size; end++) {
      if (not context.nonspecific_cleavage) {
	const int missed_cleavage_count = end - begin - 1;
	if (missed_cleavage_count > context.maximum_missed_cleavage_sites)
	  break;
      }

      const int end_index = (*cleavage_points)[end];
      assert(end_index > begin_index);
      const unsigned int peptide_size = end_index - begin_index;
      if (peptide_size < min_peptide_length)
	continue;

      assert(end_index >= 0);
      if (end_index >= p_end)
	for (int i=p_end; i<end_index; i++)
	  p_mass += fixed_parent_mass[run_sequence[i]];
      else
	for (int i=end_index; i<p_end; i++)
	  p_mass -= fixed_parent_mass[run_sequence[i]];
      p_end = end_index;

//       std::cerr << "peptide: "
// 		<< std::string(run_sequence, begin_index, peptide_size)
// 		<< " p_mass: " << p_mass << std::endl;

      int status = search_peptide(context, p_mass, m, run_sequence,
				  begin_index, peptide_size, stats,
				  db_remaining);
      if (status == -1)		// peptide was too long
	break;
      if (status == 1) {	// peptide was too short
	// so next start with a peptide at least this long
	if (end == cleavage_points_size - 1)
	  return;
	next_end = end;
      }
    }
  }
}


// Search all sequence runs for matches according to the context against the
// spectra.  Updates score_stats and the number of candidate spectra found.
void
spectrum::search_runs(const search_context &context, score_stats &stats) {
  const int no_runs = context.sequence_runs.size();
  for (int i=0; i<no_runs; i++) {
    search_run(context, context.sequence_runs[i], stats);

    if (CP.show_progress)
      std::cerr << i+1 << " of " << no_runs << " sequences, "
		<< stats.candidate_spectrum_count << " candidates\r"
		<< std::flush;
  }

  if (CP.show_progress)
    std::cerr << std::endl;
}
