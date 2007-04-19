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
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>



#ifdef __GNU__
#define NOTHROW __attribute__ ((nothrow))
// Use faster, non-threadsafe versions, since we don't use threads.  This is
// maybe 10% faster?
#define fgetsX fgets_unlocked
#define fputsX fputs_unlocked
#define fprintfX __builtin_fprintf_unlocked
#define ferrorX ferror_unlocked
#define getcX getc_unlocked
#else
#define NOTHROW
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


// FIX: does this actually help inlining?
// Return ln of n_C_k.
static inline double
ln_combination_(unsigned int n, unsigned int k) NOTHROW {
  // Occasionally happens due to problems in scoring function. (FIX)
  if (n < k)
    return 0;

  assert(0 <= k and k <= n and n < CP.ln_factorial.size());
  return CP.ln_factorial[n] - CP.ln_factorial[k] - CP.ln_factorial[n-k];
}
double
ln_combination(unsigned int n, unsigned int k) {
  return ln_combination(n, k);
}



// This is an error exit for read_spectra*.
static inline void
io_error(FILE *f, const char *message="") {
  if (ferrorX(f))
    message = "I/O error while reading (or writing) spectrum file";
  throw std::ios_base::failure(message);
}


// true iff s consists entirely of whitespace
static bool
at_eol(const char *s) {
  return s[strspn(s, " \t\r\n\f\v")] == 0;
}


// Check whether we've past the end offset, throwing exception on I/O error.
static inline bool
check_past_end(FILE *f, const long offset_end) {
  if (offset_end == -1)
    return false;
  long pos = std::ftell(f);
  if (pos == -1)
    io_error(f);
  return pos >= offset_end;
}


// Read spectra from file in ms2 format, tagging them with file_id.  Before
// reading, seek to absolute position offset_begin.  If offset_end != -1, any
// spectra that begin after position offset_end in the file are not read.

// Multiply charge spectra (e.g., +2/+3) are split into separate spectra
// having the same physical id.  Note that depending on
// offset_begin/offset_end, we may end up with the first charge and not the
// second, or vice versa.

// The peak list is initially sorted by mz.  Throws an exception on invalid
// input.  Error checking is the most stringent in this function.  Other
// readers can be a little looser since all ms2 files eventually get read here
// anyway.

// This is pretty hideous.  Is there a simpler way to read, reasonably
// efficiently, catching any input error?

std::vector<spectrum>
spectrum::read_spectra_from_ms2(FILE *f, const int file_id,
				const long offset_begin,
				const long offset_end) {
  std::vector<spectrum> spectra;

  const int bufsiz = 1024;
  char buf[bufsiz];
  char *endp;

  // first seek to offset_begin and synchronize at next "\n:"
  if (offset_begin > 0) {
    if (std::fseek(f, offset_begin-1, SEEK_SET) == -1)
      io_error(f);
    int c = getcX(f);
    while (c != '\n' and c != EOF)
      c = getcX(f);
    if (ferrorX(f))
      io_error(f);
    if (c == '\n') {
      do {
	c = getcX(f);
      } while (c != ':' and c != EOF);
      if (ferrorX(f))
	io_error(f);
      if (c == ':')
	std::ungetc(c, f);
    }
  }
  // at this point we expect to read the first ':' of a header, or EOF

  bool past_end = check_past_end(f, offset_end);
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
      if (not past_end)
	names.push_back(std::string(buf+1, buf+std::strlen(buf)-1));
      if (not fgetsX(buf, bufsiz, f))
	io_error(f, "bad ms2 format: mass/charge line expected");
      errno = 0;
      double mass = std::strtod(buf, &endp); // need double accuracy here
      if (not past_end)
	masses.push_back(mass);
      if (errno or endp == buf or mass <= 0)
	io_error(f, "bad ms2 format: bad mass");
      const char *endp0 = endp;
      int charge = std::strtol(endp0, &endp, 10);
      if (not past_end)
	charges.push_back(charge);
      if (errno or endp == endp0 or charge <= 0)
	io_error(f, "bad ms2 format: bad charge");
      if (not at_eol(endp))
	io_error(f, "bad ms2 format: junk at end of mass/charge line");
      past_end = past_end or check_past_end(f, offset_end);
      result = fgetsX(buf, bufsiz, f);
    }
    if (names.empty() and past_end)
      break;
    if (not result) {
      if (not names.empty())
	io_error(f, "bad ms2 format: spectrum has no peaks (file truncated?)");
      break;
    }
    // read peaks
    std::vector<peak> peaks;
    while (true) {
      peak p;
      errno = 0;
      p.mz = strtod(buf, &endp);
      if (errno or endp == buf or p.mz <= 0)
	io_error(f, "bad ms2 format: bad peak mz");
      const char *endp0 = endp;
      p.intensity = strtod(endp0, &endp);
      if (errno or endp == endp0 or p.intensity <= 0)
	io_error(f, "bad ms2 format: bad peak intensity");
      if (not at_eol(endp))
	io_error(f, "bad ms2 format: junk at end of peak line");
      peaks.push_back(p);

      past_end = past_end or check_past_end(f, offset_end);
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
      spectra.push_back(sp);
    }
    if (file_id != -1)
      spectrum::next_physical_id += 1;
  }
  return spectra;
}


// Filter peaks to limit their number according to TIC_cutoff_proportion, and
// to remove those too large to be fragment products.  (Peaks are assumed to
// be ordered by increasing mz.)
void
spectrum::filter_peaks(double TIC_cutoff_proportion,
		       double parent_mass_tolerance, int charge_limit) {
  if (not (0 <= TIC_cutoff_proportion and TIC_cutoff_proportion <= 1)
      or parent_mass_tolerance < 0 or charge_limit < 1)
    throw std::invalid_argument("invalid argument value");

  // remove peaks too large to be fragment products
  // (i.e., peak m/z > parent [M+H+] + parent_mass_tolerance * charge_limit)
  const double peak_mz_limit = mass + parent_mass_tolerance * charge_limit;

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


// Store the list of spectra that search_peptide will search against, and also
// build spectrum_mass_index.
void
spectrum::set_searchable_spectra(const std::vector<spectrum> &spectra) {
  searchable_spectra = spectra;
  for (std::vector<spectrum>::size_type i=0; i<spectra.size(); i++)
    spectrum_mass_index.insert(std::make_pair(spectra[i].mass, i));
}


// Calculate peaks for a synthesized mass (not mz) ladder.
static inline void
get_synthetic_Y_mass_ladder(std::vector<peak> &mass_ladder,
			    const std::vector<double> &mass_list,
			    const double fragment_C_fixed_mass) NOTHROW {
  double m = fragment_C_fixed_mass;

  const int ladder_size = mass_ladder.size();
  for (int i=ladder_size-1; i>=0; i--) {
    m += mass_list[i+1];
    mass_ladder[ladder_size-1-i] = peak(m, 1.0);
  }
}


// Calculate peaks for a synthesized mass (not mz) ladder.
static inline void
get_synthetic_B_mass_ladder(std::vector<peak> &mass_ladder,
			    const std::vector<double> &mass_list,
			    const double fragment_N_fixed_mass) NOTHROW {
  double m = fragment_N_fixed_mass;

  const int ladder_size = mass_ladder.size();
  for (int i=0; i<=ladder_size-1; i++) {
    m += mass_list[i];
    mass_ladder[i] = peak(m, 1.0);
  }
}


// Generate synthetic spectra for a set of charges.  Only the charge and peaks
// of these spectra are initialized.
static inline void
synthetic_spectra(spectrum synth_sp[/* max_fragment_charge+1 */],
		  const std::vector<double> &mass_list,
		  const double fragment_N_fixed_mass,
		  const double fragment_C_fixed_mass,
		  const double max_fragment_charge) NOTHROW {
  std::vector<peak> mass_ladder(mass_list.size()-1);

  for (int charge=1; charge<=max_fragment_charge; charge++) {
    spectrum &sp = synth_sp[charge];
    //sp.mass = fragment_peptide_mass;
    sp.charge = charge;
    sp.peaks.resize(mass_ladder.size() * (ION_MAX-ION_MIN));
    std::vector<peak>::size_type pi=0;

    for (ion ion_type=ION_MIN; ion_type<ION_MAX; ion_type++) {
      switch (ion_type) {
      case ION_Y:
	get_synthetic_Y_mass_ladder(mass_ladder, mass_list,
				    fragment_C_fixed_mass);
	break;
      case ION_B:
	get_synthetic_B_mass_ladder(mass_ladder, mass_list,
				    fragment_N_fixed_mass);
	break;
      default:
	assert(false);
      }

      for (std::vector<peak>::size_type i=0; i<mass_ladder.size(); i++) {
	sp.peaks[pi++].mz = peak::get_mz(mass_ladder[i].mz, charge);
      }
    }
    std::sort(sp.peaks.begin(), sp.peaks.end(), peak::less_mz);
  }
}


// Return the MyriMatch-style score of a (real) spectrum (x) versus a
// theoretical spectrum (y) generated from a peptide candidate.
//
// FIX: This code seems complicated.  Is the MyriMatch idea of using maps (or
// better yet, sets) for peak lists simpler?  As fast or faster?
//
// This is the innermost loop, so it's worthwhile to optimize this some.
// FIX: Should some of these vectors be arrays?
static inline double
score_spectrum(const spectrum &x, const spectrum &y) NOTHROW {
  assert (not x.peaks.empty() and not y.peaks.empty());	// FIX

  std::vector<double> peak_best_delta(y.peaks.size(),
				      CP.fragment_mass_tolerance);
  std::vector<int> peak_best_class(y.peaks.size(), -1);

  std::vector<peak>::const_iterator x_it = x.peaks.begin();
  std::vector<peak>::const_iterator y_it = y.peaks.begin();

  // Iterate over all closest pairs of peaks, collecting info on the class of
  // the best (meaning closest mz) real peak matching each theoretical peak
  // (for real/theoretical peak pairs that are no farther apart than
  // fragment_mass_tolerance).
  int y_index = 0;
  while (x_it != x.peaks.end() and y_it != y.peaks.end()) {
    const double delta = y_it->mz - x_it->mz;
    if (std::abs(delta) <= peak_best_delta[y_index]) {
      peak_best_class[y_index] = x_it->intensity_class;
      peak_best_delta[y_index] = std::abs(delta);
    }
    if (delta > 0)
      x_it++;
    else {
      y_it++;
      y_index++;
    }
  }

  assert(CP.intensity_class_count < INT_MAX);
  std::vector<unsigned> peak_hit_histogram(CP.intensity_class_count);
  std::vector<int>::const_iterator b_it;
  for (b_it = peak_best_class.begin(); b_it != peak_best_class.end(); b_it++)
    if (*b_it >= 0) {
      assert(*b_it < static_cast<int>(peak_hit_histogram.size()));
      peak_hit_histogram[*b_it] += 1;
    }

  // How many theoretical peaks overlap the real peak range?
  int valid_theoretical_peaks = 0;
  for (y_it = y.peaks.begin(); y_it != y.peaks.end(); y_it++)
    if (x.min_peak_mz <= y_it->mz and y_it->mz <= x.max_peak_mz)
      valid_theoretical_peaks++;

  // How many valid theoretical peaks were misses?
  const int peak_misses = (valid_theoretical_peaks
			   - accumulate(peak_hit_histogram.begin(),
					peak_hit_histogram.end(), 0));
#if 0
  for (unsigned int i=0; i<peak_hit_histogram.size(); i++)
    if (x.intensity_class_counts[i] < peak_hit_histogram[i])
      std::cerr << i << " C(" << x.intensity_class_counts[i] << ", "
		<< peak_hit_histogram[i] << ")" << std::endl;
  if (x.empty_peak_bins < peak_misses)
    std::cerr << "M C(" << x.empty_peak_bins << ", " << peak_misses << ")"
	      << std::endl;
  if (x.total_peak_bins < valid_theoretical_peaks)
    std::cerr << "T C(" << x.total_peak_bins << ", " << valid_theoretical_peaks
	      << ")" << std::endl;
#endif

#if 0
  for (unsigned int i=0; i<peak_hit_histogram.size(); i++)
    std::cerr << i << " C(" << x.intensity_class_counts[i] << ", "
	      << peak_hit_histogram[i] << ")" << std::endl;
  std::cerr << "C(" << x.empty_peak_bins << ", " << peak_misses << ")"
	    << std::endl;
  std::cerr << "C(" << x.total_peak_bins << ", " << valid_theoretical_peaks
	    << ")" << std::endl;
  std::cerr << std::endl;
#endif

  assert(peak_misses >= 0);

  double score = 0.0;
  for (unsigned int i=0; i<peak_hit_histogram.size(); i++)
    score += ln_combination_(x.intensity_class_counts[i],
			     peak_hit_histogram[i]);
  score += ln_combination_(x.empty_peak_bins, peak_misses);
  score -= ln_combination_(x.total_peak_bins, valid_theoretical_peaks);

  return score;
}


// Return the score of the specified spectrum against a group of synthetic
// fragment spectra.

// For +1 and +2 precursors, this does what MyriMatch does.  For +3 (and
// above), this currently takes the sum of +1b/+1y or +2b/+2y, which may or
// may not make much sense.  (Wouldn't it be better to take the max of +1b/+2y
// or +2b/+1y?)
// FIX: implement MM smart +3 model, and/or come up with something better.
static inline double
score_spectrum_over_charges(spectrum synth_sp[/* max_fragment_charge+1 */],
			    const spectrum &sp,
			    const int max_fragment_charge) NOTHROW {
  double best_score = score_spectrum(sp, synth_sp[1]);

  const int charge_limit = std::max<int>(1, sp.charge-1);
  for (int charge=2; charge<=charge_limit; charge++) {
    assert(charge <= max_fragment_charge);
    best_score += score_spectrum(sp, synth_sp[charge]);
  }
  return best_score;
}


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
  assert(m.peptide_sequence.size() == mass_list.size());

  int max_candidate_charge = 0;
  for (spmi_c_it it=candidate_spectra_info_begin;
       it != candidate_spectra_info_end; it++)
    max_candidate_charge
      = std::max<int>(max_candidate_charge,
		      spectrum::searchable_spectra[it->second].charge);
  assert(max_candidate_charge >= 1);

  const int max_fragment_charge = std::max<int>(1, max_candidate_charge-1);
  assert(max_fragment_charge <= spectrum::MAX_SUPPORTED_CHARGE);
  // e.g. +2 -> 'B' -> spectrum
  // FIX: should this be a std::vector??
  static spectrum synth_sp[spectrum::MAX_SUPPORTED_CHARGE+1];
  synthetic_spectra(synth_sp, mass_list, context.fragment_N_fixed_mass,
		    context.fragment_C_fixed_mass, max_fragment_charge);

  for (spmi_c_it candidate_it = candidate_spectra_info_begin;
       candidate_it != candidate_spectra_info_end; candidate_it++) {
    m.spectrum_index = candidate_it->second;
    spectrum &sp = spectrum::searchable_spectra[m.spectrum_index];

    sp.comparisons += 1;
    stats.candidate_spectrum_count += 1;
    if (CP.estimate_only)
      continue;

    //std::cerr << "sp " << m.spectrum_index << std::endl;
    m.score = score_spectrum_over_charges(synth_sp, sp, max_fragment_charge);
    //std::cerr << "score " << m.score << std::endl;

    // Is this score good enough to be in the best score list?  (Better scores
    // are more negative.)
    if (m.score >= stats.best_matches[m.spectrum_index].back().score)
      continue;

    // We're saving this match, so remember the mass trace info, too.
    m.mass_trace.clear();
    for (const mass_trace_list *p=mtlp; p; p=p->next)
      m.mass_trace.push_back(p->item);

    // Is this match the first candidate for this spectrum?
    if (stats.best_matches[m.spectrum_index].front().score == 0)
      stats.spectra_with_candidates += 1;

    // Insert this match in the correct position.
    int i;
    for (i=stats.best_matches[m.spectrum_index].size() - 2; i>=0; i--)
      if (stats.best_matches[m.spectrum_index][i].score <= m.score)
	break;
      else
	stats.best_matches[m.spectrum_index][i+1]
	  = stats.best_matches[m.spectrum_index][i];
    stats.best_matches[m.spectrum_index][i+1] = m;
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
  assert(remaining_positions_to_choose
	 <= m.peptide_sequence.size() - next_position_to_consider);

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
  const int min_peptide_length = std::max<int>(CP.minimum_peptide_length,
					       context.mod_count);
  const std::vector<double> &fixed_parent_mass \
    = CP.parent_mass_regime[context.mass_regime_index].fixed_residue_mass;

  const std::string &run_sequence = sequence_run.sequence;
  const std::vector<int> &cleavage_points = sequence_run.cleavage_points;
  assert(cleavage_points.size() >= 2); // endpoints assumed always present

  // This match will be passed inward and used to record information that we
  // need to remember about a match when we finally see one.  At that point, a
  // copy of this match will be saved.
  match m;			// FIX: move inward?
  m.sequence_name = sequence_run.name;

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
    const int begin_index = cleavage_points[begin];
    m.peptide_begin = begin_index + sequence_run.sequence_offset;
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

      //std::cerr << "lb..ub: " << sp_mass_lb << ".." << sp_mass_ub << std::endl;

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

      m.predicted_parent_mass = p_mass;

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
  const int no_runs = context.sequence_runs.size();
  for (int i=0; i<no_runs; i++) {
    search_run(context, context.sequence_runs[i], stats);

    if (CP.show_progress)
      std::cerr << i+1 << " of " << no_runs << " sequences, "
		<< stats.candidate_spectrum_count << " cand for "
		<< stats.spectra_with_candidates << " sp\r" << std::flush;
  }

  if (CP.show_progress)
    std::cerr << std::endl;
}
