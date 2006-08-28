

//	$Id$



#include "spectrum.h"


#include <cerrno>
#include <cstdlib>
#include <cmath>
#include <stdexcept>


// speed up reading, at some expense to portability
// fast version assumes no threads and single-precision good enough
// (mass and intensity values currently under 7 digits)
// XXX: even faster: mmap + strtol only
#ifdef SLOW
#define strtodX std::strtod
#define fgetsX std::fgets
#define getcX std::getc
#define ungetcX std::ungetc
#define ferrorX std::ferror
#else
#define strtodX std::strtof
#define fgetsX fgets_unlocked
#define getcX getc_unlocked
#define ungetcX std::ungetc
#define ferrorX ferror_unlocked
#endif


parameters parameters::the;


char *
peak::__repr__() const {
  static char temp[128];
  sprintf(temp, "<peak mz=%.4f intensity=%.4f>", mz, intensity);
  return &temp[0];
}


int spectrum::next_id = 0;


char *
spectrum::__repr__() const {
  static char temp[1024];
  if (secondary_charge)
    sprintf(temp, "<spectrum #%d '%s' %.4f/%+d (%.4f/%+d) [%d peaks]>",
	    id, name.substr(0, 800).c_str(), mass, charge, secondary_mass,
	    secondary_charge, peaks.size());
  else
    sprintf(temp, "<spectrum #%d '%s' %.4f/%+d [%d peaks]>",
	    id, name.substr(0, 800).c_str(), mass, charge, peaks.size());
  return &temp[0];
}




// Read a spectrum, in ms2 format, returning true on success.
// (Currently our spectra either have charge 1, or 2 and 3 both--in the latter
// case, we set the secondary values.)
bool
spectrum::read(FILE *f) {
  // FIX: This seems pretty ugly.  Is there a better way?
  // Could we eliminate some of the (non-critical) error checks?
  // limitations: lines truncated to bufsiz, extra chars at eol ignored

  peaks.clear();

  // read the primary header line pair
  const int bufsiz = 1024;
  char buf[bufsiz];
  
  if (not fgetsX(buf, bufsiz, f) or buf[0] != ':')
    return false;
  name = std::string(buf+1, buf+std::strlen(buf)-1);
  if (not fgetsX(buf, bufsiz, f))
    return false;
  char *endp;
  errno = 0;
  mass = std::strtof(buf, &endp);
  if (errno or endp == buf)
    return false;
  charge = std::strtol(endp, &endp, 10);
  if (errno or endp == buf)
    return false;

  // read the possible secondary header line pair
  secondary_mass = secondary_charge = 0;
  if (not fgetsX(buf, bufsiz, f))
    return false;
  if (buf[0] == ':') {
    if (not fgetsX(buf, bufsiz, f))
      return false;
    errno = 0;
    secondary_mass = std::strtof(buf, &endp);
    if (errno or endp == buf)
      return false;
    secondary_charge = std::strtol(endp, &endp, 10);
    if (errno or endp == buf)
      return false;
    if (not fgetsX(buf, bufsiz, f))
      return false;
  }
  // at this point, we've read the first peak line
  if (buf[0] == ':')
    return false;

  // read the peaks
  while (true) {
    peak p;
    errno = 0;
    p.mz = std::strtof(buf, &endp);
    if (errno or endp == buf)
      return false;
    p.intensity = std::strtof(endp, &endp);
    if (errno or endp == buf)
      return false;
    peaks.push_back(p);

    int c = getcX(f);
    if (c == EOF)
      return not ferrorX(f);
    char c2 = ungetcX(c, f);
    if (c != c2)
      return false;		// ungetc failed
    if (c == ':')
      return true;
    if (not fgetsX(buf, bufsiz, f))
      return false;
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

  if (minimum_fragment_mz < 0 or dynamic_range <= 0 or minimum_peaks < 0 or total_peaks < 0)
    throw std::invalid_argument("invalid argument value");

  // remove_isotopes (XXX: redundant w/below?)
  // >>> removes multiple entries within 0.95 Da of each other, retaining the
  // highest value. this is necessary because of the behavior of some peak
  // finding routines in commercial software

  std::sort(peaks.begin(), peaks.end(), peak::less_mz);

  // FIX: maybe better just to accumulate result in a temporary?
  for (std::vector<peak>::size_type i=0; i+1<peaks.size();)
    if (peaks[i+1].mz - peaks[i].mz < 0.95)
      if (peaks[i+1].intensity > peaks[i].intensity)
	peaks.erase(peaks.begin() + i);
      else
	peaks.erase(peaks.begin() + i+1);
    else
      i++;

  // remove parent
  // >>> set up m/z regions to ignore: those immediately below the m/z of the
  // parent ion which will contain uninformative neutral loss ions, and those
  // immediately above the parent ion m/z, which will contain the parent ion
  // and its isotope pattern

  const double parent_mz = 1.00727 + (mass - 1.00727) / charge;
  for (std::vector<peak>::size_type i=0; i<peaks.size(); i++)
    if (not (parent_mz >= peaks[i].mz
	     and parent_mz - peaks[i].mz >= 50.0 / charge
	     or parent_mz < peaks[i].mz
	     and peaks[i].mz - parent_mz >= 5.0 / charge))
    if (parent_mz >= peaks[i].mz
	? parent_mz - peaks[i].mz < 50.0 / charge
	: peaks[i].mz - parent_mz < 5.0 / charge)
      peaks.erase(peaks.begin() + i);

  // remove_low_masses
  std::vector<peak>::iterator pi;
  for (pi=peaks.begin(); pi != peaks.end() and pi->mz <= minimum_fragment_mz; pi++);
  peaks.erase(peaks.begin(), pi);


  // normalize spectrum
  // >>> use the dynamic range parameter to set the maximum intensity value
  // for the spectrum.  then remove all peaks with a normalized intensity < 1
  if (peaks.empty())
    return false;
  std::vector<peak>::iterator most_intense_peak
    = std::max_element(peaks.begin(), peaks.end(), peak::less_intensity);
  const double normalization_factor
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
  //     # * is_noise attempts to determine if the spectrum is simply noise. if the spectrum
  //     # * does not have any peaks within a window near the parent ion mass, it is considered
  //     # * noise.
  //     #limit = sp.mass / sp.charge
  //     #if sp.charge < 3:
  //     #    limit = sp.mass - 600.0
  //     #if not [ True for mz, i in sp.peaks if mz > limit ]:
  //     #    return "looks like noise"

  // clean_isotopes removes peaks that are probably C13 isotopes
  // FIX: this seems almost identical to remove_isotopes above (kill one off?)
  for (std::vector<peak>::size_type i=0; i+1<peaks.size();)
    if (peaks[i+1].mz - peaks[i].mz < 1.5)
      if (peaks[i+1].intensity > peaks[i].intensity)
	peaks.erase(peaks.begin() + i);
      else
	peaks.erase(peaks.begin() + i+1);
    else
      i++;

  // keep only the most intense peaks
  if (peaks.size() > (unsigned int) total_peaks) {
    std::sort(peaks.begin(), peaks.end(), peak::less_intensity);
    peaks.erase(peaks.begin(), peaks.end() - total_peaks);
    std::sort(peaks.begin(), peaks.end(), peak::less_mz);
  }

  return true;
}
