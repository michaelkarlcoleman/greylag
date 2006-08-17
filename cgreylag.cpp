

//	$Id$



#include "spectrum.h"


#include <cassert>
#include <cstdlib>
#include <cmath>



// what happens when an arbitrary exception is thrown?


int spectrum::next_id = 0;


// Read a spectrum, in ms2 format, returning true on success.
// (Currently our spectra either have charge 1, or 2 and 3 both--in the latter
// case, we set the secondary values.)
bool
spectrum::read(FILE *f) {
  // FIX: This seems pretty bad--there must a better way.

  peaks.clear();

  // read the primary header line pair
  char s[256];
  int r = std::fscanf(f, ":%256s%*[ \r\n]%lf%*[ \t]%d%*[ \t\r\n]",
		      s, &mass, &charge);
  if (r != 3 or std::strlen(s) >= 255)
    return false;
  name = std::string(s);

  // read the possible secondary header line pair
  int c = std::getc(f);
  if (c == EOF)
    return false;
  if (std::ungetc(c, f) != c)
    return false;
  secondary_mass = secondary_charge = 0;
  if (c == ':') {
    r = std::fscanf(f, ":%256s%*[ \r\n]%lf%*[ \t]%d%*[ \t\r\n]",
		    s, &secondary_mass, &secondary_charge);
    if (r != 3 or std::strlen(s) >= 255)
      return false;
    assert(secondary_charge == charge+1); // FIX
    //std::fprintf(stderr, "mass: %lf %lf\n", secondary_mass, predicted_secondary_mass());
    //assert(std::fabs(mass2 - ambiguous_mass()) < 0.1); // FIX
  }

  // read the peaks
  while (true) {
    int c = std::getc(f);
    if (c == EOF)
      return not std::ferror(f);
    if (c == ':')
      return std::ungetc(c, f) == ':';
    peak p;
    r = std::fscanf(f, "%lf%*[ \t]%lf%*[ \t\r\n]", &p.mz, &p.intensity);
    if (r != 2)
      return false;
    peaks.push_back(p);
  }
}
