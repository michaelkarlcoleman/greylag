



#ifndef SPECTRUM_H
#define SPECTRUM_H


#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>


class peak {
 public:
  double mz;
  double intensity;

  peak() : mz(0), intensity(0) { }
};


class spectrum {
public:
  double mass;
  int charge;
  bool ambiguous_charge;	// means charge could be charge+1, too
  std::vector<peak> peaks;
  int id;
  std::string name;

  spectrum() : mass(0), charge(0), ambiguous_charge(false) {
    id = next_id++;
    assert(next_id > 0);
  }

  double ambiguous_mass() const {
    return (mass-1.008) / charge * (charge+1); // check this
  }

protected:
  int next_id = 0;

};


// Read a spectrum, in ms2 format, setting failbit on invalid format.
// (Currently our spectra either have charge 1, or 2 and 3 both--in the latter
// case, we set 'ambiguous_charge'.)
template<typename T, typename charT, typename traits>
std::basic_istream<charT, traits>&
operator>>(std::basic_istream<charT, traits>& in, spectrum& sp) {
  // FIX: This is hideous--there must a better way.
  // read the primary header line pair
  string s;
  if (!std::getline(in, s))
    return in;
  if (s.size() < 1 or s[0] != ':') {
    in.setstate(std::ios_base::failbit);
    return in;
  }
  sp.name = s;
  if (!std::getline(in, s))
    return in;
  const char *cs = s.c_str();
  char *endp;
  std::errno = 0;
  sp.mass = strtod(cs, &endp);
  if (endp == cs or std::errno) {
    in.setstate(std::ios_base::failbit);
    return in;
  }
  cs = endp;
  std::errno = 0;
  sp.charge = strtol(cs, &endp);
  if (endp == cs or std::errno) {
    in.setstate(std::ios_base::failbit);
    return in;
  }
  // read the secondary (ambiguous) header line pair, if present
  if (in.peek() == ':') {
    if (!std::getline(in, s))	// discard name
      return in;
    if (!std::getline(in, s))
      return in;
    const char *cs = s.c_str();
    char *endp;
    std::errno = 0;
    double mass2 = strtod(cs, &endp);
    if (endp == cs or std::errno) {
      in.setstate(std::ios_base::failbit);
      return in;
    }
    cs = endp;
    std::errno = 0;
    int charge2 = strtol(cs, &endp);
    if (endp == cs or std::errno) {
      in.setstate(std::ios_base::failbit);
      return in;
    }
    if (charge2 != sp.charge+1 or fabs(mass2 - sp.ambiguous_mass()) > 0.1) {
      in.setstate(std::ios_base::failbit);
      return in;
    }      
    sp.ambiguous_charge = true;
  }
  // now read peaks
  while (in and in.peek() != ':') {
    peak p;
    if (!std::getline(in, s))
      return in;
    const char *cs = s.c_str();
    char *endp;
    std::errno = 0;
    p.mz = strtod(cs, &endp);
    if (endp == cs or std::errno) {
      in.setstate(std::ios_base::failbit);
      return in;
    }
    cs = endp;
    std::errno = 0;
    p.intensity = strtod(cs, &endp);
    if (endp == cs or std::errno) {
      in.setstate(std::ios_base::failbit);
      return in;
    }
    sp.peaks.push_back(p);
  }
  return in;
}


#endif
