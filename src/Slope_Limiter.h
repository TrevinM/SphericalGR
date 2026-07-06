// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Slope limiters for hydro
//
// Reference: Font et.al., arXiv:gr-qc/0110047
//
//================================================

#ifndef SLOPE_LIMITER_H
#define SLOPE_LIMITER_H

enum { minmod_lim, MC_lim };
//
//================================================
// base class - doesn't do much... 
//================================================
//
class Slope {
protected:
  int limiter_type;
public:
  Slope() {};
  ~Slope() {};
  // for optional argument see MC below
  virtual double limiter(double a, double b, bool minmod = false) = 0;
  inline double min(double a, double b) { return a < b ? a : b; };
  inline double min(double a, double b, double c) {
    const double ab = min(a,b);
    return ab < c ? ab : c; 
  };

  const char * Name() {
    if (limiter_type == minmod_lim) 
      return "minmod";
    else if (limiter_type == MC_lim) 
      return "MC";
    else
      return "unknown slope limiter";
  };
};
//
//================================================
// minmod limiter (see eq. (14) in Font et.al.)
//================================================
//
class minmod : public Slope {
public:
  minmod() : Slope() {
    limiter_type = minmod_lim;
  };
  ~minmod() {};
  inline double limiter(double a, double b, bool minmod = false) {
    if (a*b <= 0.0)
      return 0.0;
    else {
      if (abs(a) < abs(b))
	return a;
      else
	return b;
    }
  };
};
//
//================================================
// monotonized central-difference limiter
//    (see eq. (15) in Font et.al.)
//================================================
//
class MC : public Slope {
public:
  MC() : Slope() {
    limiter_type = MC_lim;
  };
  ~MC() {};
  inline double limiter(double a, double b, bool minmod = false) {
    //
    // use minmod
    // 
    if (minmod)
      if (a*b <= 0.0)
	return 0.0;
      else {
	if (abs(a) < abs(b))
	  return a;
	else
	  return b;
      }
    else 
      if (a*b <= 0.0)
	return 0.0;
      else {
	const double c = 0.5*(a + b);
	return min(2.0*abs(a),2.0*abs(b),abs(c))*a/abs(a);
      }
  };
};

#endif  /* SLOPE_LIMITER_H */
