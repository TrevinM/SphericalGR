// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Approximate Riemann Solvers
//
//================================================

#ifndef APP_RIE_SOLVER_H
#define APP_RIE_SOLVER_H

enum { llf_solver, hlle_solver };
//
//================================================
// base class - doesn't do much... 
//================================================
//
class Riemann_Solver {
public:
  Riemann_Solver() {};
  ~Riemann_Solver() {};
  //
  // allows as input variables left and right fluxes, conservative 
  // variables, and maxima of eigenvalues
  //
  virtual double flux(double f_L, double f_R, 
		      double q_L, double q_R,
		      double lam_0_L, double lam_p_L, double lam_m_L,
		      double lam_0_R, double lam_p_R, double lam_m_R,
		      bool verbose = false) = 0;
  virtual const char * Name() = 0; 
  inline double max(double a, double b) { return a > b ? a : b; }; 
  inline double min(double a, double b) { return a < b ? a : b; };
  inline double max(double a, double b, double c) {
    const double ab = max(a,b);
    return ab > c ? ab : c; 
  };
  inline double min(double a, double b, double c) {
    const double ab = min(a,b);
    return ab < c ? ab : c; 
  };
};
//
//================================================
// local Lax-Friedrichs (LLF)
//
// see, e.g., Thierfelder, Bernuzzi & Bruegmann, arXiv:1104.4751,
//  eq. (50) 
//================================================
//
class LLF : public Riemann_Solver {
public:
  LLF() : Riemann_Solver() {};
  ~LLF() {};
  const char * Name() { return "local Lax-Friedrichs (LLF)"; }; 
  inline double flux(double f_L, double f_R, 
		     double q_L, double q_R,
		     double lam_0_L, double lam_p_L, double lam_m_L,
		     double lam_0_R, double lam_p_R, double lam_m_R,
		     bool verbose = false) {
    const double a_L_max = max(lam_0_L,lam_p_L,lam_m_L);
    const double a_R_max = max(lam_0_R,lam_p_R,lam_m_R);
    const double a = max(a_L_max,a_R_max);
    // cout << a << endl;
    return 0.5 * (f_L + f_R - a * ( q_R - q_L ));
  };
};
//
//================================================
// Harten, Lax, van Leer and Einfeldt (HLLE), see
//     eq. (48) in DLSS
//================================================
//
class HLLE : public Riemann_Solver {
public:
  HLLE() : Riemann_Solver() {};
  ~HLLE() {};
  const char * Name() { return "HLLE"; }; 
  inline double flux(double f_L, double f_R, 
		     double q_L, double q_R,
		     double lam_0_L, double lam_p_L, double lam_m_L,
		     double lam_0_R, double lam_p_R, double lam_m_R,
		     bool verbose = false) {
    const double lam_L_max = max(lam_0_L,lam_p_L,lam_m_L);
    const double lam_R_max = max(lam_0_R,lam_p_R,lam_m_R);
    const double lam_max   = max(0.0,lam_L_max,lam_R_max);
    const double lam_L_min = min(lam_0_L,lam_p_L,lam_m_L);
    const double lam_R_min = min(lam_0_R,lam_p_R,lam_m_R);
    const double lam_min   = - min(0.0,lam_L_min,lam_R_min);
    //    return (lam_max*f_L - lam_min*f_R + lam_min*lam_max*(q_R - q_L))
    //      /(lam_max - lam_min);
    //
    if (verbose)
      cout << " lam_L_min = " << lam_L_min
	   << " lam_R_max = " << lam_R_max 
	   << " q_L = " << q_L 
	   << " q_R = " << q_R 
	   << endl;
    // see eq. (48) in DLSS, also (10.21) in Toro
    return (lam_max*f_L + lam_min*f_R - lam_min*lam_max*(q_R - q_L))
      /(lam_max + lam_min);
  };
};

#endif  /* APP_RIE_SOLVER_H */
