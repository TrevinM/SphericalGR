// Tell emacs that this is -*-c++-*- mode
//================================================
// Classes for Equation of State
//================================================
#ifndef EOS_H
#define EOS_H

#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

using namespace std;
enum { poly, gammalaw, piecepoly, idealgas, gasrad };

//
//================================================
// 
// Base class for equation of state - doesn't do much...
//
//================================================
//
class EOS {
  
protected:
  int eos_type;
  int cs_error;
  double tau_last;
  double last_rho_0, rho_current;   // used for iteration in rho_0(rho)
  ostringstream name;
  typedef double (EOS::*FCT_PTR)(double);
public:
  EOS(double tau_last_i) : tau_last(tau_last_i), last_rho_0(0.0), 
			   rho_current(0.0), cs_error(0) {};
  virtual ~EOS() {};
  double P(double rho_0, double epsilon, double temp) {
    return P(rho_0, epsilon);
  }
  virtual double P(double rho_0, double epsilon) = 0;
  virtual double rho_0(double P) = 0;     // should be used for initial data only.
  virtual double h(double rho_0, double epsilon) = 0;
  virtual double cold_eps(double rho_0) = 0;
  virtual double temperature() { return tau_last; };
  virtual double temperature(double rho_0, double epsilon) {
    P(rho_0, epsilon); return tau_last;
  };
  // 
  virtual double sound_speed(double rho_0, double epsilon) = 0;
  inline double dPdrho_0(double rho_0, double epsilon) {
    const double delta = 1.0e-3 * rho_0;
    const double P_plus = P(rho_0, epsilon);
    const double P_minus = P(rho_0 - delta, epsilon);
    return (P_plus - P_minus)/(delta);
      };
  inline double dPdepsilon(double rho_0, double epsilon) {
    const double delta = 1.0e-3 * epsilon;
    const double P_plus = P(rho_0, epsilon + delta);
    const double P_minus = P(rho_0, epsilon);
    return (P_plus - P_minus)/(delta);
  };
  //================================================
  // Return name and type of EOS
  //================================================
  string Name() {
    return name.str();
  };
  int type() { return eos_type; }
  //================================================
  // Routines needed for initial data
  //================================================
  double rho_0_of_rho(double rho) {       // needed in TOV_BH
    if (rho < 0.0) {
      cerr << "EOS: found negative rho in EOS::rho_0_of_rho!" << endl;
      return 0.0;
    }
    if (rho == 0.0) return 0.0;
    rho_current = rho;
    if (last_rho_0 == 0) last_rho_0 = rho;
    double rho_0_1 = last_rho_0;
    double rho_0_2 = 0.9 * rho_0_1;
    bool success = zbrac(&EOS::rho_0_root, rho_0_1, rho_0_2);   // make sure that roots are bracketed!
    if (!success) cout << "EOS: could not bracket roots in EOS::rho_0_of_rho!" << endl;
    const double tol = 1.e-6;
    return last_rho_0 = zbrent(&EOS::rho_0_root, rho_0_1, rho_0_2, tol*rho_0_2);
  }
  double rho_0_root(double rho_0) {       // used for rootfinding in rho_0_of_rho
    double eps = cold_eps(rho_0);
    if (!isfinite(eps)) cout << "EOS: in EOS::rho_0_root: eps = " << eps << " rho_0 = " << rho_0 << endl;
    double rho = rho_0 * ( 1.0 + eps );
    return rho - rho_current;
  }
  virtual double epsilon_root(double tau) = 0;
  //================================================
  // Numerical Recipes routines for rootfinding
  //================================================
  Bool zbrac(FCT_PTR func, Doub &x1, Doub &x2) {
    const Int NTRY=50;
    const Doub FACTOR=1.6;
    if (x1 == x2) cerr << "EOS: Bad initial range in EOS::zbrac" << endl;
    // make sure that x1 < x2
    if (x2 < x1) {
      const double temp = x2;
      x2 = x1;
      x1 = temp;
    }
    Doub f1=(this->*func)(x1);
    Doub f2=(this->*func)(x2);
    for (Int j=0;j<NTRY;j++) {
      if (f1*f2 <= 0.0) return true;
      if (abs(f1) < abs(f2) && x1 > 0.0) {
	x1 -= FACTOR*(x2-x1);
	if (x1 < 0.0) x1 = 0.0;
	f1=(this->*func)(x1);
	//cout << " f1 = " << f1 << endl;
      } else {
	f2=(this->*func)(x2 += FACTOR*(x2-x1));
	//cout << " f2 = " << f2 << endl;
      }   
    }
    return false;
  }

  Bool zbrac_org(FCT_PTR func, Doub &x1, Doub &x2) {
    const Int NTRY=50;
    const Doub FACTOR=1.6;
      if (x1 == x2) cerr << "EOS: Bad initial range in EOS::zbrac" << endl;
    Doub f1=(this->*func)(x1);
    Doub f2=(this->*func)(x2);
    for (Int j=0;j<NTRY;j++) {
      if (f1*f2 < 0.0) return true;
      if (abs(f1) < abs(f2))
	f1=(this->*func)(x1 += FACTOR*(x1-x2));
      else
	f2=(this->*func)(x2 += FACTOR*(x2-x1));
    }
    return false;
  }
  Doub zbrent(FCT_PTR func, const Doub x1, const Doub x2, const Doub tol) {
    const Int ITMAX=100;
    const Doub EPS=numeric_limits<Doub>::epsilon();
    Doub a=x1,b=x2,c=x2,d,e,fa=(this->*func)(a),fb=(this->*func)(b),fc,p,q,r,s,tol1,xm;
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
       cerr << "Root must be bracketed in zbrent" << endl;
    fc=fb;
    for (Int iter=0;iter<ITMAX;iter++) {
      if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
	c=a;
	fc=fa;
	e=d=b-a;
      }
      if (abs(fc) < abs(fb)) {
	a=b;
	b=c;
	c=a;
	fa=fb;
	fb=fc;
	fc=fa;
      }
      tol1=2.0*EPS*abs(b)+0.5*tol;
      xm=0.5*(c-b);
      if (abs(xm) <= tol1 || fb == 0.0) return b;
      if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
	s=fb/fa;
	if (a == c) {
	  p=2.0*xm*s;
	  q=1.0-s;
	} else {
	  q=fa/fc;
	  r=fb/fc;
	  p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	  q=(q-1.0)*(r-1.0)*(s-1.0);
	}
	if (p > 0.0) q = -q;
	p=abs(p);
	Doub min1=3.0*xm*q-abs(tol1*q);
	Doub min2=abs(e*q);
	if (2.0*p < (min1 < min2 ? min1 : min2)) {
	  e=d;
	  d=p/q;
	} else {
	  d=xm;
	  e=d;
	}
      } else {
	d=xm;
	e=d;
      }
      a=b;
      fa=fb;
      if (abs(d) > tol1)
	b += d;
      else
	b += SIGN(tol1,xm);
      fb=(this->*func)(b);
    }
    //cerr << "EOS: Maximum number of iterations exceeded in zbrent" << endl;
    return b;
  }


};
//
//================================================
// derived class for Polytrope
//
//    P = K * rho_0^Gamma
//
//================================================
//
class polytrope : public EOS {
private:
  double K;
  double Gamma;
  double n;
public:
  // Constructor
  polytrope(char * input_file) : EOS(0.0) { 
    eos_type = poly;
    ifstream infile;
    infile.open(input_file);
    if (!infile) {
      cerr << " EOS: Can't open " << input_file
	   << " for input. This is bad. " << endl;
    } else 
      cout << " EOS: Reading EOS polytropic parameters from file " 
	   << input_file << endl;
    char buf[100], c;
    infile.get(buf,100,'='); infile.get(c); infile >> K;
    infile.get(buf,100,'='); infile.get(c); infile >> Gamma;
    infile.close();
    cout << " EOS: Setting up polytropic EOS with Kappa = " << K 
	 << " and Gamma = " << Gamma << endl; 
    n = 1.0/(Gamma - 1.0);
    name << "polytropic EOS with K = " << K 
	 << " and Gamma = " << Gamma ;
  };
  // Destructor
  ~polytrope() {};
  //
  // Routines
  //
  inline double P(double rho_0, double epsilon) {
    if (Gamma == 2.0)
      return K * rho_0 * rho_0;
    else 
      return K * pow(rho_0,Gamma);
  };
  inline double rho_0(double P) {;
    return pow( P/K, 1.0/Gamma );
  }
  inline double dPdrho_0(double rho_0, double epsilon) {
    return epsilon * ( Gamma - 1.0 );
  };
  inline double dPdepsilon(double rho_0, double epsilon) {
    return rho_0 * ( Gamma - 1.0 );
  };
  inline double sound_speed(double rho_0, double epsilon) {
    // FIX!!!
    return 0.0;
  }
  inline double h(double rho_0, double epsilon) {
    if (Gamma == 2.0)
      return 1.0 + 2.0 * K * rho_0;
    else
      return 1.0 + (n + 1.0) * K * pow(rho_0,Gamma-1.0);
  };
  inline double cold_eps(double rho_0) {
    if (Gamma == 2.0)
      return K * rho_0;
    else
      return K * pow(rho_0,Gamma-1.0)/(Gamma - 1.0);
  };
  double epsilon_root(double tau){return 0.0;};
};


//
//================================================
// derived class for Gamma-law EOS
//
//    P = rho_0 * epsilon ( \Gamma - 1 )
//
//================================================
//
class gamma_law : public EOS {
private:
  double K;  // needed only for cold part
  double Gamma;
  double n;
public:
  // Constructor
  gamma_law(char * input_file) : EOS(0.0) { 
    eos_type = gammalaw;
    ifstream infile;
    infile.open(input_file);
    if (!infile) {
      cerr << " EOS: Can't open " << input_file
	   << " for input. This is bad. " << endl;
    } else 
      cout << " EOS: Reading EOS polytropic parameters from file " 
	   << input_file << endl;
    char buf[100], c;
    infile.get(buf,100,'='); infile.get(c); infile >> K;
    infile.get(buf,100,'='); infile.get(c); infile >> Gamma;
    infile.close();
    n = 1.0/(Gamma - 1.0);
    cout << " EOS: Setting up gamma-law EOS with Gamma = " << Gamma << endl; 
    name << "Gamma-law EOS with Gamma = " << Gamma ;
  };
  // Destructor
  ~gamma_law() {};
  //
  // Routines
  //
  inline double P(double rho_0, double epsilon) {
    return rho_0 * epsilon * ( Gamma - 1.0 );
  };
  inline double rho_0(double P) {;
    return pow( P/K, 1.0/Gamma );
  }
  inline double sound_speed(double rho_0, double epsilon) {
    // FIX!!!
    return sqrt( Gamma * P(rho_0, epsilon) / ( h(rho_0, epsilon) * rho_0 ) );
  }
  inline double dPdrho_0(double rho_0, double epsilon) {
    return epsilon * ( Gamma - 1.0 );
  };
  inline double dPdepsilon(double rho_0, double epsilon) {
    return rho_0 * ( Gamma - 1.0 );
  };
  inline double h(double rho_0, double epsilon) {
    return 1.0 + Gamma * epsilon;
  };
  inline double cold_eps(double rho_0) {
    if (rho_0 < 0.0)
      return 0.0;
    if (Gamma == 2.0)
      return K * rho_0;
    else
      return K * pow(rho_0,Gamma-1.0)/(Gamma - 1.0);
  };
  double epsilon_root(double tau){return 0.0;};
};

//
//================================================
// derived class for Piecewise Polytrope EOS
//
//    P = K_i * rho_0 ** ( \Gamma_i )
//
//================================================
//
class piece_polytrope : public EOS {
private:
  double K_0, K_1, K_2, K_3;
  double Gamma_0, Gamma_1, Gamma_2, Gamma_3;
  double Gamma_thermal;
  double epsilon_thermal, rho_0_current; // need for rootfinding
  double rho_01, rho_02, rho_03;
  double a_0, a_1, a_2, a_3;
  double P_1, P_2, P_3;
  int Piece_polytrope_type;
  double G, c, M_sun, M_sun_g, M_sun_cm;
  double mu, alpha;   // mean molecular weight; radiation constant
  double u_F;         // constant for energy-density of semi-degenerate gas
  double e_F;         // constant so that Fermi energy is e_F \rho_0^{2/3}
  double eta;         // factor for radiative internal energy
  double PI;
public:
  // Constructor
  piece_polytrope(char * input_file) : EOS(0.0), mu(1.0), eta(11./4.) { 
    eos_type = piecepoly;
    Gamma_0 = 1.35692;
    double rho_01_cgs = 1.46220e14;
    double rho_02_cgs = 5.01187e14;
    double rho_03_cgs = 1.00000e15;
    double P_1_cgs, P_2_cgs, P_3_cgs;
    double K_0_cgs, K_1_cgs, K_2_cgs, K_3_cgs;
    PI = acos(-1.0);
    double h = 6.62607e-27;
    // what should I take as Baryon mass?
    double m_B = 1.67e-24;
    // the data from J. Read's paper is given in cgs units and the code uses geometric units (i.e. c=G=1) hence we must convert units:
    G = 6.67408e-8;
    c = 2.99792e10;
    M_sun_g = 1.98847e33;
    M_sun_cm = (G * M_sun_g)/(c * c);
    const double c2 = c * c;
    const double c4 = c2 * c2;
    const double M_sun2 = M_sun_cm * M_sun_cm;
    // cgs units for rho_01
    // There are four regions for the densities given by Gamma_0 through 3 
    // (or K_0 through 3 respectively)
    // The regions are  divided at rho_01,02 and 03 (or pressures P_1, 2, 3 respectively) 
    //   K_0   |   K_1   |   K_2   |   K_3
    // Gamma_0 | Gamma_1 | Gamma_2 | Gamma_3
    //   a_0   |   a_1   |   a_2   | a_3 
    //       rho_01     rho_02    rho_03
    //        P_1        P_2       P_3
    ifstream infile;
    infile.open(input_file);
    if (!infile) {
      cerr << " EOS: Can't open " << input_file
	   << " for input. This is bad. " << endl;
    } else 
      cout << " EOS: Reading EOS polytropic parameters from file " 
	   << input_file << endl;
    char buf[200], k;
    infile.get(buf,200,'='); infile.get(k); infile >> Piece_polytrope_type;
    infile.get(buf,200,'='); infile.get(k); infile >> Gamma_thermal;
    infile.close();
    if(Piece_polytrope_type == 1){
      Gamma_1 = 3.005;
      Gamma_2 = 2.988;
      Gamma_3 = 2.851;
      P_2_cgs = pow(10., 34.384);
      name << "Piecewise Polytrope for SLy EOS";
	}
    else if(Piece_polytrope_type == 2){
      Gamma_1 = 2.442;
      Gamma_2 = 3.256;
      Gamma_3 = 2.908;
      P_2_cgs = pow(10., 33.943);
      name << "Piecewise Polytrope for AP1  EOS";
    	}
    else if(Piece_polytrope_type == 3){
      Gamma_1 = 3.224;
      Gamma_2 = 3.033;
      Gamma_3 = 1.325;
      P_2_cgs = pow(10., 34.858);
      name << "Piecewise Polytrope for MS1 EOS";
    	}
    else if(Piece_polytrope_type == 4){
      Gamma_1 = 2.595;
      Gamma_2 = 1.845;
      Gamma_3 = 1.897;
      P_2_cgs = pow(10., 34.564);
      name << "Piecewise Polytrope for H1 EOS"; 
    } 
    else if(Piece_polytrope_type == 5){
      Gamma_1 = 2.787;
      Gamma_2 = 1.951;
      Gamma_3 = 1.901;
      P_2_cgs = pow(10., 34.646);
      name << "Piecewise Polytrope for H3 EOS"; 
    } 
    else if(Piece_polytrope_type == 6){
      Gamma_0 = Gamma_thermal;
      Gamma_1 = Gamma_thermal;
      Gamma_2 = Gamma_thermal;
      Gamma_3 = Gamma_thermal;
      P_2_cgs = pow(10., 34);
      name << "Piecewise Polytrope for Gamma-law EOS "; 
    }
    else {
      cout << " Cannot find Piecewise-polytrope type " << Piece_polytrope_type << endl;
    }
    // Compute polytropic constants
    K_2_cgs = P_2_cgs / pow(rho_02_cgs, Gamma_2);
    K_1_cgs = K_2_cgs * pow(rho_02_cgs, Gamma_2 - Gamma_1);
    K_3_cgs = K_2_cgs * pow(rho_03_cgs, Gamma_2 - Gamma_3);
    K_0_cgs = K_1_cgs * pow(rho_01_cgs, Gamma_1 - Gamma_0);
    // there's no P_0_cgs (see image/graph above)
    P_1_cgs = K_1_cgs * pow(rho_01_cgs, Gamma_1);
    P_3_cgs = K_3_cgs * pow(rho_03_cgs, Gamma_3);
    //convert to code units
    rho_01 = rho_01_cgs * G/c2 * M_sun2;
    rho_02 = rho_02_cgs * G/c2 * M_sun2;
    rho_03 = rho_03_cgs * G/c2 * M_sun2;
    //
    K_0 = K_0_cgs * pow(G, 1.0 - Gamma_0) * pow(c, 2.0 * Gamma_0 - 4.0) * pow(M_sun_cm, 2.0 - 2.0 * Gamma_0);
    K_1 = K_1_cgs * pow(G, 1.0 - Gamma_1) * pow(c, 2.0 * Gamma_1 - 4.0) * pow(M_sun_cm, 2.0 - 2.0 * Gamma_1);
    K_2 = K_2_cgs * pow(G, 1.0 - Gamma_2) * pow(c, 2.0 * Gamma_2 - 4.0) * pow(M_sun_cm, 2.0 - 2.0 * Gamma_2);
    K_3 = K_3_cgs * pow(G, 1.0 - Gamma_3) * pow(c, 2.0 * Gamma_3 - 4.0) * pow(M_sun_cm, 2.0 - 2.0 * Gamma_3);
    //
    P_1 = P_1_cgs * G/c4 * M_sun2;
    P_2 = P_2_cgs * G/c4 * M_sun2;
    P_3 = P_3_cgs * G/c4 * M_sun2;
    //
    a_0 = 0.0;
    a_1 = (K_0 * pow(rho_01, Gamma_0 - 1.0))/(Gamma_0 - 1.0) 
      -   (K_1 * pow(rho_01, Gamma_1 - 1.0))/(Gamma_1 - 1.0);
    a_2 = (K_1 * pow(rho_02, Gamma_1 - 1.0))/(Gamma_1 - 1.0) 
      -   (K_2 * pow(rho_02, Gamma_2 - 1.0))/(Gamma_2 - 1.0) + a_1;
    a_3 = (K_2 * pow(rho_03, Gamma_2 - 1.0))/(Gamma_2 - 1.0) 
      -   (K_3 * pow(rho_03, Gamma_3 - 1.0))/(Gamma_3 - 1.0) + a_2;
    //
    // define radiation constant \bar alpha
    //
    alpha = (8.0 * PI * PI * PI * PI * PI)/15.0 * (G * m_B * m_B * m_B * m_B * c)/(h * h * h) * M_sun2;
    //
    // define factor for u_nucl
    //
    const double h_bar = h / (2. * PI);
    u_F = pow(3.*PI*PI,1./3.) * pow(m_B,8./3.) * pow(G,2./3.) * pow(c,2./3.) 
      / (6 * h_bar * h_bar) * pow(M_sun_cm,4./3.);
    // 
    // define factor for Fermi energy
    //
    e_F = pow(3.0*PI*PI*(h_bar/m_B)*(h_bar/c)*(h_bar/G)/2.,2./3.)/(2.0*m_B*m_B) /  pow(M_sun_cm,4./3.);
    cout << " EOS: u_F = " << u_F << " e_F = " << e_F << endl;    
    cout << " EOS: Type " << Name()  << endl;
    cout << " EOS: Setting up PWP EOS Gamma_0 = " << Gamma_0 
	 << ", Gamma_1 = " << Gamma_1 << ", Gamma_2 = " << Gamma_2 
	 << ", Gamma_3 = " << Gamma_3 << endl;
    cout << " EOS: Using Gamma_thermal = " << Gamma_thermal << endl;
    // cout << " EOS: K_0_cgs = " << K_0_cgs << ", K_1_cgs = " << K_1_cgs << ", K_2_cgs = " << K_2_cgs << ", K_3_cgs = " << K_3_cgs << endl;
    cout << " EOS: K_0 = " << K_0 << ", K_1 = " << K_1 << ", K_2 = " << K_2 << ", K_3 = " << K_3 << endl;
    // cout << " EOS: P_1_cgs = " << P_1_cgs << ", P_2_cgs = " << P_2_cgs << ", P_3_cgs = " << P_3_cgs << endl;
    // cout << " EOS: P_1 = " << P_1 << ", P_2 = " << P_2 << ", P_3 = " << P_3 << endl;
    cout << " EOS: a_0 = " << a_0 << ", a_1 = " << a_1 << ", a_2 = " << a_2 << ", a_3 = " << a_3 << endl;
    // cout << " EOS: rho_01_cgs = " << rho_01_cgs << ", rho_02_cgs = " << rho_02_cgs << ", rho_03_cgs = " << rho_03_cgs << endl;
    // cout << " EOS: rho_01 = " << rho_01 << ", rho_02 = " << rho_02 << ", rho_03 = " << rho_03 << endl;
    infile.close(); 
  };
  // Destructor
  ~piece_polytrope() {};
  //
  // Routines
  //
  inline double eps_nucl(double tau) {
    //    cout << " eps_nucl = " << u_F * pow(rho_0_current,-2./3.) * tau * tau << endl;
    return u_F * pow(rho_0_current,-2./3.) * tau * tau;
  }
  inline double eps_rad(double tau) {
    //    cout << " eps_rad = " << eta * alpha * tau * tau * tau * tau / rho_0_current << endl;
    return eta * alpha * tau * tau * tau * tau / rho_0_current;
  }
  inline double s_nucl(double rho_0, double tau) {
    const double E_F = e_F * pow(rho_0,2./3.);
    const double s = PI*PI * tau / (2.0 * E_F);
    //    cout << " s_nucl = " << s << endl;
    return s;
  }
  inline double s_rad(double rho_0, double tau) {
    const double s = 4.0*eta/3.0 * alpha * tau * tau * tau / rho_0;
    //    cout << " s_rad = " << s << endl;
    return s;
  }
  //
  //
  //
  inline double P(double rho_0, double epsilon) {
    double epsilon_cold = cold_eps(rho_0);
    epsilon_thermal = epsilon - epsilon_cold;
    if (epsilon_thermal < 0.0) {
    //   cout << " epsilon " << setprecision(16) << epsilon << " eps_cold = " << epsilon_cold << " rho_0 = "  << rho_0 << endl;
      epsilon_thermal = 0.0;
    }
    double P_cold = 0.0;
    if (rho_0 < rho_01)
      P_cold = K_0 * pow( rho_0, Gamma_0 );
    else if (rho_0 < rho_02)
      P_cold = K_1 * pow( rho_0, Gamma_1 );
    else if (rho_0 < rho_03)
      P_cold = K_2 * pow( rho_0, Gamma_2 );
    else
      P_cold = K_3 * pow( rho_0, Gamma_3 );
    //
    // Use rootfinding to find temperature tau
    //
    rho_0_current = rho_0;
    double tau_1 = 0.9*tau_last;
    double tau_2 = tau_last + 1.0e-13;
    // make sure that roots are bracketed!
    bool success = zbrac(&EOS::epsilon_root, tau_1, tau_2);
    if (!success) {
      cout << "Could not bracket roots in piece_polytrope::P for tau1 = " << tau_1 
	   << " root: " << epsilon_root(tau_1)
	   << ", tau_2 = " << tau_2 << " root: " << epsilon_root(tau_2) 
	   << " rho_0_current = " << rho_0_current << endl;
      exit(0);
    }
    const double tol = 1.e-12 * tau_2;
    // CHECK!
    tau_last = zbrent(&EOS::epsilon_root, tau_1, tau_2, tol*tau_2);    
    //
    // now compute thermal pressure
    //
    // if (tau_last != 0.0) 
    //   cout << " eps_thermal = " << epsilon_thermal 
    // 	   << " eps_thermal - eps_nucl - eps_rad = " 
    // 	   << epsilon_thermal - eps_nucl(tau_last) - eps_rad(tau_last) 
    // 	   << " tau = " << tau_last << endl;
    double tau = tau_last;
    double P_gas = (Gamma_thermal - 1.0) * eps_nucl(tau) * rho_0;
    double P_rad = eps_rad(tau) * rho_0 / 3.0;
    //
    double P_total = P_cold + P_gas + P_rad;
    return P_total;
  };
  // only for cold EOS
  inline double rho_0(double P) {
    if (P < P_1)
      return pow( P/K_0, 1.0/Gamma_0 );
    else if (P < P_2)
      return pow( P/K_1, 1.0/Gamma_1 );
    else if (P < P_3)
      return pow( P/K_2, 1.0/Gamma_2 );
    else
      return pow( P/K_3, 1.0/Gamma_3 );
  };
  double sound_speed(double rho_0, double epsilon) {
    // CHECK!!!
    //    return sqrt( Gamma_0 * P(rho_0, epsilon) / ( h(rho_0, epsilon) * rho_0 ) );
    double p = P(rho_0,epsilon);
    const double h = 1.0 + epsilon + p / rho_0;
    const double tau = tau_last;   // computed in P(rho_0, epsilon)
    double dp_drho_0_cold = 0.0;
    if (rho_0 < rho_01)
      dp_drho_0_cold =  Gamma_0 * K_0 * pow(rho_0, Gamma_0 - 1.0);
    else if (rho_0 < rho_02)
      dp_drho_0_cold =  Gamma_1 * K_1 * pow(rho_0, Gamma_1 - 1.0);
    else if (rho_0 < rho_03)
      dp_drho_0_cold =  Gamma_2 * K_2 * pow(rho_0, Gamma_2 - 1.0);
    else
      dp_drho_0_cold =  Gamma_3 * K_3 * pow(rho_0, Gamma_3 - 1.0);
    double s_r = 0.0;
    double s_n = 0.0;
    double s_ratio = 0.0;
    if (rho_0 > 0.0 && tau > 0.0) {
      s_r = s_rad(rho_0, tau);
      s_n = s_nucl(rho_0, tau);
      s_ratio = s_r / s_n;
    } 
    double sigma_s = (3.0 + 2.0 * s_ratio)/(2.0 + 6.0 * s_ratio);
    const double fac_nucl = (1./3. + 2.0 * sigma_s) * (Gamma_thermal - 1.0);    
    rho_0_current = rho_0;    // global variable used in eps_nucl and eps_rad
    const double dp_drho_0_nucl = fac_nucl * eps_nucl(tau);
    const double dp_drho_0_rad = 4./3. * sigma_s * eps_rad(tau);
    const double dp_drho_0 = dp_drho_0_cold + dp_drho_0_nucl + dp_drho_0_rad;
    const double cs = sqrt( dp_drho_0 / h );
    if (cs >= 1.0) cs_error++;
    return cs;
  }

  inline double h(double rho_0, double epsilon) {
    return 1.0 + epsilon + P(rho_0,epsilon)/rho_0;    
    // if (rho_0 < rho_01)
    //   return 1.0 + a_0 + Gamma_0/(Gamma_0 - 1.0) * K_0 * pow( rho_0, Gamma_0 - 1.0 );
    // else if (rho_0 < rho_02)
    //   return 1.0 + a_1 + Gamma_1/(Gamma_1 - 1.0) * K_1 * pow( rho_0, Gamma_1 - 1.0 );
    // else if (rho_0 < rho_03)
    //   return 1.0 + a_2 + Gamma_2/(Gamma_2 - 1.0) * K_2 * pow( rho_0, Gamma_2 - 1.0 );
    // else
    //   return 1.0 + a_3 + Gamma_3/(Gamma_3 - 1.0) * K_3 * pow( rho_0, Gamma_3 - 1.0 );;    
  };
  inline double cold_eps(double rho_0) {
    // CHECK!!!
    // return K_0 * pow(rho_0,Gamma_0-1.0)/(Gamma_0 - 1.0);
    if (rho_0 < rho_01)
      return a_0 + K_0/(Gamma_0 - 1.0) * pow( rho_0, Gamma_0 - 1.0 );
    else if (rho_0 < rho_02)
      return a_1 + K_1/(Gamma_1 - 1.0) * pow( rho_0, Gamma_1 - 1.0 );
    else if (rho_0 < rho_03)
      return a_2 + K_2/(Gamma_2 - 1.0) * pow( rho_0, Gamma_2 - 1.0 );
    else
      return a_3 + K_3/(Gamma_3 - 1.0) * pow( rho_0, Gamma_3 - 1.0 );
  };
  //======================================
  // function to be used for rootfinding
  //======================================
  double epsilon_root(double tau) {
    const double epsilon_gas = eps_nucl(tau);
    const double epsilon_rad = eps_rad(tau);
    return epsilon_thermal - epsilon_gas - epsilon_rad;
  }
};


//
//================================================
//
// Derived class for ideal gas law.  Similar to gamma-law
// with Gamma = 5/3, but with extra functionality for temperature.
//
//================================================
//
class ideal_gas : public EOS {
public:
  // Constructor
  ideal_gas() : EOS(0.0) { 
    eos_type = idealgas;
    cout << " EOS: Setting up ideal-gas EOS " << endl;
    name << "ideal-gas EOS" ;
  };
  // Destructor
  ~ideal_gas() {};
  //
  // Routines
  //
  inline double P(double rho_0, double epsilon) {
    tau_last = epsilon / 3.0;
    return (2.0 / 3.0) * rho_0 * epsilon;
  };
  inline double rho_0(double P) {;
    // just here to make compiler happy - should never use this!
    cout << " Do not call ideal_gas::rho_0(P)... " << endl;
    exit(0);
    return 0.0;
  }
  inline double sound_speed(double rho_0, double epsilon) {
    // CHECK!!!
    return sqrt( (5.0/3.0) * P(rho_0, epsilon) / ( h(rho_0, epsilon) * rho_0 ) );
  }
  inline double dPdrho_0(double rho_0, double epsilon) {
    return (2.0 / 3.0) * epsilon;
  };
  inline double dPdepsilon(double rho_0, double epsilon) {
    return (2.0 / 3.0) * rho_0;
  };
  inline double h(double rho_0, double epsilon) {
    return 1.0 + (5.0 / 3.0 ) * epsilon;
  };
  inline double cold_eps(double rho_0) {
    // just here to make compiler happy - should never use this!
    // cout << " Do not call ideal_gas::cold_eps(rho_0)... " << endl;
    // exit(0);
    return 0.0;
  };
  double epsilon_root(double tau){return 0.0;};
};

#include "GasRad.h"



#endif   /* EOS_H */
