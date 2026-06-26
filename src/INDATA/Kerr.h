// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Kerr initial data (see Liu et.al., arxiv:1001.4077
//================================================
//
class Kerr : public InData {
private:
  double M, a, a2;   
  double r_plus, r_minus;
public:
  //================================================
  // Constructor
  //================================================
  Kerr(char * indata_input, Grid * grid_i, Cosmology *cosmology) :
    InData(grid_i, cosmology) {
    indata_type = kerr;
    analytical = false;
    ifstream infile;
    infile.open(indata_input);
    if (!infile)
      cerr << "Can't open " << indata_input 
	   << " for input. This is bad. " << endl;
    else 
      cout << " Reading initial data parameters from file " 
	   << indata_input << endl;
    char buf[100], c;
    infile.get(buf,100,'='); infile.get(c); infile >> M;
    infile.get(buf,100,'='); infile.get(c); infile >> a;
    cout << " Will set up Kerr black hole with mass M = " << M << endl;
    cout << "    and angular momentum parameter a = " << a << endl; 
    cout << "===================================================" << endl;
    a2 = a*a;
    r_plus  = M + sqrt(M*M - a2);
    r_minus = M - sqrt(M*M - a2);
  };
  //================================================
  // Destructor
  //================================================
  ~Kerr();
  string Name() { return "Kerr initial data in Liu et.al. coordinates"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) { return true; } 
  //================================================
  // Analytical solution for h_{ij}
  //================================================
  double h_rr_analytical(double r, double theta, double phi, double t) {
    const double aux = r + r_plus/4.0;
    return Sigma(r,theta)*aux*aux/(r*r*r*(r_BL(r) - r_minus))/
      psi4(r,theta) - 1.0;
  };
  double h_rt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_rp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_tt_analytical(double r, double theta, double phi, double t) {
    return Sigma(r,theta)/(r*r)/psi4(r,theta) - 1.0;
  };
  double h_tp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_pp_analytical(double r, double theta, double phi, double t) {
    return A(r,theta)/Sigma(r,theta)/(r*r)/psi4(r,theta) - 1.0;
  };
  double phi_analytical(double r, double theta, double phi, double t) {
    return log(psi4(r,theta))/4.0;
  };
  //================================================
  // Analytical solution for connection coefficients
  //================================================
  double lam_r_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = false; return 0;
  }
  double lam_t_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = false; return 0;
  }
  double lam_p_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = false; return 0;
  }
  //================================================
  // Analytical solution for extrinsic curvature
  //================================================
  double a_rr_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_rt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_rp_analytical(double r, double theta, double phi, double t) {
    const double sintheta = sin(theta);
    const double sin2theta = sintheta*sintheta;
    const double aux1 = M*a*sin2theta/
      (Sigma(r,theta) * sqrt( A(r,theta) * Sigma(r,theta) ));
    const double rBL = r_BL(r);
    const double rBL2 = rBL*rBL;
    const double aux2 = 3.0*rBL2*rBL2 + 2.0*a2*rBL2 - a2*a2 - 
      a2*(rBL2 - a2)*sin2theta;
    const double aux3 = (1.0 + r_plus/(4.0*r)) / sqrt(r*(rBL - r_minus));
    return aux1 * aux2 * aux3 / (r*sintheta) / psi4(r,theta);
  };
  double a_tt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_tp_analytical(double r, double theta, double phi, double t) {
    const double costheta = cos(theta);
    const double sintheta = sin(theta);
    const double sin2theta = sintheta*sintheta;
    const double rBL = r_BL(r);
    const double aux1 = - 2.0*M*a*a2*rBL*costheta*sintheta*sin2theta/
      (Sigma(r,theta) * sqrt( A(r,theta) * Sigma(r,theta) ));
    const double aux2 = (r - r_plus/4.0) * sqrt((rBL - r_minus)/r);
    return aux1 * aux2 / (r*r*sintheta) / psi4(r,theta);
  };
  double a_pp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double K_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  //================================================
  // Analytical solution for gauge
  //================================================
  double lapse_analytical(double r, double theta, double phi, double t) {
    return 1.0/sqrt(psi4(r,theta));
  };
  double shift_r_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double shift_t_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double shift_p_analytical(double r, double theta, double phi, double t) {
    return - 2.0*M*a*r_BL(r)/A(r,theta) * r *sin(theta) ;
    //    return - 2.0*M*a*r_BL(r)/A(r,theta) ;
  };
  //================================================
  // Analytical solution for auxiliary functions
  //================================================
  double Theta_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  }
  double B_r_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  }
  double B_t_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  }
  double B_p_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  }
  //================================================
  // Analytical solution for hydro
  //================================================
  double rho_0_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double P_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double v_r_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double v_t_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double v_p_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  //================================================
  // Analytical solution for scalar field
  //================================================
  double sf_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double pi_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  //================================================
  // Analytical solution for Maxwell OR for Dual Maxwell:
  // For Maxwell need e_p and a_p only, for Dual Maxwell all a_i and as_i
  //================================================
  double e_p_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_r_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_t_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_p_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double as_r_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double as_t_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double as_p_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  //================================================
  // Analytical solution for radiation
  //================================================
  double E_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double F_0_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  }
  double F_r_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double F_t_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double F_p_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  //================================================
  // Auxiliary functions for Kerr
  //================================================
private:
  inline double r_BL(double r) {
    const double aux = ( 1.0 + r_plus/(4.0 * r) );
    return r * aux * aux;
  }  
  inline double Sigma(double r, double theta) {
    const double rBL = r_BL(r);
    const double costheta = cos(theta);
    return rBL*rBL + a2*costheta*costheta;
  }
  inline double Delta(double r) {
    const double rBL = r_BL(r);
    return rBL*rBL - 2.0*M*rBL + a2;
  }
  inline double A(double r, double theta) {
    const double rBL = r_BL(r);
    const double sintheta = sin(theta);
    const double aux = rBL*rBL + a2;
    return aux*aux - Delta(r)*a2*sintheta*sintheta;
  }
  // returns determinant/(r^4 sin^2 theta)
  inline double det(double r, double theta) {
    const double aux = 1.0 + r_plus/(4.0*r);
    return Sigma(r,theta)*aux*aux*A(r,theta)/(r*r*r*r*r*(r_BL(r) - r_minus));
  }
  inline double psi4(double r, double theta) {
    return pow(det(r,theta),1.0/3.0);
  }
};
