// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Ken's trumpet solution
//================================================
//
class KenTrumpet : public InData {
private:
  double M;   
  double x_C, y_C, z_C;
public:
  //================================================
  // Constructor
  //================================================
  KenTrumpet(char * indata_input, Grid * grid_i, Cosmology *cosmology) :
    InData(grid_i, cosmology) {
    indata_type = schwarzschild;
    analytical = true;
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
    infile.get(buf,100,'='); infile.get(c); infile >> x_C;
    infile.get(buf,100,'='); infile.get(c); infile >> y_C;
    infile.get(buf,100,'='); infile.get(c); infile >> z_C;
    cout << " Will set up Ken's trumpet solution with mass M = " << M << endl;
    cout << "    centered at x = " << x_C << ", y = " << y_C << ", z = " << z_C << endl;
    cout << "===================================================" << endl;
  };
  //================================================
  // Destructor
  //================================================
  ~KenTrumpet();
  string Name() { return "Ken's trumpet"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) { return true; } 
  //================================================
  // Analytical solution for h_{ij}
  //================================================
  double h_rr_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_rt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_rp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_tt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_tp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_pp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double phi_analytical(double r, double theta, double phi, double t) {
    return log(psi(r,theta,phi));
  };
  //================================================
  // Analytical solution for connection coefficients
  //================================================
  double lam_r_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = true; return 0;
  }
  double lam_t_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = true; return 0;
  }
  double lam_p_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = true; return 0;
  }
  //================================================
  // Analytical solution for extrinsic curvature
  //================================================
  double a_rr_analytical(double r, double theta, double phi, double t) {
    const double r_C_l = r_C(r,theta,phi);
    const double M_plus_r = M + r_C_l;
    return - 4.0/3.0 * M/(M_plus_r * M_plus_r);
  };
  double a_rt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_rp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_tt_analytical(double r, double theta, double phi, double t) {
    const double r_C_l = r_C(r,theta,phi);
    const double M_plus_r = M + r_C_l;
    return 2.0/3.0 * M/(M_plus_r * M_plus_r);
  };
  double a_tp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_pp_analytical(double r, double theta, double phi, double t) {
    const double r_C_l = r_C(r,theta,phi);
    const double M_plus_r = M + r_C_l;
    return 2.0/3.0 * M/(M_plus_r * M_plus_r);
  };
  double K_analytical(double r, double theta, double phi, double t) {
    const double r_C_l = r_C(r,theta,phi);
    const double M_plus_r = M + r_C_l;
    return M/(M_plus_r * M_plus_r);
  };
  //================================================
  // Analytical solution for gauge
  //================================================
  double lapse_analytical(double r, double theta, double phi, double t) {
    const double r_C_l = r_C(r,theta,phi);
    return r_C_l/(M + r_C_l);
  };
  double shift_r_analytical(double r, double theta, double phi, double t) {
    const double r_C_l = r_C(r,theta,phi);
    const double M_plus_r = M + r_C_l;
    return M * r_C_l/(M_plus_r * M_plus_r);
  };
  double shift_t_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double shift_p_analytical(double r, double theta, double phi, double t) {
    return 0.0;
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
  // Analytical solution for Maxwell
  //================================================
  double e_p_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_p_analytical(double r, double theta, double phi, double t) {
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
  // Auxiliary functions for Ken's trumpet
  //================================================
private:
  inline double r_C(double r, double theta, double phi) {
    const double x_l = r * sin(theta) * cos(phi);
    const double y_l = r * sin(theta) * sin(phi);
    const double z_l = r * cos(theta);
    return sqrt( (x_l - x_C)*(x_l - x_C) + (y_l - y_C)*(y_l - y_C) + (z_l - z_C)*(z_l - z_C) ); 
  }
  inline double psi(double r, double theta, double phi) {
    return sqrt( 1.0 + M / r_C(r,theta,phi) );
  }
  inline double psi4(double r, double theta, double phi) {
    double psil = psi(r,theta,phi);
    return psil*psil*psil*psil;
  }
};
