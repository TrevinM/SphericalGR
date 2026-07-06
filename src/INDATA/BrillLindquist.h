// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Initial data for head-on collision of two black holes: Brill-Lindquist data
//================================================
//
class Brill_Lindquist : public InData {
private:
  double M_1, M_2;   
  double x_1, y_1, z_1;
  double x_2, y_2, z_2;
public:
  //================================================
  // Constructor
  //================================================
  Brill_Lindquist(char * indata_input, Grid * grid_i, Cosmology * cosmology) :
    InData(grid_i, cosmology) {
    indata_type = brill_lindquist;
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
    infile.get(buf,100,'='); infile.get(c); infile >> M_1;
    infile.get(buf,100,'='); infile.get(c); infile >> x_1;
    infile.get(buf,100,'='); infile.get(c); infile >> y_1;
    infile.get(buf,100,'='); infile.get(c); infile >> z_1;
    infile.get(buf,100,'='); infile.get(c); infile >> M_2;
    infile.get(buf,100,'='); infile.get(c); infile >> x_2;
    infile.get(buf,100,'='); infile.get(c); infile >> y_2;
    infile.get(buf,100,'='); infile.get(c); infile >> z_2;
    cout << " Will set up Brill-Lindquist initial data for: " << endl;
    cout << "    black hole 1:  mass = " << M_1 << ", centered at x = " << x_1 << ", y = " << y_1 << ", z = " << z_1 << endl;
    cout << "    black hole 2:  mass = " << M_2 << ", centered at x = " << x_2 << ", y = " << y_2 << ", z = " << z_2 << endl;
    cout << "===================================================" << endl;
  };
  //================================================
  // Destructor
  //================================================
  ~Brill_Lindquist();
  string Name() { return "Brill-Lindquist initial data"; };
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
    return 0.0;
  };
  double a_rt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_rp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_tt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_tp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
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
    const double psi_l = psi(r,theta,phi);
    return 1.0/(psi_l*psi_l);
  };
  double shift_r_analytical(double r, double theta, double phi, double t) {
    return 0.0;
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
  // Auxiliary functions for Schwarzschild
  //================================================
private:
  inline double psi(double r, double theta, double phi) {
    const double x_l = r * sin(theta) * cos(phi);
    const double y_l = r * sin(theta) * sin(phi);
    const double z_l = r * cos(theta);
    const double r_1 = sqrt( (x_l - x_1)*(x_l - x_1) + (y_l - y_1)*(y_l - y_1) + (z_l - z_1)*(z_l - z_1) ); 
    const double r_2 = sqrt( (x_l - x_2)*(x_l - x_2) + (y_l - y_2)*(y_l - y_2) + (z_l - z_2)*(z_l - z_2) ); 
    return 1.0 + M_1/(2.0*r_1) + M_2/(2.0*r_2);
  }
  inline double psi4(double r, double theta, double phi) {
    double psil = psi(r,theta,phi);
    return psil*psil*psil*psil;
  }
};
