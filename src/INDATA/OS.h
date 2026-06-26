// Tell emacs that this is -*-c++-*- mode
//
//================================================
// OS Initial data 
//
// Reference: Staley et.al., CQG 29, 015003 (2012)
//================================================
//

class OS : public InData {
private:
  double M, R0, rho, r_iso;
  double x_C, y_C, z_C;
  double P_over_rho;
  double E_over_rho;   // for heated OS collapse
public:
  //================================================
  // Constructor
  //================================================
  OS(char * indata_input, Grid * grid_i, Cosmology * cosmology) : 
    InData(grid_i, cosmology) {
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
    // double rho_c;
    infile.get(buf,100,'='); infile.get(c); infile >> M;
    infile.get(buf,100,'='); infile.get(c); infile >> R0;
    infile.get(buf,100,'='); infile.get(c); infile >> x_C;
    infile.get(buf,100,'='); infile.get(c); infile >> y_C;
    infile.get(buf,100,'='); infile.get(c); infile >> z_C;
    infile.get(buf,100,'='); infile.get(c); infile >> P_over_rho;
    infile.get(buf,100,'='); infile.get(c); infile >> E_over_rho;
    // compute rest mass density
    rho = 3.0*M/(4.0*PI*R0*R0*R0);
    // compute isotropic radius
    r_iso = R0 * (1.0 - M/R0 + sqrt(1.0 - 2.0*M/R0))/2.0;
    cout << " Will set up OS initial data with mass M = " << M << ", " << endl;
    cout << "    areal radius R0 = " << R0 << ", isotropic radius r_iso = " << r_iso << " and density rho_0 = " << rho << endl;
    cout << "    centered on x = " << x_C << ", y = " << y_C << ", z = " << z_C << endl;
    cout << "    initial P / rho_0 = " << P_over_rho << endl;
    cout << "    initial E / rho_0 = " << E_over_rho << endl;
    cout << "===================================================" << endl;
  };
  //================================================
  // Destructor
  //================================================
  ~OS() { };
  string Name() { return "OS initial data"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) { return true; }; 
  //================================================
  // Analytical solution for h_{ij}
  //================================================
  double h_rr_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  }
  double h_rt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  }
  double h_rp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  }
  double h_tt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  }
  double h_tp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  }
  double h_pp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  }
  double phi_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C));
    double psi;
    if (r_C > r_iso) {  // see eq. (42) in Staley et.al.
      psi = 1.0 + M/(2.0*r_C);
    } else {
      const double temp1 = (1.0 + sqrt(1.0 - 2.0*M/R0)) * r_iso*R0*R0;
      const double temp2 = 2.0*r_iso*r_iso*r_iso + M*r_C*r_C;
      psi = sqrt(temp1/temp2);
    }
    return log(psi);
  }
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
    return 1.0;
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
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C));
    if (r_C > r_iso)
      return 0.0;
    else
      return rho;
  };
  double P_analytical(double r, double theta, double phi, double t) {
    return P_over_rho * rho_0_analytical(r, theta, phi, t);
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
    return E_over_rho * rho_0_analytical(r, theta, phi, t);
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
};
