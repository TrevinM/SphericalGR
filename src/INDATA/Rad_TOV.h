// Tell emacs that this is -*-c++-*- mode
//
//================================================
// TOV Initial data (constructed by integrating TOV equations)
// for fluid dominated by radiation pressure 
//================================================
//
#include "Rad_TOV_Solution.h"

class Rad_TOV : public InData {
private:
  Rad_TOV_Solution *tov_sol;
  double x_C, y_C, z_C;
public:
  //================================================
  // Constructor
  //================================================
  Rad_TOV(char * indata_input, Grid * grid_i, Cosmology * cosmology) : 
    InData(grid_i, cosmology) {
    indata_type = rad_tov;
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
    double rho_c;
    infile.get(buf,100,'='); infile.get(c); infile >> rho_c;
    infile.get(buf,100,'='); infile.get(c); infile >> Gamma; // irrelevant
    infile.get(buf,100,'='); infile.get(c); infile >> Kappa; // irrelevant
    infile.get(buf,100,'='); infile.get(c); infile >> x_C;
    infile.get(buf,100,'='); infile.get(c); infile >> y_C;
    infile.get(buf,100,'='); infile.get(c); infile >> z_C;
    cout << " Will set up Rad_TOV initial data with central energy density = " << rho_c << endl;
    cout << "    centered on x = " << x_C << ", y = " << y_C << ", z = " << z_C << endl;
    cout << "===================================================" << endl;
    int N_array = 1000;
    double dr_init = 1.e-2;
    tov_sol = new Rad_TOV_Solution(rho_c, dr_init, N_array); 
  };
  //================================================
  // Destructor
  //================================================
  ~Rad_TOV() { delete tov_sol; };
  string Name() { return "Rad_TOV initial data"; };
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
    return log(tov_sol->psi(r_C));
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
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C));
    return tov_sol->lapse(r_C);
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
    return tov_sol->rho_i(r_C);
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
};
