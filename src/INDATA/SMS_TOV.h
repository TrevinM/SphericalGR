// Tell emacs that this is -*-c++-*- mode
//
//================================================
// TOV Initial data for SMS stars
//
// solves TOV equations for combinations of radiation and gas 
// pressure (and internal energy densities), but then assigns them
// separately to hydro and radiation.
//
//================================================
//
#include "SMS_TOV_Solution.h"

class SMS_TOV : public InData {
private:
  SMS_TOV_Solution *tov_sol;
  EOS * eos;
  double x_C, y_C, z_C;
  bool all_clear;
public:
  //================================================
  // Constructor
  //================================================
  SMS_TOV(char * indata_input, EOS * eos_i, Grid * grid_i, Cosmology * cosmology) 
    : InData(grid_i, cosmology), eos(eos_i) {
    indata_type = tov;
    analytical = false;
    all_clear = true;
    // 
    // sanity check
    //
    if (eos->type() != 3) {
      cout << " SMS_TOV: always use ideal gas EOS for SMS_TOV! " << endl;
      exit(0);
      all_clear = false;
    }
    ifstream infile;
    infile.open(indata_input);
    if (!infile) {
      cerr << "Can't open " << indata_input 
	   << " for input. This is bad. " << endl;
      all_clear = false;
    } else 
      cout << " SMS_TOV: Reading initial data parameters from file " 
	   << indata_input << endl;
    char buf[100], c;
    double rho_c;  
    double s;     // total entropy per baryon
    double dr_init = 0.0;   // initial step size in TOV solver
    infile.get(buf,100,'='); infile.get(c); infile >> rho_c;
    infile.get(buf,100,'='); infile.get(c); infile >> s;
    infile.get(buf,100,'='); infile.get(c); infile >> x_C;
    infile.get(buf,100,'='); infile.get(c); infile >> y_C;
    infile.get(buf,100,'='); infile.get(c); infile >> z_C;
    infile.get(buf,100,'='); infile.get(c); infile >> dr_init;
    cout << " SMS_TOV:  Will set up SMS-TOV initial data with central density = " << rho_c << " M_sun^{-2} " << endl;
    cout << " SMS_TOV:  for entropy s = " << s << " k_B per baryon," << endl;
    cout << " SMS_TOV:  using initial step-size " << dr_init << endl;
    cout << " SMS_TOV:  star centered on x = " << x_C << ", y = " << y_C << ", z = " << z_C << endl;
    cout << "===================================================" << endl;
    int N_array = 10000;
    tov_sol = new SMS_TOV_Solution(rho_c, s, dr_init, N_array); 
  };
  //================================================
  // Destructor
  //================================================
  ~SMS_TOV() { delete tov_sol; };
  string Name() { return "SMS_TOV initial data"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) { return all_clear; }; 
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
    return tov_sol->rho_0(r_C);
  };
  double P_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C));
    return tov_sol->p_gas(r_C);
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
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C));
    return tov_sol->E(r_C);
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
