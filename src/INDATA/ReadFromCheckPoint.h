// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Read data from check point files
//================================================
//
#include <cmath>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <iterator>

class ReadFromCheckPoint : public InData {
  bool all_clear;
  int nr, nt, np;
  int timestep; 
  gf3d h_rr, h_rt, h_rp, h_tt, h_tp, h_pp, phi;
  gf3d a_rr, a_rt, a_rp, a_tt, a_tp, a_pp, K;
  gf3d lam_r, lam_t, lam_p;
  gf3d lapse, shift_r, shift_t, shift_p;
  gf3d Theta, B_r, B_t, B_p;
  gf3d rho_0, p, v_r, v_t, v_p;
  gf3d sf, pi, e_p, a_p;
  gf3d E, F_0, F_r, F_t, F_p;
  gf3d D, eps;
  bool h_rr_found, h_rt_found, h_rp_found, h_tt_found, h_tp_found, h_pp_found;
  bool a_rr_found, a_rt_found, a_rp_found, a_tt_found, a_tp_found, a_pp_found;
  bool phi_found, K_found;
  bool lam_r_found, lam_t_found, lam_p_found;
  bool lapse_found, shift_r_found, shift_t_found, shift_p_found;
  bool Theta_found, B_r_found, B_t_found, B_p_found;
  bool rho_0_found, p_found, v_r_found, v_t_found, v_p_found;
  bool sf_found, pi_found, e_p_found, a_p_found;
  bool E_found, F_0_found, F_r_found, F_t_found, F_p_found;
public:
  //================================================
  // Constructor
  //================================================
  ReadFromCheckPoint(int timestep_i, Grid * grid_i, Cosmology * cosmology)
    : InData(grid_i, cosmology), timestep(timestep_i) {
    cout << " INDATA: reading data from checkpoint files at timestep "
	 << timestep << endl;
    N_g = grid->N_ghosts();
    nr = grid->N_r_tot();
    nt = grid->N_theta_tot();
    np = grid->N_phi_tot();
    indata_type = chkpt;
    all_clear = true;
    PI = acos(-1.0);
  };
  //================================================
  // Destructory
  //================================================
  ~ReadFromCheckPoint() {};
  string Name() { return "checkpoint data"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) {
    // look for checkpoint files
    // filesystem::path p(".");
    // filesystem::directory_iterator start(p);
    // filesystem::directory_iterator end;
    
    int gf_counter = 2000;
    h_rr.setup(grid, 1, "h_rr", gf_counter++, +1, +1, +1);
    h_rr_found = read_gf(&h_rr, timestep);
    h_rt.setup(grid, 1, "h_rt", gf_counter++, -1, -1, -1);
    h_rt_found = read_gf(&h_rt, timestep);
    h_rp.setup(grid, 1, "h_rp", gf_counter++, +1, -1, +1);
    h_rp_found = read_gf(&h_rp, timestep);
    h_tt.setup(grid, 1, "h_tt", gf_counter++, +1, +1, +1);
    h_tt_found = read_gf(&h_tt, timestep);
    h_tp.setup(grid, 1, "h_tp", gf_counter++, -1, +1, -1);
    h_tp_found = read_gf(&h_tp, timestep);
    h_pp.setup(grid, 1, "h_pp", gf_counter++, +1, +1, +1);
    h_pp_found = read_gf(&h_pp, timestep);
    phi.setup(grid, 1, "phi", gf_counter++, +1, +1, +1);
    phi_found = read_gf(&phi, timestep);
    //
    a_rr.setup(grid, 1, "a_rr", gf_counter++, +1, +1, +1);
    a_rr_found = read_gf(&a_rr, timestep);
    a_rt.setup(grid, 1, "a_rt", gf_counter++, -1, -1, -1);
    a_rt_found = read_gf(&a_rt, timestep);
    a_rp.setup(grid, 1, "a_rp", gf_counter++, +1, -1, +1);
    a_rp_found = read_gf(&a_rp, timestep);
    a_tt.setup(grid, 1, "a_tt", gf_counter++, +1, +1, +1);
    a_tt_found = read_gf(&a_tt, timestep);
    a_tp.setup(grid, 1, "a_tp", gf_counter++, -1, +1, -1);
    a_tp_found = read_gf(&a_tp, timestep);
    a_pp.setup(grid, 1, "a_pp", gf_counter++, +1, +1, +1);
    a_pp_found = read_gf(&a_pp, timestep);
    K.setup(grid, 1, "K", gf_counter++, +1, +1, +1);
    K_found = read_gf(&K, timestep);
    //
    lam_r.setup(grid, 1, "lam_r", gf_counter++, -1, +1, +1);
    lam_r_found = read_gf(&lam_r, timestep);
    lam_t.setup(grid, 1, "lam_t", gf_counter++, +1, -1, -1);
    lam_t_found = read_gf(&lam_t, timestep);
    lam_p.setup(grid, 1, "lam_p", gf_counter++, -1, -1, +1);
    lam_p_found = read_gf(&lam_p, timestep);
    //
    lapse.setup(grid, 1, "lapse", gf_counter++, +1, +1, +1);
    lapse_found = read_gf(&lapse, timestep);
    shift_r.setup(grid, 1, "shift_r", gf_counter++, -1, +1, +1);
    shift_r_found = read_gf(&shift_r, timestep);
    shift_t.setup(grid, 1, "shift_t", gf_counter++, +1, -1, -1);
    shift_t_found = read_gf(&shift_t, timestep);
    shift_p.setup(grid, 1, "shift_p", gf_counter++, -1, -1, +1);
    shift_p_found = read_gf(&shift_p, timestep);
    //
    Theta.setup(grid, 1, "Theta", gf_counter++, +1, +1, +1);
    Theta_found = read_gf(&Theta, timestep);
    B_r.setup(grid, 1, "B_r", gf_counter++, -1, +1, +1);
    B_r_found = read_gf(&B_r, timestep);
    B_t.setup(grid, 1, "B_t", gf_counter++, +1, -1, -1);
    B_t_found = read_gf(&B_t, timestep);
    B_p.setup(grid, 1, "B_p", gf_counter++, -1, -1, +1);
    B_p_found = read_gf(&B_p, timestep);
    //
    rho_0.setup(grid, 1, "rho_0", gf_counter++, +1, +1, +1);
    rho_0_found = read_gf(&rho_0, timestep);
    p.setup(grid, 1, "p", gf_counter++, +1, +1, +1);
    p_found = read_gf(&p, timestep);
    v_r.setup(grid, 1, "v_r", gf_counter++, -1, +1, +1);
    v_r_found = read_gf(&v_r, timestep);
    v_t.setup(grid, 1, "v_t", gf_counter++, +1, -1, -1);
    v_t_found = read_gf(&v_t, timestep);
    v_p.setup(grid, 1, "v_p", gf_counter++, -1, -1, +1);
    v_p_found = read_gf(&v_p, timestep);
    // 
    sf.setup(grid, 1, "sf", gf_counter++, +1, +1, +1);
    sf_found = read_gf(&sf, timestep);
    pi.setup(grid, 1, "pi", gf_counter++, +1, +1, +1);
    pi_found = read_gf(&pi, timestep);
    e_p.setup(grid, 1, "e_p", gf_counter++, -1, -1, +1);
    e_p_found = read_gf(&e_p, timestep);
    a_p.setup(grid, 1, "a_p", gf_counter++, -1, -1, +1);
    a_p_found = read_gf(&a_p, timestep);
    //
    E.setup(grid, 1, "E", gf_counter++, +1, +1, +1);
    E_found = read_gf(&E, timestep);
    F_0.setup(grid, 1, "F_0", gf_counter++, +1, +1, +1);
    F_0_found = read_gf(&F_0, timestep);
    F_r.setup(grid, 1, "F_r", gf_counter++, -1, +1, +1);
    F_r_found = read_gf(&F_r, timestep);
    F_t.setup(grid, 1, "F_t", gf_counter++, +1, -1, -1);
    F_t_found = read_gf(&F_t, timestep);
    F_p.setup(grid, 1, "F_p", gf_counter++, -1, -1, +1);
    F_p_found = read_gf(&F_p, timestep);
    //
    // for checking purposes:
    // 
    D.setup(grid, 1, "D", gf_counter++, +1, +1, +1);
    bool found = read_gf(&D, timestep);
    if (found) dump(&D);
    eps.setup(grid, 1, "eps", gf_counter++, +1, +1, +1);
    found = read_gf(&eps, timestep);
    if (found) dump(&eps);
    if (p_found) dump(&p);
    return true;
  };
  //
  // 
  //================================================
  double h_rr_analytical(double r, double theta, double phi, double t) {
    if (h_rr_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return h_rr(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for h_rr -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  }
  double h_rt_analytical(double r, double theta, double phi, double t) {
    if (h_rt_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return h_rt(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for h_rt -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  }
  double h_rp_analytical(double r, double theta, double phi, double t) {
    if (h_rp_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return h_rp(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for h_rp -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  }
  double h_tt_analytical(double r, double theta, double phi, double t) {
    if (h_tt_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return h_tt(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for h_tt -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  }
  double h_tp_analytical(double r, double theta, double phi, double t) {
    if (h_tp_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return h_tp(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for h_tp -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  }
  double h_pp_analytical(double r, double theta, double phi, double t) {
    if (h_pp_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return h_pp(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for h_pp -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  }
  double phi_analytical(double r, double theta, double phi_angle, double t) {
    if (phi_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi_angle);
      return phi(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for phi -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  }
  //================================================
  // Analytical solution for connection coefficients
  //================================================
  double lam_r_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    if (lam_r_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      done = true;
      return lam_r(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for lam_r -- aborting..."
	   << endl;
      exit(0);
    }
    return 0;
  }
  double lam_t_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    if (lam_t_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      done = true;
      return lam_t(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for lam_t -- aborting..."
	   << endl;
      exit(0);
    }
    return 0;
  }
  double lam_p_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    if (lam_p_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      done = true;
      return lam_p(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for lam_p -- aborting..."
	   << endl;
      exit(0);
    }
    return 0;
  }
  //================================================
  // Analytical solution for extrinsic curvature
  //================================================
  double a_rr_analytical(double r, double theta, double phi, double t) {
    if (a_rr_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return a_rr(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for a_rr -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double a_rt_analytical(double r, double theta, double phi, double t) {
    if (a_rt_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return a_rt(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for a_rt -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double a_rp_analytical(double r, double theta, double phi, double t) {
    if (a_rp_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return a_rp(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for a_rp -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double a_tt_analytical(double r, double theta, double phi, double t) {
    if (a_tt_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return a_tt(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for a_tt -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double a_tp_analytical(double r, double theta, double phi, double t) {
    if (a_tp_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return a_tp(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for a_tp -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double a_pp_analytical(double r, double theta, double phi, double t) {
    if (a_pp_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return a_pp(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for a_pp -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double K_analytical(double r, double theta, double phi, double t) {
    if (K_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return K(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for K -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  //================================================
  // Analytical solution for gauge
  //================================================
  double lapse_analytical(double r, double theta, double phi, double t) {
    if (lapse_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return lapse(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for lapse -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double shift_r_analytical(double r, double theta, double phi, double t) {
    if (shift_r_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return shift_r(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for shift_r -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double shift_t_analytical(double r, double theta, double phi, double t) {
    if (shift_t_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return shift_t(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for shift_t -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double shift_p_analytical(double r, double theta, double phi, double t) {
    if (shift_p_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return shift_p(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for shift_p -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  //================================================
  // Analytical solution for auxiliary functions
  //================================================
  double Theta_analytical(double r, double theta, double phi, double t) {
    if (Theta_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return Theta(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for Theta -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double B_r_analytical(double r, double theta, double phi, double t) {
    if (B_r_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return B_r(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for B_r -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double B_t_analytical(double r, double theta, double phi, double t) {
    if (B_t_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return B_t(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for B_t -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double B_p_analytical(double r, double theta, double phi, double t) {
    if (B_p_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return B_p(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for B_p -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  //================================================
  // Analytical solution for hydro
  //================================================
  double rho_0_analytical(double r, double theta, double phi, double t) {
    if (rho_0_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return rho_0(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for rho_0 -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double P_analytical(double r, double theta, double phi, double t) {
    if (p_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return p(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for p -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double v_r_analytical(double r, double theta, double phi, double t) {
    std::cout << "In v_r_analytical!" << std::endl;
    if (v_r_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      
      return v_r(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for v_r -- aborting..."
	   << endl;
      exit(0);
    }
     return 0.0;
  };
  double v_t_analytical(double r, double theta, double phi, double t) {
    if (v_t_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return v_t(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for v_t -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double v_p_analytical(double r, double theta, double phi, double t) {
    if (v_p_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return v_p(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for v_p -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  //================================================
  // Data for scalar field
  //================================================
  double sf_analytical(double r, double theta, double phi, double t) {
    if (sf_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return sf(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for sf -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double pi_analytical(double r, double theta, double phi, double t) {
    if (pi_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return pi(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for pi -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  //================================================
  // Data for Maxwell
  //================================================
  double e_p_analytical(double r, double theta, double phi, double t) {
    if (e_p_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return e_p(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for e_p -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double a_p_analytical(double r, double theta, double phi, double t) {
    if (a_p_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return a_p(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for a_p -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  //================================================
  // Data for radiation
  //================================================
  double E_analytical(double r, double theta, double phi, double t) {
    if (E_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return E(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for E -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double F_0_analytical(double r, double theta, double phi, double t) {
    if (F_0_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return F_0(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for F_0 -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  }
  double F_r_analytical(double r, double theta, double phi, double t) {
    if (F_r_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return F_r(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for F_r -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double F_t_analytical(double r, double theta, double phi, double t) {
    if (F_t_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return F_t(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for F_t -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  double F_p_analytical(double r, double theta, double phi, double t) {
    if (F_p_found) {
      int i = grid->i_ind(r);
      int j = grid->j_ind(theta);
      int k = grid->k_ind(phi);
      return F_p(i,j,k);
    } else {
      cerr << " INDATA: Cannot find checkpoint data for F_p -- aborting..."
	   << endl;
      exit(0);
    }
    return 0.0;
  };
  //===================================================
  // set up grid functions if check point file exists
  //===================================================
  bool read_gf(gf3d * gf, int time_step) {
    bool found = false;
    ifstream infile;
    stringstream filename;
    filename << (*gf).Name() << "_" << setfill('0') << setw(8)
    	     << timestep << ".cpt" << ends;
    infile.open(filename.str().c_str());
    cout << " INDATA: looking for checkpoint file "
    	 << filename.str().c_str() << endl;
    if (infile) {
      cout << " INDATA: reading data from checkpoint file "
	   << filename.str().c_str() << endl;
      //
      // NOTE: this logic has to match that in gridfunction.checkpoint!
      //
      double value = 0.0;
      for (int i = N_g; i < nr; i++)     
	for (int j = N_g; j < nt-N_g; j++)
	  for (int k = N_g; k < np-N_g; k++) {
	    infile >> value;
	    (*gf)[i][j][k] = value;
	  }
      found = true;
      infile.close();
    }
    return found;
  }
  //===================================================
  // dump a grid function
  //===================================================
  void dump(gf3d * fct) {
  //================================================
  // find number of gridpoints
  //================================================
  int N_r = fct->dim1();
  int N_theta = fct->dim2();
  int N_phi = fct->dim3();
  int N_g = fct->N_ghosts();    fct->fill_ghosts();
    //================================================
    // first write out rays
    //================================================
    ofstream outfile;
    ostringstream rayfilename;
    rayfilename << fct->Name() << "_rays_" << N_r - 2*N_g 
		<< "_" << N_theta - 2*N_g 
		<< "_" << setfill('0') << setw(8)
		<< timestep << "_chkpt" << ends;
    outfile.open(rayfilename.str().c_str());
    if (!outfile) cerr << " Could not open file " 
		       << rayfilename.str().c_str() << endl;
    outfile.setf(ios::right);
    double tl = 0.0;
    double ta = PI;  // CHECK: adjust according to preference...  
    int np = N_phi/2;
    double pl = fct->phi(np);
    outfile << "# data 1 at theta = " << tl << " and phi = " << pl << endl;
    outfile << "# data 2 at theta = " << ta << " and phi = " << pl << endl;
    outfile << "# " << setw(22) << "r" << setw(24) << fct->Name() 
	    << setw(24) << fct->Name() << endl;
    outfile << "# " << setw(31) << " (" << setw(6) << tl << "," << pl << ")"
	    << "        (" << setw(6) << ta << "," << pl << ")" << endl;
    outfile << "#========================================================================" << endl;
    outfile.setf(ios::right);
    for (int i = 0; i < N_r; i++) {
      const double rl = fct->r(i);
      outfile << setprecision(16) << setw(24) << rl
	      << setw(24) << (*fct)(i,tl,np) 
	      << setw(24) << (*fct)(i,ta,np) << endl; 
    }
    outfile.close();
    //================================================
    // now write out slice
    //================================================
    // ofstream outfile;
    ostringstream slicefilename;
    slicefilename << fct->Name() << "_slice_" << N_r - 2*N_g 
		<< "_" << N_theta - 2*N_g 
		<< "_" << setfill('0') << setw(8)
		<< timestep << "_chkpt" << ends;
    outfile.open(slicefilename.str().c_str());
    if (!outfile) cerr << " Could not open file " 
		       << slicefilename.str().c_str() << endl;
    outfile.setf(ios::right);
    outfile << "# data for phi = " << pl << endl;
    outfile << "# " << setw(14) << "r" << setw(16) << "theta"
	    << setw(20) << fct->Name() << endl;
    outfile << "#===============================================================" << endl;
    outfile.setf(ios::right);
    for (int i = N_g; i < N_r; i++) {
      for (int j = N_g; j < N_theta - N_g; j++) { 
	double rl = fct->r(i);
	double thetal = fct->theta(j);
	outfile << setprecision(8) << setw(16) << rl 
		<< setw(16) << thetal 
		<< setprecision(16)
		<< setw(24) << (*fct)[i][j][np] << endl;
      }
    outfile << endl;
    }
  }
};

                                                                         
