// Tell emacs that this is -*-c++-*- mode
//
//================================================
// TOV with BH Initial data (constructed by solving Hamiltonian constraint and integrating TOV equations) 
//================================================
//
#include "General_TOV_Solution.h"

class TOV_BH : public InData {
private:
  Gen_TOV_Solution *tov_sol;
  EOS * eos;
  double x_C, y_C, z_C;
  double m_BH;
  double P_BH;
  double n_scale;
  bool all_clear;
# ifndef NoEllSolver
  FlatEllSolver3D * ellsolver;
# endif
  gf3d psi_NS, psi, u, delta_u, rho_bar, Oo_alpha, f, res_gr;
  int n_r, n_theta, n_phi;
  int N_g;
  double tol, tol_res, tol_tri;
  int max_it, max_step, max_iter;
  double eps;
  double PI;
public:
  //================================================
  // Constructor
  //================================================
  TOV_BH(char * indata_input, EOS * eos_i, Grid * grid_i, Cosmology * cosmology) 
    : InData(grid_i, cosmology), eos(eos_i) {
    N_g = grid->N_ghosts();
    indata_type = tov_bh;
    analytical = false;
    all_clear = true;
    ifstream infile;
    infile.open(indata_input);
    if (!infile) {
      cerr << "Can't open " << indata_input 
	   << " for input. This is bad. " << endl;
      all_clear = false;
    } else 
      cout << " Reading initial data parameters from file " 
	   << indata_input << endl;
    char buf[100], c;
    double rho_c;
    infile.get(buf,100,'='); infile.get(c); infile >> rho_c;
    infile.get(buf,100,'='); infile.get(c); infile >> x_C;
    infile.get(buf,100,'='); infile.get(c); infile >> y_C;
    infile.get(buf,100,'='); infile.get(c); infile >> z_C;
    infile.get(buf,100,'='); infile.get(c); infile >> m_BH;
    infile.get(buf,100,'='); infile.get(c); infile >> P_BH;
    infile.get(buf,100,'='); infile.get(c); infile >> n_scale; // scaling exponent of conformal factor
    infile.get(buf,100,'='); infile.get(c); infile >> tol;
    infile.get(buf,100,'='); infile.get(c); infile >> max_it;
    if (!eos) {
      cerr << " Can't set up TOV star without an eos..." << endl;
      exit(1);
    }
    cout << " Will set up TOV star with BH initial data with central density = " << rho_c << endl;
    cout << "    using EOS " << eos->Name() << endl;
    cout << "    centered on x = " << x_C << ", y = " << y_C << ", z = " << z_C << endl;
    cout << "    with BH of mass =  "  << m_BH  << " and momentum P_BH = " << P_BH << endl;
    cout << "    and conformal factor scaled by n = " << n_scale << endl;
    cout << "    accurate within tolerance = " << tol << endl;
    cout << "    completing max iterations = " << max_it << endl;
    cout << "===================================================" << endl;
    int N_array = 10000;
    double dr_init = 1.e-4;
    tov_sol = new Gen_TOV_Solution(eos, rho_c, dr_init, N_array);
    PI = acos(-1.0);
  };
  //================================================
  // Destructory
  //================================================
  ~TOV_BH() { delete tov_sol; };
  string Name() { return "TOV with BH initial data"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) {
    if (!all_clear) return false;
    cout << " INDATA: Initializing TOV star with BH initial data... " << endl;
    n_r = fct.dim1();
    n_theta = fct.dim2();
    n_phi = fct.dim3();
    //
    // set up grid functions
    //
    int gf_counter = 2000;
    psi_NS.setup(grid, 1, "psi_NS", gf_counter++, +1, +1, +1);
    psi.setup(grid, 1, "psi", gf_counter++, +1, +1, +1);
    rho_bar.setup(grid, 1, "rho_bar", gf_counter++, +1, +1, +1);
    u.setup(grid, 1, "u", gf_counter++, +1, +1, +1); 
    delta_u.setup(grid, 1, "delta_u", gf_counter++, +1, +1, +1);
    f.setup(grid, 1, "function", gf_counter++, +1, +1, +1);
    res_gr.setup(grid, 1, "residual_grid", gf_counter++, +1, +1, +1);
    Oo_alpha.setup(grid, 1, "one_over_alpha", gf_counter++, +1, +1, +1);
    u.equals(0.0); // initialize u
    tol_res = tol; // residual tolerance
    tol_tri = tol; // trilinos tolerance 
    max_step = max_it; // iterations for updating u 
    max_iter = max_it; // iterations for trilinos
    // 
    // set up TOV variables
    //
    bool success = setup_TOV(); // find psi_NS, rho_bar, one over alpha
    if (!success) return success;
    //
    // initialize --> solve Hamiltonian constraint
    //
    success = Solve_Hamiltonian();
    return success;
  };
  //
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
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    return log(psi(i, j, k));
    /*
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C));
    return log(tov_sol->psi(r_C));
    */ 
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
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double rl = grid->r(i);
    const double ctl = grid->costheta(j);
    const double psil = psi(i,j,k);
    const double psi6 = psil*psil*psil*psil*psil*psil;
    return 3.0*P_BH*ctl/(rl*rl*psi6);
  };
  double a_rt_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double rl = grid->r(i);
    const double stl = grid->sintheta(j);
    const double psil = psi(i,j,k);
    const double psi6 = psil*psil*psil*psil*psil*psil;
    return - 3.0*P_BH*stl/(2.0*rl*rl*psi6);
  };
  double a_rp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_tt_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double rl = grid->r(i);
    const double ctl = grid->costheta(j);
    const double psil = psi(i,j,k);
    const double psi6 = psil*psil*psil*psil*psil*psil;
    return - 3.0*P_BH*ctl/(2.0*rl*rl*psi6);
  };
  double a_tp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_pp_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double rl = grid->r(i);
    // const double stl = grid->sintheta(j);
    // const double st2 = stl*stl;
    const double ctl = grid->costheta(j);
    const double psil = psi(i,j,k);
    const double psi6 = psil*psil*psil*psil*psil*psil;
    return - 3.0*P_BH*ctl/(2.0*rl*rl*psi6);
  };
  double K_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  //================================================
  // Analytical solution for gauge
  //================================================
  double lapse_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    return 1.0/(psi(i, j, k)*psi(i, j, k));
    /*
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C));
    return tov_sol->lapse(r_C);
    */
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
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double r_l = psi_NS.r(i);
    const double ctl = psi_NS.costheta(j);
    const double stl = psi_NS.sintheta(j);
    const double z_l = r_l * ctl;
    const double s_l = r_l * stl;
    const double r_NS = sqrt((z_l - z_C)*(z_l - z_C) + s_l*s_l); 
    const double rho_0_NS = tov_sol->rho_0(r_NS);
    const double eps_NS = eos->cold_eps(rho_0_NS);
    const double rho_NS = rho_0_NS * (1.0 + eps_NS);
    const double rho_NS_BH = pow((psi(i, j, k)/psi_NS(i, j ,k)), n_scale) 
      * rho_NS;
    // double rho_0 = -0.5 + sqrt(0.25 + rho_NS);
    return eos->rho_0_of_rho(rho_NS_BH);
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

private:
  //====================================================
  // Initialize TOV star
  // ====================================================     
  bool setup_TOV() {
    if (x_C != 0.0 || y_C != 0.0) {
      cout << " TOV_BH: Center of neutron star has to be placed on z-axis! " << endl;
      return false;
    }
    for (int i = 0; i < n_r; i++)
      for (int j = 0; j < n_theta; j++)
        for (int k = 0; k < n_phi; k++) {
	  const double r_l = psi_NS.r(i);
	  const double ctl = psi_NS.costheta(j);
	  const double stl = psi_NS.sintheta(j);
	  const double z_l = r_l * ctl;
	  const double s_l = r_l * stl;
	  const double r_NS = sqrt((z_l - z_C)*(z_l - z_C) + s_l*s_l); 
	  const double psi_l = tov_sol->psi(r_NS);
	  const double rho_0 = tov_sol->rho_0(r_NS); 
	  const double rho_i = compute_rho_i(rho_0);
	  const double rho_l = rho_0 + rho_i;
	  const double rho_bar_l = rho_l * pow(psi_l,-n_scale);
	  const double Oo_alpha_l = m_BH/(2.0*r_l);
	  if (!finite(rho_bar_l)) {
	    cout << " OOOOOOOOOOOOps in Setup_TOV" << endl;
	    cout << "rho_bar_l = " << rho_bar_l << endl;
	    cout << "rhos = " << rho_0 << "  " << rho_i << endl;
	    exit(0);
	  }
	  psi_NS[i][j][k] = psi_l;
	  rho_bar[i][j][k] = rho_bar_l;
	  Oo_alpha[i][j][k] = Oo_alpha_l;
	}
    return true;
  };  

  double compute_rho_i(double rho_0) {
    //    double P = Kappa * pow(rho_0, Gamma); 
    //    double rho_i = P/(Gamma-1.0);
    double eps = eos->cold_eps(rho_0);
    return rho_0 * eps;
  };
  //====================================================                                        
  // Solve Hamiltonian constraint                                                               
  //====================================================                                        
  bool Solve_Hamiltonian() {
#ifndef NoEllSolver
    ellsolver = new FlatEllSolver3D(grid);
    //                                                                         
    // Allocate and set up elliptic solver                                    
    //   
    int step = 0;
    compute_psi();
    double res_norm = Hamiltonian_Residual();
    cout << " TOV star with BH: initial residual: " << res_norm << endl;
    while (res_norm > tol_res && step < max_step) {
      step++;
      compute_f(); // compute f
      ellsolver->SetupSolver(1.0,f);
      ellsolver->SetRHS(-1.0,res_gr);
      //                                                                      
      // Solve and get solution                                                
      //                                                                       
      int num_it;
      const double tol_tri = tol;
      ellsolver->Solve(max_iter,num_it,tol_tri);
      cout << " Solved Hamiltonian constraint to residual " << res_norm  <<
	" in " << num_it << " iterations. " << endl;             
      ellsolver->GetSolution(delta_u); //compute delta_u
      update_u(); //update u
      compute_psi(); //update psi
      res_norm = Hamiltonian_Residual(); //update residual
      cout << " TOV star with BH: residual after " << step << " steps: " << res_norm << endl;
    }
    delete ellsolver;
    bool success = false;
    if (res_norm < tol) success = true;
    return success; 
#else
    cout << " Can't construct TOV with BH initial data without an Elliptic Solver!! " << endl;
    return false;
#endif
  };
   
  double Hamiltonian_Residual() {
    for (int i = N_g; i < n_r-N_g; i++)
      for (int j = N_g; j < n_theta-N_g; j++)
	for (int k = N_g; k < n_phi-N_g; k++) {
	  const double rl = grid->r(i);
	  const double ctl = grid->costheta(j);
	  const double psi_l = psi(i, j, k);
	  const double psi_NS_l = psi_NS(i, j, k);
	  const double psi5pn = pow(psi_l, 5+n_scale); // 5pn: 5 + n
	  const double psi7 = psi_l*psi_l*psi_l*psi_l*psi_l*psi_l*psi_l;
	  const double psi_NS5pn = pow(psi_NS_l, 5+n_scale);
	  const double A2 = 9.0*P_BH*P_BH/(2.0*rl*rl*rl*rl) * (1.0 + 2.0 * ctl*ctl); 
	  const double h_of_u = 2.0*PI*rho_bar[i][j][k]*(psi_NS5pn-psi5pn) - A2 / psi7 / 8.0;
	  res_gr[i][j][k] = u.Laplace(i,j,k) - h_of_u;
	  if (!finite(res_gr[i][j][k])) {
	    cout << "Oooops!  " << i << " " << j << " " << k << endl;
	    cout << "Laplace = " << u.Laplace(i,j,k) << ", h_of_u = " << h_of_u << endl; 
	    cout << "psi5pn = " << psi5pn << ", psi_NS5pn = " << psi_NS5pn << endl;
	    cout << "rho_bar = " << rho_bar[i][j][k] << " psi7 = " << psi7 << endl;
	    exit(0);
	  }
	}
    return res_gr.L2_norm();
  };

  void compute_f() {
    for (int i = 0; i < n_r; i++)
      for (int j = 0; j < n_theta; j++)
        for (int k = 0; k < n_phi; k++) {
	  const double rl = grid->r(i);
	  const double ctl = grid->costheta(j);
	  const double psi_l = psi(i, j, k);
	  const double psi8 = psi_l*psi_l*psi_l*psi_l*psi_l*psi_l*psi_l*psi_l;
	  const double A2 = 9.0*P_BH*P_BH/(2.0*rl*rl*rl*rl) * (1.0 + 2.0 * ctl*ctl); 
	  const double rho_bar_l = rho_bar(i, j, k);
	  f[i][j][k] = (5.0+n_scale)*2.0*PI*rho_bar_l*pow(psi_l,4+n_scale) - 7.0*A2/psi8/8.0;
	}
  };

  void compute_psi() {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++)
	for (int k = 0; k < n_phi; k++) {
          const double psi_NS_l = psi_NS[i][j][k];
          const double Oo_alpha_l = Oo_alpha[i][j][k];
          const double u_l = u[i][j][k];
	  psi[i][j][k] = psi_NS_l + Oo_alpha_l + u_l;
	}
  };

  void update_u() {
     for (int i = 0; i < n_r; i++)
       for (int j = 0; j < n_theta; j++)
         for (int k = 0; k < n_phi; k++) {
           u[i][j][k] += delta_u(i,j,k);
         }
  };
};

                                                                         
