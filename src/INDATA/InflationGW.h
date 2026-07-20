//
// Tell emacs that this is -*-c++-*- mode
//================================================
// Initial for Inflaton field with a gravitational wave distrubance
//================================================
//
#include "../scalarpotential.h"

class InflationGW : public InData {
private:
#ifndef NoEllSolver
	FlatEllSolver3D* laplace;
	VecLaplacian* veclaplacian;
#endif
	gf3d psi, rho, res, u, delta_psi, sf, K; 
	// all vector and tensor components rescaled...
    // gf3d V;
	gf3d s_r, s_t, s_p;     // components of momentum densities (upstairs)
	gf3d W_r, W_t, W_p;     // vector potential for ext.  curv. (indices up)
	gf3d del_W_r, del_W_t, del_W_p;  // corrections...
	gf3d res_r, res_t, res_p;  // (upstairs...)
	// now extrinsic curvature: use *initial data* rescaling here,
	// and *upper* (geometrically) rescaled indices, i.e.
	// \bar A^{ij} = \psi^10 A^{ij}
	// need to (a) lower indices, and (b) convert to BSSN rescaling in
	// functions a_ij_analytical...
	gf3d A_rr, A_rt, A_rp, A_tt, A_tp, A_pp, A2;
	int n_r, n_theta, n_phi;
	int N_g;
	int max_it;   // parameters for elliptic solver
    int nakamura_type; //0 for A, 1 for B. See BGH 2026
	int branch;
	double tol;
	double phi_0, sigma, epsilon, m, GW_amp, step_factor, mom_step_factor;   // parameters for initial data
	double PI;
	bool all_clear;
	QuadraticPotential* potential;
	ostringstream indata_name;
	VecDoub psi_r;
public:
	//================================================
	// Constructor
	//================================================
	InflationGW(char* indata_input, Grid* grid_i, Cosmology* cosmology) :
		InData(grid_i, cosmology) {
		N_g = grid->N_ghosts();
		indata_type = em_wave;
		all_clear = true;
		ifstream infile;
		infile.open(indata_input);
		if (!infile) {
			cerr << " INFLAIONGW Can't open " << indata_input
				<< " for input. This is bad. " << endl;
			all_clear = false;
		}
		else
			cout << " INFLATIONGW: Reading initial data parameters from file "
			<< indata_input << endl;
		char buf[100], c;
		infile.get(buf, 100, '='); infile.get(c); infile >> phi_0;
		infile.get(buf, 100, '='); infile.get(c); infile >> sigma;
		infile.get(buf, 100, '='); infile.get(c); infile >> epsilon;
		infile.get(buf, 100, '='); infile.get(c); infile >> branch;

		infile.get(buf, 100, '='); infile.get(c); infile >> nakamura_type;
		infile.get(buf, 100, '='); infile.get(c); infile >> GW_amp;

		infile.get(buf, 100, '='); infile.get(c); infile >> max_it;
		infile.get(buf, 100, '='); infile.get(c); infile >> tol;

		infile.get(buf, 100, '='); infile.get(c); infile >> step_factor;
		infile.get(buf, 100, '='); infile.get(c); infile >> mom_step_factor;


		cout << " INFLATIONGW: Will set up Inflation GW initial data with" << endl;
		cout << "      inflaton parameters phi_0 = " << phi_0
			<< ", sigma = " << sigma
			<< ", epsilon = " << epsilon
			<< ", branch = " << branch << endl;

		cout << "      and Nakamura type = " << nakamura_type << " with amplitude " 
            << GW_amp << endl;
		cout << " INFLATIONGW: Will run elliptic solver with max_it = " << max_it
			<< " and tol = " << tol << endl;
		
		cout << " INFLATIONGW: Solving with step factor = " << step_factor 
		<< " and momentum step factor = " << mom_step_factor << endl;
		analytical = false;
		PI = acos(-1.0);
		indata_name << "Inflation GW initial data";
		potential = new QuadraticPotential();
	};
	//================================================
	// Destructor
	//================================================
	~InflationGW() {};
	string Name() {
		return indata_name.str();
	};
	//================================================
	// Initializer
	//================================================
	bool Initialize(gf3d& fct) {
		if (!all_clear) return false;  // don't even try...
		cout << " INFLATIONGW: Initializing Inflation GW initial data... " << endl;
		n_r = fct.dim1();
		n_theta = fct.dim2();
		n_phi = fct.dim3();
		//
		// set up grid functions
		//
		int gf_counter = 2000;
		psi.setup(grid, 1, "psi", gf_counter++, +1, +1, +1);
		delta_psi.setup(grid, 1, "delta_psi", gf_counter++, +1, +1, +1);
		rho.setup(grid, 1, "rho", gf_counter++, +1, +1, +1);
		res.setup(grid, 1, "res", gf_counter++, +1, +1, +1);
		u.setup(grid, 1, "u", gf_counter++, +1, +1, +1);
		sf.setup(grid, 1, "sf", gf_counter++, +1, +1, +1);
		K.setup(grid, 2, "K", gf_counter++, +1, +1, +1);  
		s_r.setup(grid, 1, "s_r", gf_counter++, -1, +1, +1);
		s_t.setup(grid, 1, "s_t", gf_counter++, +1, -1, -1);
		s_p.setup(grid, 1, "s_p", gf_counter++, -1, -1, +1);
		// rescaled vector potential for A^{ij}, indices upstairs
		W_r.setup(grid, 1, "W_r", gf_counter++, -1, +1, +1);
		W_t.setup(grid, 1, "W_t", gf_counter++, +1, -1, -1);
		W_p.setup(grid, 1, "W_p", gf_counter++, -1, -1, +1);
		del_W_r.setup(grid, 1, "del_W_r", gf_counter++, -1, +1, +1);
		del_W_t.setup(grid, 1, "del_W_t", gf_counter++, +1, -1, -1);
		del_W_p.setup(grid, 1, "del_W_p", gf_counter++, -1, -1, +1);
		res_r.setup(grid, 1, "res_r", gf_counter++, -1, +1, +1);
		res_t.setup(grid, 1, "res_t", gf_counter++, +1, -1, -1);
		res_p.setup(grid, 1, "res_p", gf_counter++, -1, -1, +1);
		A_rr.setup(grid, 1, "A_rr", gf_counter++, +1, +1, +1);
		A_rt.setup(grid, 1, "A_rt", gf_counter++, -1, -1, -1);
		A_rp.setup(grid, 1, "A_rp", gf_counter++, +1, -1, +1);
		A_tt.setup(grid, 1, "A_tt", gf_counter++, +1, +1, +1);
		A_tp.setup(grid, 1, "A_tp", gf_counter++, -1, +1, -1);
		A_pp.setup(grid, 1, "A_pp", gf_counter++, +1, +1, +1);
		A2.setup(grid, 1, "A2", gf_counter++, +1, +1, +1);
		// AE_p.setup(grid, 1, "AE_p", gf_counter++, -1, -1, +1); 
		// AE_dot_p.setup(grid, 1, "AE_p", gf_counter++, -1, -1, +1); 
		//
		int i;
		psi_r.resize(n_r);
		for (i = 0; i < n_r; i++) {
			psi_r[i] = 0.0;
		}
		//================================================
		// compute source functions and initialize psi etc...
		//================================================
		//
		psi.equals(1.0);
		W_r.equals(0.0);
		W_t.equals(0.0);
		W_p.equals(0.0);

		//
#ifdef SR
		cout << "INFLATIONGW: Running without GR" << endl;
		Compute_Fields();
		dump(&psi);

		K.equals(0.0);
		s_r.equals(0.0);
		s_t.equals(0.0);
		s_p.equals(0.0);
		rho.equals(0.0);

		W_r.equals(0.0);
		W_t.equals(0.0);
		W_p.equals(0.0);

		res_r.equals(0.0);
		res_p.equals(0.0);
		res_t.equals(0.0);

		A_rr.equals(0.0);
		A_rt.equals(0.0);
		A_rp.equals(0.0);
		A_tt.equals(0.0);
		A_tp.equals(0.0);
		A_pp.equals(0.0);
		cout << " INFLATIONGW: done with initialization!" << endl;
		return true;
#else
	#ifdef NoEllSolver
			cout << " INFLATIONGW: Can't construct Inflation Gw initial data without an Elliptic Solver!! " << endl;
			return false;
	#else
			veclaplacian = new VecLaplacian(grid, true);
			veclaplacian->SetupSolver();
	#endif
			double res = Solve_Constraints();
			cout << " INFLATIONGW: Solved constraints momentum residual = " << Momentum_Residual() << endl;
			Compute_Aij();
			dump(&psi);
			if (res < tol) {
				cout << " INFLATIONGW: done with initialization! " << endl;
				return true;
			}
			else {
				cout << " INFLATIONGW: initialization did not converge... " << endl;
				return false;
			}
#endif
		};

	//================================================
	// Solve constraints
	//================================================
	double Solve_Constraints(bool verbose = false) {
		int step = 0;
		const double tol_tri = 1.e-12;
		cout << " INFLATIONGW: Computing scalar field and sources" << endl;
		Compute_Fields();
		Compute_Sources();
		cout << " INFLATIONGW: Solving for psi" << endl;
		Solve_Psi(tol_tri); 
		cout << " INFLATIONGW: Finding initial A_ij" << endl;
		Compute_Aij();
		cout << " INFLATIONGW: psi part of hamiltonian has residual: " << Hamiltonian_Psi_Residual() << endl;
		double res = Residual();
		cout << " INFLATIONGW: initial constraint residual = " << res << endl;
		while (abs(res) > tol && step < max_it) {
			step++;
			// call individual constraint solver with slightly smaller tolerances
			double current_tol = res / 1.e3;
			cout << " INFLATIONGW: current tolerance = " << current_tol << endl;
			Solve_K();
			Solve_Momentum(tol_tri, current_tol, verbose);
			res = Residual();
			cout << " INFLATIONGW: after " << step
				<< " steps constraint residual = " << res << endl;
			if (verbose)
				cout << " ============================================================"
				<< endl;
		}
		return Residual();
	}
	//================================================
	// Total residual
	//================================================
	double Residual() {
		const double mom_res = Momentum_Residual();
		cout << " INFLATIONGW: Momentum residual = " << mom_res << endl;
		Compute_Aij();
		const double ham_res = Hamiltonian_K_Residual();
		cout << " INFLATIONGW: Hamiltonian K residual = " << ham_res << endl;
		return sqrt(mom_res * mom_res + ham_res * ham_res);
	}

	//================================================
	// Analytical solution for h_{ij}
	//================================================
	double h_rr_analytical(double rl, double thetal, double phil, double t) {
		return 0.0;
	};
	double h_rt_analytical(double r, double theta, double phi, double t) {
		return 0.0;
	};
	double h_rp_analytical(double r, double theta, double phi, double t) {
		return 0.0;
	};
	double h_tt_analytical(double rl, double thetal, double phil, double t) {
		return 0.0;
	};
	double h_tp_analytical(double r, double theta, double phi, double t) {
		return 0.0;
	};
	double h_pp_analytical(double rl, double thetal, double phil, double t) {
		return 0.0;
	};
	double phi_analytical(double rl, double thetal, double phil, double tl) {
		// reality check:
		if (thetal == 0.0) cout << " ouch - found theta = 0 in phi_analytical " << endl;
		VecDoub* r = grid->r();
		VecDoub* theta = grid->theta();
		VecDoub* phi = grid->phi();
		if (rl != (*r)[grid->i_ind(rl)])
			cout << " i index not right - rl = " << rl << " i = " << grid->i_ind(rl) << " r(i) = "
			<< (*r)[grid->i_ind(rl)] << endl;
		if (!((thetal == (*theta)[grid->j_ind(thetal)]) || ((*theta)[grid->j_ind(thetal)] == PI - thetal)))
			//      cout << " j index not right - thetal = " << thetal << " PI - thetal = " << PI - thetal <<  " j = " << grid->j_ind(thetal) << " theta(j) = " << (*theta)[grid->j_ind(thetal)] << endl;
			if (phil != (*phi)[grid->k_ind(phil)])
				cout << " k index not right - phil = " << phil << " k = " << grid->k_ind(phil) << " phi(k) = "
				<< (*phi)[grid->k_ind(phil)] << endl;
		int i = grid->i_ind(rl);
		int j = grid->j_ind(thetal);
		int k = grid->k_ind(phil);
		double psil = psi(i, j, k);
		return log(psil);
	};
	//================================================
	// Analytical solution for connection coefficients
	//================================================
	double lam_r_analytical(double r, double theta, double phi, double t,
													bool& done) {
		done = true; return 0;
	}
	double lam_t_analytical(double r, double theta, double phi, double t,
													bool& done) {
		done = true; return 0;
	}
	double lam_p_analytical(double r, double theta, double phi, double t,
													bool& done) {
		done = true; return 0;
	}
	//================================================
	// Analytical solution for extrinsic curvature
	// NOTE: in initialize we compute \bar A_{ij} = psi^6 \tilde A_{ij},
	// now need BSSN rescaling of extrinsic curvature...
	//================================================
	double a_rr_analytical(double r, double theta, double phi, double t) {
		int i = grid->i_ind(r);
		int j = grid->j_ind(theta);
		int k = grid->k_ind(phi);
		double psi6 = pow(psi(i, j, k), 6);
		return A_rr(i, j, k) / psi6;
	};
	double a_rt_analytical(double r, double theta, double phi, double t) {
		int i = grid->i_ind(r);
		int j = grid->j_ind(theta);
		int k = grid->k_ind(phi);
		double psi6 = pow(psi(i, j, k), 6);
		return A_rt(i, j, k) / psi6;
	};
	double a_rp_analytical(double r, double theta, double phi, double t) {
		int i = grid->i_ind(r);
		int j = grid->j_ind(theta);
		int k = grid->k_ind(phi);
		double psi6 = pow(psi(i, j, k), 6);
		return A_rp(i, j, k) / psi6;
	};
	double a_tt_analytical(double r, double theta, double phi, double t) {
		int i = grid->i_ind(r);
		int j = grid->j_ind(theta);
		int k = grid->k_ind(phi);
		double psi6 = pow(psi(i, j, k), 6);
		return A_tt(i, j, k) / psi6;
	};
	double a_tp_analytical(double r, double theta, double phi, double t) {
		int i = grid->i_ind(r);
		int j = grid->j_ind(theta);
		int k = grid->k_ind(phi);
		double psi6 = pow(psi(i, j, k), 6);
		return A_tp(i, j, k) / psi6;
	};
	double a_pp_analytical(double r, double theta, double phi, double t) {
		int i = grid->i_ind(r);
		int j = grid->j_ind(theta);
		int k = grid->k_ind(phi);
		double psi6 = pow(psi(i, j, k), 6);
		return A_pp(i, j, k) / psi6;
	};
	double K_analytical(double r, double theta, double phi, double t) {
		int i = grid->i_ind(r);
		int j = grid->j_ind(theta);
		int k = grid->k_ind(phi);
		return K(i,j,k);
	};
	//================================================
	// Analytical solution for gauge
	//================================================
	double lapse_analytical(double rl, double thetal, double phil, double t) {
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
	double rho_0_analytical(double rl, double thetal, double phil, double tl) {
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
		int i = grid->i_ind(r);
		int j = grid->j_ind(theta);
		int k = grid->k_ind(phi);
		return sf(i,j,k);
	};
	double pi_analytical(double r, double theta, double phi, double t) {
		return 0.0;
	};
	//================================================
	// Analytical solution for Maxwell OR for Dual Maxwell:
	// For Maxwell need e_p and a_p only, for Dual Maxwell all a_i and as_i
	//================================================
	// returns correct function at t=0 only
	double e_p_analytical(double r, double theta, double phi, double t) {
		return 0.0;
	};
	//
	// provides flat-space analytical result at any time,
	// but only for r_0 = 0; see
	// eq. (16) in Knapp, Walker & Baumgarte, 2002
	//
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
private:
	//===============================================================
	// Compute fields (including rho_ADM and S^i)
	//===============================================================
#include "INFLATIONGW/Compute_Fields.h"
	//===============================================================
	// Solve Hamiltonian constraint
	//===============================================================
#include "INFLATIONGW/Solve_Hamiltonian.h"
	//===============================================================
	// Solve momentum constraints
	//===============================================================
#include "INFLATIONGW/Solve_Momentum.h"
	//=============================================================== 
	// Compute \bar A^{ij} from vector potential W^i
	//=============================================================== 
#include "INFLATIONGW/Solve_Aij.h"
	//================================================
	// Dump function
	//================================================
	void dump(gf3d* fct, int step = 0) {
		ofstream outfile_slice;
		ostringstream filename_slice;
		// slice files 
		filename_slice << "output/" << fct->Name() << "_indata_slice_"
			<< step << "_" << n_r - 2 * N_g
			<< "_" << n_theta - 2 * N_g << ends;
		outfile_slice.open(filename_slice.str().c_str());
		for (int i = N_g; i < n_r; i++) {
			for (int j = N_g; j < n_theta - N_g; j++)
				outfile_slice << setprecision(16) << setw(24) << fct->r(i)
				<< setprecision(16) << setw(24) << fct->theta(j)
				<< setprecision(16) << setw(24) << (*fct)(i, j, N_g) << endl;
			outfile_slice << endl;
		}
		outfile_slice.close();
		// ray files 
		ofstream outfile_rays;
		ostringstream filename_rays;
		filename_rays << "output/" << fct->Name() << "_indata_rays_"
			<< step << "_" << n_r - 2 * N_g
			<< "_" << n_theta - 2 * N_g << ends;
		outfile_rays.open(filename_rays.str().c_str());
		for (int i = N_g; i < n_r; i++) {
			outfile_rays << setprecision(16) << setw(24) << fct->r(i)
				<< setprecision(16) << setw(24) << (*fct)(i, 0.0, N_g)
				<< setprecision(16) << setw(24) << (*fct)(i, PI / 2., N_g) << endl;
		}
		outfile_rays.close();
	}
	//================================================
	// K Error function
	//================================================
	void dump_error_monitor(double K_variable) {
    	ofstream outfile("output/error_monitor", ios::app);

    	outfile << setprecision(16)
            	<< K_variable << endl;

    	outfile.close();
}
};

