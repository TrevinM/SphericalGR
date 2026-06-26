// Tell emacs that this is -*-c++-*- mode
//
//================================================
// 
// Initial for electromagnetic waves in electrovacuum, expressed in terms of two
// vector potentials A and *A with
//
//       B = curl A   and E = curl *A
//
// Will generate divergence-free A and A* by writing
//       A = A_tor   + curl V_A
//       *A = *A_tor + curl *V_A
// where only non-zero components of A_tor, A_tor, V_A and *V_A
// are the phi components, which
// depend on r and theta only (i.e. all purely toroidal).  Then:
//       B_pol = curl A_tor         B_tor = curl curl V_A
//       E_pol = curl *A_tor        E_tor = curl curl *V_A
//================================================
//
class DualEMWave : public InData {
private:
#ifndef NoEllSolver
	FlatEllSolver3D* laplace;
	VecLaplacian* veclaplacian;
#endif
	gf3d psi, rho, res, u, delta_psi;
	// all vector and tensor components rescaled...
	gf3d a_r, a_t, a_p;   // components of vector potential A (indices down)  
	gf3d as_r, as_t, as_p;  // components of vector potential *A (indices down)
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
	// gf3d v_a_p, v_as_p;     // phi component of tor. vectors that generate pol. parts of A and *A
	int n_r, n_theta, n_phi;
	int N_g;
	int n_psi;  // power of psi in vectors A and V_A, ie. A \propto psi^n etc
	int max_it;   // parameters for elliptic solver
	double tol;
	double a1_amp, a2_amp, a3_amp, as1_amp, as2_amp, as3_amp;   // parameters for initial data
	int n_a, n_v_a, n_as, n_v_as;       // powers of r in respective seed functions
	double r0;
	int flat;
	double PI;
	bool all_clear;
	ostringstream indata_name;
public:
	//================================================
	// Constructor
	//================================================
	DualEMWave(char* indata_input, Grid* grid_i, Cosmology* cosmology) :
		InData(grid_i, cosmology) {
		N_g = grid->N_ghosts();
		indata_type = em_wave;
		all_clear = true;
		ifstream infile;
		infile.open(indata_input);
		if (!infile) {
			cerr << " DUALEMWAVE Can't open " << indata_input
				<< " for input. This is bad. " << endl;
			all_clear = false;
		}
		else
			cout << " DUALEMWAVE: Reading initial data parameters from file "
			<< indata_input << endl;
		char buf[100], c;
		infile.get(buf, 100, '='); infile.get(c); infile >> a1_amp;
		infile.get(buf, 100, '='); infile.get(c); infile >> a2_amp;
		infile.get(buf, 100, '='); infile.get(c); infile >> a3_amp;
		infile.get(buf, 100, '='); infile.get(c); infile >> as1_amp;
		infile.get(buf, 100, '='); infile.get(c); infile >> as2_amp;
		infile.get(buf, 100, '='); infile.get(c); infile >> as3_amp;

		infile.get(buf, 100, '='); infile.get(c); infile >> n_psi;
		infile.get(buf, 100, '='); infile.get(c); infile >> r0;
		infile.get(buf, 100, '='); infile.get(c); infile >> flat;
		infile.get(buf, 100, '='); infile.get(c); infile >> max_it;
		infile.get(buf, 100, '='); infile.get(c); infile >> tol;
		infile.get(buf, 100, '='); infile.get(c); infile >> tau_star;
		infile.get(buf, 100, '='); infile.get(c); infile >> xi;

		cout << " DUALEMWAVE: Will set up E&M wave initial data with" << endl;
		cout << "      amplitude parameters a1_amp = " << a1_amp
			<< ", a2_amp = " << a2_amp
			<< ", a3_amp = " << a3_amp
			<< ", as1_amp = " << as1_amp
			<< ", as2_amp = " << as2_amp
			<< ", as3_amp = " << as3_amp << endl;

		cout << "      and exponents n_a = " << n_a << " n_as = " << n_as
			<< " n_v_a = " << n_v_a << " n_v_as = " << n_v_as << endl;
		cout << "      power of psi in seed functions: n_psi = " << n_psi << endl;
		cout << "      and off-center parameter r_0 = " << r0 << endl;
		cout << " DUALEMWAVE: Will run elliptic solver with max_it = " << max_it
			<< " and tol = " << tol << endl;
		analytical = false;
		PI = acos(-1.0);
		indata_name << "Dual E&M wave initial data";
	};
	//================================================
	// Destructor
	//================================================
	~DualEMWave() {};
	string Name() {
		return indata_name.str();
	};
	//================================================
	// Initializer
	//================================================
	bool Initialize(gf3d& fct) {
		if (!all_clear) return false;  // don't even try...
		cout << " DUALEMWAVE: Initializing E&M wave initial data... " << endl;
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
		//
		// Careful: weird parity across equator:
		//     A: tor (phi-component): behaves like vector
		//        pol (r,theta): behaves like pseudo-vector
		//    *A: tor (r-theta): behaves like pseudo-vector
		//        pol (phi): behaves like vector
		//
		// toroidal, poloidal here refers to direction of electric field in either polarization
		//
		a_r.setup(grid, 1, "a_r", gf_counter++, -1, +1, -1);   // pseudo-vector
		a_t.setup(grid, 1, "a_t", gf_counter++, +1, -1, +1);   // pseudo-vector
		a_p.setup(grid, 1, "a_p", gf_counter++, -1, -1, +1);   // vector
		as_r.setup(grid, 1, "as_r", gf_counter++, -1, +1, -1); // pseudo-vector
		as_t.setup(grid, 1, "as_t", gf_counter++, +1, -1, +1); // pseudo-vector
		as_p.setup(grid, 1, "as_p", gf_counter++, -1, -1, +1); // vector
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
		//================================================
		// compute source functions and initialize psi etc...
		//================================================
		//
		psi.equals(1.0);
		//
#ifdef NoEllSolver
		cout << " DUALEMWAVE: Can't construct EM initial data without an Elliptic Solver!! " << endl;
		return false;
#else
		veclaplacian = new VecLaplacian(grid, true);
		veclaplacian->SetupSolver();
#endif
		double res = Solve_Constraints();
		dump(&psi);
		if (res < tol) {
			cout << " DUALEMWAVE: done with initialization! " << endl;
			return true;
		}
		else {
			cout << " DUALEMWAVE: initialization did not converge... " << endl;
			return false;
		}
	}
	//================================================
	// Solve constraints
	//================================================
	double Solve_Constraints(bool verbose = false) {
		Compute_Fields();
		Compute_Sources();
		double res = Residual();
		cout << " DUALEMWAVE: initial constraint residual = " << res << endl;
		int step = 0;
		const double tol_tri = 1.e-12;
		while (res > tol && step < max_it) {
			step++;
			// call individual constraint solver with slightly smaller tolerances
			double current_tol = max(tol / 3., res / 1.e3);
			cout << " DUALEMWAVE: current tolerance = " << current_tol << endl;
			Compute_Fields();
			Compute_Sources();
			Solve_Momentum(tol_tri, current_tol, verbose);
			Compute_Fields();
			Compute_Sources();
			Compute_Aij();
			Solve_Hamiltonian(tol_tri, current_tol, verbose);
			res = Residual();
			cout << " DUALEMWAVE: after " << step
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
		const double ham_res = Hamiltonian_Residual();
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
		return 0.0;
	};
	//================================================
	// Analytical solution for gauge
	//================================================
	double lapse_analytical(double rl, double thetal, double phil, double t) {
		int i = grid->i_ind(rl);
		int j = grid->j_ind(thetal);
		int k = grid->k_ind(phil);
		double psil = psi(i, j, k);
		// return 0.2 + 0.8*rl*rl/(1.0 + rl*rl);
		// return 1.0;
		return 1.0 / (psil * psil);
	};
	double shift_r_analytical(double r, double theta, double phi, double t) {
		// return 2.0 * r / (1.0 + r*r);
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
		return 0.0;
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
		return compute_a_r(r, theta);
	};
	double a_t_analytical(double r, double theta, double phi, double t) {
		return compute_a_t(r, theta);
	};
	double a_p_analytical(double r, double theta, double phi, double t) {
		return compute_a_p(r, theta);
	};
	double as_r_analytical(double r, double theta, double phi, double t) {
		return compute_as_r(r, theta);
	};
	double as_t_analytical(double r, double theta, double phi, double t) {
		return compute_as_t(r, theta);
	};
	double as_p_analytical(double r, double theta, double phi, double t) {
		return compute_as_p(r, theta);
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
#include "DUALMAXWELL/Compute_Fields.h"
	//===============================================================
	// Solve Hamiltonian constraint
	//===============================================================
#include "DUALMAXWELL/Solve_Hamiltonian.h"
	//===============================================================
	// Solve momentum constraints
	//===============================================================
#include "DUALMAXWELL/Solve_Momentum.h"
	//=============================================================== 
	// Compute \bar A^{ij} from vector potential W^i
	//=============================================================== 
	void Compute_Aij() {
		for (int i = N_g; i < n_r - N_g; i++) {
			const double rl = rho.r(i);
			const double r2 = rl * rl;
			for (int j = N_g; j < n_theta - N_g; j++) {
				const double stl = rho.sintheta(j);
				const double st2 = stl * stl;
				const double csl = rho.costheta(j);
				const double ctl = csl / stl;
				for (int k = N_g; k < n_phi - N_g; k++) {
					const double div = W_r.dr(i, j, k) + 2. * W_r(i, j, k) / rl
						+ W_t.dtheta(i, j, k) / rl + W_t(i, j, k) * ctl / rl;
					// compute covariant derivatives D^i W^j
					const double DrWr = W_r.dr(i, j, k);
					const double DtWr = (W_r.dtheta(i, j, k) - W_t(i, j, k)) / rl;
					const double DpWr = -W_p(i, j, k) / rl;
					const double DrWt = W_t.dr(i, j, k);
					const double DtWt = (W_t.dtheta(i, j, k) + W_r(i, j, k)) / rl;
					const double DpWt = -ctl * W_p(i, j, k) / rl;
					const double DrWp = W_p.dr(i, j, k);
					const double DtWp = W_p.dtheta(i, j, k) / rl;
					const double DpWp = (W_r(i, j, k) + ctl * W_t(i, j, k)) / rl;
					A_rr[i][j][k] = 2.0 * DrWr - 2. / 3. * div;
					A_rt[i][j][k] = DrWt + DtWr;
					A_rp[i][j][k] = DrWp + DpWr;
					A_tt[i][j][k] = 2.0 * DtWt - 2. / 3. * div;
					A_tp[i][j][k] = DtWp + DpWt;
					A_pp[i][j][k] = 2.0 * DpWp - 2. / 3. * div;
					// if (i == N_g && j == N_g && k == N_g) 
					//   cout << " reality check: "
					//        << A_rr(i,j,k) + A_tt(i,j,k) + A_pp(i,j,k) << endl;
					A2[i][j][k] = A_rr(i, j, k) * A_rr(i, j, k)
						+ 2.0 * A_rt(i, j, k) * A_rt(i, j, k)
						+ 2.0 * A_rp(i, j, k) * A_rp(i, j, k)
						+ A_tt(i, j, k) * A_tt(i, j, k)
						+ 2.0 * A_tp(i, j, k) * A_tp(i, j, k)
						+ A_pp(i, j, k) * A_pp(i, j, k);
				}
			}
		}

	}
	//===============================================================
	// Compute physical curl: returns scaled eps^ijk ( D_j A_k - D_k A_j )
	// where D_i is covariant derivative with respect to reference metric,
	// but epsilon is physical (rescaled) 3D epsilon: psi^{-6}[ijk]
	//===============================================================
	void curl(gf3d& a_r, gf3d& a_t, gf3d& a_p,
						double& curl_a_r, double& curl_a_t, double& curl_a_p,
						int i, int j, int k) {
		const double rl = a_r.r(i);
		const double stl = a_r.sintheta(j);
		const double ctl = a_r.costheta(j);
		const double cot = ctl / stl;
		const double psil = psi(i, j, k);
		const double psim6 = 1. / (psil * psil * psil * psil * psil * psil);
		// const double D_r_A_t = rl*a_t.dr(i,j,k);
		// const double D_r_A_p = rl*stl*a_p.dr(i,j,k);
		// const double D_t_A_r = a_r.dtheta(i,j,k) - a_t(i,j,k);
		// const double D_t_A_p = rl*stl*a_p.dtheta(i,j,k);
		// const double D_p_A_r = - stl*a_p(i,j,k);
		// const double D_p_A_t = - rl*ctl*a_p(i,j,k);
		curl_a_r = psim6 * (a_p.dtheta(i, j, k) + cot * a_p(i, j, k)) / rl;
		curl_a_t = -psim6 * (a_p.dr(i, j, k) + a_p(i, j, k) / rl);
		curl_a_p = psim6 * (a_t.dr(i, j, k) - (a_r.dtheta(i, j, k) - a_t(i, j, k)) / rl);
	}

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

};

