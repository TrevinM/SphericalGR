//================================================
//
// NOTE: In this version all hydro variables are properly rescaled.
//
// REFERENCE: Thierfelder, Bernuzzi & Bruegmann, arXiv:1104.4751 (TBB)
//      Duez, Liu, Shapiro & Stephens, arXiv:astro-ph/0503420 (DLSS)
//      Farris, Li, Liu & Shapiro, PRD 78, 024023 (2008) (FLLS)
//================================================
//
// NOTE: use state vectors for conservative variables and fluxes 
// as defined in (21) and (24) of TBB - i.e. *with* factors of determinant of 
// metric and lapse
//
// i.e.: gridfunction D, e.g., is D-component of q_D = gamma * D
//
//================================================
//
// File that includes routines for class RadHydro
//
//================================================
#include "Grid.h"
#include "Matter.h"
#include "tensors.h"
//================================================
//
// Read input from file "RadHydro_Input"
//
//================================================
int RadHydro::ReadInput(int & N_particles, double & R_max) {
  int error = 0;
  ifstream infile;
  infile.open("RadHydro_Input");
  if (!infile) {
    cerr << " RADHYDRO: Can't open input file RadHydro_Input " << endl;
    cerr << " RADHYDRO: Will use default values " << endl;
    error = 1;
    return error;
  }
  char buf[500],c;
  infile.get(buf,500,'='); infile.get(c); infile >> slope_limiter_type;
  infile.get(buf,500,'='); infile.get(c); infile >> riemann_type;
  infile.get(buf,500,'='); infile.get(c); infile >> partial;
  infile.get(buf,500,'='); infile.get(c); infile >> f_atm;
  infile.get(buf,500,'='); infile.get(c); infile >> f_thr;
  infile.get(buf,500,'='); infile.get(c); infile >> N_particles;
  infile.get(buf,500,'='); infile.get(c); infile >> R_max;
  if(infile.eof()) {
    cerr << " RADHYDRO: Error reading input file RADHYDRO_Input " << endl;
    cerr << " RADHYDRO: Will use default values " << endl;
    error = 2;
  }
  return error;
}
//================================================
//
// Initialize fluid
//
//================================================
void RadHydro::Initialize(state *s, curvature *c, diagnostics *d) {
  indata->Initialize_RadHydro(aux->rho_0,aux->p,aux->v_r,aux->v_t,aux->v_p,
			      aux->E, aux->f, aux->f_r, aux->f_t, aux->f_p);
  cout << " eps = " << aux->eps(3,3,3) << " p = "  << aux->p(3,3,3)
       << " D = " << last->D(3,3,3) << " tau = " << last->tau(3,3,3)
       << " E = " << aux->E(3,3,3) << " rho_0 = " << aux->rho_0(3,3,3) << endl;
  rho_atm = f_atm * aux->rho_0.max();
  rho_thr = f_thr * rho_atm;
  for (int i = N_g; i < N_r; i++)    
    for (int j = N_g; j < N_t - N_g; j++)
      for (int k = N_g; k < N_p - N_g; k++) {
	if (aux->rho_0(i,j,k) < rho_thr) {
	  aux->rho_0[i][j][k] = rho_atm;
	  aux->eps[i][j][k] = eos->cold_eps(aux->rho_0(i,j,k));
	}
	aux->eps[i][j][k] = (3.0 / 2.0 ) * aux->p(i,j,k) / aux->rho_0(i,j,k);
      }	
  //
  // Now compute conserved quantities on initial slice
  //
  Compute_Conserved_Vars(last,s,c);
  //
  // "Recover" primitive variables - deals with atmosphere for t=0
  //
  cout << " eps = " << aux->eps(3,3,3) << " p = "  << aux->p(3,3,3) 
       << " D = " << last->D(3,3,3) << " tau = " << last->tau(3,3,3) 
       << " E = " << aux->E(3,3,3) << " rho_0 = " << aux->rho_0(3,3,3) << endl;
  Recovery(last,s,c);
  cout << " eps = " << aux->eps(3,3,3) << " p = "  << aux->p(3,3,3) 
       << " D = " << last->D(3,3,3) << " tau = " << last->tau(3,3,3) 
       << " E = " << aux->E(3,3,3) << " rho_0 = " << aux->rho_0(3,3,3) << endl; 
  //
  //   Test_Recovery(s, c);
  //
  // Finally: add field variables to particle dump functions:
  //
  if (particles != NULL) {
    particles->add_to_dump_list("Particle_Fct_List", s);
    particles->add_to_dump_list("Particle_Fct_List", c);
    particles->add_to_dump_list("Particle_Fct_List", d);
  }
};

//================================================
//
// Compute right-hand sides...
//
//================================================
void RadHydro::Compute_RHS(state *s, curvature *c, double time) {
  //
  // inverse Cowling approximation: radhydro-without-radhydro...
  //
  if (cowling == 2)
    derivs->equals(0.0);   
  else {
    //
    // Otherwise: compute primitives from matter state inter
    //
    Recovery(inter, s, c);
    //
    // Construct fluxes; involves
    // - reconstruction at cell interfaces (in Reconstruct_i)
    // - computation of fluxes at interfaces (in Compute_Fluxes)
    // - solution of local Riemann problem (in riemann-solver->flux)
    //
    Construct_Fluxes(s, c);
    //
    // then compute time derivatives of conserved quantities
    //
    dot_conserved_vars(inter, s, c);
    //
    // Outer boundaries (CHECK: may have to fix...)
    //
    for (int i = 0; i < derivs->N_fcts; i++)
      derivs->fct_list[i]->derivs_outerboundary(inter->fct_list[i]);
    //
    // finally: compute right-hand sides for particles
    // 
    if (particles != NULL) 
      particles->Compute_RHS(&s->lapse, &aux->W,
			     &s->shift_r, &s->shift_t, &s->shift_p, 
			     &aux->v_r, &aux->v_t, &aux->v_p); 
  }
};

//================================================
// Compute conserved variables (multiplied with determinant of metric, 
//    i.e. members of q vector!) from primitive ones
//
// Either on entire grid, or locally.  
// First global version which calls local version
//
//================================================
int RadHydro::Compute_Conserved_Vars(radhydro_state *m, state *s, curvature *c) {
  int error = 0;
  for (int i = N_g; i < N_r; i++ ) {
    double r_l = grid->r(i);
    for (int j = N_g; j < N_t-N_g; j++ ) {
      double sintheta_l = grid->sintheta(j);
      for (int k = N_g; k < N_p-N_g; k++ ) {
  	double h_rr_l = s->h_rr(i,j,k);
  	double h_rt_l = s->h_rt(i,j,k);
  	double h_rp_l = s->h_rp(i,j,k);
  	double h_tt_l = s->h_tt(i,j,k);
  	double h_tp_l = s->h_tp(i,j,k);
  	double h_pp_l = s->h_pp(i,j,k);
  	double phi_l  = s->phi(i,j,k);
  	double lapse_l = s->lapse(i,j,k);
  	double det_l = c->det(i,j,k);
	double D_l, S_r_l, S_t_l, S_p_l, tau_l, W_l, v2_l;
	double tau_rad_l, S_rad_r_l, S_rad_t_l, S_rad_p_l;
	const double rho_l = aux->rho_0(i,j,k);
	const double eps_l = aux->eps(i,j,k);
	double p_l = eos->P(rho_l,eps_l);
	error += Compute_Conserved_Vars(rho_l, eps_l, p_l,
					aux->v_r(i,j,k), aux->v_t(i,j,k), aux->v_p(i,j,k),
					aux->E(i,j,k), aux->f(i,j,k), 
					aux->f_r(i,j,k), aux->f_t(i,j,k), aux->f_p(i,j,k),
					h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, 
					lapse_l, phi_l, det_l, r_l, sintheta_l, 
					D_l, S_r_l, S_t_l, S_p_l, tau_l, W_l, v2_l,
					tau_rad_l, S_rad_r_l, S_rad_t_l, S_rad_p_l);
	m->D[i][j][k] = D_l;
	m->S_r[i][j][k] = S_r_l;
	m->S_t[i][j][k] = S_t_l;
	m->S_p[i][j][k] = S_p_l;
	m->tau[i][j][k] = tau_l;
	//
	m->tau_rad[i][j][k] = tau_rad_l;
	m->S_rad_r[i][j][k] = S_rad_r_l;
	m->S_rad_t[i][j][k] = S_rad_t_l;
	m->S_rad_p[i][j][k] = S_rad_p_l;      }
    }
  }
  m->D.fill_ghosts();
  m->S_r.fill_ghosts();
  m->S_t.fill_ghosts();
  m->S_p.fill_ghosts();
  m->tau.fill_ghosts();
  m->S_rad_r.fill_ghosts();
  m->S_rad_t.fill_ghosts();
  m->S_rad_p.fill_ghosts();
  m->tau_rad.fill_ghosts();
  return error;
};
//
//================================================
// Now local version
//================================================
//
int RadHydro::Compute_Conserved_Vars(double rho_0, double eps, double p,
				     double v_r, double v_t, double v_p,
				     double E, double f, double f_r, double f_t, double f_p,
				     double h_rr, double h_rt, double h_rp,
				     double h_tt, double h_tp, double h_pp,
				     double lapse, double phi, double det,
				     double r, double sintheta,
				     double & D, 
				     double & S_r, double & S_t, double & S_p, 
				     double & tau, double & W, double & v2,
				     double & tau_rad, 
				     double & S_rad_r, double & S_rad_t, double & S_rad_p,
				     bool verbose) {
  int error = 0;
  // CHECK: assuming that v^i *is* rescaled!
  const double det_3D = exp(6.0*phi)*sqrt(det);
  const double e4p = exp(4.0 * phi);
  tensor g_conf(1.0 + h_rr, h_rt, h_rp, 1.0 + h_tt, h_tp, 1.0 + h_pp);
  //
  vect v_up(v_r, v_t, v_p);
  vect f_up(f_r, f_t, f_p);
  vect v_low;
  vect f_low;
  double vdotv = 0.0;
  // lower index on v and compute dot product
  for (int a = 0; a < 3; a++) {
    v_low[a] = 0.0;
    f_low[a] = 0.0;
    for (int b = 0; b < 3; b++) {
      v_low[a] += e4p * g_conf[a][b] * v_up[b];
      f_low[a] += e4p * g_conf[a][b] * f_up[b];
      vdotv += e4p * g_conf[a][b] * v_up[a] * v_up[b];
    }
  }
  v2 = vdotv;
  if (v2 >= 1.0) {
    cout << " Trouble in Compute_Conserved_Vars: v2 = " << v2 << " r = " << r << endl;
    error = 1;
  }
  //
  // compute Lorentz factor W
  //
  W = 1.0 / sqrt( 1.0 - vdotv );   
  const double ut = W / lapse;
  //
  // compute matter density D
  //
  if (partial)
    D = det_3D * r * r * sintheta * W * rho_0;
  else
    D = det_3D * W * rho_0;
  //
  // compute momentum density S_i (indices down-stairs!) (again: *rescaled!*)
  //
  const double W2 = W*W;
  const double hl = eos->h(rho_0,eps);
  S_r = det_3D * W2 * rho_0 * hl * v_low[0];
  S_t = det_3D * W2 * rho_0 * hl * v_low[1];
  S_p = det_3D * W2 * rho_0 * hl * v_low[2];
  //
  // compute energy density tau
  //
  //  tau = W2 * rho_0 * hl - eos->P(rho_0,eps) - D;
  if (partial) 
    tau = det_3D * r * r * sintheta * ( W2 * rho_0 * hl - p ) - D;
  else
    tau = det_3D * ( W2 * rho_0 * hl - p ) - D;
  //
  if (verbose) 
    cout << " in Compute_Conserved_Vars: sintheta = " << sintheta 
	 << setprecision(12) << " tau = " << tau
	 << " v_r = " << v_r
	 << " v_t = " << v_t
	 << " v_p = " << v_p
	 << " W2 = " << W 
	 << " det_3D = " << det_3D << endl;
  //
  // now radiation quantities (Note: u_i = W v_i)
  //
  if (partial)
    tau_rad = det_3D * r * r * sintheta * (4.0/3.0*W2*E + 2.0*lapse*W*f - E/3.0 );
  else
    tau_rad = det_3D * (4.0/3.0*W2*E + 2.0*W*f - E/3.0 );
  S_rad_r =   det_3D * (4.0/3.0*W2*E*v_low[0] + f*W*v_low[0] + W*f_low[0] );
  S_rad_t =   det_3D * (4.0/3.0*W2*E*v_low[1] + f*W*v_low[1] + W*f_low[1] );
  S_rad_p =   det_3D * (4.0/3.0*W2*E*v_low[2] + f*W*v_low[2] + W*f_low[2] );
  return error;
};
//================================================
// Compute primitive variables from conserved ones:
//
// pass in current time-level of conserved variables, computes
// primitive variables
//================================================
int RadHydro::Recovery(radhydro_state *m, state *s, curvature *c) {
  int error = 0;
  //
  // set up atmosphere treatment...
  // 
  const double rho_atm = f_atm * aux->rho_0.max();
  //  const double rho_atm = f_atm * 0.4;
  rho_thr = f_thr * rho_atm;
  const double eps_atm = eos->cold_eps(rho_atm);
  //  cout << " eps_atm = " << eps_atm << endl;
  const double p_atm = eos->P(rho_atm,eps_atm); 
  //
  // find radius of outermost "particle", if they are in use
  // 
  double r_surface = grid->r(N_r-N_g);
  if (particles != NULL) r_surface = particles->RadiusOuterParticle();
  //
  // now go to each grid-point
  //
  int i_test = -65;
  int j_test = 3;
  int k_test = 3;
  // NOTE TOLERANCE! - consider adjusting...
  const double tol = 1.e-8;
  for (int i = N_g; i < N_r; i++ ) {
    double r_l = grid->r(i);
    for (int j = N_g; j < N_t - N_g; j++ ) {
      double sintheta_l = grid->sintheta(j);
      for (int k = N_g; k < N_p - N_g; k++ ) {
	const double det_3D = exp(6.0*s->phi(i,j,k))*sqrt( c->det(i,j,k) );
	//
	// define local conserved quantities
	// RECALL: gridfunctions are components for "q" vector; see eq. (50) in
	//         Montero, Baumgarte & Mueller
	//
	double Dl = m->D[i][j][k] / det_3D;
	double taul = m->tau[i][j][k] / det_3D;
	// CHECK: do we need to divide by det_3D??
	double tau_radl = m->tau_rad[i][j][k];
	if (partial) {
	  Dl /= r_l * r_l * sintheta_l;
	  taul /= r_l * r_l * sintheta_l;
	  tau_radl /= r_l * r_l * sintheta_l;
	}
	//
	// store S_i and metric in tensors to compute S2
	//
	vect S(m->S_r[i][j][k]/det_3D, m->S_t[i][j][k]/det_3D, m->S_p[i][j][k]/det_3D);
	vect S_rad(m->S_rad_r(i,j,k), m->S_rad_t(i,j,k), m->S_rad_p(i,j,k));	
	vect S_up;
	vect S_rad_up;
	//
	const double em4p = exp(- 4.0 * s->phi(i,j,k));
	tensor g_conf(1.0 + s->h_rr(i,j,k), s->h_rt(i,j,k), s->h_rp(i,j,k), 
		      1.0 + s->h_tt(i,j,k), s->h_tp(i,j,k), 1.0 + s->h_pp(i,j,k));
	tensor gup = g_conf.inverse();
	//
	double S2 = 0.0;
	for (int a = 0; a < 3; a++) {
	  S_up[a] = 0.0;
	  S_rad_up[a] = 0.0;
	  for (int b = 0; b < 3; b++) {
	    S_up[a] += em4p * gup[a][b] * S[b];
	    S_rad_up[a] += em4p * gup[a][b] * S_rad[b];
	    S2 += em4p * gup[a][b] * S[a] * S[b];
	  }
	}
	//
	// Now check condition (34) in Etienne et.al., PRD 77, 084002 (2008) - 
	// CHECK: put in a "safety" factor of 0.99 or so??
	//
	const double S2_max = taul*(taul + 2.0*Dl);
	if (S2 > 0.98 * S2_max) {
	  const double factor = sqrt(0.98 * S2_max / S2);
	  m->S_r[i][j][k] *= factor;
	  m->S_t[i][j][k] *= factor;
	  m->S_p[i][j][k] *= factor;
	  for (int a = 0; a < 3; a++) {
	    S[a] *= factor;
	    S_up[a] *= factor;
	  }
	  S2 *= factor*factor;
	  //  rho_0_guess = 0.0;  // will send this into atmosphere...
	}
	double p_guess = aux->p[i][j][k];  // take old pressure value as starting value for iteration
	//
	// compute "guess" value for rho_0 to decide whether it's atmosphere or not
	//
	double rho_0_guess = 0.0;
	double W_guess = 1.0;
	const double tauDp_guess = taul + p_guess + Dl;
	if (tauDp_guess > 0.0) {
	  W_guess = tauDp_guess/sqrt( fabs( tauDp_guess*tauDp_guess - S2 ) );	
	  rho_0_guess = Dl / W_guess;
	} 
	//
	// Now, if rho_0_guess too small...
	//
	if (rho_0_guess <= rho_thr || taul < 0.0 || Dl < 0.0 ) {
	// CHECK: follow outermost particle for OS collapse:
	//    || r_l > r_surface ) {
	  // cout << " Found atmosphere! tau = " << taul << " det = " << det_3D << endl;
	  //===========================================================
	  // ATMOSPHERE
	  //===========================================================
	  aux->rho_0[i][j][k] = rho_atm;
	  aux->v_r[i][j][k] = 0.0;
	  aux->v_t[i][j][k] = 0.0;
	  aux->v_p[i][j][k] = 0.0;
	  aux->eps[i][j][k] = eps_atm;
	  aux->p[i][j][k] = p_atm;
	  // CHECK atmosphere...
	  aux->E[i][j][k] = 0.0;
	  aux->f[i][j][k] = 0.0;
	  aux->f_r[i][j][k] = 0.0;
	  aux->f_t[i][j][k] = 0.0;
	  aux->f_p[i][j][k] = 0.0;
	  //
	  // update conserved variables accordingly
	  //
	  double h_rr_l = s->h_rr(i,j,k);
	  double h_rt_l = s->h_rt(i,j,k);
	  double h_rp_l = s->h_rp(i,j,k);
	  double h_tt_l = s->h_tt(i,j,k);
	  double h_tp_l = s->h_tp(i,j,k);
	  double h_pp_l = s->h_pp(i,j,k);
	  double phi_l  = s->phi(i,j,k);
	  double lapse_l = s->lapse(i,j,k);
	  double det_l = c->det(i,j,k);
	  //
	  double D_l, S_r_l, S_t_l, S_p_l, tau_l, W_l, v2_l;
	  double tau_rad_l, S_rad_r_l, S_rad_t_l, S_rad_p_l;
	  Compute_Conserved_Vars(rho_atm, eps_atm, p_atm, 0.0, 0.0, 0.0,
				 aux->E(i,j,k), aux->f(i,j,k), 
				 aux->f_r(i,j,k), aux->f_t(i,j,k), aux->f_p(i,j,k),
				 h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, 
				 lapse_l, phi_l, det_l,
				 r_l, sintheta_l, 
				 D_l, S_r_l, S_t_l, S_p_l, tau_l, W_l, v2_l,
				 tau_rad_l, S_rad_r_l, S_rad_t_l, S_rad_p_l);
	  m->D[i][j][k] = D_l;
	  m->S_r[i][j][k] = S_r_l;
	  m->S_t[i][j][k] = S_t_l;
	  m->S_p[i][j][k] = S_p_l;
	  m->tau[i][j][k] = tau_l;
	  aux->W[i][j][k] = W_l;
	  //
	  m->tau_rad[i][j][k] = tau_rad_l;
	  m->S_rad_r[i][j][k] = S_rad_r_l;
	  m->S_rad_t[i][j][k] = S_rad_t_l;
	  m->S_rad_p[i][j][k] = S_rad_p_l;	  
	} else {	  
	  //===========================================================
	  // Not in atmosphere - now use Newton-Raphson to iterate...
	  //===========================================================
	  double dfdp = 0.0;
	  double f_guess = f(p_guess,taul,Dl,S2,dfdp);
	  int it = 0;
	  if (i == i_test && j == j_test && k == k_test)
	    cout << " Recovery: p_guess = " << setprecision(16) << p_guess << " after " 
		 << it << " iterations, f_guess = " << f_guess << endl;
	  int it_max = 50;
	  // NOTE: tolerance for relative error in p - check same criterion below
	  while (it < it_max && abs(f_guess) > tol*p_guess) {
	    it++;
	    // 
	    // improved value of p:
	    //
	    p_guess -= f_guess / dfdp;
	    f_guess = f(p_guess,taul,Dl,S2,dfdp);
	    if (i == i_test && j == j_test && k == k_test)
	      cout << " Recovery: p_guess = " << p_guess << " after " 
	    	   << it << " iterations, f_guess = " << f_guess 
	    	   << " tau*(tau + 2 D) - S2 = " << taul*(taul + 2.0*Dl) - S2 << endl;
	  }
	  // NOTE: again...
	  if (abs(f_guess) < tol*p_guess) {
	    //
	    // iteration converged... compute primitive variables
	    //
	    const double taupD = taul + p_guess + Dl;
	    const double taupD2 = taupD * taupD;
	    const double Wl = taupD / sqrt( taupD2 - S2 );
	    // if (Wl != 1.0) cout << " at r = " << r_l << " Wl = " << setprecision(16)<< Wl 
	    // 			<< " taul = " << taul << " p_guess = " << p_guess 
	    // 			<< " Dl = " << Dl << " i = " << i << endl;
	    aux->W[i][j][k] = Wl;
	    aux->rho_0[i][j][k] = Dl/Wl;
	    aux->eps[i][j][k] = ( sqrt(taupD2 - S2) - Wl * p_guess - Dl) / Dl;
	    // if (aux->eps(i,j,k) < 0.0) 
	    if (i == i_test && j == j_test && k == k_test)
	      cout << " In Recovery : taupD = " << taupD << " p_guess = " 
		   << p_guess << " P " << eos->P(aux->rho_0(i,j,k),aux->eps(i,j,k)) << " Dl = " << Dl 
		   << " S2 = " << S2 << " i " << i << " rho_0 "
		   << aux->rho_0(i,j,k) << " eps " << aux->eps(i,j,k) << endl;
	    aux->v_r[i][j][k] = S_up[0] / taupD;
	    aux->v_t[i][j][k] = S_up[1] / taupD;
	    aux->v_p[i][j][k] = S_up[2] / taupD;
	    aux->p[i][j][k] = eos->P(aux->rho_0(i,j,k),aux->eps(i,j,k));
	    //
	    // now compute radiation variables...
	    // (see Eqs. (65) and (66) in FLLS)
	    // 
	    const double lapse_l = s->lapse(i,j,k);
	    const double shift_rl = s->shift_r(i,j,k);
	    const double shift_tl = s->shift_t(i,j,k);
	    const double shift_pl = s->shift_p(i,j,k);
	    //	    const double tau_rad_l = tau_rad(i,j,k);
	    const double v_dot_S = aux->v_r(i,j,k) * m->S_rad_r(i,j,k)
	      + aux->v_t(i,j,k) * m->S_rad_t(i,j,k) 
	      + aux->v_p(i,j,k) * m->S_rad_p(i,j,k);
	    const double s_0 = tau_radl;
	    const double s_1 = Wl * ( tau_radl - v_dot_S );
	    const double A_00 = det_3D * (4./3.*Wl*Wl - 1./3.);
	    const double A_01 = 2.*det_3D*Wl;
	    const double A_10 = det_3D*Wl;
	    const double A_11 = det_3D;
	    aux->f[i][j][k] = (A_10*s_0 - A_00*s_1)/(A_01*A_10 - A_11*A_00);
	    aux->E[i][j][k] = s_0/A_00 - A_01 * aux->f(i,j,k) / A_00;
	    if (!isfinite(aux->E(i,j,k))) 
	      cout << i << "   " << j << " v_dot_S = " << v_dot_S 
		   << " S_rad_r = " << m->S_rad_r(i,j,k) << " S_r = " << m->S_r(i,j,k) 
		   << " A_00 = " << A_00 << endl; 

	    // see Eq. (67) in FLLS, 
	    //   except note different definition of v^i used here
	    aux->f_r[i][j][k] = S_rad_up[0] / (Wl * det_3D) 
	      - 4./3.*aux->E(i,j,k)*Wl * aux->v_r(i,j,k)
	      - aux->f(i,j,k) * aux->v_r(i,j,k);
	    aux->f_t[i][j][k] = S_rad_up[1] / (Wl * det_3D) 
	      - 4./3.*aux->E(i,j,k)*Wl * aux->v_t(i,j,k)
	      - aux->f(i,j,k) * aux->v_t(i,j,k);
	    aux->f_p[i][j][k] = S_rad_up[2] / (Wl * det_3D) 
	      - 4./3.*aux->E(i,j,k)*Wl * aux->v_p(i,j,k)
	      - aux->f(i,j,k) * aux->v_p(i,j,k);
	    //
	    // in atmosphere after all...
	    //
	    if (aux->rho_0[i][j][k] <= rho_thr || aux->eps[i][j][k] <= eps_atm) {
	      //===========================================================
	      // ATMOSPHERE
	      //===========================================================
	      aux->rho_0[i][j][k] = rho_atm;
	      aux->v_r[i][j][k] = 0.0;
	      aux->v_t[i][j][k] = 0.0;
	      aux->v_p[i][j][k] = 0.0;
	      aux->eps[i][j][k] = eps_atm;
	      aux->p[i][j][k] = p_atm;
	      // CHECK atmosphere...
	      aux->E[i][j][k] = 0.0;
	      aux->f[i][j][k] = 0.0;
	      aux->f_r[i][j][k] = 0.0;
	      aux->f_t[i][j][k] = 0.0;
	      aux->f_p[i][j][k] = 0.0;
	      //
	      // update conserved variables accordingly
	      //
	      double h_rr_l = s->h_rr(i,j,k);
	      double h_rt_l = s->h_rt(i,j,k);
	      double h_rp_l = s->h_rp(i,j,k);
	      double h_tt_l = s->h_tt(i,j,k);
	      double h_tp_l = s->h_tp(i,j,k);
	      double h_pp_l = s->h_pp(i,j,k);
	      double phi_l  = s->phi(i,j,k);
	      double lapse_l = s->lapse(i,j,k);
	      double det_l = c->det(i,j,k);
	      //
	      double D_l, S_r_l, S_t_l, S_p_l, tau_l, W_l, v2_l;
	      double tau_rad_l, S_rad_r_l, S_rad_t_l, S_rad_p_l;
	      Compute_Conserved_Vars(rho_atm, eps_atm, p_atm, 0.0, 0.0, 0.0,
				     aux->E(i,j,k), aux->f(i,j,k), 
				     aux->f_r(i,j,k), aux->f_t(i,j,k), aux->f_p(i,j,k),
				     h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, 
				     lapse_l, phi_l, det_l, r_l, sintheta_l, 
				     D_l, S_r_l, S_t_l, S_p_l, tau_l, W_l, v2_l,
				     tau_rad_l, S_rad_r_l, S_rad_t_l, S_rad_p_l);
	      m->D[i][j][k] = D_l;
	      m->S_r[i][j][k] = S_r_l;
	      m->S_t[i][j][k] = S_t_l;
	      m->S_p[i][j][k] = S_p_l;
	      m->tau[i][j][k] = tau_l;
	      aux->W[i][j][k] = W_l;
	      //
	      m->tau_rad[i][j][k] = tau_rad_l;
	      m->S_rad_r[i][j][k] = S_rad_r_l;
	      m->S_rad_t[i][j][k] = S_rad_t_l;
	      m->S_rad_p[i][j][k] = S_rad_p_l;	  
	    }
	  } else { 
	    //
	    // iteration did not converge
	    //
	    // const double tauDp = taul + p_guess + Dl;
	    // const double W_guess = tauDp/sqrt(tauDp * tauDp - S2);
	    cout << " In Recovery: iteration did not converge at r = " << r_l << " i,j,k = "  
	    	 << i << ", " << j << ", " << k << " S2 = " << S2 << " tau(tau+2D) " 
		 << taul*(taul+ 2.0*Dl) << " rho_0_guess = " << rho_0_guess << " rho_thr = " << rho_thr << endl;
	    //
	    // take some emergency measures...
	    //
	    aux->rho_0[i][j][k] = rho_0_guess;
	    aux->eps[i][j][k] = eos->cold_eps(rho_0_guess);
	    aux->p[i][j][k] = eos->P(rho_0_guess,aux->eps(i,j,k));
	    aux->W[i][j][k] = 1.0;
	    aux->v_r[i][j][k] = 0.0;
	    aux->v_t[i][j][k] = 0.0;
	    aux->v_p[i][j][k] = 0.0;
	    error++;
	  }
	}
      }
    }
  };
  aux->rho_0.fill_ghosts();
  aux->eps.fill_ghosts();
  aux->v_r.fill_ghosts();
  aux->v_t.fill_ghosts();
  aux->v_p.fill_ghosts();
  aux->p.fill_ghosts();
  aux->W.fill_ghosts();
  m->D.fill_ghosts();
  m->S_r.fill_ghosts();
  m->S_t.fill_ghosts();
  m->S_p.fill_ghosts();
  m->tau.fill_ghosts();
  m->S_rad_r.fill_ghosts();
  m->S_rad_t.fill_ghosts();
  m->S_rad_p.fill_ghosts();
  m->tau_rad.fill_ghosts();
  aux->E.fill_ghosts();
  aux->f.fill_ghosts();
  aux->f_r.fill_ghosts();
  aux->f_t.fill_ghosts();
  aux->f_p.fill_ghosts();
  recovery_error_counter += error;
  return error;
};
//============================================================
// compute function f(p) (see eq. (A.5) in TBB)
//============================================================
double RadHydro::f(double p, double tau, double D, double S2, double & dfdp) {
  const double taupD = tau + p + D;
  const double taupD2 = taupD * taupD;
  //
  // compute W(p) (eq. (A.2))
  //
  const double W = taupD / sqrt( taupD2 - S2 );
  //
  // compute rho_0(p) and D(p) (eqs. (A.3) and (A.4))
  // 
  const double rho_0 = D/W;
  const double eps = ( sqrt(taupD2 - S2) - W * p - D) / D;
  //
  // compute partial rho / partial p (eq. (A.8))
  //
  const double drhodp = D * S2 / ( taupD2 * sqrt( taupD2 - S2 ) );
  //
  // compute partial eps / partial p (eq. (A.9))
  //
  const double depsdp = p * S2 / ( D * ( taupD2 - S2 ) * sqrt( taupD2 - S2 ) );
  //
  // compute P and its derivatives from rho_0 and eps using eos 
  //
  const double P = eos->P(rho_0,eps);
  const double chi = eos->dPdrho_0(rho_0,eps);
  const double kappa = eos->dPdepsilon(rho_0,eps);
  //
  // compute dfdp  (eq. (A.7))
  //
  dfdp = 1.0 - chi*drhodp - kappa*depsdp;
  //
  //return f = p - P
  //
  return p - P;
};

//
//================================================
// Reconstruct primitive variables at grid interfaces
//
// Fluxes are stored at cell interfaces,
//
//      p_L_i = p_i-1/2-eps
//      p_R_i = p_i-1/2+eps
//
//  Will loop over interfaces, and compute
// left and right fluxes at each interface
//
//     |   x   |   x   |   x   |   x   |   x   |
//                         i
//                   i-1/2
//           df_m    df_c    df_p
//              sigma_m sigma_p     
//                    L R
//
// For logarithmic grid: assume that cell interfaces are half-way between 
// grid-points
//
//================================================
//
int RadHydro::Reconstruct_r(gf3d & fct, gf3d & fct_L, gf3d & fct_R) {
  int error = 0;
  //
  // loop over interior grid interfaces
  //
  for (int i = N_g; i < N_r - N_g + 1; i++ ) {
    const double rmt = grid->r(i-2);
    const double rmo = grid->r(i-1);
    const double rl  = grid->r(i);
    const double rpo = grid->r(i+1); 
    for (int j = N_g; j < N_t - N_g; j++ ) 
      for (int k = N_g; k < N_p - N_g; k++ ) {
  	// recall: index i refers to interface i-1/2...
  	const double df_m = ( fct[i-1][j][k] - fct[i-2][j][k] ) / (rmo - rmt);
  	const double df_c = ( fct[i][j][k]   - fct[i-1][j][k] ) / (rl  - rmo);
  	const double df_p = ( fct[i+1][j][k] - fct[i][j][k]   ) / (rpo - rl);
	bool minmod;
	// use minmod no matter what in innermost 2 gridcells? 
	if (i < N_g + 2)
	  minmod = true;
	else
	  minmod = false;
	minmod = false;   // no...
  	const double sigma_m = slope->limiter(df_m,df_c,minmod);
  	const double sigma_p = slope->limiter(df_c,df_p,minmod);
  	// reconstruct at i-1/2
  	fct_L[i][j][k] = fct[i-1][j][k] + 0.5 * sigma_m * ( rl - rmo );
  	fct_R[i][j][k] = fct[i][j][k]   - 0.5 * sigma_p * ( rl - rmo );
      }
  }
  return error;
};
int RadHydro::Reconstruct_t(gf3d & fct, gf3d & fct_L, gf3d & fct_R) {
  int error = 0;
  //
  // loop over interior grid interfaces
  //
  for (int i = N_g; i < N_r - N_g; i++ )
    for (int j = N_g; j < N_t - N_g + 1; j++ ) {
      const double thetamt = grid->theta(j-2);
      const double thetamo = grid->theta(j-1);
      const double thetal = grid->theta(j);
      const double thetapo = grid->theta(j+1);
      for (int k = N_g; k < N_p - N_g; k++ ) {
	const double df_m = ( fct[i][j-1][k] - fct[i][j-2][k] ) / ( thetamo - thetamt );
	const double df_c = ( fct[i][j][k]   - fct[i][j-1][k] ) / ( thetal  - thetamo );
	const double df_p = ( fct[i][j+1][k] - fct[i][j][k]   ) / ( thetapo - thetal  );
	bool minmod = false;
	const double sigma_m = slope->limiter(df_m,df_c,minmod);
	const double sigma_p = slope->limiter(df_c,df_p,minmod);
	fct_L[i][j][k] = fct[i][j-1][k] + 0.5 * sigma_m * ( thetal - thetamo );
	fct_R[i][j][k] = fct[i][j][k]   - 0.5 * sigma_p * ( thetal - thetamo );
      }
    }
  return error;
};
int RadHydro::Reconstruct_p(gf3d &fct, gf3d & fct_L, gf3d & fct_R) {
  int error = 0;
  //
  // loop over interior grid interfaces
  //
  const double dphi = grid->delta_phi();    // so far grid in phi is uniform...
  for (int i = N_g; i < N_r - N_g; i++ )
    for (int j = N_g; j < N_t - N_g; j++ )
      for (int k = N_g; k < N_p - N_g + 1; k++ ) {
	const double df_m = ( fct[i][j][k-1] - fct[i][j][k-2] ) / dphi;
	const double df_c = ( fct[i][j][k] -   fct[i][j][k-1] ) / dphi;
	const double df_p = ( fct[i][j][k+1] - fct[i][j][k]   ) / dphi;
	bool minmod = false;
	const double sigma_m = slope->limiter(df_m,df_c,minmod);
	const double sigma_p = slope->limiter(df_c,df_p,minmod);
	fct_L[i][j][k] = fct[i][j][k-1] + 0.5 * sigma_m * dphi;
	fct_R[i][j][k] = fct[i][j][k]   - 0.5 * sigma_p * dphi;
      }
  return error;
};
//
//================================================
// Compute eigenvalues (see eqs. (35) and (36) in TBB)
//================================================
//
int RadHydro::lambda(double lapse, double v, double v2, 
		  double shift, double gamma, double cs,
		  double & lambda_0, double & lambda_p, double & lambda_m ) {
  int error = 0;
  if (v2 >= 1.0) { 
    cout << " Trouble in max_lambda: v2 = " << v2 << endl;
    error++;
  }
  if (cs >= 1.0) {
    cout << " Trouble in max_lambda: cs = " << cs << endl;
    error++;
  }
  const double cs2 = cs*cs;
  const double omv2cs2 = 1.0 - v2*cs2; 
  const double omcs2   = 1.0 - cs2; 
  const double help = cs * sqrt( (1.0 - v2) * ( gamma * omv2cs2 - v*v * omcs2 ) );
  lambda_0 = lapse*v - shift;
  lambda_p = lapse / omv2cs2 * ( v * omcs2 + help ) - shift;
  lambda_m = lapse / omv2cs2 * ( v * omcs2 - help ) - shift;
  return error;
};

//
//================================================
// Construct fluxes at interfaces
//================================================
//
int RadHydro::Construct_Fluxes(state *s, curvature *c) {
  int error = 0;
  //
  //==============================================
  // deal with r-direction
  //==============================================
  //
  // first reconstruct primitive variables at cell interfaces
  //
  Reconstruct_r(aux->rho_0,aux->rho_L,aux->rho_R);
  Reconstruct_r(aux->v_r,aux->v_r_L,aux->v_r_R);
  Reconstruct_r(aux->v_t,aux->v_t_L,aux->v_t_R);
  Reconstruct_r(aux->v_p,aux->v_p_L,aux->v_p_R);
  Reconstruct_r(aux->eps,aux->eps_L,aux->eps_R);
  Reconstruct_r(aux->E,aux->E_L,aux->E_R);
  Reconstruct_r(aux->f,aux->f_L,aux->f_R);
  Reconstruct_r(aux->f_r,aux->f_r_L,aux->f_r_R);
  Reconstruct_r(aux->f_t,aux->f_t_L,aux->f_t_R);
  Reconstruct_r(aux->f_p,aux->f_p_L,aux->f_p_R);
  //
  // now go to each interior cell interface
  //
  int i_test = -66;
  int j_test = 3;
  int k_test = 3;
  for (int i = N_g; i < N_r - N_g + 1; i++ ) {
    // NOTE: interface half-way between grid points:
    double r_l = 0.5 * ( grid->r(i-1) + grid->r(i) );
    for (int j = N_g; j < N_t - N_g; j++ ) {
      double sintheta_l = grid->sintheta(j);
      for (int k = N_g; k < N_p - N_g; k++ ) {
  	//
  	// interpolate grid functions to left interfaces
  	//
  	double h_rr_l = s->h_rr(r_l,j,k);
  	double h_rt_l = s->h_rt(r_l,j,k);
  	double h_rp_l = s->h_rp(r_l,j,k);
  	double h_tt_l = s->h_tt(r_l,j,k);
  	double h_tp_l = s->h_tp(r_l,j,k);
  	double h_pp_l = s->h_pp(r_l,j,k);
  	double phi_l  = s->phi(r_l,j,k);
  	double lapse_l = s->lapse(r_l,j,k);
  	double shift_r_l = s->shift_r(r_l,j,k);
  	double det_l = c->det(r_l,j,k);
  	double g_up_rr_l = ((1.0 + h_tt_l)*(1.0 + h_pp_l) - h_tp_l*h_tp_l)/det_l;
  	g_up_rr_l *= exp(- 4.0*phi_l);
  	//
  	// compute conserved variables, fluxes and max eigenvalues from left states
  	//
  	double D_L, S_r_L, S_t_L, S_p_L, tau_L, W_L, v2_L;
	double tau_rad_L, S_rad_r_L, S_rad_t_L, S_rad_p_L;
  	const double rho_L_l = aux->rho_L(i,j,k);
  	const double eps_L_l = aux->eps_L(i,j,k);
  	double p_L = eos->P(rho_L_l,eps_L_l);
  	Compute_Conserved_Vars(rho_L_l, eps_L_l, p_L,
			       aux->v_r_L(i,j,k), aux->v_t_L(i,j,k), aux->v_p_L(i,j,k),
			       aux->E_L(i,j,k), aux->f_L(i,j,k), 
			       aux->f_r_L(i,j,k), aux->f_t_L(i,j,k), aux->f_p_L(i,j,k),
			       h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, lapse_l, phi_l, det_l,
			       r_l, sintheta_l, 
			       D_L, S_r_L, S_t_L, S_p_L, tau_L, W_L, v2_L,
			       tau_rad_L, S_rad_r_L, S_rad_t_L, S_rad_p_L);
  	double f_D_L, f_S_r_L, f_S_t_L, f_S_p_L, f_tau_L;
	double f_tau_rad_L, f_S_rad_r_L, f_S_rad_t_L, f_S_rad_p_L;
	const double h_L = eos->h(rho_L_l,eps_L_l);
  	Compute_Fluxes(0, r_l, sintheta_l, lapse_l, phi_l, det_l, 
  		       D_L, aux->v_r_L(i,j,k), shift_r_l, S_r_L, S_t_L, S_p_L, p_L, h_L, W_L,
		       tau_L, aux->E_L(i,j,k), aux->f_L(i,j,k), aux->f_r_L(i,j,k), 
		       tau_rad_L, S_rad_r_L, S_rad_t_L, S_rad_p_L,
  		       f_D_L, f_S_r_L, f_S_t_L, f_S_p_L, f_tau_L,
		       f_tau_rad_L, f_S_rad_r_L, f_S_rad_t_L, f_S_rad_p_L);
	double cs = eos->sound_speed(rho_L_l,eps_L_l);
	// //
	// // Now compute sound speed in grid frame
	// //
	// const double W2 = W_L*W_L;
	// const double c2 = cs*cs;
	// const double aa = (W2 * ( 1.0 - c2)  + cs*cs) / lapse_l;
	// const double bb = - 2.0 *( W2 * ( 1.0 -  
	double lam_0_L, lam_p_L, lam_m_L;
 	error += lambda(lapse_l, aux->v_r_L(i,j,k), v2_L, shift_r_l, g_up_rr_l, cs, 
			lam_0_L, lam_p_L, lam_m_L);
  	//
  	// compute conserved variables, fluxes and max eigenvalues from right states
  	//
  	double D_R, S_r_R, S_t_R, S_p_R, tau_R, W_R, v2_R;
	double tau_rad_R, S_rad_r_R, S_rad_t_R, S_rad_p_R;
  	const double rho_R_l = aux->rho_R(i,j,k);
  	const double eps_R_l = aux->eps_R(i,j,k);
  	double p_R = eos->P(rho_R_l,eps_R_l);
  	Compute_Conserved_Vars(rho_R_l, eps_R_l, p_R,
			       aux->v_r_R(i,j,k), aux->v_t_R(i,j,k), aux->v_p_R(i,j,k), 
			       aux->E_R(i,j,k), aux->f_R(i,j,k), 
			       aux->f_r_R(i,j,k), aux->f_t_R(i,j,k), aux->f_p_R(i,j,k),
			       h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, lapse_l, phi_l, det_l,
			       r_l, sintheta_l, 
			       D_R, S_r_R, S_t_R, S_p_R, tau_R, W_R, v2_R,
			       tau_rad_R, S_rad_r_R, S_rad_t_R, S_rad_p_R);
  	double f_D_R, f_S_r_R, f_S_t_R, f_S_p_R, f_tau_R;
	double f_tau_rad_R, f_S_rad_r_R, f_S_rad_t_R, f_S_rad_p_R;
	const double h_R = eos->h(rho_R_l,eps_R_l);
  	Compute_Fluxes(0, r_l, sintheta_l, lapse_l, phi_l, det_l, 
  		       D_R, aux->v_r_R(i,j,k), shift_r_l, S_r_R, S_t_R, S_p_R, p_R, h_R, W_R,
		       tau_R, aux->E_R(i,j,k), aux->f_R(i,j,k), aux->f_r_R(i,j,k), 
		       tau_rad_R, S_rad_r_R, S_rad_t_R, S_rad_p_R,
  		       f_D_R, f_S_r_R, f_S_t_R, f_S_p_R, f_tau_R,
		       f_tau_rad_R, f_S_rad_r_R, f_S_rad_t_R, f_S_rad_p_R);
	cs = eos->sound_speed(rho_R_l,eps_R_l);
	double lam_0_R, lam_p_R, lam_m_R;
 	error += lambda(lapse_l, aux->v_r_R(i,j,k), v2_R, shift_r_l, g_up_rr_l, cs, 
			lam_0_R, lam_p_R, lam_m_R);
  	//
  	// now compute fluxes using approximate Riemann solver
  	//
  	aux->f_D_r[i][j][k]   = riemann_solver->flux(f_D_L,   f_D_R,   D_L,   D_R,   
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
  	aux->f_S_r_r[i][j][k] = riemann_solver->flux(f_S_r_L, f_S_r_R, S_r_L, S_r_R, 
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
  	aux->f_S_t_r[i][j][k] = riemann_solver->flux(f_S_t_L, f_S_t_R, S_t_L, S_t_R, 
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
  	aux->f_S_p_r[i][j][k] = riemann_solver->flux(f_S_p_L, f_S_p_R, S_p_L, S_p_R, 
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
  	aux->f_tau_r[i][j][k] = riemann_solver->flux(f_tau_L, f_tau_R, tau_L, tau_R, 
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
	// CHECK: is this right??
	cs = sqrt(1.0/3.0);
 	error += lambda(lapse_l, aux->v_r_L(i,j,k), v2_L, shift_r_l, g_up_rr_l, cs, 
			lam_0_L, lam_p_L, lam_m_L);
 	error += lambda(lapse_l, aux->v_r_R(i,j,k), v2_R, shift_r_l, g_up_rr_l, cs, 
			lam_0_R, lam_p_R, lam_m_R);

  	aux->f_S_rad_r_r[i][j][k] = riemann_solver->flux(f_S_rad_r_L, f_S_rad_r_R, S_rad_r_L, S_rad_r_R, 
						    lam_0_L, lam_p_L, lam_m_L, lam_0_R, lam_p_R, lam_m_R);
  	aux->f_S_rad_t_r[i][j][k] = riemann_solver->flux(f_S_rad_t_L, f_S_rad_t_R, S_rad_t_L, S_rad_t_R, 
						    lam_0_L, lam_p_L, lam_m_L, lam_0_R, lam_p_R, lam_m_R);
  	aux->f_S_rad_p_r[i][j][k] = riemann_solver->flux(f_S_rad_p_L, f_S_rad_p_R, S_rad_p_L, S_rad_p_R, 
						    lam_0_L, lam_p_L, lam_m_L, lam_0_R, lam_p_R, lam_m_R);
  	aux->f_tau_rad_r[i][j][k] = riemann_solver->flux(f_tau_rad_L, f_tau_rad_R, tau_rad_L, tau_rad_R, 
						    lam_0_L, lam_p_L, lam_m_L, lam_0_R, lam_p_R, lam_m_R);
	if (i == i_test && j == j_test && k == k_test) 
	  cout << " in Construct_Flux : " << setprecision(16) << aux->f_tau_r(i,j,k) 
	       << " tau = " << inter->tau(i,j,k) 
	       << " tau_R = " << tau_R << " tau_L = " << tau_L
	       << endl;
     }
    }
  }	
  //
  //==============================================
  // deal with theta-direction
  //==============================================
  //
  // first reconstruct primite variables at cell interfaces
  //
  Reconstruct_t(aux->rho_0,aux->rho_L,aux->rho_R);
  Reconstruct_t(aux->v_r,aux->v_r_L,aux->v_r_R);
  Reconstruct_t(aux->v_t,aux->v_t_L,aux->v_t_R);
  Reconstruct_t(aux->v_p,aux->v_p_L,aux->v_p_R);
  Reconstruct_t(aux->eps,aux->eps_L,aux->eps_R);
  Reconstruct_t(aux->E,aux->E_L,aux->E_R);
  Reconstruct_t(aux->f,aux->f_L,aux->f_R);
  Reconstruct_t(aux->f_r,aux->f_r_L,aux->f_r_R);
  Reconstruct_t(aux->f_t,aux->f_t_L,aux->f_t_R);
  Reconstruct_t(aux->f_p,aux->f_p_L,aux->f_p_R);
  //
  // now go to each interior cell interface
  //
  for (int i = N_g; i < N_r - N_g; i++ ) {
    double r_l = grid->r(i);
    for (int j = N_g; j < N_t - N_g + 1; j++ ) {
      double theta_l = 0.5* ( grid->theta(j-1) + grid->theta(j) );
      double sintheta_l = sin(theta_l);
      for (int k = N_g; k < N_p - N_g; k++ ) {
	//
	// interpolate grid functions to left interfaces
	//
	double h_rr_l = s->h_rr(i,theta_l,k);
	double h_rt_l = s->h_rt(i,theta_l,k);
	double h_rp_l = s->h_rp(i,theta_l,k);
	double h_tt_l = s->h_tt(i,theta_l,k);
	double h_tp_l = s->h_tp(i,theta_l,k);
	double h_pp_l = s->h_pp(i,theta_l,k);
	double phi_l  = s->phi(i,theta_l,k);
	double lapse_l = s->lapse(i,theta_l,k);
	double shift_t_l = s->shift_t(i,theta_l,k);
	double det_l = c->det(i,theta_l,k);
	double g_up_tt_l = ((1.0 + h_rr_l)*(1.0 + h_pp_l) - h_rp_l*h_rp_l)/det_l;
	g_up_tt_l *= exp(- 4.0*phi_l);
	//
	// compute conserved variables, fluxes and max eigenvalues from left states
	//
	double D_L, S_r_L, S_t_L, S_p_L, tau_L, W_L, v2_L;
	double tau_rad_L, S_rad_r_L, S_rad_t_L, S_rad_p_L;
	const double rho_L_l = aux->rho_L(i,j,k);
	const double eps_L_l = aux->eps_L(i,j,k);
	double p_L = eos->P(rho_L_l,eps_L_l);
	Compute_Conserved_Vars(rho_L_l, eps_L_l, p_L,
			       aux->v_r_L(i,j,k), aux->v_t_L(i,j,k), aux->v_p_L(i,j,k), 
			       aux->E_L(i,j,k), aux->f_L(i,j,k), 
			       aux->f_r_L(i,j,k), aux->f_t_L(i,j,k), aux->f_p_L(i,j,k),
			       h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, lapse_l, phi_l, det_l,
			       r_l, sintheta_l, 
			       D_L, S_r_L, S_t_L, S_p_L, tau_L, W_L, v2_L,
			       tau_rad_L, S_rad_r_L, S_rad_t_L, S_rad_p_L);
	double f_D_L, f_S_r_L, f_S_t_L, f_S_p_L, f_tau_L;
	double f_tau_rad_L, f_S_rad_r_L, f_S_rad_t_L, f_S_rad_p_L;
	const double h_L = eos->h(rho_L_l,eps_L_l);
	Compute_Fluxes(1, r_l, sintheta_l, lapse_l, phi_l, det_l, 
		       D_L, aux->v_t_L(i,j,k), shift_t_l, S_r_L, S_t_L, S_p_L, p_L, h_L, W_L, 
		       tau_L, aux->E_L(i,j,k), aux->f_L(i,j,k), aux->f_t_L(i,j,k), 
		       tau_rad_L, S_rad_r_L, S_rad_t_L, S_rad_p_L,
		       f_D_L, f_S_r_L, f_S_t_L, f_S_p_L, f_tau_L,
		       f_tau_rad_L, f_S_rad_r_L, f_S_rad_t_L, f_S_rad_p_L);
	double cs = eos->sound_speed(rho_L_l,eps_L_l);
	double lam_0_L, lam_p_L, lam_m_L;
 	error += lambda(lapse_l, aux->v_t_L(i,j,k), v2_L, shift_t_l, g_up_tt_l, cs, lam_0_L, lam_p_L, lam_m_L);
	//
	// compute conserved variables, fluxes and max eigenvalues from right states
	//
	double D_R, S_r_R, S_t_R, S_p_R, tau_R, W_R, v2_R;
	double tau_rad_R, S_rad_r_R, S_rad_t_R, S_rad_p_R;
	const double rho_R_l = aux->rho_R(i,j,k);
	const double eps_R_l = aux->eps_R(i,j,k);
	double p_R = eos->P(rho_R_l,eps_R_l);
	Compute_Conserved_Vars(rho_R_l, eps_R_l, p_R,
			       aux->v_r_R(i,j,k), aux->v_t_R(i,j,k), aux->v_p_R(i,j,k), 
			       aux->E_R(i,j,k), aux->f_R(i,j,k), 
			       aux->f_r_R(i,j,k), aux->f_t_R(i,j,k), aux->f_p_R(i,j,k),
			       h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, lapse_l, phi_l, det_l,
			       r_l, sintheta_l, 
			       D_R, S_r_R, S_t_R, S_p_R, tau_R, W_R, v2_R,
			       tau_rad_R, S_rad_r_R, S_rad_t_R, S_rad_p_R);
	double f_D_R, f_S_r_R, f_S_t_R, f_S_p_R, f_tau_R;
	double f_tau_rad_R, f_S_rad_r_R, f_S_rad_t_R, f_S_rad_p_R;
	const double h_R = eos->h(rho_R_l,eps_R_l);
	Compute_Fluxes(1, r_l, sintheta_l, lapse_l, phi_l, det_l, 
		       D_R, aux->v_t_R(i,j,k), shift_t_l, S_r_R, S_t_R, S_p_R, p_R, h_R, W_R, 
		       tau_R, aux->E_R(i,j,k), aux->f_R(i,j,k), aux->f_t_R(i,j,k), 
		       tau_rad_R, S_rad_r_R, S_rad_t_R, S_rad_p_R,
		       f_D_R, f_S_r_R, f_S_t_R, f_S_p_R, f_tau_R,
		       f_tau_rad_R, f_S_rad_r_R, f_S_rad_t_R, f_S_rad_p_R);
	cs = eos->sound_speed(rho_R_l,eps_R_l);
	double lam_0_R, lam_p_R, lam_m_R;
 	error += lambda(lapse_l, aux->v_t_R(i,j,k), v2_R, shift_t_l, g_up_tt_l, cs, lam_0_R, lam_p_R, lam_m_R);
	//
	// now compute fluxes using approximate Riemann solver
	//
	aux->f_D_t[i][j][k]   = riemann_solver->flux(f_D_L,   f_D_R,   D_L,   D_R,   
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
	aux->f_S_r_t[i][j][k] = riemann_solver->flux(f_S_r_L, f_S_r_R, S_r_L, S_r_R, 
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
	aux->f_S_t_t[i][j][k] = riemann_solver->flux(f_S_t_L, f_S_t_R, S_t_L, S_t_R, 
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
	aux->f_S_p_t[i][j][k] = riemann_solver->flux(f_S_p_L, f_S_p_R, S_p_L, S_p_R, 
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
	aux->f_tau_t[i][j][k] = riemann_solver->flux(f_tau_L, f_tau_R, tau_L, tau_R, 
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);

	// CHECK!!
	cs = sqrt(1.0/3.0);
 	error += lambda(lapse_l, aux->v_t_L(i,j,k), v2_L, shift_t_l, g_up_tt_l, cs, 
			lam_0_L, lam_p_L, lam_m_L);
 	error += lambda(lapse_l, aux->v_t_R(i,j,k), v2_R, shift_t_l, g_up_tt_l, cs, 
			lam_0_R, lam_p_R, lam_m_R);
  	aux->f_S_rad_r_t[i][j][k] = riemann_solver->flux(f_S_rad_r_L, f_S_rad_r_R, S_rad_r_L, S_rad_r_R, 
							 lam_0_L, lam_p_L, lam_m_L, 
							 lam_0_R, lam_p_R, lam_m_R);
  	aux->f_S_rad_t_t[i][j][k] = riemann_solver->flux(f_S_rad_t_L, f_S_rad_t_R, S_rad_t_L, S_rad_t_R, 
							 lam_0_L, lam_p_L, lam_m_L, 
							 lam_0_R, lam_p_R, lam_m_R);
  	aux->f_S_rad_p_t[i][j][k] = riemann_solver->flux(f_S_rad_p_L, f_S_rad_p_R, S_rad_p_L, S_rad_p_R, 
							 lam_0_L, lam_p_L, lam_m_L, 
							 lam_0_R, lam_p_R, lam_m_R);
  	aux->f_tau_rad_t[i][j][k] = riemann_solver->flux(f_tau_rad_L, f_tau_rad_R, tau_rad_L, tau_rad_R, 
							 lam_0_L, lam_p_L, lam_m_L, 
							 lam_0_R, lam_p_R, lam_m_R);

      }
    }
  }	
  //
  //==============================================
  // deal with phi-direction
  //==============================================
  //
  // first reconstruct primite variables at cell interfaces
  //
  Reconstruct_p(aux->rho_0,aux->rho_L,aux->rho_R);
  Reconstruct_p(aux->v_r,aux->v_r_L,aux->v_r_R);
  Reconstruct_p(aux->v_t,aux->v_t_L,aux->v_t_R);
  Reconstruct_p(aux->v_p,aux->v_p_L,aux->v_p_R);
  Reconstruct_p(aux->eps,aux->eps_L,aux->eps_R);
  Reconstruct_p(aux->E,aux->E_L,aux->E_R);
  Reconstruct_p(aux->f,aux->f_L,aux->f_R);
  Reconstruct_p(aux->f_r,aux->f_r_L,aux->f_r_R);
  Reconstruct_p(aux->f_t,aux->f_t_L,aux->f_t_R);
  Reconstruct_p(aux->f_p,aux->f_p_L,aux->f_p_R);
  //
  // now go to each interior cell interface
  //
  for (int i = N_g; i < N_r - N_g; i++ ) {
    double r_l = grid->r(i);
    for (int j = N_g; j < N_t - N_g; j++ ) {
      double sintheta_l = grid->sintheta(j);
      for (int k = N_g; k < N_p - N_g + 1; k++ ) {
	double p_l = 0.5 * ( grid->phi(k-1) + grid->phi(k) );
	//
	// interpolate grid functions to left interfaces
	//
	double h_rr_l = s->h_rr(i,j,p_l);
	double h_rt_l = s->h_rt(i,j,p_l);
	double h_rp_l = s->h_rp(i,j,p_l);
	double h_tt_l = s->h_tt(i,j,p_l);
	double h_tp_l = s->h_tp(i,j,p_l);
	double h_pp_l = s->h_pp(i,j,p_l);
	double phi_l  = s->phi(i,j,p_l);
	double lapse_l = s->lapse(i,j,p_l);
	double shift_p_l = s->shift_p(i,j,p_l);
	double det_l = c->det(i,j,p_l);
	double g_up_pp_l = ((1.0 + h_rr_l)*(1.0 + h_tt_l) - h_rt_l*h_rt_l)/det_l;
	g_up_pp_l *= exp(- 4.0*phi_l);
	//
	// compute conserved variables, fluxes and max eigenvalues from left states
	//
	double D_L, S_r_L, S_t_L, S_p_L, tau_L, W_L, v2_L;
	double tau_rad_L, S_rad_r_L, S_rad_t_L, S_rad_p_L;
	const double rho_L_l = aux->rho_L(i,j,k);
	const double eps_L_l = aux->eps_L(i,j,k);
	double p_L = eos->P(rho_L_l,eps_L_l);
	Compute_Conserved_Vars(rho_L_l, eps_L_l, p_L,
			       aux->v_r_L(i,j,k), aux->v_t_L(i,j,k), aux->v_p_L(i,j,k), 
			       aux->E_L(i,j,k), aux->f_L(i,j,k), 
			       aux->f_r_L(i,j,k), aux->f_t_L(i,j,k), aux->f_p_L(i,j,k),
			       h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, lapse_l, phi_l, det_l,
			       r_l, sintheta_l, 
			       D_L, S_r_L, S_t_L, S_p_L, tau_L, W_L, v2_L,
			       tau_rad_L, S_rad_r_L, S_rad_t_L, S_rad_p_L);
	double f_D_L, f_S_r_L, f_S_t_L, f_S_p_L, f_tau_L;
	double f_tau_rad_L, f_S_rad_r_L, f_S_rad_t_L, f_S_rad_p_L;
	const double h_L = eos->h(rho_L_l,eps_L_l);
	Compute_Fluxes(2, r_l, sintheta_l, lapse_l, phi_l, det_l, 
		       D_L, aux->v_p_L(i,j,k), shift_p_l, S_r_L, S_t_L, S_p_L, p_L, h_L, W_L,
		       tau_L, aux->E_L(i,j,k), aux->f_L(i,j,k), aux->f_p_L(i,j,k), 
		       tau_rad_L, S_rad_r_L, S_rad_t_L, S_rad_p_L,
		       f_D_L, f_S_r_L, f_S_t_L, f_S_p_L, f_tau_L,
		       f_tau_rad_L, f_S_rad_r_L, f_S_rad_t_L, f_S_rad_p_L);
	double cs = eos->sound_speed(rho_L_l,eps_L_l);
	double lam_0_L, lam_p_L, lam_m_L;
 	error += lambda(lapse_l, aux->v_p_L(i,j,k), v2_L, shift_p_l, g_up_pp_l, cs, 
			lam_0_L, lam_p_L, lam_m_L);
	//
	// compute conserved variables, fluxes and max eigenvalues from right states
	//
	double D_R, S_r_R, S_t_R, S_p_R, tau_R, W_R, v2_R;
	double tau_rad_R, S_rad_r_R, S_rad_t_R, S_rad_p_R;
	const double rho_R_l = aux->rho_R(i,j,k);
	const double eps_R_l = aux->eps_R(i,j,k);
	double p_R = eos->P(rho_R_l,eps_R_l);
	Compute_Conserved_Vars(rho_R_l, eps_R_l, p_R,
			       aux->v_r_R(i,j,k), aux->v_t_R(i,j,k), aux->v_p_R(i,j,k), 
			       aux->E_R(i,j,k), aux->f_R(i,j,k), 
			       aux->f_r_R(i,j,k), aux->f_t_R(i,j,k), aux->f_p_R(i,j,k),
			       h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, lapse_l, phi_l, det_l,
			       r_l, sintheta_l, 
			       D_R, S_r_R, S_t_R, S_p_R, tau_R, W_R, v2_R,
			       tau_rad_R, S_rad_r_R, S_rad_t_R, S_rad_p_R);
	double f_D_R, f_S_r_R, f_S_t_R, f_S_p_R, f_tau_R;
	double f_tau_rad_R, f_S_rad_r_R, f_S_rad_t_R, f_S_rad_p_R;
	const double h_R = eos->h(rho_R_l,eps_R_l);
	Compute_Fluxes(2, r_l, sintheta_l, lapse_l, phi_l, det_l, 
		       D_R, aux->v_p_R(i,j,k), shift_p_l, S_r_R, S_t_R, S_p_R, p_R, h_R, W_R,
		       tau_R, aux->E_R(i,j,k), aux->f_R(i,j,k), aux->f_p_R(i,j,k), 
		       tau_rad_R, S_rad_r_R, S_rad_t_R, S_rad_p_R,
		       f_D_R, f_S_r_R, f_S_t_R, f_S_p_R, f_tau_R,
		       f_tau_rad_R, f_S_rad_r_R, f_S_rad_t_R, f_S_rad_p_R);
	cs = eos->sound_speed(rho_R_l,eps_R_l);
	double lam_0_R, lam_p_R, lam_m_R;
 	error += lambda(lapse_l, aux->v_p_R(i,j,k), v2_R, shift_p_l, g_up_pp_l, cs, 
			lam_0_R, lam_p_R, lam_m_R);
	//
	// now compute fluxes using approximate Riemann solver
	//
	aux->f_D_p[i][j][k]   = riemann_solver->flux(f_D_L,   f_D_R,   D_L,   D_R,   
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
	aux->f_S_r_p[i][j][k] = riemann_solver->flux(f_S_r_L, f_S_r_R, S_r_L, S_r_R, 
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
	aux->f_S_t_p[i][j][k] = riemann_solver->flux(f_S_t_L, f_S_t_R, S_t_L, S_t_R, 
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
	aux->f_S_p_p[i][j][k] = riemann_solver->flux(f_S_p_L, f_S_p_R, S_p_L, S_p_R, 
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
	aux->f_tau_p[i][j][k] = riemann_solver->flux(f_tau_L, f_tau_R, tau_L, tau_R, 
						     lam_0_L, lam_p_L, lam_m_L, 
						     lam_0_R, lam_p_R, lam_m_R);
	// CHECK!!
	cs = sqrt(1.0/3.0);
 	error += lambda(lapse_l, aux->v_p_L(i,j,k), v2_L, shift_p_l, g_up_pp_l, cs, 
			lam_0_L, lam_p_L, lam_m_L);
 	error += lambda(lapse_l, aux->v_p_R(i,j,k), v2_R, shift_p_l, g_up_pp_l, cs, 
			lam_0_R, lam_p_R, lam_m_R);
  	aux->f_S_rad_r_p[i][j][k] = riemann_solver->flux(f_S_rad_r_L, f_S_rad_r_R, S_rad_r_L, S_rad_r_R, 
							 lam_0_L, lam_p_L, lam_m_L, 
							 lam_0_R, lam_p_R, lam_m_R);
  	aux->f_S_rad_t_p[i][j][k] = riemann_solver->flux(f_S_rad_t_L, f_S_rad_t_R, S_rad_t_L, S_rad_t_R, 
							 lam_0_L, lam_p_L, lam_m_L, 
							 lam_0_R, lam_p_R, lam_m_R);
  	aux->f_S_rad_p_p[i][j][k] = riemann_solver->flux(f_S_rad_p_L, f_S_rad_p_R, S_rad_p_L, S_rad_p_R, 
							 lam_0_L, lam_p_L, lam_m_L, 
							 lam_0_R, lam_p_R, lam_m_R);
  	aux->f_tau_rad_p[i][j][k] = riemann_solver->flux(f_tau_rad_L, f_tau_rad_R, tau_rad_L, tau_rad_R, 
							 lam_0_L, lam_p_L, lam_m_L, 
							 lam_0_R, lam_p_R, lam_m_R);
      }
    }
  }
  //
  return error;
};
//
//================================================
// Compute Fluxes from conserved radhydro variables 
//
// Get correct component by 
// (1) passing in corresponding component of v and shift  -- i.e., pass 
//     in v^r and \beta^r, get r components of fluxes)
// (2) pass in correct index (0 for r, 1 for theta, 2 for phi)
//
//================================================
//
int RadHydro::Compute_Fluxes(int index, double r, double sintheta,
			     double lapse, double phi, double det, double D,
			     double v, double shift, double S_r, double S_t, double S_p, 
			     double p, double h, double W,
			     double tau, double E, double f, double F,  // note: f = f, F = f^i...
			     double tau_rad,
			     double S_rad_r, double S_rad_t, double S_rad_p, 
			     double & f_D, 
			     double & f_S_r, double & f_S_t, double & f_S_p,
			     double & f_tau, double & f_tau_rad,
			     double & f_S_rad_r, double & f_S_rad_t, double & f_S_rad_p, bool verbose) {
  int error = 0;
#ifdef SR
  double det_3D = 1.0;
  const double vel = v;
#else
  double det_3D = exp(6.0*phi) * sqrt(det);
  const double vel = v - shift/lapse;
#endif
  f_D = lapse * D * vel;
  f_S_r = lapse * ( S_r * vel + det_3D * (index == 0) * p); 
  f_S_t = lapse * ( S_t * vel + det_3D * (index == 1) * p); 
  f_S_p = lapse * ( S_p * vel + det_3D * (index == 2) * p); 
  if (partial) {
    f_tau = lapse * ( tau * vel + det_3D * r*r*sintheta * p * v );
    f_tau_rad = lapse * (tau_rad * vel + det_3D * r*r*sintheta * (W * F - W * f * v + E * vel/3.0) );
  } else { 
    f_tau = lapse * ( tau * vel + det_3D * p * v );
    f_tau_rad = lapse * (tau_rad * vel + det_3D * (W * F - W * f * v + E * v/3.0) );
  }
  if (verbose) {
    cout << " In Compute_Fluxes for r = " << r << " sintheta = " << sintheta << endl;
    cout << "tau = " << tau << " vel = " << vel << " p = " 
	 << p << " D = " << D << " f_tau = " << f_tau << endl;
  }
  //
  // compute lower components of v_i from S_i:  v_i = S_i / (W D h)
  //
  // FIX this for partial approach - need to compute u_i differently at r = 0...
  //
  //  double Dl =  D / det_3D;
  double Dl =  D;
  if (partial)
    Dl /= r*r*sintheta;
  const double WDh = W * Dl * h;
  const double v_low_r = S_r / WDh;
  const double v_low_t = S_t / WDh;
  const double v_low_p = S_p / WDh;
  f_S_rad_r = lapse * ( S_rad_r * vel + det_3D* ( W*v_low_r*(F - f*v) + (index == 0) * E / 3.0 ));  
  f_S_rad_t = lapse * ( S_rad_t * vel + det_3D* ( W*v_low_t*(F - f*v) + (index == 1) * E / 3.0 ));  
  f_S_rad_p = lapse * ( S_rad_p * vel + det_3D* ( W*v_low_p*(F - f*v) + (index == 2) * E / 3.0 ));  
  if (!isfinite(f_S_rad_r)) cout << " GOT IT! " 
			       << " r = " << r << " sintheta = " << sintheta << " D = " << D
			       << endl;
  return error;
};
//
//================================================
// 
// Compute right-hand sides...
//
// Three parts:
//   - divergence of fluxes
//   - flat connection terms
//   - source terms
//
//================================================
//
int RadHydro::dot_conserved_vars(radhydro_state *m, state *s, curvature *c) {
  int error = 0;
  int i_test = -66;
  int j_test = 3;
  int k_test = 3;
  //
  // loop over interior grid
  //
  for (int i = N_g; i < N_r - N_g; i++ ) {
    const double rl = grid->r(i);
    //
    // compute dr : distance between cell-interfaces (different
    // from distances between grid points!
    //
    const double dr = 0.5 * ( grid->r(i+1) - grid->r(i-1) );
    for (int j = N_g; j < N_t - N_g; j++ ) {
      const double sintheta = grid->sintheta(j);
      const double rst = rl * sintheta;
      const double sin2theta = sintheta * sintheta;
      const double costheta = grid->costheta(j);
      const double cottheta = costheta / sintheta;
      const double dtheta = 0.5 * ( grid->theta(j+1) - grid->theta(j-1) );
      for (int k = N_g; k < N_p - N_g; k++ ) {
	//
	const double dphi = grid->delta_phi();
	const double lapsel = s->lapse(i,j,k);
	// const double phil = s->phi(i,j,k);
	// const double detl = c->det(i,j,k);
	// const double shift_rl = s->shift_r(i,j,k);
	// const double shift_tl = s->shift_t(i,j,k);
	// const double shift_pl = s->shift_p(i,j,k);
	//================================================
	// first compute divergence of fluxes - note extra factors 
	// and terms due to rescaling
	//================================================
	derivs->D[i][j][k] = -
	  (aux->f_D_r[i+1][j][k] - aux->f_D_r[i][j][k])/dr -
	  (aux->f_D_t[i][j+1][k] - aux->f_D_t[i][j][k])/dtheta/rl -
	  (aux->f_D_p[i][j][k+1] - aux->f_D_p[i][j][k])/dphi/rst;
	derivs->S_r[i][j][k] = -
	  (aux->f_S_r_r[i+1][j][k] - aux->f_S_r_r[i][j][k])/dr -
	  (aux->f_S_r_t[i][j+1][k] - aux->f_S_r_t[i][j][k])/dtheta/rl -
	  (aux->f_S_r_p[i][j][k+1] - aux->f_S_r_p[i][j][k])/dphi/rst;
	derivs->S_t[i][j][k] = -
	  (aux->f_S_t_r[i+1][j][k] - aux->f_S_t_r[i][j][k])/dr -
	  (aux->f_S_t_t[i][j+1][k] - aux->f_S_t_t[i][j][k])/dtheta/rl -
	  (aux->f_S_t_p[i][j][k+1] - aux->f_S_t_p[i][j][k])/dphi/rst;
	derivs->S_p[i][j][k] = -
	  (aux->f_S_p_r[i+1][j][k] - aux->f_S_p_r[i][j][k])/dr -
	  (aux->f_S_p_t[i][j+1][k] - aux->f_S_p_t[i][j][k])/dtheta/rl -
	  (aux->f_S_p_p[i][j][k+1] - aux->f_S_p_p[i][j][k])/dphi/rst;
	derivs->tau[i][j][k] = -
	  (aux->f_tau_r[i+1][j][k] - aux->f_tau_r[i][j][k])/dr -
	  (aux->f_tau_t[i][j+1][k] - aux->f_tau_t[i][j][k])/dtheta/rl -
	  (aux->f_tau_p[i][j][k+1] - aux->f_tau_p[i][j][k])/dphi/rst;
	if (i == i_test && j == j_test && k == k_test) cout << " 1 " << setprecision(16) << derivs->tau(i,j,k) << " f_tau_r = " << aux->f_tau_r(i,j,k) << "   " << aux->f_tau_r(i+1,j,k) << endl;
	// if (i == 2 && j == 2 && k == 2) cout << " 1 : tau-dot = " << derivs->tau[i][j][k] << endl;
	derivs->tau_rad[i][j][k] = -
	  (aux->f_tau_rad_r[i+1][j][k] - aux->f_tau_rad_r[i][j][k])/dr -
	  (aux->f_tau_rad_t[i][j+1][k] - aux->f_tau_rad_t[i][j][k])/dtheta/rl -
	  (aux->f_tau_rad_p[i][j][k+1] - aux->f_tau_rad_p[i][j][k])/dphi/rst;
	derivs->S_rad_r[i][j][k] = -
	  (aux->f_S_rad_r_r[i+1][j][k] - aux->f_S_rad_r_r[i][j][k])/dr -
	  (aux->f_S_rad_r_t[i][j+1][k] - aux->f_S_rad_r_t[i][j][k])/dtheta/rl -
	  (aux->f_S_rad_r_p[i][j][k+1] - aux->f_S_rad_r_p[i][j][k])/dphi/rst;
	derivs->S_rad_t[i][j][k] = -
	  (aux->f_S_rad_t_r[i+1][j][k] - aux->f_S_rad_t_r[i][j][k])/dr -
	  (aux->f_S_rad_t_t[i][j+1][k] - aux->f_S_rad_t_t[i][j][k])/dtheta/rl -
	  (aux->f_S_rad_t_p[i][j][k+1] - aux->f_S_rad_t_p[i][j][k])/dphi/rst;
	derivs->S_rad_p[i][j][k] = -
	  (aux->f_S_rad_p_r[i+1][j][k] - aux->f_S_rad_p_r[i][j][k])/dr -
	  (aux->f_S_rad_p_t[i][j+1][k] - aux->f_S_rad_p_t[i][j][k])/dtheta/rl -
	  (aux->f_S_rad_p_p[i][j][k+1] - aux->f_S_rad_p_p[i][j][k])/dphi/rst;

	//================================================
	// compute fluxes at cell centers for flat connection terms
	// as well as extra terms due to rescaling
	//================================================
	double f_D_r_c, f_D_t_c, f_D_p_c;
	double f_S_r_r_c, f_S_r_t_c, f_S_r_p_c;
	double f_S_t_r_c, f_S_t_t_c, f_S_t_p_c;
	double f_S_p_r_c, f_S_p_t_c, f_S_p_p_c;
	double f_tau_r_c, f_tau_t_c, f_tau_p_c;
	double f_tau_rad_r_c, f_tau_rad_t_c, f_tau_rad_p_c;
	double f_S_rad_r_r_c, f_S_rad_r_t_c, f_S_rad_r_p_c;
	double f_S_rad_t_r_c, f_S_rad_t_t_c, f_S_rad_t_p_c;
	double f_S_rad_p_r_c, f_S_rad_p_t_c, f_S_rad_p_p_c;
	const double hl = eos->h(aux->rho_0(i,j,k),aux->eps(i,j,k));
	bool verbose = false;
	// radial fluxes
	Compute_Fluxes(0, rl, sintheta, s->lapse[i][j][k], s->phi[i][j][k], c->det[i][j][k],
		       m->D[i][j][k], aux->v_r[i][j][k], s->shift_r[i][j][k],
		       m->S_r[i][j][k], m->S_t[i][j][k], m->S_p[i][j][k],
		       aux->p[i][j][k], hl, aux->W[i][j][k], 
		       m->tau[i][j][k], aux->E[i][j][k], aux->f[i][j][k], 
		       aux->f_r[i][j][k], m->tau_rad[i][j][k],
		       m->S_rad_r[i][j][k], m->S_rad_t[i][j][k], m->S_rad_p[i][j][k],
		       f_D_r_c, f_S_r_r_c, f_S_t_r_c, f_S_p_r_c, f_tau_r_c,
		       f_tau_rad_r_c, f_S_rad_r_r_c, f_S_rad_t_r_c, f_S_rad_p_r_c, verbose);
	// theta fluxes
	Compute_Fluxes(1, rl, sintheta, s->lapse[i][j][k], s->phi[i][j][k], c->det[i][j][k],
		       m->D[i][j][k], aux->v_t[i][j][k], s->shift_t[i][j][k],
		       m->S_r[i][j][k], m->S_t[i][j][k], m->S_p[i][j][k],
		       aux->p[i][j][k], hl, aux->W[i][j][k], 
		       m->tau[i][j][k], aux->E[i][j][k], aux->f[i][j][k], 
		       aux->f_t[i][j][k], m->tau_rad[i][j][k],
		       m->S_rad_r[i][j][k], m->S_rad_t[i][j][k], m->S_rad_p[i][j][k],
		       f_D_t_c, f_S_r_t_c, f_S_t_t_c, f_S_p_t_c, f_tau_t_c,
		       f_tau_rad_t_c, f_S_rad_r_t_c, f_S_rad_t_t_c, f_S_rad_p_t_c);

	// // phi fluxes
	Compute_Fluxes(2, rl, sintheta, s->lapse[i][j][k], s->phi[i][j][k], c->det[i][j][k],
		       m->D[i][j][k], aux->v_p[i][j][k], s->shift_p[i][j][k],
		       m->S_r[i][j][k], m->S_t[i][j][k], m->S_p[i][j][k],
		       aux->p[i][j][k], hl, aux->W[i][j][k], 
		       m->tau[i][j][k], aux->E[i][j][k], aux->f[i][j][k], 
		       aux->f_p[i][j][k], m->tau_rad[i][j][k],
		       m->S_rad_r[i][j][k], m->S_rad_t[i][j][k], m->S_rad_p[i][j][k],
		       f_D_p_c, f_S_r_p_c, f_S_t_p_c, f_S_p_p_c, f_tau_p_c,
		       f_tau_rad_p_c, f_S_rad_r_p_c, f_S_rad_t_p_c, f_S_rad_p_p_c);
	//
	if (!partial) {
	  derivs->D[i][j][k] += - f_D_r_c * 2.0 / rl - f_D_t_c * cottheta / rl; 
	  derivs->tau[i][j][k] += - f_tau_r_c * 2.0 / rl - f_tau_t_c * cottheta / rl; 
	  // if (i == 2 && j == 2 && k == 2) cout << " 1 : tau-dot = " << derivs->tau[i][j][k] << endl;
	  derivs->tau_rad[i][j][k] += ( - f_tau_rad_r_c * 2.0 - f_tau_rad_t_c * cottheta ) / rl; 
	}
	if (i == i_test && j == j_test && k == k_test) cout << " 2 " << derivs->tau(i,j,k) << " D = " << derivs->D(i,j,k) << endl;
	derivs->S_r[i][j][k] += (f_S_t_t_c + f_S_p_p_c) / rl 
	  - f_S_r_r_c * 2.0 / rl - f_S_r_t_c * cottheta / rl;  // term due to rescaling
	derivs->S_t[i][j][k] += - f_S_r_t_c / rl + f_S_t_r_c / rl + f_S_p_p_c * cottheta / rl  
	  - f_S_t_r_c * 2.0 / rl - f_S_t_t_c * cottheta / rl  // term due to rescaling 
	  - f_S_t_r_c / rl;   // term due to rescaling
	derivs->S_p[i][j][k] += - f_S_r_p_c / rl - f_S_t_p_c * cottheta / rl 
	  + f_S_p_r_c / rl + f_S_p_t_c * cottheta / rl
	  - f_S_p_r_c * 2.0 / rl - f_S_p_t_c * cottheta / rl
	  - f_S_p_r_c / rl - cottheta / rl * f_S_p_t_c;  // terms due to rescaling
	derivs->S_rad_r[i][j][k] += ( f_S_rad_t_t_c + f_S_rad_p_p_c 
					  - f_S_rad_r_r_c * 2.0  - f_S_rad_r_t_c * cottheta ) / rl;  
	derivs->S_rad_t[i][j][k] += ( - f_S_rad_r_t_c + f_S_rad_t_r_c + f_S_rad_p_p_c * cottheta   
					  - f_S_rad_t_r_c * 2.0 - f_S_rad_t_t_c * cottheta 
					  - f_S_rad_t_r_c ) / rl;   // term due to rescaling
	derivs->S_rad_p[i][j][k] += ( - f_S_rad_r_p_c - f_S_rad_t_p_c * cottheta 
					  + f_S_rad_p_r_c + f_S_rad_p_t_c * cottheta
					  - f_S_rad_p_r_c * 2.0 - f_S_rad_p_t_c * cottheta
					  // terms due to rescaling:
					  - f_S_rad_p_r_c - cottheta * f_S_rad_p_t_c  ) / rl; 
	//================================================
	// Finally compute source terms on right-hand side
	//================================================
	//
	// first compute components of stress-energy tensor (see eqs. (27) - (30) in TBB)
	//
	const double rho_0l = aux->rho_0(i,j,k);
	const double enth = eos->h(rho_0l, aux->eps(i,j,k));
	const double p = eos->P(rho_0l, aux->eps(i,j,k));   // CHECK: do we need this call?
	const double Wl = aux->W(i,j,k);
	const double rhohW2 = rho_0l * enth * Wl*Wl;
	const double u0 = Wl / lapsel;
	//
	// NOTE: We will do the following calculation in terms of physical, i.e. *not* rescaled
	// variables.  Why?  D_i \gamma_jk in rest of code is not rescaled,
	// and it is easier to use that variable here, too.
	//
	//
	// physical shift, velocity and flux
	vect shift(s->shift_r(i,j,k), s->shift_t(i,j,k)/rl, s->shift_p(i,j,k)/rst);
	vect v_up(aux->v_r(i,j,k), aux->v_t(i,j,k)/rl, aux->v_p(i,j,k)/rst);   
	vect f_up(aux->f_r(i,j,k), aux->f_t(i,j,k)/rl, aux->f_p(i,j,k)/rst);
	vect v_down;
	vect f_down;
 	tensor g_conf(1.0 + s->h_rr(i,j,k), 
		      rl * s->h_rt(i,j,k), 
		      rst * s->h_rp(i,j,k),
 		      rl*rl * (1.0 + s->h_tt(i,j,k)), 
		      rl*rst * s->h_tp(i,j,k),  
		      rst*rst*(1.0 + s->h_pp(i,j,k)));
	tensor g_up = g_conf.inverse();
	const double e4p = exp( 4.0 * s->phi(i,j,k) );
	//
	// lower index of v^i and f^i 
	// (result is physical velocity, not rescaled)
	//
	for (int a = 0; a < 3; a++) {
	  v_down[a] = 0.0;
	  f_down[a] = 0.0;
	  for (int b = 0; b < 3; b++) {
	    v_down[a] += e4p * g_conf[a][b] * v_up[b];
	    f_down[a] += e4p * g_conf[a][b] * f_up[b];
	  }
	}
	//
	// compute spatial components of four-velocity
	//
	vect u_up;
	vect u_down;
	for (int a = 0; a < 3; a++) {
	  u_up[a] = Wl * ( v_up[a] -  shift[a] / lapsel );
	  u_down[a] = Wl * v_down[a];
	}
	//
	// compute components of hydro and radiation stress-energy tensor (use *unrescaled* shift)
	//
	const double El = aux->E(i,j,k);
	const double fl = aux->f(i,j,k);
	const double T00 = ( rhohW2 - p ) / (lapsel * lapsel);   // upper components 
	//	const double R00 = 4./3. * El * u0 * u0 + 2.0 * aux->f(i,j,k) * u0 / lapsel 
	//	  - El / (3.0 * lapsel * lapsel);
	const double rho_bar = ( 4.0*Wl*Wl - 1.0) * El / 3.0 + 2.0 * Wl * fl; 
	vect T0i_up;  // none of these are rescaled!
	vect T0i_down;  // also: index 0 is up in both cases...
	tensor Tij_up;
	// vect R0i_up;
	// vect R0i_down;
	vect j_bar_up;
	vect j_bar_down;
	// tensor Rij_up;
	tensor S_bar_up;
	for (int a = 0; a < 3; a++) {
	  T0i_up[a] = ( rhohW2 * ( v_up[a] - shift[a]/lapsel ) + p * shift[a] / lapsel ) / lapsel;
	  // check shift term in next line:
	  // R0i_up[a] = 4./3. * El * u0 * u_up[a] + aux->F_0(i,j,k) * u_up[a] + u0 * F_up[a] 
	  //   - El * shift[a] / (3.0 * lapsel * lapsel);  
	  T0i_down[a] = rhohW2 * v_down[a] / lapsel;
	  // R0i_down[a] = 4./3. * El * u0 * u_down[a] + aux->F_0(i,j,k) * u_down[a] + u0 * F_down[a];
	  j_bar_up[a] =   4./3.*El*Wl*Wl*v_up[a]   + Wl*f_up[a]   + fl*Wl*v_up[a];
	  j_bar_down[a] = 4./3.*El*Wl*Wl*v_down[a] + Wl*f_down[a] + fl*Wl*v_down[a];
	  for (int b = 0; b < 3; b++) {
	    Tij_up[a][b] = rhohW2 * (v_up[a] - shift[a]/lapsel) * (v_up[b] - shift[b]/lapsel)
	      + p * ( g_up[a][b]/e4p - shift[a]*shift[b]/(lapsel*lapsel) );
	    // Rij_up[a][b] = 4./3. * El * u_up[a]*u_up[b] + F_up[a]*u_up[b] + u_up[a]*F_up[b]
	      // CHECK shift terms...
	    //  + El * ( g_up[a][b]/e4p - shift[a]*shift[b]/(lapsel*lapsel) ) / 3.0; 
	    S_bar_up[a][b] = ( 4.*Wl*Wl*v_up[a]*v_up[b] + g_up[a][b]/e4p ) * El / 3.0 +
	      Wl*f_up[a]*v_up[b] + Wl*v_up[a]*f_up[b];
	  }
	}
	//
	// Assign derivatives of physical, spatial metric to rank-3 tensor D_gamma
	//    example: D_gamma[0][1][2] = D_r \gamma_{tp} 
	//                              = e^{4 \phi} ( D_r \bar \gamma_{tp} + 4 \bar \gamma_{tp} D_r \phi ) 
	//
	rank3tens D_gamma(e4p * ( c->Dr_e_rr(i,j,k) + 4.0 * g_conf[0][0] * s->phi.dr(i,j,k) ),
			  e4p * ( c->Dr_e_rt(i,j,k) + 4.0 * g_conf[0][1] * s->phi.dr(i,j,k) ),
			  e4p * ( c->Dr_e_rp(i,j,k) + 4.0 * g_conf[0][2] * s->phi.dr(i,j,k) ),
			  e4p * ( c->Dr_e_tt(i,j,k) + 4.0 * g_conf[1][1] * s->phi.dr(i,j,k) ),
			  e4p * ( c->Dr_e_tp(i,j,k) + 4.0 * g_conf[1][2] * s->phi.dr(i,j,k) ),
			  e4p * ( c->Dr_e_pp(i,j,k) + 4.0 * g_conf[2][2] * s->phi.dr(i,j,k) ),
			  e4p * ( c->Dt_e_rr(i,j,k) + 4.0 * g_conf[0][0] * s->phi.dtheta(i,j,k) ),
			  e4p * ( c->Dt_e_rt(i,j,k) + 4.0 * g_conf[0][1] * s->phi.dtheta(i,j,k) ),
			  e4p * ( c->Dt_e_rp(i,j,k) + 4.0 * g_conf[0][2] * s->phi.dtheta(i,j,k) ),
			  e4p * ( c->Dt_e_tt(i,j,k) + 4.0 * g_conf[1][1] * s->phi.dtheta(i,j,k) ),
			  e4p * ( c->Dt_e_tp(i,j,k) + 4.0 * g_conf[1][2] * s->phi.dtheta(i,j,k) ),
			  e4p * ( c->Dt_e_pp(i,j,k) + 4.0 * g_conf[2][2] * s->phi.dtheta(i,j,k) ),
			  e4p * ( c->Dp_e_rr(i,j,k) + 4.0 * g_conf[0][0] * s->phi.dphi(i,j,k) ),
			  e4p * ( c->Dp_e_rt(i,j,k) + 4.0 * g_conf[0][1] * s->phi.dphi(i,j,k) ),
			  e4p * ( c->Dp_e_rp(i,j,k) + 4.0 * g_conf[0][2] * s->phi.dphi(i,j,k) ),
			  e4p * ( c->Dp_e_tt(i,j,k) + 4.0 * g_conf[1][1] * s->phi.dphi(i,j,k) ),
			  e4p * ( c->Dp_e_tp(i,j,k) + 4.0 * g_conf[1][2] * s->phi.dphi(i,j,k) ),
			  e4p * ( c->Dp_e_pp(i,j,k) + 4.0 * g_conf[2][2] * s->phi.dphi(i,j,k) ) );
	//
	// assign components of physical extrinsic cuvature to tensor K 
	//
	tensor Aij(                   s->a_rr(i,j,k), 
		    rl              * s->a_rt(i,j,k), 
		    rl*sintheta     * s->a_rp(i,j,k),
		    rl*rl           * s->a_tt(i,j,k), 
		    rl*rl*sintheta  * s->a_tp(i,j,k), 
		    rl*rl*sin2theta * s->a_pp(i,j,k));
	const double Kl = s->K(i,j,k);
	tensor Kij;
	for (int a = 0; a < 3; a++) 
	  for (int b = 0; b < 3; b++) 
	    Kij[a][b] = e4p * ( Aij[a][b] + g_conf[a][b] * Kl / 3.0 );
	//
	// compute covariant derivative of shift with respect to background metric: \hat D_i \beta^j
	// result is physical scaling (i.e. not rescaled)
	//
	const double shift_r_r =   s->shift_r.dr(i,j,k);
	const double shift_r_t =   s->shift_r.dtheta(i,j,k) - s->shift_t(i,j,k);
	const double shift_r_p =   s->shift_r.dphi(i,j,k) - sintheta * s->shift_p(i,j,k);
	const double shift_t_r = ( s->shift_t.dr(i,j,k) ) / rl;
	const double shift_t_t = ( s->shift_t.dtheta(i,j,k) + s->shift_r(i,j,k) ) / rl;
	const double shift_t_p = ( s->shift_t.dphi(i,j,k) - costheta * s->shift_p(i,j,k) ) / rl;
	const double shift_p_r = ( s->shift_p.dr(i,j,k) ) / (rst);
	const double shift_p_t = ( s->shift_p.dtheta(i,j,k) ) / (rst);
	const double shift_p_p = ( s->shift_p.dphi(i,j,k) + sintheta * s->shift_r(i,j,k) 
				  + costheta * s->shift_t(i,j,k) ) / (rst);
	// CHECK that indices are right...
	tensor D_shift( shift_r_r, shift_r_t, shift_r_p,
			shift_t_r, shift_t_t, shift_t_p,
			shift_p_r, shift_p_t, shift_p_p );
	vect D_lapse(s->lapse.dr(i,j,k), s->lapse.dtheta(i,j,k), s->lapse.dphi(i,j,k));
	//
	// square root of spatial and spacetime metric determinant
	//
	const double det_3D = exp( 6.0 * s->phi(i,j,k) ) * sqrt( c->det(i,j,k) );
	const double det_4D = s->lapse(i,j,k) * det_3D;
	double det_4D_tau = det_4D;
	double det_3D_tau = det_3D;
	if (partial) { 
	  det_4D_tau *= rl*rl*sintheta;
	  det_3D_tau *= rl*rl*sintheta;
	}
	//
	// now compute source terms, and take care of rescaling
	//
	derivs->S_r[i][j][k] += - det_4D * T00 * lapsel * D_lapse[0];
	derivs->S_t[i][j][k] += - det_4D * T00 * lapsel * D_lapse[1] / rl;
	derivs->S_p[i][j][k] += - det_4D * T00 * lapsel * D_lapse[2] / rst;
	derivs->S_rad_r[i][j][k] += - det_3D * rho_bar * D_lapse[0];
	derivs->S_rad_t[i][j][k] += - det_3D * rho_bar * D_lapse[1] / rl;
	derivs->S_rad_p[i][j][k] += - det_3D * rho_bar * D_lapse[2] / rst;
	for (int a = 0; a < 3; a++) {
	  derivs->S_r[i][j][k] += det_4D * T0i_down[a] * D_shift[a][0];  // CHECK ordering of indices
	  derivs->S_t[i][j][k] += det_4D * T0i_down[a] * D_shift[a][1] / rl;
	  derivs->S_p[i][j][k] += det_4D * T0i_down[a] * D_shift[a][2] / rst;
	  derivs->tau[i][j][k] += det_4D_tau * ( - T00*shift[a] - T0i_up[a] ) * D_lapse[a];
	  derivs->tau_rad[i][j][k] += - det_3D_tau * j_bar_up[a] * D_lapse[a];
	  derivs->S_rad_r[i][j][k] += det_3D * j_bar_down[a] * D_shift[a][0];  // CHECK ordering of indices
	  derivs->S_rad_t[i][j][k] += det_3D * j_bar_down[a] * D_shift[a][1] / rl;
	  derivs->S_rad_p[i][j][k] += det_3D * j_bar_down[a] * D_shift[a][2] / rst;
	  for (int b = 0; b < 3; b++) {
	    const double term_ab = T00*shift[a]*shift[b]/2.0 + T0i_up[a]*shift[b] + Tij_up[a][b]/2.0;
	    // const double rad_ab = R00*shift[a]*shift[b]/2.0 + R0i_up[a]*shift[b] + Rij_up[a][b]/2.0;
	    derivs->S_r[i][j][k] += det_4D * term_ab * D_gamma[0][a][b];
	    derivs->S_t[i][j][k] += det_4D * term_ab * D_gamma[1][a][b] / rl;
	    derivs->S_p[i][j][k] += det_4D * term_ab * D_gamma[2][a][b] / rst;
	    derivs->tau[i][j][k] += det_4D_tau * 2.0 * term_ab * Kij[a][b];
	    derivs->tau_rad[i][j][k] += det_4D_tau * S_bar_up[a][b] * Kij[a][b];
	    derivs->S_rad_r[i][j][k] += det_4D * 0.5 * S_bar_up[a][b] * D_gamma[0][a][b];
	    derivs->S_rad_t[i][j][k] += det_4D * 0.5 * S_bar_up[a][b] * D_gamma[1][a][b] / rl;
	    derivs->S_rad_p[i][j][k] += det_4D * 0.5 * S_bar_up[a][b] * D_gamma[2][a][b] / rst;
	  } 
	}
	if (i == i_test && j == j_test && k == k_test) cout << " 3 " << setprecision(16) << derivs->tau(i,j,k) << " D = " << derivs->D(i,j,k) << endl;

	//	if (i == 2 && j == 2 && k == 2) cout << " 3 : tau-dot = " << derivs->tau[i][j][k] << endl;
	//
	//============================
	// add radiation force:
	//============================
	//
	// CHECK: KLUDGE!!!  Assume that T = m ( P / rho_0 ) and 
	//                   4 pi B = a_R T^4 = a_R m^4 ( P / rho_0)^4
	// Start by defining value for a_R m^4
	// value for flat-space test 1:
	//	const double am4 = 1.2345e10;
	// value for flat-space test 4:
	const double am4 = 1.38889e8;
	// value for OS test (with R_i = 5 M), E_i / rho_0_i = 1.e-4, P_i/rho_0_i = 1.e-5
	// const double am4 = 1.9099e13;
	//	const double am4 = 1.958e17;
	// const double am4 = 0.0;
	const double PoverRho = p / rho_0l;
	// const double FourPiB = am4 * PoverRho * PoverRho * PoverRho * PoverRho;
	const double FourPiB = El;
	// now define value for absorbtion and scattering opacities
	// CHECK: again, KLUDGE - should be provided by EOS
	// computed from tau^a = kappa^a \rho R = 50 initially (for OS)
	double kappa_a = 5238.3;
	//	double kappa_a = 0.24;
	// double kappa_a = 0.08;
	// if (rho_0l < rho_thr) kappa_a = 0.0;
	// const double kappa_a = 0.0;
	const double kappa_s = 0.0;
	const double gl = rho_0l * ( kappa_a * (El - FourPiB) * Wl + 
					 (kappa_a + kappa_s) * fl);
	const double g_r_down = rho_0l * ( kappa_a * (El - FourPiB) * Wl * v_down[0] + 
					   (kappa_a + kappa_s) * f_down[0] );
	const double g_t_down = rho_0l * ( kappa_a * (El - FourPiB) * Wl * v_down[1] + 
					   (kappa_a + kappa_s) * f_down[1] ) / rl;
	const double g_p_down = rho_0l * ( kappa_a * (El - FourPiB) * Wl * v_down[2] + 
					   (kappa_a + kappa_s) * f_down[2] ) / rst;
	derivs->tau[i][j][k] += det_4D_tau * gl;
	derivs->S_r[i][j][k] += det_4D * g_r_down;
	derivs->S_t[i][j][k] += det_4D * g_t_down;
	derivs->S_p[i][j][k] += det_4D * g_p_down;
	derivs->tau_rad[i][j][k] -= det_4D_tau * gl;
	derivs->S_rad_r[i][j][k] -= det_4D * g_r_down;
	derivs->S_rad_t[i][j][k] -= det_4D * g_t_down;
	derivs->S_rad_p[i][j][k] -= det_4D * g_p_down;
	//	if (i == 2 && j == 2 && k == 2) cout << " 4 : tau-dot = " << derivs->tau[i][j][k] << endl;
	if (i == i_test && j == j_test && k == k_test) { 
	  cout << " 4 " << setprecision(16) << derivs->tau(i,j,k) << " D = " << derivs->D(i,j,k) << endl;
	  cout << " gl = " << gl << " rho_0l = " << rho_0l << " p = " << p 
	       << " FourPiB = " << FourPiB << endl;
	}
      }
    }
  }
  return error;
};
//
//================================================
//
// Compute ADM sources
//
// NOTE: will compute unrescaled sources, as expected in field equations
//================================================
//
void RadHydro::ADM_Sources(state *s, curvature *c) {
  const double fourthirds = 4.0 / 3.0;
  for (int i = N_g; i < N_r - N_g; i++) {     
    const double rl = grid->r(i);
    for (int j = N_g; j < N_t - N_g; j++) {
      const double sinthetal = grid->sintheta(j);
      const double rst = rl * sinthetal;
      for (int k = N_g; k < N_p - N_g; k++) {
	const double e4p = exp(4.0 * s->phi(i,j,k));
	tensor gup(c->gup_rr(i,j,k), c->gup_rt(i,j,k), c->gup_rp(i,j,k),
		   c->gup_tt(i,j,k), c->gup_tp(i,j,k), c->gup_pp(i,j,k));
	tensor g_conf = gup.inverse();
	//
  	// adopt "Valencia" definition of v^i, use unrescaled velocity now
	//
	vect v_up(aux->v_r(i,j,k), aux->v_t(i,j,k)/rl, aux->v_p(i,j,k)/rst);
	vect f_up(aux->f_r(i,j,k), aux->f_t(i,j,k)/rl, aux->f_p(i,j,k)/rst);
  	double vdotv = 0.0;
  	for (int a = 0; a < 3; a++)
  	  for (int b = 0; b < 3; b++) 
  	    vdotv += e4p * g_conf[a][b] * v_up[a] * v_up[b];
  	//
  	// compute Lorentz factor W
  	//
	double Wl = aux->W(i,j,k);
	double W2 = Wl * Wl;
 	//
  	// compute u_i (lower indices!)
  	//
  	vect u;
	vect v_low;
	vect f_low;
  	for (int a = 0; a < 3; a++) {
  	  u[a] = 0.0;
	  v_low[a] = 0.0;
	  f_low[a] = 0.0;
  	  for (int b = 0; b < 3; b++) {
  	    u[a] += Wl * e4p * g_conf[a][b] * v_up[b];
	    v_low[a] += e4p * g_conf[a][b] * v_up[b];
	    f_low[a] += e4p * g_conf[a][b] * f_up[b];
	  }
  	}
  	//
  	// find some thermodynamic quantities
  	//	
  	const double rho_rest = aux->rho_0(i,j,k);
	const double epsilon = aux->eps(i,j,k);
  	const double enth = eos->h(rho_rest,epsilon);  
  	const double P = eos->P(rho_rest,epsilon); 
	const double E = aux->E(i,j,k);
	const double f = aux->f(i,j,k);
  	//
  	// Now compute ADM density...
  	// 
	adm_sources->rho_ADM[i][j][k] = rho_rest * enth * W2 - P + fourthirds * W2 * E - E/3.0 + 2.0 * Wl * f;
  	//
  	// ... fluxes ...
  	//
	adm_sources->S_r[i][j][k] = rho_rest*enth*Wl*u[0] + fourthirds*E*W2*v_low[0] + Wl*f_low[0] + f*Wl*v_low[0];
  	adm_sources->S_t[i][j][k] = rho_rest*enth*Wl*u[1] + fourthirds*E*W2*v_low[1] + Wl*f_low[1] + f*Wl*v_low[1];
  	adm_sources->S_p[i][j][k] = rho_rest*enth*Wl*u[2] + fourthirds*E*W2*v_low[2] + Wl*f_low[2] + f*Wl*v_low[2];
  	//
  	// ... stresses: hydro first...
  	//
  	adm_sources->S_rr[i][j][k] = P*e4p*g_conf[0][0] + rho_rest*enth*u[0]*u[0];
  	adm_sources->S_rt[i][j][k] = P*e4p*g_conf[0][1] + rho_rest*enth*u[0]*u[1];
  	adm_sources->S_rp[i][j][k] = P*e4p*g_conf[0][2] + rho_rest*enth*u[0]*u[2];
  	adm_sources->S_tt[i][j][k] = P*e4p*g_conf[1][1] + rho_rest*enth*u[1]*u[1];
  	adm_sources->S_tp[i][j][k] = P*e4p*g_conf[1][2] + rho_rest*enth*u[1]*u[2];
  	adm_sources->S_pp[i][j][k] = P*e4p*g_conf[2][2] + rho_rest*enth*u[2]*u[2];
	//
	// ... then add radiation stuff ...
	//
  	adm_sources->S_rr[i][j][k] += fourthirds*W2*E*v_low[0]*v_low[0] + E*g_conf[0][0]/3.0 + Wl*(v_low[0]*f_low[0] + f_low[0]*v_low[0]);
  	adm_sources->S_rt[i][j][k] += fourthirds*W2*E*v_low[0]*v_low[1] + E*g_conf[0][1]/3.0 + Wl*(v_low[0]*f_low[1] + f_low[0]*v_low[1]);
  	adm_sources->S_rp[i][j][k] += fourthirds*W2*E*v_low[0]*v_low[2] + E*g_conf[0][2]/3.0 + Wl*(v_low[0]*f_low[2] + f_low[0]*v_low[2]);
  	adm_sources->S_tt[i][j][k] += fourthirds*W2*E*v_low[1]*v_low[1] + E*g_conf[1][1]/3.0 + Wl*(v_low[1]*f_low[1] + f_low[1]*v_low[1]);
  	adm_sources->S_tp[i][j][k] += fourthirds*W2*E*v_low[1]*v_low[2] + E*g_conf[1][2]/3.0 + Wl*(v_low[1]*f_low[2] + f_low[1]*v_low[2]);
  	adm_sources->S_pp[i][j][k] += fourthirds*W2*E*v_low[2]*v_low[2] + E*g_conf[2][2]/3.0 + Wl*(v_low[2]*f_low[2] + f_low[2]*v_low[2]);
	//
  	// ... and trace of stress:
  	//
	// CHECK radiation contribution!
  	adm_sources->trace_S[i][j][k] = 3.0 * P + rho_rest*enth*(Wl*Wl - 1.0) + E*(W2 - 1./3.) + 2.*Wl*f;
      } 
    }
  }
};
//
//================================================
// 
// Rest Mass
//
//================================================
//
double RadHydro::Rest_Mass(radhydro_state *m) {
  double rest_mass = 0.0;
  if (N_t == 6) {
    // assuming spherical symmetry...
    const double PI = acos(-1.0);
    for (int i = N_g; i < N_r - N_g; i++) {
      const double rl = grid->r(i);
      const double dr = grid->delta_r(i);
      rest_mass += m->D(i,2,2)*rl*rl*4.0*PI*dr;
    } 
  } else {
    for (int i = N_g; i < N_r - N_g; i++ ) {
      const double rl = grid->r(i);
      const double dr = grid->delta_r(i);
      for (int j = N_g; j < N_t - N_g; j++ ) {
	const double dtheta = grid->delta_theta(j);
	const double sintheta = grid->sintheta(j);
	for (int k = N_g; k < N_p - N_g; k++ ) {
	  const double dphi = grid->delta_phi(k);
	  if (partial)
	    rest_mass += m->D(i,j,k) * dr*dtheta*dphi;
	  else 
	    rest_mass += m->D(i,j,k) * rl*rl*sintheta*dr*dtheta*dphi;
	}
      }
    }
  }
#ifdef EQSYMMETRY
  return 2.0 * rest_mass;
#else
  return rest_mass;
#endif
}
//
//================================================
// 
// Test recovery
//
//================================================
//
int RadHydro::Test_Recovery(state *s, curvature *c) {
  int i_test = N_g;
  int j_test = N_g;
  int k_test = N_g;
  aux->f_r[i_test][j_test][k_test] = 3.e-5;
  cout << " RADHYDRO: rho_0 initially = " << aux->rho_0(i_test,j_test,k_test) 
       << ", v_r = " << aux->v_r(i_test,j_test,k_test) 
       << ", f_r = " << aux->f_r(i_test,j_test,k_test) 
       << endl;
  cout << " RADHYDRO: Computing conserveds from primitives..." << endl;
  Compute_Conserved_Vars(last, s, c);
  // const double Dl =  4.04157e-11;
  // const double factor = 0.01;
  // const double taul = 4.69192e-12;
  // const double Srl = - factor * 3.3e-9;
  // const double pl_max = 1.e-11;
  // double dfdp = 0.0;
  // for (double pl = - pl_max; pl <= pl_max; pl += pl_max/20) {
  //   cout << " pl = " << setw(16) << pl << " f = " << f(pl, taul, Dl, 1.e-21, dfdp) << " dfdp = " << dfdp << endl;
  // }
  for (int i = 0; i < N_r; i++ )
    for (int j = 0; j < N_t; j++ )
      for (int k = 0; k < N_p; k++ ) {
	aux->rho_0[i][j][k] = 0.0;
	aux->v_r[i][j][k] = 0.0;
	aux->f_r[i][j][k] = 0.0;
      }
  // 	const double det_3D = exp(6.0*s->phi(i,j,k))*sqrt( c->det(i,j,k) );
  // 	double r_l = grid->r(i);
  // 	double sintheta_l = grid->sintheta(j);
  // 	last->D[i][j][k] = Dl * det_3D;
  // 	last->tau[i][j][k] = taul * det_3D;
  // 	if (partial) {
  // 	  last->D[i][j][k] *= r_l * r_l * sintheta_l;
  // 	  last->tau[i][j][k] *= r_l * r_l * sintheta_l;
  // 	}	
  // 	last->S_r[i][j][k] = - Srl * det_3D;
  // 	last->S_t[i][j][k] = 0.0;
  //     }
  int error = Recovery(last, s, c);
  cout << " RADHYDRO: recovered primitives with " << error << " errors " << endl;
  cout << " RADHYDRO: rho_0 after recovery = " << aux->rho_0(i_test,j_test,k_test) << endl;
  cout << " RADHYDRO: v_r after recovery = " << aux->v_r(i_test,j_test,k_test) << endl;
  cout << " RADHYDRO: f_r after recovery = " << aux->f_r(i_test,j_test,k_test) << endl;
  return error;
};

