// Tell emacs that this is -*-c++-*- mode
//====================================================
//
// Manages allocation of grid and regridding
//
// Uses uniform grids in x and y; r and theta are given as functions
// of x and y.  These functions are specified in r_fct and theta_fct.
//    
//====================================================
//
#include "nr3.h"
//
#ifndef GRID_H
#define GRID_H
//
class Grid {
private:
  int N_theta, N_phi, N_r;
  Doub r_max_init, r_max_fin;
  int regrids, regrid_counter;
  int regrid_type, grid_type;
  Doub tau_star;    // for self-similar regridding
  Doub last_selfsim_ratio;
  bool selfsimregrid;
  Doub courant; 
  Doub cutoff;
  Doub s_param, t_param, t_amp, theta_param;
  double r_focus, x_focus;
  bool tracking;
  Doub r_max_current, r_max_old, r_max_new;
  int N_ghost;
  Doub r_max_factor;
  bool printed_warning;
  Doub PI;
  Doub dx, dy, dphi;
  VecDoub x_v, y_v;
  VecDoub r_v, theta_v, phi_v;
  VecDoub dxdr_v, dydtheta_v;
  VecDoub ddxdr_v, ddydtheta_v;
  VecDoub costheta_v, sintheta_v;
  Doub  eps;  // helps with conversion from coordinates to indices works
  int i_ind_sav;  // needed in "hunt" routine for i_ind
  //===================================================
  // Constructor
  //===================================================
public:
  Grid() : regrid_counter(0), N_ghost(4), printed_warning(false), eps(1.e-12)
  {
#ifdef EIGHTHORDER
    N_ghost = 4;
    cout << " GRID: Using N_ghost = " << N_ghost << " for eighth-oder spatial differencing." << endl;
#elif SIXTHORDER
    N_ghost = 4;
    cout << " GRID: Using N_ghost = " << N_ghost << " for sixth-oder spatial differencing." << endl;
#else
    N_ghost = 3;
    cout << " GRID: Using N_ghost = " << N_ghost << " for fourth-oder spatial differencing." << endl;
#endif 
    // 
    // read input from file "Grid_Input"
    // 
    PI = acos(-1.0);
    ifstream infile;
    infile.open("Grid_Input");
    if (!infile) {
      cerr << " GRID: Can't open input file Grid_Input - this is bad! " << endl;
    }
    char buf[500], c;
    infile.get(buf,500,'='); infile.get(c); infile >> N_r;
    infile.get(buf,500,'='); infile.get(c); infile >> N_theta;
    infile.get(buf,500,'='); infile.get(c); infile >> N_phi;
    infile.get(buf,500,'='); infile.get(c); infile >> r_max_init;
    infile.get(buf,500,'='); infile.get(c); infile >> r_max_fin;
    infile.get(buf,500,'='); infile.get(c); infile >> courant;    
    infile.get(buf,500,'='); infile.get(c); infile >> regrid_type;    
    infile.get(buf,500,'='); infile.get(c); infile >> regrids;
    infile.get(buf,500,'='); infile.get(c); infile >> cutoff;
    infile.get(buf,500,'='); infile.get(c); infile >> tau_star;
    infile.get(buf,500,'='); infile.get(c); infile >> grid_type;
    infile.get(buf,500,'='); infile.get(c); infile >> s_param;
    infile.get(buf,500,'='); infile.get(c); infile >> theta_param;
    infile.get(buf,500,'='); infile.get(c); infile >> t_param;
    infile.get(buf,500,'='); infile.get(c); infile >> t_amp;
    infile.get(buf,500,'='); infile.get(c); infile >> tracking;
    infile.get(buf,500,'='); infile.get(c); infile >> r_focus;
    infile.get(buf,500,'='); infile.get(c); infile >> x_focus;
     //
    if(infile.eof()) {
      cerr << " GRID: Error reading input file Grid_Input " << endl;
    }
    infile.close();
    //
#ifdef AXISYMMETRY
    cout << " GRID: setting N_phi = 1 in Axisymmetry " << endl;
    N_phi = 1;
#endif
    //
    // echo input...
    // 
    cout << " GRID:  N_r      = " << N_r << endl;
    cout << " GRID:  N_theta  = " << N_theta << endl;
    cout << " GRID:  N_phi     = " << N_phi << endl;
    cout << " GRID:  r_max_init = " << r_max_init 
	 << ", r_max_fin = " << r_max_fin << endl;
    cout << " GRID:  using Courant factor of " << courant << endl;
    cout << " GRID:  using " << N_ghost << " ghost zones " << endl;
    //
    if (regrids == 0) {
      cout << " GRID: no regridding. " << endl;
    } else {
      if (regrid_type == 0) {
        selfsimregrid = false;
        cout << " GRID: using accuracy regridding, maximally " << regrids 
	          << " regriddings with cutoff = " << cutoff << endl;
      } else if (regrid_type == 1) {
        selfsimregrid = true;
        cout << " GRID: using self-similar regridding with tau_star = "
            << tau_star << " and cutoff = " << cutoff << endl;
      }
    }
    cout << " GRID: Using grid parameters s_param = " << s_param 
	 << " theta_param = " << theta_param << endl;
   if (grid_type == 0) {
    cout << " GRID: Using tanh grid parameters t_param = " << t_param 
        << " t_amp = " << t_amp << endl;
   } else if (grid_type == 1) {
    cout << " GRID: Using offset grid parameters r_focus = " << r_focus
        << " x_focus = " << x_focus << endl;
   }

    // 
    r_v = VecDoub(N_r + 2*N_ghost);
    x_v = VecDoub(N_r + 2*N_ghost);
    dxdr_v = VecDoub(N_r + 2*N_ghost);
    ddxdr_v = VecDoub(N_r + 2*N_ghost);
    theta_v = VecDoub(N_theta + 2*N_ghost);
    y_v = VecDoub(N_theta + 2*N_ghost);
    dydtheta_v = VecDoub(N_theta + 2*N_ghost);
    ddydtheta_v = VecDoub(N_theta + 2*N_ghost);
    costheta_v = VecDoub(N_theta + 2*N_ghost);
    sintheta_v = VecDoub(N_theta + 2*N_ghost);
    phi_v = VecDoub(N_phi + 2*N_ghost);
    //
    r_max_current = r_max_new = r_max_old = r_max_init;
    //
    // change r_max in equal factor
    //
    if (regrid_type == 0) {
      r_max_factor = pow(r_max_fin/r_max_init,1.0/regrids);
    }
    else if (regrid_type == 1) { 
      r_max_factor = 1. / cutoff;
      last_selfsim_ratio = r_max_current / tau_star;
    } 
    //
    // needed in i_ind (for "hunt" routine):
    // 
    //    dj = MIN(1,(int)pow((Doub)(N_r + 2*N_ghost),0.25));
  }
  ~Grid() {
    cout << " GRID: Regridded " << regrid_counter << " times. " << endl; 
  }
  //=================================================
  // Return grid specifications
  //=================================================
  int N_r_int() { return N_r; }
  int N_r_tot() { return N_r + 2 * N_ghost; }
  int N_theta_int() { return N_theta; }
  int N_theta_tot() { return N_theta + 2 * N_ghost; }
  int N_phi_int() { return N_phi; }
  int N_phi_tot() { return N_phi + 2 * N_ghost; }
  int N_ghosts() { return N_ghost; }
  Doub r_max() { return r_max_current; }
  Doub r_max_final() {return r_max_fin; }
  Doub r_parameter() { return s_param; }
  Doub theta_parameter() { return theta_param; }
  Doub courant_factor() { return courant; }
  Doub delta_x() { return dx; }
  Doub delta_y() { return dy; }
  Doub delta_phi() { return dphi; }
  Doub delta_r(int i) { return dx / dxdr_v[i]; }
  Doub delta_theta(int j) { return dy / dydtheta_v[j]; }
  Doub delta_phi(int k) { return dphi; }
  VecDoub *r() { return &r_v; }
  VecDoub *theta() { return &theta_v; }
  VecDoub *phi() { return &phi_v; }
  VecDoub *dxdr() { return &dxdr_v; }
  VecDoub *ddxdr() { return &ddxdr_v; }
  VecDoub *dydtheta() { return &dydtheta_v; }
  VecDoub *ddydtheta() { return &ddydtheta_v; }
  VecDoub *costheta() { return &costheta_v; }
  VecDoub *sintheta() { return &sintheta_v; }
  Doub r(int i) { return r_v[i]; }
  Doub theta(int j) { return theta_v[j]; }
  Doub phi(int k) { return phi_v[k]; }
  Doub dxdr(int i) { return dxdr_v[i]; }
  Doub ddxdr(int i) { return ddxdr_v[i]; }
  Doub sintheta(int j) { return sintheta_v[j]; }
  Doub costheta(int j) { return costheta_v[j]; }
  Doub dydtheta(int j) { return dydtheta_v[j]; } 
  Doub ddydtheta(int j) { return ddydtheta_v[j]; } 
  Doub RegridCriterion(Doub tau_c) {
    return r_max_current / (tau_star - tau_c) / last_selfsim_ratio;
  }
  int Regrid_Type() { return regrid_type; }
  
  //=================================================
  // Set up grid
  //=================================================
  int Setup_Grid(VecDoub & r, VecDoub & r2, VecDoub & theta, 
		 VecDoub & sintheta, VecDoub & sin2theta, 
		 VecDoub & costheta, VecDoub & phi) {
    Setup_Radial_Grid(r, r2);
    Setup_Angular_Grid(theta, sintheta, sin2theta, costheta, phi);
    return 1;
  }
  //=================================================
  // regrid
  //=================================================
  bool TimeToRegrid(Doub criterion) {
    if (criterion > cutoff) {
      if (regrid_counter < regrids) 
    	return true;
      else {
    	if (!printed_warning) {
    	  printed_warning = true;
    	  cout << " GRID: exceeded maximum number of regrids. " << endl;	
    	}
    	return false;
      }
    } else
      return false;
  }
  int Regrid(VecDoub & r_new) {
    //
    // function returns vector with new radial gridpoints, but
    // doesn't do anything else yet -- need to complete regridding by
    // calling Setup_Radial_Grid.
    //
    regrid_counter++;
    //
    // kind of a hack: want to compute new r_max, use it temporarily in
    // r_fcct, but then want to restore old r_max so that it can be used
    // in i_ind during regridding (where old grid is needed)...
    r_max_old = r_max_current;
    if (regrid_type == 2) {
      r_max_new = r_max_current - cutoff;
    }
    else {
      r_max_new = r_max_current * r_max_factor;
    }
    // ... therefore temporarily set r_max_current to r_max_new...
    r_max_current = r_max_new;
    cout << " GRID: regridding with r_max = " << r_max_current << endl;
    Doub temp1, temp2;
    for (int i = 0; i < N_r_tot(); i++) {
      r_new[i] = r_fct(x_v[i], temp1, temp2);
    }
    // ... but then restore it
    r_max_current = r_max_old;
    return regrid_counter;
  }
  //=================================================
  // Setup radial grid
  //=================================================
  int Setup_Radial_Grid(VecDoub & r, VecDoub & r2) {
    //
    // now set up radial grid for good, to complete regridding,
    // and set up r_max_current to r_max_new
    // 
    r_max_current = r_max_new;
    cout << " GRID: Setting up radial grid with N_r = " << N_r 
	 << ", r_max = " << r_max_current << ", and A = " 
	 << s_param << endl;
    // 
    // build radial grid from uniform grid in x...
    //
    dx = 1.0 / Doub(N_r);  // N_r = N_r_int!
    for (int i = 0; i < N_r_tot(); i++) {
      x_v[i] = (i - (N_ghost - 0.5))*dx;
      r_v[i] = r_fct(x_v[i], dxdr_v[i], ddxdr_v[i]);
      r[i] = r_v[i];
      r2[i] = r_v[i] * r_v[i];
      // cout << i << " r " << r[i] << " dxdr " 
      // 	   << dxdr_v[i] << "ddxdr " << ddxdr_v[i] << endl;
    }
    cout << " GRID: dr across origin:     " << dx/dxdr_v[N_ghost] << endl;
    cout << " GRID: dr at outer boundary: " << dx/dxdr_v[N_r-2] << endl;
    return 1;
  }
  //=================================================
  // Setup angular grid
  //=================================================
  int Setup_Angular_Grid(VecDoub & theta, VecDoub & sintheta, 
			 VecDoub & sin2theta, VecDoub & costheta,
			 VecDoub & phi) {
    cout << " GRID: Setting up angular grid with N_theta = " << N_theta 
	 << ", N_phi = " << N_phi << " and B = " << theta_param << endl;
    //
    // build theta grid from uniform grid in y
    // NOTE: in equatorial symmetry y covers range [0, 0.5] only
    //
#ifdef EQSYMMETRY
    dy = 0.5 / Doub(N_theta);  // N_theta = N_theta_int!
#else
    dy = 1.0 / Doub(N_theta);
#endif
    for (int j = 0; j < N_theta_tot(); j++) {
      y_v[j] = (j - (N_ghost - 0.5))*dy;
      theta_v[j] = theta_fct(y_v[j], dydtheta_v[j], ddydtheta_v[j]);
      costheta_v[j] = cos(theta_v[j]);
      //      sintheta_v[j] = sqrt(1.0 - costheta_v[j]*costheta_v[j]);
      sintheta_v[j] = sin(theta_v[j]);
      theta[j] = theta_v[j];
      sintheta[j] = sintheta_v[j];
      sin2theta[j] = sintheta_v[j] * sintheta_v[j];
      costheta[j] = costheta_v[j];
    }
#ifndef EQSYMMETRY
    //
    // if there is no equatorial symmetry make sin's and cos's
    // symmetric about equator
    //
    const double N_theta_middle = N_theta_tot()/2;
    for (int j = 0; j < N_theta_middle; j++) {
      sintheta[N_theta_middle + j] =   sintheta[N_theta_middle - 1 - j];
      costheta[N_theta_middle + j] = - costheta[N_theta_middle - 1 - j];
    }
    // for (int j = 0; j < N_theta_tot(); j++) 
    //   cout << " j = " << j << " theta = " << setprecision(16)
    // 	   << theta[j] << " sin = " << sintheta[j]
    // 	   << " cos = " << costheta[j] << endl;
#endif
    //
    // phi
    //
    dphi = 2.0*PI/Doub(N_phi);
    for (int k = 0; k < N_phi_tot(); k++) {
      phi_v[k] = (k - (N_ghost - 0.5))*dphi;
      phi[k] = phi_v[k];
    }
    return 1;
  }
  //================================================
  // Compute r as function of x -- specify function here!
  // Note:
  //    dx / d r = 1 / r'
  //    d^2 x/ d r^2 = - r'' / ( r' )^3
  //
  // Here: r = r_max sinh(r_param x) / sinh(r_param)
  //
  //================================================
  // Doub r_fct(Doub x, Doub & x_prime, Doub & x_dprime) {
  //   Doub r = 0.0;
  //   if (s_param > 0.) {
  //     const Doub sinhA = sinh(s_param);
  //     const Doub sinhAx = sinh(s_param * x);
  //     const Doub coshAx = cosh(s_param * x);
  //     r = r_max_current * sinhAx / sinhA;
  //     const Doub r_prime = r_max_current * s_param * coshAx / sinhA;
  //     const Doub r_dprime = s_param * s_param * r;
  //     x_prime = 1.0 / r_prime;
  //     x_dprime = - r_dprime / (r_prime * r_prime * r_prime); 
  //   } else {
  //     r = r_max_current * x;
  //     x_prime = 1.0 / r_max_current;
  //     x_dprime = 0.0;
  //   }
  //   return r;
  // }  
  //================================================
  // Compute r as function of x -- specify function here!
  // Note:
  //    dx / d r = 1 / r'
  //    d^2 x/ d r^2 = - r'' / ( r' )^3
  //
  // Here: r = r_max / (1 + A) ( A tanh(t_param x) / tanh(t_param)
  //                             + sinh(s_param x) / sinh(s_param) )
  //
  //================================================
  Doub r_fct_tanh(Doub x, Doub & x_prime, Doub & x_dprime) {
    Doub r = 0.0;
    if (s_param > 0.) {
      const Doub A = t_amp;
      const Doub amp = r_max_current / ( 1.0 + A);
      const Doub sinhS = sinh(s_param);
      const Doub sinhSx = sinh(s_param * x);
      const Doub coshSx = cosh(s_param * x);
      //
      const Doub tanhT = tanh(t_param);
      const Doub tanhTx = tanh(t_param * x);
      const Doub coshTx = cosh(t_param * x);
      const Doub cosh2Tx = coshTx * coshTx;
      // 
      r = amp * ( A * tanhTx / tanhT + sinhSx / sinhS);
      //
      const Doub r_prime = amp * (A * t_param / (tanhT * cosh2Tx) +
				  s_param * coshSx / sinhS );
      const Doub r_dprime = amp * (-2.0 * A * t_param * t_param / cosh2Tx *
				   tanhTx / tanhT +
				   s_param * s_param * sinhSx / sinhS );
      x_prime = 1.0 / r_prime;
      x_dprime = - r_dprime / (r_prime * r_prime * r_prime); 
    } else {
      r = r_max_current * x;
      x_prime = 1.0 / r_max_current;
      x_dprime = 0.0;
    }
    return r;
  }
  //================================================
  // Compute r as function of x -- specify function here!
  // Note:
  //    dx / d r = 1 / r'
  //    d^2 x/ d r^2 = - r'' / ( r' )^3
  //
  // Here: r = r_max (sinh(sx-sa)+sinh(sa))/(sinh(s-sa)+sinh(sa))
  //
  //  s_param -> flatness
  //  t_param -> r location of flatness 
  //  t_amp   -> x location of flatness (used if t_param=0)
  //
  //================================================
  Doub r_fct_offset(Doub x, Doub & x_prime, Doub & x_dprime) {
    double r = 0.0;

    if (s_param > 0.) {
      int sign = 1;
      if (x < 0.) {
        x = -x;
        sign = -1;
      }
      if (r_focus > 0.) {
        const double arg = (r_max_current/r_focus - 1.) * (1./sinh(s_param)) + (1./tanh(s_param));
        const double arccosh = log((arg + 1.) / (arg - 1.)) / 2.;
        x_focus = (1./s_param)*arccosh;
      }

      const double sinhSxSa = sinh(s_param * (x -  x_focus));
      const double sinhSa = sinh(s_param * x_focus);
      const double sinhSSa = sinh(s_param * (1. - x_focus));
      const double coshSxSa = cosh(s_param * (x -  x_focus));

      r = sign * r_max_current * (sinhSxSa + sinhSa)/(sinhSSa + sinhSa);

      const double r_prime = sign * r_max_current * s_param * 
          coshSxSa / (sinhSSa + sinhSa);
      const double r_dprime = sign * r_max_current * s_param * s_param * 
          sinhSxSa / (sinhSSa + sinhSa); 
      x_prime = 1.0 / r_prime;
      x_dprime = - r_dprime / (r_prime * r_prime * r_prime); 
    } else {
      r = r_max_current * x;
      x_prime = 1.0 / r_max_current;
      x_dprime = 0.0;
    }
    return r;
  }  

  double r_fct(Doub x, Doub & x_prime, Doub & x_dprime) {
    if (grid_type == 0)
      return r_fct_tanh(x, x_prime, x_dprime);
    else if (grid_type == 1)
      return r_fct_offset(x, x_prime, x_dprime);
    else {
      cout << " GRID: Unknown grid type" << endl;
      exit(1);
    }
  }
  //================================================
  // Compute theta as function of y -- specify function here!
  // Note:
  //    dy / d theta = 1 / theta'
  //    d^2 y/ d theta^2 = - theta'' / ( theta' )^3
  //================================================
  Doub theta_fct(Doub y, Doub & y_prime, Doub & y_dprime) {
    Doub siny = sin(2.0 * PI * y);
    Doub cosy = cos(2.0 * PI * y);
    Doub theta = PI * y + theta_param*siny;
    Doub dthetady = PI + 2.0 * PI * theta_param * cosy;
    y_prime = 1.0 / dthetady;
    y_dprime = 4.0 * PI * PI * theta_param * siny / 
      ( dthetady * dthetady * dthetady );
    if (!(y_prime > 0.)) {
      cerr << " Poor choice for theta(y) -- abort now!! " << endl;
      exit(0);
    }
    return theta;
  }
  //
  //==================================================
  // finally routines that find indices for coordinates, eg., return
  //    index i with r[i] <= rl < r[i+1]
  //==================================================
  //
public:
  // int i_ind(const double rl) {
  //   if (s_param > 0.0) {
  //     // first find x
  //     const double temp = sinh(s_param) * rl / r_max_current;
  //     const double xl = log(temp + sqrt(temp*temp + 1.0)) / s_param;
  //     const double x_ind = xl + (N_ghost - 0.5)*dx;
  //     return int(x_ind/dx + eps);
  //   } else {
  //     const double xl = rl / r_max_current;
  //     const double x_ind = xl + (N_ghost - 0.5)*dx;
  //     return int(x_ind/dx + eps);
  //   }
  // }

  int i_ind(const Doub rl)
  {
	int jl=i_ind_sav, jm, ju, inc=1;
	int mm = 2;
	int n = N_r_tot();
	// if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
	Bool ascnd=true;
	if (rl >= r_v[n-1]) return n-1;
	if (jl < 0 || jl > n-1) {
		jl=0;
		ju=n-1;
	} else {
		if (rl >= r_v[jl] == ascnd) {
			for (;;) {
				ju = jl + inc;
				if (ju >= n-1) { ju = n-1; break;}
				else if (rl < r_v[ju] == ascnd) break;
				else {
					jl = ju;
					inc += inc;
				}
			}
		} else {
			ju = jl;
			for (;;) {
				jl = jl - inc;
				if (jl <= 0) { jl = 0; break;}
				else if (rl >= r_v[jl] == ascnd) break;
				else {
					ju = jl;
					inc += inc;
				}
			}
		}
	}
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1;
		if (rl >= r_v[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	i_ind_sav = jl;
	return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}

  int j_ind(const double thetal) {
    if (thetal <= theta_v[0]) return 0;
    int max_ind = N_theta + 2*N_ghost - 1;
    if (thetal >= theta_v[max_ind]) return max_ind;
    int jlo = 0;
    while ( (theta_v[jlo] <= thetal) && (jlo < max_ind) ) {
      jlo += 1;
    }
    return jlo - 1;
  }
  int k_ind(double phil) { return int(phil/dphi + (N_ghost - 0.5) + eps); }
};

#endif /* GRID_H */
