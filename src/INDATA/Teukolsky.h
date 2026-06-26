// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Linear Wave Initial Data (Teukolsky wave)
// (see pages 320 ff in Numerical Relativity)
//================================================
//
class LinWave : public InData {
private:
  double Amp;   
  double lambda;
  double r0;
  int l;
  int m;
  double r_exp;
  double term0, A_term2, B_term2, C_term2;  // terms for expansions of 
  double A_term4, B_term4, C_term4, expo;   // A, B, and C
  double delta_r; // for numerical derivatives...
public:
  //================================================
  // Constructor
  //================================================
  LinWave(char * indata_input, Grid * grid_i, Cosmology *cosmology) :
    InData(grid_i, cosmology) {
    indata_type = linwave;
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
    infile.get(buf,100,'='); infile.get(c); infile >> Amp;
    infile.get(buf,100,'='); infile.get(c); infile >> lambda;
    infile.get(buf,100,'='); infile.get(c); infile >> r0;
    infile.get(buf,100,'='); infile.get(c); infile >> l;
    infile.get(buf,100,'='); infile.get(c); infile >> m;
    infile.get(buf,100,'='); infile.get(c); infile >> r_exp;
    cout << " Will set up l = " << l << ", m = " << m << " Teukolsky wave with parameters " << endl;
    cout << "        Amplitude = " << Amp << endl;
    cout << "           lambda = " << lambda << endl;
    cout << "               r0 = " << r0 << endl;
    cout << "  Using Expansion inside r_exp = " << r_exp << endl;
    cout << "===================================================" << endl;
    if (m != 0) {
      cout << " INDATA/LinWave.h currently implement for m = 0 only!!! " 
	   << endl;
      exit(0);
    }
    //
    // compute expansion coefficients
    // 
    const double lam2 = lambda*lambda;
    const double lam4 = lam2 * lam2;
    const double lam6 = lam4 * lam2;
    const double lam8 = lam6 * lam2;
    const double lam10 = lam8 * lam2;
    const double lam14 = lam10 * lam4;
    const double lam18 = lam14 * lam4;
    const double r02 = r0 * r0;
    const double r04 = r02 * r02;
    const double r06 = r04 * r02;
    const double r08 = r06 * r02;
    expo = Amp * exp(-(r0*r0)/(lam2));
    term0 = - 8 * expo * (20*r04*lam2 - 60*r02*lam4 + 
			  15*lam6) / (5*lam10);
    A_term2 = 8 * expo * (-56*r06*lam2 + 420*r04*lam4 - 
			  630*r02*lam6 + 105*lam8) / (35*lam14);
    A_term4 = - 2 * expo * (288*r08*lam2 - 4032*r06*lam4 + 15120*r04*lam6 -
			    15120*r02*lam8 + 1890*lam10) / (315*lam18);
    B_term2 = 8 * expo * (-56*r06*lam2 + 420*r04*lam4 - 
			  630*r02*lam6 + 105*lam8) / (21*lam14);
    B_term4 = - 2 * expo * (288*r08*lam2 - 4032*r06*lam4 + 15120*r04*lam6 -
			    15120*r02*lam8 + 1890*lam10) / (135*lam18); 
    C_term2 = 8 * expo * (-56*r06*lam2 + 420*r04*lam4 - 
			  630*r02*lam6 + 105*lam8) / (15*lam14);
    C_term4 = - 26 * expo * (288*r08*lam2 - 4032*r06*lam4 + 15120*r04*lam6 -
			     15120*r02*lam8 + 1890*lam10) / (945*lam18);
    // CHECK!!
    delta_r = 1.e-4;
  };
  //================================================
  // Destructor
  //================================================
  ~LinWave();
  const char * Name() { return "Teukolsky wave initial data"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) {
    double r = 0.001;
    double t = 0.0;
    // cout << " F(r,t) = " << F(r,t) << endl;
    // cout << " F1(r,t) = " << F1(r,t) << endl;
    // cout << " F2(r,t) = " << F2(r,t) << endl;
    // cout << " F3(r,t) = " << F3(r,t) << endl;
    // cout << " F4(r,t) = " << F4(r,t) << endl;
    return true; 
  } 
  //================================================
  // Analytical solution for h_{ij}
  //================================================
  double h_rr_analytical(double r, double theta, double phi, double t) {
    return A(r,t)*f_rr(theta,phi);
  };
  double h_rt_analytical(double r, double theta, double phi, double t) {
    return B(r,t)*f_rt(theta,phi);
  };
  double h_rp_analytical(double r, double theta, double phi, double t) {
    return B(r,t)*f_rp(theta,phi);
  };
  double h_tt_analytical(double r, double theta, double phi, double t) {
    return C(r,t)*f_tt1(theta,phi) + A(r,t)*f_tt2(theta,phi);
  };
  double h_tp_analytical(double r, double theta, double phi, double t) {
    return (A(r,t) - 2.0*C(r,t))*f_tp(theta,phi);
  };
  double h_pp_analytical(double r, double theta, double phi, double t) {
    return C(r,t)*f_pp1(theta,phi) + A(r,t)*f_pp2(theta,phi);
  };
  double phi_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  //================================================
  // Analytical solution for connection coefficients
  //================================================
  double lam_r_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    double lam_r, lam_t, lam_p;
    done = compute_lambdas(r, theta, phi, t, lam_r, lam_t, lam_p);
    return lam_r;
  }
  double lam_t_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    double lam_r, lam_t, lam_p;
    done = compute_lambdas(r, theta, phi, t, lam_r, lam_t, lam_p);
    return lam_t;
  }
  double lam_p_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    double lam_r, lam_t, lam_p;
    done = compute_lambdas(r, theta, phi, t, lam_r, lam_t, lam_p);
    return lam_p;
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
  double sf_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
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
  //================================================
  // Auxiliary functions for linear wave
  //================================================
private:
  //================================================
  // radial functions
  //================================================
  inline double A(double r, double t) {
    if (r < r_exp && t == 0.0)
      return A_exp(r);
    else
      return 3.0 * ( r * ( r * F2(r,t) + 3.0*F1(r,t) ) + 3.0*F(r,t) ) / (r*r*r*r*r);
  };
  inline double B(double r, double t) {
    if (r < r_exp && t == 0.0)
      return B_exp(r);
    else
      return - ( r * ( r * ( r * F3(r,t) + 3.0*F2(r,t) ) + 6.0*F1(r,t) ) + 6.0*F(r,t) ) / (r*r*r*r*r);
  };
  inline double C(double r, double t) {
    if (r < r_exp && t == 0.0)
      return C_exp(r);
    else
      return 0.25 * ( r * ( r * ( r * ( r * F4(r,t) + 2.0*F3(r,t) ) + 9.0*F2(r,t)) + 21.0*F1(r,t)) + 21.0*F(r,t) ) / (r*r*r*r*r);
  };
  inline double F(double r, double t) {
    const double um = (t - r - r0);
    const double vm = (t + r - r0);
    const double uml = um / lambda;
    const double vml = vm / lambda;
    const double up = (t - r + r0);
    const double vp = (t + r + r0);
    const double upl = up / lambda;
    const double vpl = vp / lambda;
    const double uc = (t - r);
    const double vc = (t + r);
    // const double ucl = uc / lambda;
    // const double vcl = vc / lambda;
    const double pterm = ( uc * exp(-upl*upl) - 
			   vc * exp(-vpl*vpl) );
    const double mterm = ( uc * exp(-uml*uml) - 
			   vc * exp(-vml*vml) );
    return 0.5 * Amp * ( pterm + mterm );
  }
  inline double F1(double r, double t) {
    const double uml = (t - r - r0) / lambda;
    const double vml = (t + r - r0) / lambda;
    const double upl = (t - r + r0) / lambda;
    const double vpl = (t + r + r0) / lambda;
    const double ucl = (t - r) / lambda;
    const double vcl = (t + r) / lambda;
    const double pterm =  ( (1.0 - 2.0*ucl*upl) * exp(-upl*upl) + 
			    (1.0 - 2.0*vcl*vpl) * exp(-vpl*vpl) );
    const double mterm =  ( (1.0 - 2.0*ucl*uml) * exp(-uml*uml) + 
			    (1.0 - 2.0*vcl*vml) * exp(-vml*vml) );
    return 0.5 * Amp * ( pterm + mterm );
  }
  inline double F2(double r, double t) {
    const double uml = (t - r - r0) / lambda;
    const double vml = (t + r - r0) / lambda;
    const double upl = (t - r + r0) / lambda;
    const double vpl = (t + r + r0) / lambda;
    const double ucl = (t - r) / lambda;
    const double vcl = (t + r) / lambda;
    const double pterm = ( (-4.0*upl - 2.0*ucl + 4.0*ucl*upl*upl) * exp(-upl*upl) - 
			   (-4.0*vpl - 2.0*vcl + 4.0*vcl*vpl*vpl) * exp(-vpl*vpl) ) / lambda;
    const double mterm = ( (-4.0*uml - 2.0*ucl + 4.0*ucl*uml*uml) * exp(-uml*uml) - 
			   (-4.0*vml - 2.0*vcl + 4.0*vcl*vml*vml) * exp(-vml*vml) ) / lambda;
    return 0.5 * Amp * ( pterm + mterm );
  }
  inline double F3(double r, double t) {
    const double uml = (t - r - r0) / lambda;
    const double vml = (t + r - r0) / lambda;
    const double upl = (t - r + r0) / lambda;
    const double vpl = (t + r + r0) / lambda;
    const double ucl = (t - r) / lambda;
    const double vcl = (t + r) / lambda;
    const double pterm = ( (3.0 - 6.0*upl*upl - 6.0*ucl*upl + 4.0*ucl*upl*upl*upl) * exp(-upl*upl) + 
			   (3.0 - 6.0*vpl*vpl - 6.0*vcl*vpl + 4.0*vcl*vpl*vpl*vpl) * exp(-vpl*vpl) ) / (lambda*lambda);
    const double mterm = ( (3.0 - 6.0*uml*uml - 6.0*ucl*uml + 4.0*ucl*uml*uml*uml) * exp(-uml*uml) + 
			   (3.0 - 6.0*vml*vml - 6.0*vcl*vml + 4.0*vcl*vml*vml*vml) * exp(-vml*vml) ) / (lambda*lambda);
    return -2.0 * 0.5 * Amp * ( pterm + mterm );
  }
  inline double F4(double r, double t) {
    const double uml = (t - r - r0) / lambda;
    const double vml = (t + r - r0) / lambda;
    const double upl = (t - r + r0) / lambda;
    const double vpl = (t + r + r0) / lambda;
    const double ucl = (t - r) / lambda;
    const double vcl = (t + r) / lambda;
    const double pterm = ( ( 12.0*upl + 3.0*ucl - 8.0*upl*upl*upl - 12.0*ucl*upl*upl + 4.0*ucl*upl*upl*upl*upl) * exp(-upl*upl) -
			   ( 12.0*vpl + 3.0*vcl - 8.0*vpl*vpl*vpl - 12.0*vcl*vpl*vpl + 4.0*vcl*vpl*vpl*vpl*vpl) * exp(-vpl*vpl) ) 
      / (lambda*lambda*lambda);
    const double mterm = ( ( 12.0*uml + 3.0*ucl - 8.0*uml*uml*uml - 12.0*ucl*uml*uml + 4.0*ucl*uml*uml*uml*uml) * exp(-uml*uml) -
			   ( 12.0*vml + 3.0*vcl - 8.0*vml*vml*vml - 12.0*vcl*vml*vml + 4.0*vcl*vml*vml*vml*vml) * exp(-vml*vml) ) 
      / (lambda*lambda*lambda);
    return 4.0 * 0.5 * Amp * ( pterm + mterm );
  }
  //================================================
  // angular functions
  //================================================
  inline double f_rr(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    const double st2 = stl * stl;
    if (l == 2 && m == 0)
      return 2.0 - 3.0*st2;
    else if (l == 2 && m == 2)
      return stl*stl*cos(2.0*phi);
   else if (l == 3 && m == 0)
     return ctl*(2.0 - 5.0*st2);
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };
  inline double f_rt(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    const double st2 = stl * stl;
    if (l == 2 && m == 0)
      return - 3.0*stl*ctl;
    else if (l == 2 && m == 2)
      return stl*ctl*cos(2.0*phi);
   else if (l == 3 && m == 0)
     return ctl*(2.0 - 5.0*st2);
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };
  inline double f_rp(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    if (l == 2 && m == 0)
      return 0.0;
    else if (l == 2 && m == 2)
      return - stl*sin(2.0*phi); 
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };  
  inline double f_tt1(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    if (l == 2 && m == 0)
      return 3.0*stl*stl;
    else if (l == 2 && m == 2)
      return (1.0 + ctl*ctl)*cos(2.0*phi);
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };
  inline double f_tt2(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    if (l == 2 && m == 0)
      return -1.0;
    else if (l == 2 && m == 2)
      return - cos(2.0*phi);
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };
  inline double f_tp(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    if (l == 2 && m == 0)
      return 0.0;
    else if (l == 2 && m == 2)
      return ctl*sin(2.0*phi);
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };
  inline double f_pp1(double theta, double phi) {
    return - f_tt1(theta,phi);
  };
  inline double f_pp2(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    if (l == 2 && m == 0)
      return 3.0*stl*stl - 1.0;
    else if (l == 2 && m == 2)
      return ctl*ctl*cos(2.0*phi);
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };
  //======================================================
  // coefficients for expansion of A, B, and C about r = 0
  //======================================================
  inline double A_exp(double r) {
    return term0 + A_term2 * r*r + A_term4 * r*r*r*r;
  };
  //======================================================
  inline double B_exp(double r) {
    return term0 + B_term2 * r*r + B_term4 * r*r*r*r;
  };
  //======================================================
  inline double C_exp(double r) {
    return term0 + C_term2 * r*r + C_term4 * r*r*r*r;
  };
  //======================================================
  // Derivatives of radial functions...
  //======================================================
  inline double dAdr(double r, double t) {
    if (r < r_exp && t == 0.0)
      return 2.0*A_term2*r + 4.0*A_term4*r*r*r;
    else
      return ( 1.0  * (A(r + 3.0*delta_r,t) - A(r - 3.0*delta_r,t)) -
	       9.0  * (A(r + 2.0*delta_r,t) - A(r - 2.0*delta_r,t)) +
	       45.0 * (A(r +     delta_r,t) - A(r -     delta_r,t)) ) / 
	( 60.0*delta_r);
  };
  inline double dBdr(double r, double t) {
    if (r < r_exp && t == 0.0)
      return 2.0*B_term2*r + 4.0*B_term4*r*r*r;
    else
      return ( 1.0  * (B(r + 3.0*delta_r,t) - B(r - 3.0*delta_r,t)) -
	       9.0  * (B(r + 2.0*delta_r,t) - B(r - 2.0*delta_r,t)) +
	       45.0 * (B(r +     delta_r,t) - B(r -     delta_r,t)) ) / 
	( 60.0*delta_r);

  };
  inline double dCdr(double r, double t) {
    if (r < r_exp && t == 0.0)
      return 2.0*C_term2*r + 4.0*C_term4*r*r*r;
    else
      return ( 1.0  * (C(r + 3.0*delta_r,t) - C(r - 3.0*delta_r,t)) -
	       9.0  * (C(r + 2.0*delta_r,t) - C(r - 2.0*delta_r,t)) +
	       45.0 * (C(r +     delta_r,t) - C(r -     delta_r,t)) ) / 
	( 60.0*delta_r);
  };
  //================================================
  // derivatives of angular functions
  //================================================
  inline double df_rrdt(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    if (l == 2 && m == 0)
      return - 6.0*stl*ctl;
    else if (l == 2 && m == 2)
      return 2.0*stl*ctl*cos(2.0*phi);
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };
  inline double df_rtdt(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    if (l == 2 && m == 0)
      return - 3.0*ctl*ctl + 3.0*stl*stl;
    else if (l == 2 && m == 2)
      return (ctl*ctl - stl*stl)*cos(2.0*phi);
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };
  inline double df_rpdt(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    if (l == 2 && m == 0)
      return 0.0;
    else if (l == 2 && m == 2)
      return - ctl*sin(2.0*phi); 
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };  
  inline double df_tt1dt(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    if (l == 2 && m == 0)
      return 6.0*stl*ctl;
    else if (l == 2 && m == 2)
      return - 2.0*stl*ctl*cos(2.0*phi);
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };
  inline double df_tt2dt(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    if (l == 2 && m == 0)
      return 0.0;
    else if (l == 2 && m == 2)
      return 0.0;
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };
  inline double df_tpdt(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    if (l == 2 && m == 0)
      return 0.0;
    else if (l == 2 && m == 2)
      return - stl*sin(2.0*phi);
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };
  inline double df_pp1dt(double theta, double phi) {
    return - df_tt1dt(theta,phi);
  };
  inline double df_pp2dt(double theta, double phi) {
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    if (l == 2 && m == 0)
      return 6.0*stl*ctl;
    else if (l == 2 && m == 2)
      return -2.0*stl*ctl*cos(2.0*phi);
    else {
      cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
      return 0.0;
    }
  };
  //================================================
  // compute lambdas...
  //================================================
  bool compute_lambdas(double r, double theta, double phi, double t,
		       double & lam_r, double & lam_t, double & lam_p) { 
    const double stl = grid->sintheta(grid->j_ind(theta));
    const double ctl = grid->costheta(grid->j_ind(theta));
    const double r2 = r * r;
    const double st2 = stl*stl;
    //
    const double AA = A(r,t);
    const double BB = B(r,t);
    const double CC = C(r,t);
    //
    const double h_rr = AA * f_rr(theta,phi);
    const double h_rt = BB * f_rt(theta,phi);
    const double h_rp = BB * f_rp(theta,phi);
    const double h_tt = CC * f_tt1(theta,phi) + AA * f_tt2(theta,phi);
    const double h_tp = (AA - 2.0*CC) * f_tp(theta,phi);
    const double h_pp = CC * f_pp1(theta,phi) + AA * f_pp2(theta,phi);
    tensor g(1.0 + h_rr,
	     r * h_rt,
	     r * stl * h_rp,
	     r2 * (1.0 + h_tt),
	     r2 * stl * h_tp,
	     r2 * st2 * (1.0 + h_pp));
    tensor gup = g.inverse();
    // see eqs. (25) in Baumgarte et.al., PRD 87, 044026 (2013) 
    const double D_r_g_rr =           dAdr(r,t)*f_rr(theta,phi);
    const double D_r_g_rt = r *       dBdr(r,t)*f_rt(theta,phi);
    const double D_r_g_rp = r * stl * dBdr(r,t)*f_rp(theta,phi);
    const double D_r_g_tt = r2 *     (dCdr(r,t)*f_tt1(theta,phi) +
				      dAdr(r,t)*f_tt2(theta,phi) );
    const double D_r_g_tp = r2 * stl*(dAdr(r,t) - 2.0 * dCdr(r,t) )
      *f_tp(theta,phi);
    const double D_r_g_pp = r2 * st2*(dCdr(r,t)*f_pp1(theta,phi) +
				      dAdr(r,t)*f_pp2(theta,phi) );
    //
    const double D_t_g_rr =            AA * df_rrdt(theta,phi) - 2.0*h_rt;
    const double D_t_g_rt = r *       (BB * df_rtdt(theta,phi) + h_rr - h_tt);
    const double D_t_g_rp = r * stl * (BB * df_rpdt(theta,phi) - h_tp);
    const double D_t_g_tt = r2    *   (CC * df_tt1dt(theta,phi) +
				       AA * df_tt2dt(theta,phi) + 2.0 * h_rt);
    const double D_t_g_tp = r2 * stl* ((AA - 2.0*CC)*df_tpdt(theta,phi) + h_rp);
    const double D_t_g_pp = r2 * st2* (CC * df_pp1dt(theta,phi) + 
				       AA * df_pp2dt(theta,phi));
    // CHECK: assuming m = 0: no phi dependence, and h_rp = h_tp = 0.0;
    const double D_p_g_rr = 0.0;
    const double D_p_g_rt = 0.0;
    const double D_p_g_rp = r * stl * (stl * h_rr + ctl * h_rt - stl * h_pp);
    const double D_p_g_tt = 0.0;
    const double D_p_g_tp = r2* stl * (stl * h_rt + ctl * h_tt - ctl * h_pp);
    const double D_p_g_pp = 0.0;
    //
    rank3tens D_g(D_r_g_rr, D_r_g_rt, D_r_g_rp, D_r_g_tt, D_r_g_tp, D_r_g_pp, 
		  D_t_g_rr, D_t_g_rt, D_t_g_rp, D_t_g_tt, D_t_g_tp, D_t_g_pp, 
		  D_p_g_rr, D_p_g_rt, D_p_g_rp, D_p_g_tt, D_p_g_tp, D_p_g_pp);
    //
    lam_r = 0.0;
    lam_t = 0.0;
    lam_p = 0.0;
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
	for (int l = 0; l < 3; l++) {
	  lam_r += gup[j][k] * gup[0][l] *
	    ( D_g[j][k][l] + D_g[k][j][l] - D_g[l][j][k] );
	  lam_t += gup[j][k] * gup[1][l] *
	    ( D_g[j][k][l] + D_g[k][j][l] - D_g[l][j][k] );
	  lam_p += gup[j][k] * gup[2][l] *
	    ( D_g[j][k][l] + D_g[k][j][l] - D_g[l][j][k] );
	}
    lam_r *= 0.5;
    lam_t *= 0.5 * r;
    lam_p *= 0.5 * r * stl;
    return true;
  };
};

