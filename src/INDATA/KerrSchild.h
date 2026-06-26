// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Kerr-Schild initial data 
//================================================
//
class KerrSchild : public InData {
private:
  double x_C, y_C, z_C;
  double M;
  double a, a2; 
  double eps;   // small number for computing derivatives
public:
  //================================================
  // Constructor
  //================================================
  KerrSchild(char * indata_input, Grid * grid_i, Cosmology * cosmology) :
    InData(grid_i, cosmology) {
    indata_type = kerrschild;
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
    infile.get(buf,100,'='); infile.get(c); infile >> M;
    infile.get(buf,100,'='); infile.get(c); infile >> a;
    cout << " Will set up Kerr black hole with mass M = " << M << endl;
    cout << "    and angular momentum parameter a = " << a << endl; 
    cout << "===================================================" << endl;
    a2 = a*a;
    eps = 1.e-6;
  };
  //================================================
  // Destructor
  //================================================
  ~KerrSchild() {};
  string Name() { return "Kerr initial data in Kerr-Schild coordinates"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) { 
    //
    // tests...
    //
    // on z-axis:
    // double x = 0.0;
    // double y = 0.0;
    // double z = 1.0;
    // double r = 1.0;
    // double theta = 0.0;
    // double phi = 0.0;
    // cout << " on z-axis (cartesian): l = (" << l_x(x,y,z) << "," << l_y(x,y,z) << "," << l_z(x,y,z) << ")" 
    // 	 << ",  H = " << H(x,y,z) << endl; 
    // cout << " on z-axis (spherical): l = (" << l_r(r,theta,phi) << "," << l_t(r,theta,phi) << "," << l_p(r,theta,phi) << ")" << endl; 
    // // on x-axis:
    // x = 1.0;
    // y = 0.0;
    // z = 0.0;
    // r = 1.0;
    // theta = PI/2.0;
    // phi = 0.0;
    // cout << " on x-axis (cartesian): l = (" << l_x(x,y,z) << "," << l_y(x,y,z) << "," << l_z(x,y,z) << ")" 
    // 	 << ",  H = " << H(x,y,z)<< endl; 
    // cout << " on x-axis (spherical): l = (" << l_r(r,theta,phi) << "," << l_t(r,theta,phi) << "," << l_p(r,theta,phi) << ")" << endl; 
    return true; 
  } 
  //================================================
  // Analytical solution for h_{ij}
  //================================================
  double h_rr_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    return 2.0*H(x,y,z)*l_r(r,theta,phi)*l_r(r,theta,phi);
  };
  double h_rt_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    return 2.0*H(x,y,z)*l_r(r,theta,phi)*l_t(r,theta,phi)/r;
  };
  double h_rp_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    return 2.0*H(x,y,z)*l_r(r,theta,phi)*l_p(r,theta,phi)/(r*sin(theta));
  };
  double h_tt_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    return 2.0*H(x,y,z)*l_t(r,theta,phi)*l_t(r,theta,phi)/(r*r);
  };
  double h_tp_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    return 2.0*H(x,y,z)*l_t(r,theta,phi)*l_p(r,theta,phi)/(r*r*sin(theta));
  };
  double h_pp_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    return 2.0*H(x,y,z)*l_p(r,theta,phi)*l_p(r,theta,phi)/(r*r*sin(theta)*sin(theta));
  };
  double phi_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  //================================================
  // Analytical solution for connection coefficients
  //================================================
  double lam_r_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = false; return 0;
  }
  double lam_t_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = false; return 0;
  }
  double lam_p_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = false; return 0;
  }
  //================================================
  // Analytical solution for extrinsic curvature
  //================================================
  double a_rr_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    // in cartesian coordinates...
    vect D_H(H_x(x,y,z), H_y(x,y,z), H_z(x,y,z));
    vect l(l_x(x,y,z), l_y(x,y,z), l_z(x,y,z));
    // notation D_l[1][2] = partial_1 l_2 = l_2,1
    tensor D_l(l_xx(x,y,z), l_yx(x,y,z), l_zx(x,y,z),
    	       l_xy(x,y,z), l_yy(x,y,z), l_zy(x,y,z),
    	       l_xz(x,y,z), l_yz(x,y,z), l_zz(x,y,z));
    tensor A_cart;
    double H_l = H(x,y,z);
    double lapse_l = lapse_analytical(r,theta,phi,t);
    for (int i = 0; i<3; i++)
      for (int j = 0; j<3; j++) {
	A_cart[i][j] = lapse_l*( l[i]*D_H[j] + l[j]*D_H[i] + H_l*(D_l[i][j] + D_l[j][i] ) );
	//	A_ana[i][j] = 2.0*M/(r*r*sqrt(1.0 + 2.0*M/r))*(eta[i][j] - (2.0 + M/r)*l[i]*l[j]); 
	for (int k = 0; k<3; k++) 
	  A_cart[i][j] += 2.0*lapse_l*H_l*l[k] * (l[i]*l[j]*D_H[k] + H_l*( l[i]*D_l[k][j] + l[j]*D_l[k][i] ));
      } 
    tensor g(1.0 + 2.0*H_l*l[0]*l[0],
	           2.0*H_l*l[0]*l[1],
	           2.0*H_l*l[0]*l[2],
	     1.0 + 2.0*H_l*l[1]*l[1],
	           2.0*H_l*l[1]*l[2],
	     1.0 + 2.0*H_l*l[2]*l[2]);
    tensor gup = g.inverse();
    A_cart.remove_trace(gup,g);
    // transform into spherical polar coordinates
    tensor A = Cartesian_to_Spherical(A_cart, r, theta, phi);
    return A[0][0];
  };
  double a_rt_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    // in cartesian coordinates...
    vect D_H(H_x(x,y,z), H_y(x,y,z), H_z(x,y,z));
    vect l(l_x(x,y,z), l_y(x,y,z), l_z(x,y,z));
    // notation D_l[1][2] = partial_1 l_2 = l_2,1
    tensor D_l(l_xx(x,y,z), l_yx(x,y,z), l_zx(x,y,z),
    	       l_xy(x,y,z), l_yy(x,y,z), l_zy(x,y,z),
    	       l_xz(x,y,z), l_yz(x,y,z), l_zz(x,y,z));
    tensor A_cart;
    double H_l = H(x,y,z);
    double lapse_l = lapse_analytical(r,theta,phi,t);
    for (int i = 0; i<3; i++)
      for (int j = 0; j<3; j++) {
	A_cart[i][j] = lapse_l*( l[i]*D_H[j] + l[j]*D_H[i] + H_l*(D_l[i][j] + D_l[j][i] ) );
	for (int k = 0; k<3; k++) 
	  A_cart[i][j] += 2.0*lapse_l*H_l*l[k] * (l[i]*l[j]*D_H[k] + H_l*( l[i]*D_l[k][j] + l[j]*D_l[k][i] ));
      } 
    tensor g(1.0 + 2.0*H_l*l[0]*l[0],
	           2.0*H_l*l[0]*l[1],
	           2.0*H_l*l[0]*l[2],
	     1.0 + 2.0*H_l*l[1]*l[1],
	           2.0*H_l*l[1]*l[2],
	     1.0 + 2.0*H_l*l[2]*l[2]);
    tensor gup = g.inverse();
    A_cart.remove_trace(gup,g);
    // transform into spherical polar coordinates
    tensor A = Cartesian_to_Spherical(A_cart, r, theta, phi);
    return A[0][1]/r;
  };
  double a_rp_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    // in cartesian coordinates...
    vect D_H(H_x(x,y,z), H_y(x,y,z), H_z(x,y,z));
    vect l(l_x(x,y,z), l_y(x,y,z), l_z(x,y,z));
    // notation D_l[1][2] = partial_1 l_2 = l_2,1
    tensor D_l(l_xx(x,y,z), l_yx(x,y,z), l_zx(x,y,z),
    	       l_xy(x,y,z), l_yy(x,y,z), l_zy(x,y,z),
    	       l_xz(x,y,z), l_yz(x,y,z), l_zz(x,y,z));
    tensor A_cart;
    double H_l = H(x,y,z);
    double lapse_l = lapse_analytical(r,theta,phi,t);
    for (int i = 0; i<3; i++)
      for (int j = 0; j<3; j++) {
	A_cart[i][j] = lapse_l*( l[i]*D_H[j] + l[j]*D_H[i] + H_l*(D_l[i][j] + D_l[j][i] ) );
	for (int k = 0; k<3; k++) 
	  A_cart[i][j] += 2.0*lapse_l*H_l*l[k] * (l[i]*l[j]*D_H[k] + H_l*( l[i]*D_l[k][j] + l[j]*D_l[k][i] ));
      } 
    tensor g(1.0 + 2.0*H_l*l[0]*l[0],
	           2.0*H_l*l[0]*l[1],
	           2.0*H_l*l[0]*l[2],
	     1.0 + 2.0*H_l*l[1]*l[1],
	           2.0*H_l*l[1]*l[2],
	     1.0 + 2.0*H_l*l[2]*l[2]);
    tensor gup = g.inverse();
    A_cart.remove_trace(gup,g);
    // transform into spherical polar coordinates
    tensor A = Cartesian_to_Spherical(A_cart, r, theta, phi);
    return A[0][2]/(r*sin(theta));
  };
  double a_tt_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    // in cartesian coordinates...
    vect D_H(H_x(x,y,z), H_y(x,y,z), H_z(x,y,z));
    vect l(l_x(x,y,z), l_y(x,y,z), l_z(x,y,z));
    // notation D_l[1][2] = partial_1 l_2 = l_2,1
    tensor D_l(l_xx(x,y,z), l_yx(x,y,z), l_zx(x,y,z),
    	       l_xy(x,y,z), l_yy(x,y,z), l_zy(x,y,z),
    	       l_xz(x,y,z), l_yz(x,y,z), l_zz(x,y,z));
    tensor A_cart;
    double H_l = H(x,y,z);
    double lapse_l = lapse_analytical(r,theta,phi,t);
    for (int i = 0; i<3; i++)
      for (int j = 0; j<3; j++) {
	A_cart[i][j] = lapse_l*( l[i]*D_H[j] + l[j]*D_H[i] + H_l*(D_l[i][j] + D_l[j][i] ) );
	for (int k = 0; k<3; k++) 
	  A_cart[i][j] += 2.0*lapse_l*H_l*l[k] * (l[i]*l[j]*D_H[k] + H_l*( l[i]*D_l[k][j] + l[j]*D_l[k][i] ));
      } 
    tensor g(1.0 + 2.0*H_l*l[0]*l[0],
	           2.0*H_l*l[0]*l[1],
	           2.0*H_l*l[0]*l[2],
	     1.0 + 2.0*H_l*l[1]*l[1],
	           2.0*H_l*l[1]*l[2],
	     1.0 + 2.0*H_l*l[2]*l[2]);
    tensor gup = g.inverse();
    A_cart.remove_trace(gup,g);
    // transform into spherical polar coordinates
    tensor A = Cartesian_to_Spherical(A_cart, r, theta, phi);
    return A[1][1]/(r*r);
  };
  double a_tp_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    // in cartesian coordinates...
    vect D_H(H_x(x,y,z), H_y(x,y,z), H_z(x,y,z));
    vect l(l_x(x,y,z), l_y(x,y,z), l_z(x,y,z));
    // notation D_l[1][2] = partial_1 l_2 = l_2,1
    tensor D_l(l_xx(x,y,z), l_yx(x,y,z), l_zx(x,y,z),
    	       l_xy(x,y,z), l_yy(x,y,z), l_zy(x,y,z),
    	       l_xz(x,y,z), l_yz(x,y,z), l_zz(x,y,z));
    tensor A_cart;
    double H_l = H(x,y,z);
    double lapse_l = lapse_analytical(r,theta,phi,t);
    for (int i = 0; i<3; i++)
      for (int j = 0; j<3; j++) {
	A_cart[i][j] = lapse_l*( l[i]*D_H[j] + l[j]*D_H[i] + H_l*(D_l[i][j] + D_l[j][i] ) );
	for (int k = 0; k<3; k++) 
	  A_cart[i][j] += 2.0*lapse_l*H_l*l[k] * (l[i]*l[j]*D_H[k] + H_l*( l[i]*D_l[k][j] + l[j]*D_l[k][i] ));
      } 
    tensor g(1.0 + 2.0*H_l*l[0]*l[0],
	           2.0*H_l*l[0]*l[1],
	           2.0*H_l*l[0]*l[2],
	     1.0 + 2.0*H_l*l[1]*l[1],
	           2.0*H_l*l[1]*l[2],
	     1.0 + 2.0*H_l*l[2]*l[2]);
    tensor gup = g.inverse();
    A_cart.remove_trace(gup,g);
    // transform into spherical polar coordinates
    tensor A = Cartesian_to_Spherical(A_cart, r, theta, phi);
    return A[1][2]/(r*r*sin(theta));
  };
  double a_pp_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    // in cartesian coordinates...
    vect D_H(H_x(x,y,z), H_y(x,y,z), H_z(x,y,z));
    vect l(l_x(x,y,z), l_y(x,y,z), l_z(x,y,z));
    // notation D_l[1][2] = partial_1 l_2 = l_2,1
    tensor D_l(l_xx(x,y,z), l_yx(x,y,z), l_zx(x,y,z),
    	       l_xy(x,y,z), l_yy(x,y,z), l_zy(x,y,z),
    	       l_xz(x,y,z), l_yz(x,y,z), l_zz(x,y,z));
    tensor A_cart;
    double H_l = H(x,y,z);
    double lapse_l = lapse_analytical(r,theta,phi,t);
    for (int i = 0; i<3; i++)
      for (int j = 0; j<3; j++) {
	A_cart[i][j] = lapse_l*( l[i]*D_H[j] + l[j]*D_H[i] + H_l*(D_l[i][j] + D_l[j][i] ) );
	for (int k = 0; k<3; k++) 
	  A_cart[i][j] += 2.0*lapse_l*H_l*l[k] * (l[i]*l[j]*D_H[k] + H_l*( l[i]*D_l[k][j] + l[j]*D_l[k][i] ));
      } 
    tensor g(1.0 + 2.0*H_l*l[0]*l[0],
	           2.0*H_l*l[0]*l[1],
	           2.0*H_l*l[0]*l[2],
	     1.0 + 2.0*H_l*l[1]*l[1],
	           2.0*H_l*l[1]*l[2],
	     1.0 + 2.0*H_l*l[2]*l[2]);
    tensor gup = g.inverse();
    A_cart.remove_trace(gup,g);
    // transform into spherical polar coordinates
    tensor A = Cartesian_to_Spherical(A_cart, r, theta, phi);
    return A[2][2]/(r*r*sin(theta)*sin(theta));
  };
  double K_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    // in cartesian coordinates...
    vect D_H(H_x(x,y,z), H_y(x,y,z), H_z(x,y,z));
    vect l(l_x(x,y,z), l_y(x,y,z), l_z(x,y,z));
    // notation D_l[1][2] = partial_1 l_2 = l_2,1
    tensor D_l(l_xx(x,y,z), l_yx(x,y,z), l_zx(x,y,z),
	       l_xy(x,y,z), l_yy(x,y,z), l_zy(x,y,z),
	       l_xz(x,y,z), l_yz(x,y,z), l_zz(x,y,z));
    tensor A_cart;
    double H_l = H(x,y,z);
    double lapse_l = lapse_analytical(r,theta,phi,t);
    for (int i = 0; i<3; i++)
      for (int j = 0; j<3; j++) {
	A_cart[i][j] = lapse_l*( l[i]*D_H[j] + l[j]*D_H[i] + H_l*(D_l[i][j] + D_l[j][i] ) );
	for (int k = 0; k<3; k++) 
	  A_cart[i][j] += 2.0*lapse_l*H_l*l[k] * (l[i]*l[j]*D_H[k] + H_l*( l[i]*D_l[k][j] + l[j]*D_l[k][i] ));
      } 
    // compute trace...
    tensor g(1.0 + 2.0*H_l*l[0]*l[0],
	           2.0*H_l*l[0]*l[1],
	           2.0*H_l*l[0]*l[2],
	     1.0 + 2.0*H_l*l[1]*l[1],
	           2.0*H_l*l[1]*l[2],
	     1.0 + 2.0*H_l*l[2]*l[2]);
    tensor gup = g.inverse();
    // if (r > 1.0 && r < 1.01 && theta > 1.0 && theta < 1.2 && phi < 3.0) {
    //   cout << "r = " << r << endl;
    //   g.print();
    //   cout << " A_cart: " << endl;
    //   A_cart.print();
    //   double K_alt = 0.0;
    //   for (int i = 0; i<3; i++) {
    // 	K_alt += 2.0 * lapse_l* lapse_l* lapse_l * (1.0 + H_l) * l[i] * D_H[i] + 2.0 * lapse_l * H_l * D_l[i][i];
    //   }
    //   cout << " K = " << A_cart.trace(gup) 
    // 	   << " or " << K_alt << endl;
    // 	//	   << ", analytical = " << 2.0*M/(r*r*(1.0 + 2.0*M/r)*sqrt(1.0 + 2.0*M/r))*(1.0 + 3.0*M/r) << endl;
    //   double ll = 0.0;
    //   tensor eta(1.0,0.0,0.0,1.0,0.0,1.0);
    //   for (int i = 0; i<3; i++)
    // 	for (int j = 0; j<3; j++) {
    // 	  ll += eta[i][j] * l[i] * l[j];
    // 	}
    //   cout << " check l's " << ll << endl;
    // }
    return A_cart.trace(gup);
  };
  //================================================
  // Analytical solution for gauge
  //================================================
  double lapse_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    return 1.0/sqrt(1.0 + 2.0*H(x,y,z));
  };
  double shift_r_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double Hl = H(x,y,z);
    return 2.0*Hl/(1.0 + 2.0*Hl) * l_r(r,theta,phi);
  };
  double shift_t_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double Hl = H(x,y,z);
    return 2.0*Hl/(1.0 + 2.0*Hl) * l_t(r,theta,phi) / r;
  };
  double shift_p_analytical(double r, double theta, double phi, double t) {
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double Hl = H(x,y,z);
    return 2.0*Hl/(1.0 + 2.0*Hl) * l_p(r,theta,phi) / (r*sin(theta));
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
  //================================================
  // lots of auxiliary functions...
  //================================================
  inline double r_KS(double x, double y, double z) {
    const double xyza = x*x + y*y + z*z - a2; 
    return sqrt(xyza/2.0 + sqrt( xyza*xyza/4.0 + a2*z*z));
  };
  //================================================
  inline double H(double x, double y, double z) {
    const double r = r_KS(x,y,z);
    const double costheta = z/r;
    return M*r/(r*r + a2*costheta*costheta);
  };
  //================================================
  inline double H_x(double x, double y, double z) {
    return (H(x+eps,y,z) - H(x-eps,y,z))/(2.0*eps);
  }
  inline double H_y(double x, double y, double z) {
    return (H(x,y+eps,z) - H(x,y-eps,z))/(2.0*eps);
  }
  inline double H_z(double x, double y, double z) {
    return (H(x,y,z+eps) - H(x,y,z-eps))/(2.0*eps);
  }
  //================================================
  inline double H_xx(double x, double y, double z) {
    return (H(x+eps,y,z) - 2.0*H(x,y,z) + H(x-eps,y,z))/(eps*eps);
  }
  inline double H_yy(double x, double y, double z) {
    return (H(x,y+eps,z) - 2.0*H(x,y,z) + H(x,y-eps,z))/(eps*eps);
  }
  inline double H_zz(double x, double y, double z) {
    return (H(x,y,z+eps) - 2.0*H(x,y,z) + H(x,y,z-eps))/(eps*eps);
  }
  inline double H_xy(double x, double y, double z) {
    return (H(x+eps,y+eps,z) - H(x+eps,y-eps,z) - H(x-eps,y+eps,z) + H(x-eps,y-eps,z))/(4.0*eps*eps);
  }
  inline double H_xz(double x, double y, double z) {
    return (H(x+eps,y,z+eps) - H(x+eps,y,z-eps) - H(x-eps,y,z+eps) + H(x-eps,y,z-eps))/(4.0*eps*eps);
  }
  inline double H_yz(double x, double y, double z) {
    return (H(x,y+eps,z+eps) - H(x,y+eps,z-eps) - H(x,y-eps,z+eps) + H(x,y-eps,z-eps))/(4.0*eps*eps);
  }
  //================================================
  inline double l_x(double x, double y, double z) {
    const double r = r_KS(x,y,z);
    return (r*x + a*y)/(r*r + a2);
  }
  inline double l_y(double x, double y, double z) {
    const double r = r_KS(x,y,z);
    return (r*y - a*x)/(r*r + a2);
  }
  inline double l_z(double x, double y, double z) {
    const double r = r_KS(x,y,z);
    return z/r;
  }
  //================================================
  inline double l_xx(double x, double y, double z) {
    return (l_x(x+eps,y,z) - l_x(x-eps,y,z))/(2.0*eps);
  }
  inline double l_xy(double x, double y, double z) {
    return (l_x(x,y+eps,z) - l_x(x,y-eps,z))/(2.0*eps);
  }
  inline double l_xz(double x, double y, double z) {
    return (l_x(x,y,z+eps) - l_x(x,y,z-eps))/(2.0*eps);
  }
  //================================================
  inline double l_yx(double x, double y, double z) {
    return (l_y(x+eps,y,z) - l_y(x-eps,y,z))/(2.0*eps);
  }
  inline double l_yy(double x, double y, double z) {
    return (l_y(x,y+eps,z) - l_y(x,y-eps,z))/(2.0*eps);
  }
  inline double l_yz(double x, double y, double z) {
    return (l_y(x,y,z+eps) - l_y(x,y,z-eps))/(2.0*eps);
  }
  //================================================
  inline double l_zx(double x, double y, double z) {
    return (l_z(x+eps,y,z) - l_z(x-eps,y,z))/(2.0*eps);
  }
  inline double l_zy(double x, double y, double z) {
    return (l_z(x,y+eps,z) - l_z(x,y-eps,z))/(2.0*eps);
  }
  inline double l_zz(double x, double y, double z) {
    return (l_z(x,y,z+eps) - l_z(x,y,z-eps))/(2.0*eps);
  }
  //================================================
  // lower components of l in spherical polar coords
  //================================================
  inline double l_r(double r, double theta, double phi) {
    const double st = sin(theta);
    const double ct = cos(theta);
    const double sp = sin(phi);
    const double cp = cos(phi);
    const double x = r * st * cp;
    const double y = r * st * sp;
    const double z = r * ct;
    return st * cp * l_x(x,y,z) + st * sp * l_y(x,y,z) + ct * l_z(x,y,z);
  }
  inline double l_t(double r, double theta, double phi) {
    const double st = sin(theta);
    const double ct = cos(theta);
    const double sp = sin(phi);
    const double cp = cos(phi);
    const double x = r * st * cp;
    const double y = r * st * sp;
    const double z = r * ct;
    return r * ct * cp * l_x(x,y,z) + r * ct * sp * l_y(x,y,z) - r * st * l_z(x,y,z);
  }
  inline double l_p(double r, double theta, double phi) {
    const double st = sin(theta);
    const double ct = cos(theta);
    const double sp = sin(phi);
    const double cp = cos(phi);
    const double x = r * st * cp;
    const double y = r * st * sp;
    const double z = r * ct;
    return - r * st * sp * l_x(x,y,z) + r * st * cp * l_y(x,y,z);
  }


};
