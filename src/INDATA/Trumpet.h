// Tell emacs that this is -*-c++-*- mode
//================================================
// Schwarzschild initial data in trumpet topology
//================================================
//
#include "SBN.h"

class Trumpet : public InData {
private:
  double M;   
  double x_C, y_C, z_C;
  SBN *sbn;
  int array_length;
  double Rmax;
  double C;
public:
  //================================================
  // Constructor
  //================================================
  Trumpet(char * indata_input, Grid * grid_i, Cosmology * cosmology) :
    InData(grid_i, cosmology) {
    indata_type = trumpet;
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
    infile.get(buf,100,'='); infile.get(c); infile >> M;
    infile.get(buf,100,'='); infile.get(c); infile >> x_C;
    infile.get(buf,100,'='); infile.get(c); infile >> y_C;
    infile.get(buf,100,'='); infile.get(c); infile >> z_C;
    cout << " Will set up trumpet black hole with mass M = " << M << endl;
    cout << "    centered on x = " << x_C << ", y = " << y_C << ", z = " << z_C << endl;
    cout << "===================================================" << endl;
    // allocate SBN class
    array_length = 10000000;
    double r_max = grid_i->r_max();
    Rmax = 1.5 * r_max;
    sbn = new SBN(Rmax,array_length);
    // finally compute "magical" value of C
    C = 3.0*sqrt(3.0)*M*M/4.0;
  };
  //================================================
  // Destructor
  //================================================
  ~Trumpet() { delete sbn; };
  string Name() { return "maximally sliced trumpet initial data"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) { return true; } 
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
    return log(sbn->conFactor(abs(r_C),M));
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
    // find cartesian coordinates
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C) );
    // set up normal vector
    const double n_x = (x - x_C)/r_C;
    const double n_y = (y - y_C)/r_C;
    const double n_z = (z - z_C)/r_C;
    const double psi = sbn->conFactor(abs(r_C),M);
    const double R_C = psi*psi*r_C;
    const double factor = - C/(R_C*R_C*R_C);
    // compute extrinsic curvature in cartesian coordinates
    tensor a_cart(factor * ( 3.0 * n_x * n_x - 1.0 ),
		  factor * ( 3.0 * n_x * n_y ),
		  factor * ( 3.0 * n_x * n_z ),
		  factor * ( 3.0 * n_y * n_y - 1.0 ),
		  factor * ( 3.0 * n_y * n_z ),
		  factor * ( 3.0 * n_z * n_z - 1.0 ) );
    // transform to spherical polar coordinates
    tensor a = Cartesian_to_Spherical(a_cart, r, theta, phi);
    // and return rr component
    return a[0][0];
  };
  double a_rt_analytical(double r, double theta, double phi, double t) {
    // find cartesian coordinates
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C) );
    // set up normal vector
    const double n_x = (x - x_C)/r_C;
    const double n_y = (y - y_C)/r_C;
    const double n_z = (z - z_C)/r_C;
    const double psi = sbn->conFactor(abs(r_C),M);
    const double R_C = psi*psi*r_C;
    const double factor = - C/(R_C*R_C*R_C);
    // compute extrinsic curvature in cartesian coordinates
    tensor a_cart(factor * ( 3.0 * n_x * n_x - 1.0 ),
		  factor * ( 3.0 * n_x * n_y ),
		  factor * ( 3.0 * n_x * n_z ),
		  factor * ( 3.0 * n_y * n_y - 1.0 ),
		  factor * ( 3.0 * n_y * n_z ),
		  factor * ( 3.0 * n_z * n_z - 1.0 ) );
    // transform to spherical polar coordinates
    tensor a = Cartesian_to_Spherical(a_cart, r, theta, phi);
    // and return rt component
    return a[0][1]/r;
  };
  double a_rp_analytical(double r, double theta, double phi, double t) {
    // find cartesian coordinates
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C) );
    // set up normal vector
    const double n_x = (x - x_C)/r_C;
    const double n_y = (y - y_C)/r_C;
    const double n_z = (z - z_C)/r_C;
    const double psi = sbn->conFactor(abs(r_C),M);
    const double R_C = psi*psi*r_C;
    const double factor = - C/(R_C*R_C*R_C);
    // compute extrinsic curvature in cartesian coordinates
    tensor a_cart(factor * ( 3.0 * n_x * n_x - 1.0 ),
		  factor * ( 3.0 * n_x * n_y ),
		  factor * ( 3.0 * n_x * n_z ),
		  factor * ( 3.0 * n_y * n_y - 1.0 ),
		  factor * ( 3.0 * n_y * n_z ),
		  factor * ( 3.0 * n_z * n_z - 1.0 ) );
    // transform to spherical polar coordinates
    tensor a = Cartesian_to_Spherical(a_cart, r, theta, phi);
    // and return rp component
    return a[0][2]/(r*sin(theta));
  };
  double a_tt_analytical(double r, double theta, double phi, double t) {
    // find cartesian coordinates
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C) );
    // set up normal vector
    const double n_x = (x - x_C)/r_C;
    const double n_y = (y - y_C)/r_C;
    const double n_z = (z - z_C)/r_C;
    const double psi = sbn->conFactor(abs(r_C),M);
    const double R_C = psi*psi*r_C;
    const double factor = - C/(R_C*R_C*R_C);
    // compute extrinsic curvature in cartesian coordinates
    tensor a_cart(factor * ( 3.0 * n_x * n_x - 1.0 ),
		  factor * ( 3.0 * n_x * n_y ),
		  factor * ( 3.0 * n_x * n_z ),
		  factor * ( 3.0 * n_y * n_y - 1.0 ),
		  factor * ( 3.0 * n_y * n_z ),
		  factor * ( 3.0 * n_z * n_z - 1.0 ) );
    // transform to spherical polar coordinates
    tensor a = Cartesian_to_Spherical(a_cart, r, theta, phi);
    // and return tt component
    return a[1][1]/(r*r);
  };
  double a_tp_analytical(double r, double theta, double phi, double t) {
    // find cartesian coordinates
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C) );
    // set up normal vector
    const double n_x = (x - x_C)/r_C;
    const double n_y = (y - y_C)/r_C;
    const double n_z = (z - z_C)/r_C;

    const double psi = sbn->conFactor(abs(r_C),M);
    const double R_C = psi*psi*r_C;
    const double factor = - C/(R_C*R_C*R_C);
    // compute extrinsic curvature in cartesian coordinates
    tensor a_cart(factor * ( 3.0 * n_x * n_x - 1.0 ),
		  factor * ( 3.0 * n_x * n_y ),
		  factor * ( 3.0 * n_x * n_z ),
		  factor * ( 3.0 * n_y * n_y - 1.0 ),
		  factor * ( 3.0 * n_y * n_z ),
		  factor * ( 3.0 * n_z * n_z - 1.0 ) );
    // transform to spherical polar coordinates
    tensor a = Cartesian_to_Spherical(a_cart, r, theta, phi);
    // and return tp component
    return a[1][2]/(r*r*sin(theta));
  };
  double a_pp_analytical(double r, double theta, double phi, double t) {
    // find cartesian coordinates
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C) );
    // set up normal vector
    const double n_x = (x - x_C)/r_C;
    const double n_y = (y - y_C)/r_C;
    const double n_z = (z - z_C)/r_C;
    const double psi = sbn->conFactor(abs(r_C),M);
    const double R_C = psi*psi*r_C;
    const double factor = - C/(R_C*R_C*R_C);
    // compute extrinsic curvature in cartesian coordinates
    tensor a_cart(factor * ( 3.0 * n_x * n_x - 1.0 ),
		  factor * ( 3.0 * n_x * n_y ),
		  factor * ( 3.0 * n_x * n_z ),
		  factor * ( 3.0 * n_y * n_y - 1.0 ),
		  factor * ( 3.0 * n_y * n_z ),
		  factor * ( 3.0 * n_z * n_z - 1.0 ) );
    // transform to spherical polar coordinates
    tensor a = Cartesian_to_Spherical(a_cart, r, theta, phi);
    // and return rr component
    return a[2][2]/(r*r*sin(theta)*sin(theta));
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
    return sbn->lapse(abs(r_C),M);
  };
  double shift_r_analytical(double r, double theta, double phi, double t) {
    // find cartesian coordinates
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C) );
    // set up normal vector
    const double n_x = (x - x_C)/r_C;
    const double n_y = (y - y_C)/r_C;
    const double n_z = (z - z_C)/r_C;
    // compute shift vector in cartesian coordinates
    const double shift_norm = sbn->shift(abs(r_C),M);
    vect beta_cart(shift_norm * n_x,
		   shift_norm * n_y,
		   shift_norm * n_z );
    // transform to spherical polar coordinates
    vect beta = Cartesian_to_Spherical_upper(beta_cart,r,theta,phi);
    // return r component
    return beta[0];
  };
  double shift_t_analytical(double r, double theta, double phi, double t) {
    // find cartesian coordinates
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C) );
    // set up normal vector
    const double n_x = (x - x_C)/r_C;
    const double n_y = (y - y_C)/r_C;
    const double n_z = (z - z_C)/r_C;
    // compute shift vector in cartesian coordinates
    const double shift_norm = sbn->shift(abs(r_C),M);
    vect beta_cart(shift_norm * n_x,
		   shift_norm * n_y,
		   shift_norm * n_z );
    // transform to spherical polar coordinates
    vect beta = Cartesian_to_Spherical_upper(beta_cart,r,theta,phi);
    // return r component
    return beta[1]*r;
  };
  double shift_p_analytical(double r, double theta, double phi, double t) {
    // find cartesian coordinates
    const double x = r * sin(theta) * cos(phi);
    const double y = r * sin(theta) * sin(phi);
    const double z = r * cos(theta);
    const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C) );
    // set up normal vector
    const double n_x = (x - x_C)/r_C;
    const double n_y = (y - y_C)/r_C;
    const double n_z = (z - z_C)/r_C;
    // compute shift vector in cartesian coordinates
    const double shift_norm = sbn->shift(abs(r_C),M);
    vect beta_cart(shift_norm * n_x,
		   shift_norm * n_y,
		   shift_norm * n_z );
    // transform to spherical polar coordinates
    vect beta = Cartesian_to_Spherical_upper(beta_cart,r,theta,phi);
    // return r component
    return beta[2]*r*sin(theta);
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
};

