// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Radiation Hydro Shock test - 
// see App. C in Farris et.al., PRD 78, 024023 (2008)
//
// Note: will assume axisymmetry; i.e. will line up x-axis of the 
// Farris et.al. solution with the symmetry axis here - i.e. the z-axis
//================================================
//

class RadHydroShockTest : public InData {
private:
  double L;  // outer boundary in supplied datafile
  double delta_z;  // spacing of data in arrays 
  int N_data;  // length of arrays read in from datafile
  VecDoub z_a, rho_a, P_a, uz_a, E_a, Fz_a;   // arrays for analytical solution
  double z_disc;  // look in datafile for directions to move discontinuity
  double beta, gamma;   // likewise for beta
  int N_g;
  bool all_clear;
public:
  //================================================
  // Constructor
  //================================================
  RadHydroShockTest(char * indata_input, Grid * grid_i, Cosmology * cosmology) :
    InData(grid_i, cosmology) {
    N_g = grid->N_ghosts();
    indata_type = radhydroshocktest;
    all_clear = true;
    analytical = true;
    ifstream infile;
    infile.open(indata_input);
    if (!infile) {
      cerr << " RADHYDROSHOCKTEST: Can't open " << indata_input 
	   << " for input. This is bad. " << endl;
      all_clear = false;
    } else 
      cout << " RADHYDROSHOCKTEST: Reading analytic solution from file " 
	   << indata_input << endl;
    //
    // read datafile once to find number of entries
    //
    string line;
    N_data = 0;
    z_disc = 0.0;
    beta = 0.0;
    while (getline(infile, line)) {
      if (line.at(0) != '#') N_data++;
      else {
	istringstream words(line);
	do { 
	  // also look for directions to move discontinuity...
	  string word;
	  words >> word;
	  if (word == "x_disc") {
	    words >> word; // read one more for = sign
	    words >> z_disc;
	  }
	  // ... or boost
	  if (word == "beta") {
	    words >> word; // read one more for = sign
	    words >> beta;
	  }
	} while (words);
      }
    }
    cout << " RADHYDROSHOCKTEST: found " << N_data 
	 << " data points in file " << indata_input << endl;
    gamma = 1.0 / sqrt(1.0 - beta*beta);
    cout << " RADHYDROSHOCKTEST: will use z_disc = " << z_disc 
	 << " and beta = " << beta << ", gamma = " << gamma << endl;
    //
    // now create arrays for data
    // 
    z_a = VecDoub(N_data);
    rho_a = VecDoub(N_data);
    P_a = VecDoub(N_data);
    uz_a = VecDoub(N_data);
    E_a = VecDoub(N_data);
    Fz_a = VecDoub(N_data);
    //
    // rewind file
    //
    infile.clear();
    infile.seekg(0, infile.beg);
    //
    // and read data into arrays
    //
    int i = 0;
    while (getline(infile, line)) {
      if (line.at(0) != '#') { 
	stringstream ss;
	ss << line;
	ss >> z_a[i] >> rho_a[i] >> P_a[i] >> uz_a[i] >> E_a[i] >> Fz_a[i];
	i++;
      }
    }   
    // outer boundary and spacing
    L = z_a[N_data - 1];
    delta_z = L/double(N_data/2 - 1);
    cout << " RADHYDROSHOCKTEST: outer boundary at " << L 
	 << ", delta_z = " << delta_z << endl; 
    infile.close();
  };
  //================================================
  // Destructor
  //================================================
  ~RadHydroShockTest() { };
  string Name() { return "radiation hydro shock initial data"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) { return all_clear; } 
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
    return 0.0;
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
  double rho_0_analytical(double r, double theta, double phi, double t) {
    double z = r * cos(theta) - z_disc;
    // boost to unprimed frame
    z = gamma*z;
    return interp_data(z, rho_a);
  };
  double P_analytical(double r, double theta, double phi, double t) {
    double z = r * cos(theta) - z_disc;
    return interp_data(z, P_a);
  };
  double v_r_analytical(double r, double theta, double phi, double t) {
    double costheta = cos(theta);
    double z = r * costheta - z_disc;
    // boost to unprimed frame
    z = gamma*z;
    double uz = interp_data(z, uz_a);
    double u0 = sqrt(1. + uz*uz);
    // now boost back to primed frame
    uz = gamma*(uz - beta*u0);
    u0 = sqrt(1. + uz*uz);
    double vz = uz / u0;
    return vz * costheta;
  };
  double v_t_analytical(double r, double theta, double phi, double t) {
    double costheta = cos(theta);
    double z = r * costheta - z_disc;
    // boost to unprimed frame
    z = gamma*z;
    double uz = interp_data(z, uz_a);
    double u0 = sqrt(1. + uz*uz);
    // now boost back to primed frame
    uz = gamma*(uz - beta*u0);
    u0 = sqrt(1. + uz*uz);
    double vz = uz / u0;
    return - sin(theta) * vz; // rescaled...
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
    double z = r * cos(theta) - z_disc;
    // boost to unprimed frame
    z = gamma*z;
    return interp_data(z, E_a);
  };
  double F_0_analytical(double r, double theta, double phi, double t) {
    double costheta = cos(theta);
    double z = r * costheta - z_disc;
    // boost to unprimed frame
    z = gamma*z;
    double Fz = interp_data(z, Fz_a);
    double uz = interp_data(z, uz_a);
    double u0 = sqrt(1.0 + uz*uz);
    double F0 = uz * Fz / u0;
    // now boost back to primed frame
    F0 = gamma*(F0 - beta*Fz);
    return F0;
  }
  double F_r_analytical(double r, double theta, double phi, double t) {
    double costheta = cos(theta);
    double z = r * costheta - z_disc;
    // boost to unprimed frame
    z = gamma*z;
    double Fz = interp_data(z, Fz_a);
    double uz = interp_data(z, uz_a);
    double u0 = sqrt(1.0 + uz*uz);
    double F0 = uz * Fz / u0;
    // now boost back to primed frame
    Fz = gamma*(Fz - beta*F0);
    return Fz * costheta;
  };
  double F_t_analytical(double r, double theta, double phi, double t) {
    double costheta = cos(theta);
    double z = r * costheta - z_disc;
    // boost to unprimed frame
    z = gamma*z;
    double Fz = interp_data(z, Fz_a);
    double uz = interp_data(z, uz_a);
    double u0 = sqrt(1.0 + uz*uz);
    double F0 = uz * Fz / u0;
    // now boost back to primed frame
    Fz = gamma*(Fz - beta*F0);
    return - Fz * sin(theta);  // rescaled...
  };
  double F_p_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
private:
  //================================================
  // Methods for interpolation of data
  //================================================
  double interp_data(double z, VecDoub f_a) {
    double f = 0;
    // first: handle points outside of available data 
    if (z < -L) {
      f = f_a[0];
    } else if (z > L) {
      f = f_a[N_data-1];
    } else {
      // otherwise find index i so that i and i+1 bracket z:
      // recall: two gridpoints at 0 in datafile - therefore treat 
      // positive and negative z separately.
      int i_low = 0;
      if (z < 0) {
	i_low = int( (z + L) / delta_z );   
	// make sure we don't extrapolate across z = 0:
	if (i_low > N_data/2 - 3) i_low = N_data/2 - 3;
	// also:
	if (i_low < 1) i_low = 1;
      } else {
	i_low = int( z / delta_z ) + N_data/2;
	// make sure we don't extrapolate across z = 0:
	if (i_low < N_data/2 + 1) i_low = N_data/2 + 1;
	// also:
	if (i_low > N_data -3) i_low = N_data - 3;
      } 
      // prepare for fourth-order interpolation
      int order = 4;
      Doub *z_int, *f_int;
      z_int = new Doub[order];
      f_int = new Doub[order];
      // Now move i_low over by one, so it becomes lowest point of 
      // four-point stencil
      i_low -= order/2 - 1;
      for (int i = 0; i < order; i++) {
	z_int[i] = z_a[i_low + i];      // ...and fill stencils with data
	f_int[i] = f_a[i_low + i];
      }
      //      cout << z_int[0] << "  " << z << "  " << z_int[order-1] << endl;
      //
      // Now follow numerical recipes routine Poly_interp
      // 
      Doub *c, *d; 
      c = new Doub[order];  // allocate arrays that store differences in tableaus
      d = new Doub[order];
      Doub dif = abs(z - z_int[0]);
      Doub dift;
      int ns=0;
      for (int i = 0; i < order; i++ ) {
	if ((dift=abs(z - z_int[i])) < dif) {
	  ns = i;
	  dif = dift;
	}
	c[i] = f_int[i];
	d[i] = f_int[i];
      }
      f = f_int[ns--];
      Doub den, ho, hp, w, dy;
      for (int m = 1; m < order; m++) {
	for (int i = 0; i < order-m; i++) {
	  ho = z_int[i] - z;
	  hp = z_int[i+m] - z;
	  w = c[i+1] - d[i];
	  den = ho - hp;
	  den = w/den;
	  d[i] = hp*den;
	  c[i] = ho*den;
	}
	z += (dy=(2*(ns+1) < (order - m) ? c[ns+1] : d[ns--]));
      }
      delete z_int;
      delete f_int;
      delete c;
      delete d;
    }
    return f;
  }
};

