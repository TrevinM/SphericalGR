// Tell emacs that this is -*-c++-*- mode
//====================================================
// Numerically Solve for \bar A^ij from W^i NOW /bar A_ij from W^i
//====================================================

void Compute_Aij() {
    for (int i = N_g; i < n_r - N_g; i++) {
        const double rl = rho.r(i);
        const double r2 = rl * rl;
        const double r4 = r2 * r2;
        for (int j = N_g; j < n_theta - N_g; j++) {
            const double thetal = rho.theta(j);
            const double stl = rho.sintheta(j);
            const double st2 = stl * stl;
            const double st4 = st2 * st2;
            const double csl = rho.costheta(j);
            const double ctl = csl / stl;
            for (int k = N_g; k < n_phi - N_g; k++) {
                const double div = W_r.dr(i, j, k) + 2. * W_r(i, j, k) / rl
                    + W_t.dtheta(i, j, k) / rl + W_t(i, j, k) * ctl / rl;
                // compute contravariant derivatives D^i W^j
                const double DrWr = W_r.dr(i, j, k);
                const double DtWr = (W_r.dtheta(i, j, k) - W_t(i, j, k)) / rl;
                const double DpWr = -W_p(i, j, k) / rl;
                const double DrWt = W_t.dr(i, j, k);
                const double DtWt = (W_t.dtheta(i, j, k) + W_r(i, j, k)) / rl;
                const double DpWt = -ctl * W_p(i, j, k) / rl;
                const double DrWp = W_p.dr(i, j, k);
                const double DtWp = (W_p.dtheta(i, j, k) + ctl * W_p(i,j,k)) / rl;
                const double DpWp = (W_r(i, j, k) + ctl * W_t(i, j, k)) / rl;
                // Contravariant, but should be the same if rescaled as conformally flat
                A_rr[i][j][k] = 2.0 * DrWr - 2. / 3. * div + An_rr(rl, thetal);
                A_rt[i][j][k] = (DrWt + DtWr) + An_rt(rl, thetal);
                A_rp[i][j][k] = (DrWp + DpWr) + An_rp(rl, thetal);
                A_tt[i][j][k] = (2.0 * DtWt) - 2. / 3. * div + An_tt(rl, thetal);
                A_tp[i][j][k] = (DtWp + DpWt) + An_tp(rl, thetal);
                A_pp[i][j][k] = (2.0 * DpWp) - 2. / 3. * div + An_pp(rl, thetal);
                //if (i == N_g && j == N_g && k == N_g) 
                //  cout << " reality check: "
                //       << A_rr(i,j,k) + A_tt(i,j,k) + A_pp(i,j,k) << endl;
                A2[i][j][k] = A_rr(i, j, k) * A_rr(i, j, k)
                    + 2.0 * A_rt(i, j, k) * A_rt(i, j, k)
                    + 2.0 * A_rp(i, j, k) * A_rp(i, j, k)
                    + A_tt(i, j, k) * A_tt(i, j, k)
                    + 2.0 * A_tp(i, j, k) * A_tp(i, j, k)
                    + A_pp(i, j, k) * A_pp(i, j, k);
            }
        }
    }
    for (int j = N_g; j < n_theta - N_g; j++) {
      for (int k = N_g; k < n_phi - N_g; k++) {
        for (int i = n_r - N_g; i < n_r; i++) {
          const double rl = rho.r(i);
          A_rr[i][j][k] = (A_rr[n_r - N_g - 1][j][k] * rho.r(n_r - N_g - 1) * rho.r(n_r - N_g - 1)) / (rl * rl);
          A_rt[i][j][k] = (A_rt[n_r - N_g - 1][j][k] * rho.r(n_r - N_g - 1) * rho.r(n_r - N_g - 1)) / (rl * rl);
          A_tt[i][j][k] = (A_tt[n_r - N_g - 1][j][k] * rho.r(n_r - N_g - 1) * rho.r(n_r - N_g - 1)) / (rl * rl);
          A_rp[i][j][k] = (A_rp[n_r - N_g - 1][j][k] * rho.r(n_r - N_g - 1) * rho.r(n_r - N_g - 1)) / (rl * rl);
          A_pp[i][j][k] = (A_pp[n_r - N_g - 1][j][k] * rho.r(n_r - N_g - 1) * rho.r(n_r - N_g - 1)) / (rl * rl);
          A2[i][j][k] = A_rr(i, j, k) * A_rr(i, j, k)
                    + 2.0 * A_rt(i, j, k) * A_rt(i, j, k)
                    + 2.0 * A_rp(i, j, k) * A_rp(i, j, k)
                    + A_tt(i, j, k) * A_tt(i, j, k)
                    + 2.0 * A_tp(i, j, k) * A_tp(i, j, k)
                    + A_pp(i, j, k) * A_pp(i, j, k);
        }
      }
    }
  A_rr.fill_ghosts();
  A_rt.fill_ghosts();
  A_tt.fill_ghosts();
  A_rp.fill_ghosts();
  A_tp.fill_ghosts();
  A_pp.fill_ghosts();
  A2.fill_ghosts();
}


//================================================
  // rescaled spherical polar components of \hat A_{ij} for m=0 version of model (A)
  // (compare (3.7) in Shibata & Nakamura) ONLY FAMILY B HAS SIGMA IMPLEMENTED
  //================================================
  inline double An_rr(double r, double theta) {
    const double r2 = r*r;
    const double sigma2 = sigma * sigma;
      if (nakamura_type == 0) {
	const double costheta = cos(theta);
	return GW_amp * exp(-r2/2.0) * (1. - 3.0 * costheta * costheta);
      } else if (nakamura_type == 1) {
	const double costheta = cos(theta);
	return GW_amp * exp(-r2/sigma2) * (5. -  2. * r2 / sigma2) *
	  (1. - 3.0 * costheta * costheta);
      }
      else {
        return 0.0;
      }
    };

  inline double An_rt(double r, double theta) {
    const double r2 = r*r;
    const double sigma2 = sigma * sigma;
      if (nakamura_type == 0) {
	return GW_amp * exp(-r2/2.0) * (3.0 - r2) * sin(theta) * cos(theta);
      } else if (nakamura_type == 1) {
	return GW_amp * exp(-r2/sigma2) * (15.0 - 20.*r2/sigma2 + 4.*r2*r2/(sigma2*sigma2)) *
	  sin(theta) * cos(theta);
      }
      else {
        return 0.0;
      }
  };

  inline double An_rp(double r, double theta) {
    return 0.0;
  }
  inline double An_tt(double r, double theta) {
    const double r2 = r*r;
    const double r4 = r2*r2;
      if (nakamura_type == 0) {
	const double costheta = cos(theta);
	const double sintheta = sin(theta);
	return - GW_amp * exp(-r2/2.0) *
	  (2 - 6.*costheta*costheta + (6 - 8*r2 + r4)*sintheta*sintheta) / 4.0;
      } else if (nakamura_type == 1) {
        const double costheta = cos(theta);
      const double cos2theta = costheta*costheta;
      const double r6 = r4*r2;
      const double sigma2 = sigma * sigma;
      const double sigma4 = sigma2 * sigma2;
      const double sigma6 = sigma4 * sigma2;
      return (GW_amp / 8.0) * exp(-r2/sigma2) *
	(-80. + 256.*r2/sigma2 -136.*r4/sigma4 + 16.*r6/sigma6 +
	 2.*(60. - 136.*r2/sigma2 + 68.*r4/sigma4 - 8.*r6/sigma6) * cos2theta);
      }
      else {
        return 0.0;
      }
  };
  inline double An_tp(double r, double theta) {
    return 0.0;
  }
  inline double An_pp(double r, double theta) {
    const double r2 = r*r;
    const double r4 = r2*r2;
      if (nakamura_type == 0) {
	const double costheta = cos(theta);
	const double sintheta = sin(theta);
	return GW_amp * exp(-r2/2.0) *
	  (- 2 + 6.*costheta*costheta +
	   (6 - 8*r2 + r4)*sintheta*sintheta) / 4.0;
      } else if (nakamura_type == 1) {
        const double costheta = cos(theta);
        const double sigma2 = sigma * sigma;
        const double sigma4 = sigma2 * sigma2;
        const double sigma6 = sigma4 * sigma2;
        const double cos2theta = costheta * costheta;
        const double r6 = r4*r2;
	return (GW_amp / 8.) * exp(-r2/sigma2) *
	  ( 40. - 240.*r2/sigma2 + 136.*r4/sigma4 - 16.*r6/sigma6 +
	   (4.*r2/sigma2)*(56 - 34.*r2/sigma2 + 4.*r4/sigma4) * cos2theta );      
      }
      else
      {
        return 0.0;
      }
  };