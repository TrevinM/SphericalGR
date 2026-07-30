// Tell emacs that this is -*-c++-*- mode
//================================================
// Classes for different choices of the Bona-Masso function f(alpha)
// (see Baumgarte & de Oliveira, 2022)
//================================================
#ifndef BONAMASSO_H
#define BONAMASSO_H

//================================================
// Base class
//================================================

class BonaMasso_f {
protected:
  ostringstream bonamasso_name;
public:
  BonaMasso_f() {};
  virtual ~BonaMasso_f() {};
  //================================================
  // virtual function f(alpha)
  //================================================
  virtual double operator()(double alpha) = 0;
  virtual string Name() = 0;
};

//================================================
// generalized 1+log
//================================================
class GenOnePlusLog : public BonaMasso_f {
private:
  double k;
public:
  GenOnePlusLog(double par) : BonaMasso_f(), k(par) {
    cout << " BONAMASSO: setting up f(alpha) for generalized one_plus_log..." << endl;
    cout << " BONAMASSO: using parameter k = " << k << endl;
    bonamasso_name << "generalized one-plus_log with k = " << k;
  }
  ~GenOnePlusLog() {};
  double operator()(double alpha) { return k / alpha; };
  string Name() { return bonamasso_name.str(); };
};
    
//================================================
// Ken's Bona-Masso function
//================================================
class Ken : public BonaMasso_f {
public:
  Ken() : BonaMasso_f() {
    cout << " BONAMASSO: setting up f(alpha) for Ken's trumpet slicing..." << endl;
    bonamasso_name << "Ken's trumpet slicing" ;
  }
  ~Ken() {};
  double operator()(double alpha) { return (1.0 - alpha) / alpha; };
  string Name() { return bonamasso_name.str(); };
};

//================================================
// Gauge-Shock avoiding (see Alcubierre, 1997)
//================================================
class GaugeShockAvoid : public BonaMasso_f {
private:
  double kappa;
public:
  GaugeShockAvoid(double par) : BonaMasso_f(), kappa(par) {
    cout << " BONAMASSO: setting up f(alpha) for gauge-shock-avoiding slices..." << endl;
    cout << " BONAMASSO: using parameter kappa = " << kappa << endl;
    bonamasso_name << "gauge-shock-avoiding f with kappa = " << kappa ;
  }
  ~GaugeShockAvoid() {};
  double operator()(double alpha) { return 1.0 + kappa / (alpha*alpha); };
  string Name() { return bonamasso_name.str(); };
};

//================================================
// linear Gauge-Shock avoiding (see Alcubierre, 2003)
//================================================
class GaugeShockAvoid_lin : public BonaMasso_f {
private:
  double a_0;
public:
  GaugeShockAvoid_lin(double par) : BonaMasso_f(), a_0(par) {
    cout << " BONAMASSO: setting up f(alpha) for linear gauge-shock-avoiding slices..." << endl;
    cout << " BONAMASSO: using parameter a_0 = " << a_0 << endl;
    bonamasso_name << "linear gauge-shock-avoiding f with a_0 = " << a_0 ;
  }
  ~GaugeShockAvoid_lin() {};
  double operator()(double alpha) { 
    return a_0*a_0 / (2.*alpha + (a_0 - 2.0)*alpha*alpha); 
  }
  string Name() { return bonamasso_name.str(); };
};

//================================================
// Cosmo - Brady
//================================================
class Cosmo : public BonaMasso_f {
private:
  double k;
public:
  Cosmo(double par) : BonaMasso_f(), k(par) {
    cout << " BONAMASSO: setting up f(alpha) for Cosmo slicing..." << endl;
    cout << " BONAMASSO: using barameter k = " << k << endl;
    bonamasso_name << "Cosmo" ;
  }
  ~Cosmo() {};
  double operator()(double alpha) { return k * exp(-alpha) / alpha; };
  string Name() { return bonamasso_name.str(); };
};

#endif /* BONAMASSO_H */
