// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Parameters for asymptotic Cosmological models
//
//================================================

#ifndef COSMOLOGY_H
#define COSMOLOGY_H

enum { minkowski, desitter, radiation };
//
//================================================
// base class - doesn't do much... 
//================================================
//
class Cosmology {
protected:
  double PI;
  double current_time;
public:
  Cosmology() { PI = acos(-1.0); current_time = 0.0; };
  virtual ~Cosmology() {};
  //
  double a() { return a(current_time); }
  virtual double a(double time) = 0;
  double Hubble() { return Hubble(current_time); }
  virtual double Hubble(double time) = 0;
  double K0() { return - 3.0*Hubble(current_time); }
  double K0(double time) { return - 3.0*Hubble(time); };
  virtual double Lambda() = 0;
  double rho0() { return rho0(current_time); }
  virtual double rho0(double time) = 0;
  virtual const char * Name() = 0; 
  double update_time(double dt) { current_time += dt; return current_time; }
  double time() { return current_time; }
  virtual int type() = 0;
};
//
//================================================
// Minkowski
//================================================
//
class Minkowski : public Cosmology {
public:
  Minkowski() : Cosmology() {};
  ~Minkowski() {};
  //
  double a(double time) { return 1.0; };
  double Hubble(double time) { return 0.0; };
  double Lambda() { return 0.0; }
  double rho0(double time) { return 0.0; }
  const char * Name() { return "Minkowski"; }; 
  int type() { return minkowski; };
};
//
//================================================
// de Sitter
//================================================
//
class DeSitter : public Cosmology {
private:
  double lambda;
  double H;
public:
  DeSitter() : Cosmology() { 
    // read input from Cosmology_Input
    ifstream infile;
    infile.open("Cosmology_Input");
    if (!infile) {
      cerr << " Can't open file Cosmology_Input for input  -- setting Lambda to zero!! " << endl;
      lambda = 0.0;
    } else {
      char buf[500],c;
      infile.get(buf,500,'='); infile.get(c); infile >> lambda;
      cout << " COSMOLOGY: Read lambda = " << lambda << " from Cosmology_Input" << endl;
      cout << "===================================================" << endl;
 
    }
    H = sqrt(lambda/3.0); 
  };
  ~DeSitter() {};
  //
  double a(double time) { return exp(H*time); };
  double Hubble(double time) { return H; };
  double Lambda() { return lambda; };
  double rho0(double time) { return 0.0; }
  const char * Name() { return "de Sitter"; }; 
  int type() { return desitter; };
};


//
//================================================
// (pure) Radiation
//================================================
//
class Radiation : public Cosmology {
private:
  double t0;
public:
  Radiation() : Cosmology() { 
    // read input from Cosmology_Input
    ifstream infile;
    infile.open("Cosmology_Input");
    if (!infile) {
      cerr << " Can't open file Cosmology_Input for input  -- setting t0 to unity!! " << endl;
      t0 = 1.0;
    } else {
      char buf[500],c;
      double lambda;  // irrelevant - won't be used 
      infile.get(buf,500,'='); infile.get(c); infile >> lambda;
      infile.get(buf,500,'='); infile.get(c); infile >> t0;
      cout << " COSMOLOGY: Read t0 = " << t0 << " from Cosmology_Input" << endl;
      cout << "===================================================" << endl;
 
    }
  };
  ~Radiation() {};
  //
  double a(double time) { return sqrt( ( time + t0)/t0 ); };
  double Hubble(double time) { return 1.0/(2.0*(time + t0)); };
  double Lambda() { return 0.0; };
  double rho0(double time) { return 3.0/(32.0 * PI * (time +t0) * (time + t0)); }
  const char * Name() { return "radiation-dominated Universe"; }; 
    int type() { return radiation; };
};



#endif  /* COSMOLOGY_H */
