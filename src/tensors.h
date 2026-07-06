// Tell emacs that this is -*-c++-*- mode
//================================================
// classes that handle 3D tensors of different ranks
//================================================
#include "nr3.h"

#ifndef TENSORS_H
#define TENSORS_H

#define CHECKBOUNDS
// #define CHECKSYMMETRY
//================================================
// First: RANK 1
//================================================
class vect {
 private:
  int n;
  double *data;
 public:
  //
  // Constructors
  //
  vect() { 
    n = 3; 
    data = new double[n];  
    for (int i=0; i<n; i++) data[i] = 0.0; 
  };
  // Copy constructor
  vect(const vect & rhs) {
    n = 3;
    data = new double[n];
    for (int i=0; i<n; i++) data[i] = rhs[i];
  };
  // Copy constructor
  vect(const double vx, const double vy, const double vz) {
    n = 3;
    data = new double[n];
    data[0] = vx;
    data[1] = vy;
    data[2] = vz;
  };
  //
  // Destructor
  //
  ~vect() { delete data; };
  //
  // subscript
  //
  inline double & operator[](const int i) {
#ifdef CHECKBOUNDS
    if (i < 0 || i >= n) 
      cerr << "Vector subscript out of bounds in vect! " << endl;
#endif
    return data[i];
  };
  inline const double & operator[](const int i) const {
#ifdef CHECKBOUNDS
    if (i < 0 || i >= n) 
      cerr << "Vector subscript out of bounds in vect! " << endl;
#endif
    return data[i];
  };

  void print() {
    cout << " (" << data[0] << "," << data[1] << "," << data[2] << ") " << endl;
  };
  //
  // Copy operator
  //
  vect & operator=(const vect & rhs) {
    for (int i=0; i<n; i++) data[i] = rhs[i];
    return *this;
  };
};

//
//================================================
// Now: RANK 2
//================================================
//
class tensor {
 private:
  int n;
  int nn;
  double **data;
 public:
  //
  // Constructors
  //
  tensor() { 
    n = 3;
    nn = n*n;
    data = new double*[n];
    data[0] = new double[nn];
    for (int i = 1; i<n; i++) data[i] = data[i-1] + n;
    for (int i=0; i<n; i++) 
      for (int j=0; j<n; j++) 
	data[i][j] = 0.0; 
  };
  // Copy constructor
  tensor(const tensor & rhs) {
    n = 3;
    nn = n*n;
    data = new double*[n];
    data[0] = new double[nn];
    for (int i = 1; i<n; i++) data[i] = data[i-1] + n;
    for (int i=0; i<n; i++) 
      for (int j=0; j<n; j++) 
	data[i][j] = rhs[i][j]; 
  };
  // constructor for symmetric tensor
  tensor(const double hxx, const double hxy, const double hxz,
	 const double hyy, const double hyz, const double hzz) {
    n = 3;
    nn = n*n;
    data = new double*[n];
    data[0] = new double[nn];
    for (int i = 1; i<n; i++) data[i] = data[i-1] + n;
    data[0][0] = hxx;
    data[0][1] = data[1][0] = hxy;
    data[0][2] = data[2][0] = hxz;
    data[1][1] = hyy;
    data[2][1] = data[1][2] = hyz;
    data[2][2] = hzz;
  };
  // constructor for general tensor
  tensor(const double hxx, const double hxy, const double hxz,
	 const double hyx, const double hyy, const double hyz,
	 const double hzx, const double hzy, const double hzz) {
    n = 3;
    nn = n*n;
    data = new double*[n];
    data[0] = new double[nn];
    for (int i = 1; i<n; i++) data[i] = data[i-1] + n;
    data[0][0] = hxx;
    data[0][1] = hxy;
    data[0][2] = hxz;
    data[1][0] = hyx;
    data[1][1] = hyy;
    data[1][2] = hyz;
    data[2][0] = hzx;
    data[2][1] = hzy;
    data[2][2] = hzz;
  };
  //
  // Destructor
  //
  ~tensor() { delete data[0]; delete data; };
  //
  // subscript
  //
  inline double * operator[](const int i) {  // returns pointer!
#ifdef CHECKBOUNDS
    if (i < 0 || i >= n) 
      cerr << "Vector subscript out of bounds in tensor! " << endl;
#endif
    return data[i];
  };
  inline const double * operator[](const int i) const {
#ifdef CHECKBOUNDS
    if (i < 0 || i >= n) 
      cerr << "Vector subscript out of bounds in tensor! " << endl;
#endif
    return data[i];
  };
  //
  // addition
  //
  inline tensor operator +(const tensor & A2) {
    return tensor(data[0][0] + A2[0][0],
		  data[0][1] + A2[0][1],
		  data[0][2] + A2[0][2],
		  data[1][0] + A2[1][0],
		  data[1][1] + A2[1][1],
		  data[1][2] + A2[1][2],
		  data[2][0] + A2[2][0],
		  data[2][1] + A2[2][1],
		  data[2][2] + A2[2][2]);
  };
  //
  // multiplication with scalar
  //
  inline tensor operator *(const double factor) {
    return tensor(factor * data[0][0],
		  factor * data[0][1],		  
		  factor * data[0][2],
		  factor * data[1][0],		  
		  factor * data[1][1],
		  factor * data[1][2],
		  factor * data[2][0],		  
		  factor * data[2][1],
		  factor * data[2][2]);
  };
  //
  // print
  // 
  void print() {
    cout << " ( " << data[0][0] << "," << data[1][0] << "," << data[2][0] << " ) " << endl;
    cout << " ( " << data[0][1] << "," << data[1][1] << "," << data[2][1] << " ) " << endl;
    cout << " ( " << data[0][2] << "," << data[1][2] << "," << data[2][2] << " ) " << endl;
  };
  //
  // Copy operator
  //
  tensor & operator=(const tensor & rhs) {
    for (int i=0; i<n; i++) 
      for (int j=0; j<n; j++) 
  	data[i][j] = rhs[i][j];
    return *this;
  };
  //
  // determinant and inverse
  //
  inline double determinant() {
    return data[0][0] * ( data[1][1] * data[2][2] - data[2][1] * data[1][2] )
      - data[0][1] * ( data[1][0] * data[2][2] - data[1][2] * data[2][0] ) 
      + data[0][2] * ( data[1][0] * data[2][1] - data[1][1] * data[2][0] );
  };
  tensor inverse() {
#ifdef CHECKSYMMETRY
    if (data[0][1] != data[1][0] || data[0][2] != data[2][0] ||
	data[1][2] != data[2][1]) 
      cerr << "Inverse defined only for symmetric tensor! " << endl;

#endif
    double det = determinant();
    double xx =    ( data[1][1] * data[2][2] - data[1][2] * data[2][1] ) / det;
    double xy =  - ( data[1][0] * data[2][2] - data[1][2] * data[2][0] ) / det;
    double xz =    ( data[1][0] * data[2][1] - data[1][1] * data[2][0] ) / det;
    double yy =    ( data[0][0] * data[2][2] - data[0][2] * data[2][0] ) / det;
    double yz =  - ( data[0][0] * data[2][1] - data[0][1] * data[2][0] ) / det; 
    double zz =    ( data[0][0] * data[1][1] - data[0][1] * data[1][0] ) / det;
    tensor inverse(xx,xy,xz,yy,yz,zz);
    return inverse;
  };
  //
  // trace and traceless part
  //
  inline double trace(tensor g) {
    double tr = 0.0;
    for (int i = 0; i<n; i++)
      for (int j = 0; j<n; j++)
	tr += g[i][j]*data[i][j];
    return tr;
  };
  double remove_trace(tensor g, tensor gup) { // assumes taking trace of tensor with UPPER indices
    double tr = trace(g);
    for (int i = 0; i<n; i++)
      for (int j = 0; j<n; j++)
	data[i][j] -= gup[i][j] * tr / n;
    return trace(g);
  };
};

//
//================================================
// Now: RANK 3
//================================================
// Note: assume symmetry on last two indices only!
class rank3tens {
 private:
  int n, nn, nnn;
  double ***data;
 public:
    //
  // Constructors
  //
  rank3tens() { 
    n = 3;
    nn = n*n;
    nnn = n*n*n;
    data = new double**[n];
    data[0] = new double*[nn];
    data[0][0] = new double[nnn];
    for (int j = 1; j<n; j++) data[0][j] = data[0][j-1] + n;
    for (int i = 1; i<n; i++) {
      data[i] = data[i-1] + n;
      data[i][0] = data[i-1][0] + nn;
      for (int j = 1; j<n; j++) 
	data[i][j] = data[i][j-1] + n; 
    }
    for (int i=0; i<n; i++) 
      for (int j=0; j<n; j++) 
	for (int k=0; k<n; k++)
	  data[i][j][k] = 0.0;
  };
  // Copy constructor
  rank3tens(const rank3tens & rhs) {
    n = 3;
    nn = n*n;
    nnn = n*n*n;
    data = new double**[n];
    data[0] = new double*[nn];
    data[0][0] = new double[nnn];
    for (int j = 1; j<n; j++) data[0][j] = data[0][j-1] + n;
    for (int i = 1; i<n; i++) {
      data[i] = data[i-1] + n;
      data[i][0] = data[i-1][0] + nn;
      for (int j = 1; j<n; j++) 
	data[i][j] = data[i][j-1]; 
    }
    for (int i=0; i<n; i++) 
      for (int j=0; j<n; j++) 
	for (int k=0; k<n; k++)
	  data[i][j][k] = rhs[i][j][k];
  };
  // Copy constructor
  rank3tens(const double Dxhxx, const double Dxhxy, const double Dxhxz,
	    const double Dxhyy, const double Dxhyz, const double Dxhzz,
	    const double Dyhxx, const double Dyhxy, const double Dyhxz,
	    const double Dyhyy, const double Dyhyz, const double Dyhzz,
	    const double Dzhxx, const double Dzhxy, const double Dzhxz,
	    const double Dzhyy, const double Dzhyz, const double Dzhzz) {
    n = 3;
    nn = n*n;
    nnn = n*n*n;
    data = new double**[n];
    data[0] = new double*[nn];
    data[0][0] = new double[nnn];
    for (int j = 1; j<n; j++) data[0][j] = data[0][j-1] + n;
    for (int i = 1; i<n; i++) {
      data[i] = data[i-1] + n;
      data[i][0] = data[i-1][0] + nn;
      for (int j = 1; j<n; j++) 
	data[i][j] = data[i][j-1] + n; 
    }

    data[0][0][0] = Dxhxx;
    data[0][0][1] = data[0][1][0] = Dxhxy;
    data[0][0][2] = data[0][2][0] = Dxhxz;
    data[0][1][1] = Dxhyy;
    data[0][2][1] = data[0][1][2] = Dxhyz;
    data[0][2][2] = Dxhzz;

    data[1][0][0] = Dyhxx;
    data[1][0][1] = data[1][1][0] = Dyhxy;
    data[1][0][2] = data[1][2][0] = Dyhxz;
    data[1][1][1] = Dyhyy;
    data[1][2][1] = data[1][1][2] = Dyhyz;
    data[1][2][2] = Dyhzz;

    data[2][0][0] = Dzhxx;
    data[2][0][1] = data[2][1][0] = Dzhxy;
    data[2][0][2] = data[2][2][0] = Dzhxz;
    data[2][1][1] = Dzhyy;
    data[2][2][1] = data[2][1][2] = Dzhyz;
    data[2][2][2] = Dzhzz;
    
  };
  //
  // subscript
  //
  inline double ** operator[](const int i) {  // returns pointer!
#ifdef CHECKBOUNDS
    if (i < 0 || i >= n) 
      cerr << "Vector subscript out of bounds in rank3tens! " << endl;
#endif
    return data[i];
  };
  inline const double * const * operator[](const int i) const {
#ifdef CHECKBOUNDS
    if (i < 0 || i >= n) 
      cerr << "Vector subscript out of bounds in rank3tens! " << endl;
#endif
    return data[i];
  };
  //
  // Destructor
  //
  ~rank3tens() { delete data[0][0]; delete data[0]; delete data; };
}; 



#endif   // TENSORS_H
