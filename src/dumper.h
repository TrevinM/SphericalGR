// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing all the dumping stuff
//================================================
//

#ifndef DUMPER_H
#define DUMPER_H

class InData;
class Slicing;
class Gauge;
#include <ctime>
#include "nr3.h"
#include "gridfunction.h"
//#include "InData.h"
//#include "Slicing.h"
//#include "Gauge.h"


// #define _DEBUG_


class dumper {
private:
  int dump_step;
  InData *indata;
  Slicing *slicing;
  Gauge *gauge;
  double PI;
public:
  //================================================
  // Constructor
  //================================================
  dumper(int step, InData *indata_i, Slicing *slicing_i, Gauge *gauge_i):
    dump_step(step), indata(indata_i), slicing(slicing_i), gauge(gauge_i)
  {
    PI = acos(-1.0);
  };
  //================================================
  // Destructor
  //================================================
  ~dumper() {};
  //================================================
  // reset dump_step
  //================================================
  int reset_dumpstep(int step) {
    dump_step = step;
    return dump_step;
  }
  //================================================
  // "Conditional dump"
  //================================================
  bool time_to_dump(int timestep) { return (timestep % dump_step == 0); };
  void cond_dump(double phys_time, double prop_time, int timestep, gf3d *fct,
		 const char * suffix = "") {
    // time to dump?      
    if (time_to_dump(timestep) || strcmp(suffix,""))  
      dump(phys_time, prop_time, timestep, fct, suffix);
  };
  //================================================
  // forced dump
  //================================================
  void dump(double phys_time, double prop_time, int timestep, gf3d *fct, 
		  const char * suffix = ""); 
  //================================================
  // "Conditional slice" dump
  //================================================
  void cond_slice(double phys_time, double prop_time, int timestep, gf3d * fct,
		  const char * suffix = "") {
    if (time_to_dump(timestep) || strcmp(suffix,""))  
      slice(phys_time,prop_time,timestep,fct,suffix);
  };
  //================================================
  // forced slice dump
  //================================================
  void slice(double phys_time, double prop_time, int timestep, gf3d * fct,
	     const char * suffix = ""); 
};
 
#endif /* DUMPER_H */

