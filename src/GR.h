#include "Manager.h"

#include "Grid.h"
#include "Cosmology.h"
#include "EOS.h"
#include "InData.h"
#include "Slicing.h"
#include "Gauge.h"
#include "dumper.h"
#include "Monitor.h"
#include "Profiles.h"
#include "WaveExtraction.h"
#include "Photons.h"
#include "CheckPoint.h"

class GR {
public:
    static Manager* manager;

    static Grid* grid;
    static Cosmology* cosmology;
    static EOS* eos;
    static InData* indata;
    static Slicing* slicing;
    static Gauge* gauge;
    static dumper* dump;
    static Monitor* monitor;
    static Profiles* profiles;
    static WaveExtraction* waves;
    static Photons* photons;
    static CheckPoint* checkpoint;
};