#ifndef CONTAINER_H
#define CONTAINER_H

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
#include "Tracker.h"

class Container {
public:
    // Evolution, times, states
    inline static Manager* manager = nullptr;

    // Resolution 
    inline static Grid* grid = nullptr;
    inline static Cosmology* cosmology = nullptr;
    inline static EOS* eos = nullptr;
    inline static InData* indata = nullptr;
    inline static Slicing* slicing = nullptr;
    inline static Gauge* gauge = nullptr;
    inline static dumper* dump = nullptr;
    inline static Monitor* monitor = nullptr;
    inline static Profiles* profiles = nullptr;
    inline static WaveExtraction* waves = nullptr;
    inline static Photons* photons = nullptr;
    inline static CheckPoint* checkpoint = nullptr;

    // Location and values of various extrema
    inline static Tracker* tracker = nullptr;

};


#endif