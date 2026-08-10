#include "Container.h"

Manager* Container::manager = nullptr;
Grid* Container::grid = nullptr;
Cosmology* Container::cosmology = nullptr;
EOS* Container::eos = nullptr;
InData* Container::indata = nullptr;
Slicing* Container::slicing = nullptr;
Gauge* Container::gauge = nullptr;
dumper* Container::dump = nullptr;
Monitor* Container::monitor = nullptr;
Profiles* Container::profiles = nullptr;
WaveExtraction* Container::waves = nullptr;
Photons* Container::photons = nullptr;
CheckPoint* Container::checkpoint = nullptr;