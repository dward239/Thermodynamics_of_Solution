#pragma once
#include "decay_network.h"
#include "lumped_thermal.h"
#include "radiolysis_sources.h"
#include <map>
#include <string>
#include <vector>
namespace thermo {
struct VesselState {
    double timeSeconds=0.0;
    double temperatureK=298.15;
    std::map<std::string,double> nuclideAtoms;
    double headspaceGasMol=0.0;
};
struct VesselModel {
    double headspaceVolumeM3=0.0;
    double sensibleHeatFraction=1.0;
    LumpedThermalProperties thermal;
    std::vector<RadiolysisYield> radiolysisYields;
};
struct VesselStepResult {
    VesselState state;
    DecayStep decay;
    std::map<std::string,double> grossRadiolysisMolByPhaseSpecies;
    double headspacePressurePa=0.0;
};
VesselStepResult advanceClosedVessel(const VesselState& initial,
    double dtSeconds,const DecayNetwork& network,const VesselModel& model,
    double externallySpecifiedNetGasAddedMol=0.0);
} // namespace thermo
