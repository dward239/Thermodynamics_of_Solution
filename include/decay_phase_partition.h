#pragma once
#include "decay_network.h"
#include "phase_partition.h"
#include <map>
#include <string>
#include <vector>
namespace thermo {
// Experimental narrow bridge: treats each element as one neutral,
// nonreactive pseudo-solute distributed by supplied phase coefficients.
// It is NOT actual isotope-specific aqueous chemical speciation.
struct DecayPartitionResult {
    DecayStep decay;
    std::map<std::string,double> elementTotalMol;
    PartitionSolution phaseDistribution;
};
DecayPartitionResult advanceDecayAndLinearRepartition(
    const DecayNetwork& network,
    const std::map<std::string,double>& initialNuclideAtoms,
    double dtSeconds,
    const std::vector<PartitionPhase>& phases,
    const std::string& referencePhase,
    const std::map<std::string,std::map<std::string,double>>& elementPartitionRatios);
} // namespace thermo
