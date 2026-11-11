#pragma once
#include <map>
#include <string>
#include <vector>
namespace thermo {
// Restriction: neutral, dilute, nonreacting components; independent linear
// partition equilibria with fixed phase volumes and temperature.
// Each partition ratio is concentration(phase)/concentration(reference phase).
// Concentration units: mol/m3 of each phase. Gas concentration is mol/m3 gas.
struct PartitionPhase {
    std::string id;
    double volumeM3=0.0;
};
struct PartitionSolute {
    std::string id;
    double totalMoles=0.0;
    std::map<std::string,double> concentrationRatioToReference;
};
struct PartitionSolution {
    // Solute -> phase -> moles and concentrations (mol/m3)
    std::map<std::string,std::map<std::string,double>> moles;
    std::map<std::string,std::map<std::string,double>> concentrationMolM3;
};
PartitionSolution solveLinearPhasePartition(
    const std::vector<PartitionPhase>& phases,
    const std::string& referencePhase,
    const std::vector<PartitionSolute>& solutes);
} // namespace thermo
