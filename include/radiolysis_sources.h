#pragma once
#include <map>
#include <string>
#include <vector>
namespace thermo {
// Gross, user-specified radiation-chemical G value. Net yield requires a
// kinetic model for consumption/recombination and radiation quality.
struct RadiolysisYield {
    std::string phaseId;
    std::string speciesId;
    double moleculesPer100eV = 0.0;
};
std::map<std::string,double> grossRadiolyticProductionMol(
    const std::map<std::string,double>& depositedEnergyJByPhase,
    const std::vector<RadiolysisYield>& yields);
} // namespace thermo
