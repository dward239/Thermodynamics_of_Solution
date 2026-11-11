#pragma once
#include <map>
#include <string>
#include <vector>
namespace thermo {
// Restricted general mass-action equilibrium: neutral IDEAL-DILUTE species.
// Species may inhabit distinct fixed-volume liquids or a gas headspace.
// No solids, charges, liquid-liquid activity coupling, or phase creation.
struct IdealComponent {
    std::string id;
    double totalMol=0.0;
};
struct IdealPhase {
    std::string id;
    double volumeM3=0.0;
    double standardConcentrationMolM3=1000.0;
};
struct IdealSpecies {
    std::string id;
    std::string phaseId;
    std::map<std::string,int> componentStoichiometry;
    double standardMuOverRT=0.0;
};
struct IdealMultiphaseModel {
    std::vector<IdealComponent> components;
    std::vector<IdealPhase> phases;
    std::vector<IdealSpecies> species;
};
struct IdealEquilibriumResult {
    std::map<std::string,double> speciesMoles;
    std::map<std::string,double> concentrationsMolM3;
    std::map<std::string,double> componentResidualMol;
    int iterations=0;
    double maxNormalizedResidual=0.0;
};
IdealEquilibriumResult solveIdealMultiphase(
    const IdealMultiphaseModel& model,
    int maxIterations=150,
    double relativeTolerance=1e-10);
} // namespace thermo
