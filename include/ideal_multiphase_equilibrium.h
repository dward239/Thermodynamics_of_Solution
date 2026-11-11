#pragma once
#include <map>
#include <string>
#include <vector>
namespace thermo {
// Restricted general mass-action equilibrium: neutral IDEAL-DILUTE species.
// Species may inhabit distinct fixed-volume liquids or a gas headspace.
// Base solve excludes solids; optional bounded pure-solid extension below.
// No charges, liquid-liquid activity coupling, or phase creation.
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
    // Component Lagrange multipliers divided by RT (natural-log standard).
    std::map<std::string,double> componentPotentialOverRT;
    int iterations=0;
    double maxNormalizedResidual=0.0;
};
IdealEquilibriumResult solveIdealMultiphase(
    const IdealMultiphaseModel& model,
    int maxIterations=150,
    double relativeTolerance=1e-10);
// Pure-solid candidates initially absent; no solid solutions or seed inventories.
// At equilibrium: present solid affinity ~0, absent solid affinity <=0.
struct IdealPureSolid {
    std::string id;
    std::map<std::string,int> componentStoichiometry;
    double standardMuOverRT=0.0;
};
struct IdealSolidEquilibriumResult {
    IdealEquilibriumResult fluidSolution;
    std::map<std::string,double> solidMoles;
    // Affinity = sum_i nu_is * lambda_i/(RT) - mu_s^0/(RT).
    // Positive affinity denotes supersaturation; equilibrium requires <=0.
    std::map<std::string,double> solidAffinityOverRT;
};
// Enumerates candidate active-solid subsets. Exponential in solids count,
// so <=8 candidates. Experimental neutral, ideal-dilute, fixed-volume model.
IdealSolidEquilibriumResult solveIdealMultiphaseWithSolids(
    const IdealMultiphaseModel& fluidModel,
    const std::vector<IdealPureSolid>& solidCandidates,
    int maxIterations=150,
    double relativeTolerance=1e-10);
} // namespace thermo
