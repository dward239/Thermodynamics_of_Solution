#pragma once
#include <map>
#include <string>
#include <vector>
namespace thermo {
// Structural domain model: a phase and a physical compartment are distinct concepts.
enum class PhaseKind { AqueousLiquid, OrganicLiquid, Gas, Solid };
struct PhaseDefinition {
    std::string id;
    PhaseKind kind = PhaseKind::AqueousLiquid;
    double volumeM3 = 0.0; // occupies physical volume; not its container volume
    double temperatureK = 298.15;
};
struct ChemicalSpecies {
    std::string id;
    std::string phaseId;
    int charge = 0;
    // Stoichiometric counts of conserved chemical components (often elements).
    std::map<std::string, int> composition;
};
struct EquilibriumReaction {
    std::string id;
    // Reaction convention: reactants negative, products positive.
    std::map<std::string, int> coefficients;
    double logKReference = 0.0;
};
struct MultiphaseSystem {
    std::vector<PhaseDefinition> phases;
    std::vector<ChemicalSpecies> species;
    std::vector<EquilibriumReaction> reactions;
    // Structural, charge and component-conservation validation only.
    // NOT a multiphase equilibrium or EOS solver.
    void validate() const;
};
} // namespace thermo
