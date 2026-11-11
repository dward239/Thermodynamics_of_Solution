#pragma once
#include "decay_network.h"
#include "multiphase_system.h"
#include <map>
#include <string>
namespace thermo {
// Convert nuclear atom inventory into chemical-element mole totals.
// Nuclides with identical chemical element labels aggregate into one total.
// Does not assign isotopes to molecular species, phases or oxidation states.
std::map<std::string,double> elementalMolesFromNuclides(
    const DecayNetwork& network,
    const std::map<std::string,double>& nuclideAtoms);
// Verify that elemental inventory names are represented by available
// chemical species compositions before attempting any chemistry coupling.
void validateElementCoverage(const MultiphaseSystem& system,
    const std::map<std::string,double>& elementalMoles);
} // namespace thermo
