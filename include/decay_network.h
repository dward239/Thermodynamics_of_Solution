#pragma once
#include <map>
#include <string>
#include <vector>
namespace thermo {
// One parent nucleus transforms to exactly one daughter nucleus per event.
// Fission, delayed particle multiplication, irradiation production and external
// source terms are NOT covered by this restricted initial decay module.
struct DaughterBranch {
    std::string daughter;
    double fraction = 0.0;
};
struct NuclideDefinition {
    std::string id;
    std::string element;
    double halfLifeSeconds = 0.0; // 0 denotes stable
    double emittedEnergyJPerDecay = 0.0; // released non-neutrino radiation energy
    std::vector<DaughterBranch> daughters;
    // Sum <=1: remaining emitted radiation energy escapes this modeled system.
    std::map<std::string,double> depositionFractions;
};
struct DecayStep {
    std::map<std::string,double> finalAtoms;
    std::map<std::string,double> decaysByNuclide;
    double emittedEnergyJ = 0.0;
    std::map<std::string,double> depositedEnergyJByPhase;
};
class DecayNetwork {
public:
    explicit DecayNetwork(std::vector<NuclideDefinition> definitions);
    const std::vector<NuclideDefinition>& nuclides() const { return definitions_; }
    // Probability-conserving continuous-time decay propagation by uniformization.
    // Large stiffness is explicitly rejected rather than producing unbounded work.
    DecayStep advance(const std::map<std::string,double>& initialAtoms,
                      double durationSeconds) const;
private:
    std::vector<NuclideDefinition> definitions_;
};
} // namespace thermo
