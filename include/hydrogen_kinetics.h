#pragma once

namespace thermo {

// Fixed-coefficient, one-species (H2) approximation. Two reservoirs:
// dissolved H2 and H2 in an ideal-gas headspace.
// Gross radiolytic H2 is added to dissolved inventory; a first-order
// consumption sink explicitly removes it from the modeled H2 inventory.
// Liquid/gas transfer has a user-supplied equilibrium *mole ratio*.
// This is NOT a chemical reaction network or Henry-law parameter database.
struct HydrogenKineticsParameters {
    double transferRatePerSecond=0.0;
    double dissolvedConsumptionRatePerSecond=0.0;
    // (moles H2 in liquid)/(moles H2 in gas) at interphase equilibrium.
    // For fixed volumes, ratio = (c_liquid/c_gas)*(V_liquid/V_gas).
    double equilibriumDissolvedToGasMoleRatio=0.0;
    unsigned substeps=32; // backward Euler, 1..100000
};

struct HydrogenKineticsState {
    double dissolvedMol=0.0;
    double headspaceMol=0.0;
    double cumulativeConsumedMol=0.0;
};

struct HydrogenKineticsStep {
    HydrogenKineticsState state;
    double grossGeneratedMol=0.0;
    double consumedDuringStepMol=0.0;
    double netTransferToHeadspaceMol=0.0;
    // initial inventories + gross - final inventories - consumed
    double closureResidualMol=0.0;
};

// Production is distributed uniformly across dt; the fraction produced
// and transferred within the step depends on temporal resolution.
// Positivity-preserving implicit substeps and a conservation audit.
HydrogenKineticsStep advanceHydrogenKinetics(
    const HydrogenKineticsState& initial,
    double dtSeconds,
    double grossH2GeneratedMol,
    const HydrogenKineticsParameters& parameters);

} // namespace thermo
