#include "hydrogen_kinetics.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace thermo {
namespace {
void finiteNonnegative(double value,const char* name) {
    if(!std::isfinite(value)||value<0.0)
        throw std::invalid_argument(name);
}
}
HydrogenKineticsStep advanceHydrogenKinetics(
    const HydrogenKineticsState& initial,
    double dtSeconds,
    double grossH2GeneratedMol,
    const HydrogenKineticsParameters& p) {
    finiteNonnegative(initial.dissolvedMol,"Invalid initial dissolved hydrogen");
    finiteNonnegative(initial.headspaceMol,"Invalid initial headspace hydrogen");
    finiteNonnegative(initial.cumulativeConsumedMol,"Invalid cumulative consumption");
    finiteNonnegative(dtSeconds,"Invalid hydrogen kinetics time step");
    finiteNonnegative(grossH2GeneratedMol,"Invalid gross H2 generation");
    finiteNonnegative(p.transferRatePerSecond,"Invalid interphase transfer rate");
    finiteNonnegative(p.dissolvedConsumptionRatePerSecond,"Invalid H2 consumption rate");
    finiteNonnegative(p.equilibriumDissolvedToGasMoleRatio,"Invalid equilibrium mole ratio");
    if(p.substeps<1||p.substeps>100000)
        throw std::invalid_argument("Hydrogen kinetics substeps must be 1 to 100000");
    if(dtSeconds==0 && grossH2GeneratedMol!=0)
        throw std::invalid_argument("H2 generation during zero time");

    HydrogenKineticsStep out;
    out.state=initial;
    out.grossGeneratedMol=grossH2GeneratedMol;
    const double initialTotal=initial.dissolvedMol+initial.headspaceMol;
    if(!std::isfinite(initialTotal)||!std::isfinite(initialTotal+grossH2GeneratedMol))
        throw std::overflow_error("Hydrogen inventory overflow");
    if(dtSeconds==0)return out;

    const double step=dtSeconds/static_cast<double>(p.substeps);
    const double generatedStep=grossH2GeneratedMol/p.substeps;
    // Backward-Euler implicit linear 2x2 system:
    // dL/dt = q - k_c*L - k_t*(L-r*G)
    // dG/dt = k_t*(L-r*G)
    // with q constant across the time step.
    const double a=1.0+step*(p.dissolvedConsumptionRatePerSecond+
                             p.transferRatePerSecond);
    const double b=-step*p.transferRatePerSecond*
                    p.equilibriumDissolvedToGasMoleRatio;
    const double c=-step*p.transferRatePerSecond;
    const double d=1.0+step*p.transferRatePerSecond*
                    p.equilibriumDissolvedToGasMoleRatio;
    // Expand determinant analytically to avoid cancellation when transfer
    // rates are much faster than the time step.
    const double hkc=step*p.dissolvedConsumptionRatePerSecond;
    const double hkt=step*p.transferRatePerSecond;
    const double r=p.equilibriumDissolvedToGasMoleRatio;
    const double det=1.0+hkc+hkt+hkt*r+hkc*hkt*r;
    if(!std::isfinite(a)||!std::isfinite(b)||!std::isfinite(c)||
       !std::isfinite(d)||!std::isfinite(det)||det<=0.0)
        throw std::invalid_argument("Hydrogen kinetics time step/rates overflow");

    for(unsigned i=0;i<p.substeps;++i) {
        const double rhsL=out.state.dissolvedMol+generatedStep;
        const double rhsG=out.state.headspaceMol;
        const double liquid=(rhsL*d-b*rhsG)/det;
        const double gas=(a*rhsG-c*rhsL)/det;
        const double consumed=step*p.dissolvedConsumptionRatePerSecond*liquid;
        if(!std::isfinite(liquid)||!std::isfinite(gas)||
           !std::isfinite(consumed)||liquid<0||gas<0||consumed<0)
            throw std::runtime_error("Nonphysical hydrogen kinetics result");
        out.state.dissolvedMol=liquid;
        out.state.headspaceMol=gas;
        out.state.cumulativeConsumedMol+=consumed;
        out.consumedDuringStepMol+=consumed;
        if(!std::isfinite(out.state.cumulativeConsumedMol)||
           !std::isfinite(out.consumedDuringStepMol))
            throw std::overflow_error("Hydrogen consumption integral overflow");
    }
    out.netTransferToHeadspaceMol=
        out.state.headspaceMol-initial.headspaceMol;
    out.closureResidualMol=initialTotal+grossH2GeneratedMol-
        out.state.dissolvedMol-out.state.headspaceMol-out.consumedDuringStepMol;
    if(!std::isfinite(out.closureResidualMol)||
       std::abs(out.closureResidualMol)>
       1e-9*std::max(1.0,initialTotal+grossH2GeneratedMol))
        throw std::runtime_error("Hydrogen molecular-inventory audit failed");
    return out;
}
} // namespace thermo
