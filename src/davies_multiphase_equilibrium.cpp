#include "davies_multiphase_equilibrium.h"
#include "aqueous_engine.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

namespace thermo {
namespace {
constexpr double kLn10=2.3025850929940456840179914546844;
constexpr double kMaximumDaviesIonicStrengthMolL=0.1;
double ionicStrength(const IdealMultiphaseModel& model,
                     const IdealEquilibriumResult& result) {
    std::map<std::string,IdealPhase> phases;
    for(const auto& p:model.phases)phases.emplace(p.id,p);
    double value=0.0;
    for(const auto& s:model.species) {
        if(s.charge==0)continue;
        const auto it=phases.find(s.phaseId);
        if(it==phases.end() ||
           it->second.kind!=IdealFluidPhaseKind::AqueousLiquid)
            throw std::invalid_argument("Davies activity requires an aqueous ion phase");
        const double concentrationMolL=
            result.speciesMoles.at(s.id)/(it->second.volumeM3*1000.0);
        value+=0.5*static_cast<double>(s.charge)*s.charge*concentrationMolL;
    }
    if(!std::isfinite(value)||value<0.0 ||
       value>kMaximumDaviesIonicStrengthMolL+1e-10)
        throw std::invalid_argument(
            "Davies activity restricted to ionic strength <= 0.1 mol/L");
    return value;
}
double gammaForCharge(int charge,double ionicStrengthMolL,double temperatureK) {
    if(charge==0)return 1.0;
    const double gamma=std::exp(kLn10*
        aqueous::daviesLogGamma(charge,ionicStrengthMolL,temperatureK));
    if(!std::isfinite(gamma)||gamma<=0.0)
        throw std::runtime_error("Nonphysical Davies coefficient");
    return gamma;
}
}
DaviesMultiphaseResult solveDaviesMultiphase(
    const IdealMultiphaseModel& model,
    const DaviesMultiphaseOptions& options) {
    if(model.electroneutralPhaseId.empty())
        throw std::invalid_argument(
            "Davies coupling requires a designated charge-balanced aqueous phase");
    if(options.maxOuterIterations<1||options.maxOuterIterations>100000 ||
       options.maxInnerIterations<1 ||options.maxInnerIterations>100000 ||
       !std::isfinite(options.ionicStrengthRelativeTolerance)||
       options.ionicStrengthRelativeTolerance<=0||
       options.ionicStrengthRelativeTolerance>=0.01 ||
       !std::isfinite(options.innerRelativeTolerance)||
       options.innerRelativeTolerance<=0||options.innerRelativeTolerance>=0.01||
       !std::isfinite(options.relaxation)||options.relaxation<=0||
       options.relaxation>1.0)
        throw std::invalid_argument("Invalid Davies multiphase convergence options");
    // This function validates the temperature and retains the conservative
    // 273.15-323.15 K domain of the original aqueous Davies implementation.
    (void)aqueous::daviesLogGamma(1,0.0,options.temperatureK);
    // Also validates phase, charge, free component and species definitions.
    auto ideal=solveIdealMultiphase(
        model,options.maxInnerIterations,options.innerRelativeTolerance);
    double currentI=ionicStrength(model,ideal);
    std::map<std::string,IdealPhase> phases;
    for(const auto& p:model.phases)phases.emplace(p.id,p);
    for(int iteration=1;iteration<=options.maxOuterIterations;++iteration) {
        IdealMultiphaseModel corrected=model;
        std::map<std::string,double> gamma;
        for(auto& s:corrected.species) {
            const double coefficient=gammaForCharge(
                s.charge,currentI,options.temperatureK);
            gamma[s.id]=coefficient;
            if(s.charge!=0)
                s.standardMuOverRT+=std::log(coefficient);
        }
        auto state=solveIdealMultiphase(
            corrected,options.maxInnerIterations,options.innerRelativeTolerance);
        const double resultI=ionicStrength(model,state);
        const double denom=std::max(1e-12,std::max(resultI,currentI));
        const double relativeResidual=std::abs(resultI-currentI)/denom;
        if(relativeResidual<=options.ionicStrengthRelativeTolerance) {
            DaviesMultiphaseResult output;
            output.equilibrium=std::move(state);
            output.ionicStrengthMolL=resultI;
            output.relativeIonicStrengthResidual=relativeResidual;
            output.outerIterations=iteration;
            output.activityCoefficients=std::move(gamma);
            for(const auto& s:model.species) {
                const auto& phase=phases.at(s.phaseId);
                const double concentration=output.equilibrium.concentrationsMolM3.at(s.id);
                const double activity=output.activityCoefficients.at(s.id)*
                    concentration/phase.standardConcentrationMolM3;
                if(!std::isfinite(activity)||activity<=0.0)
                    throw std::runtime_error("Nonphysical species activity");
                output.activities[s.id]=activity;
            }
            return output;
        }
        currentI+=(resultI-currentI)*options.relaxation;
        if(!std::isfinite(currentI)||currentI<0.0)
            throw std::runtime_error("Invalid Davies outer ionic strength iterate");
    }
    throw std::runtime_error(
        "Davies multiphase ionic-strength fixed point did not converge");
}
} // namespace thermo
