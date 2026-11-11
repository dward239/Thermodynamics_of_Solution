#include "adaptive_vessel_integration.h"
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

namespace thermo {
namespace {
constexpr double kGasConstant = 8.31446261815324;

void positiveFinite(double value, const char* what) {
    if(!std::isfinite(value) || value <= 0.0)
        throw std::invalid_argument(what);
}

void checkOptions(const AdaptiveVesselOptions& o) {
    positiveFinite(o.initialStepSeconds, "Invalid initial adaptive step");
    positiveFinite(o.minimumStepSeconds, "Invalid minimum adaptive step");
    positiveFinite(o.maximumStepSeconds, "Invalid maximum adaptive step");
    positiveFinite(o.relativeTolerance, "Invalid adaptive relative tolerance");
    positiveFinite(o.absoluteTemperatureToleranceK, "Invalid absolute temperature tolerance");
    positiveFinite(o.absolutePressureTolerancePa, "Invalid absolute pressure tolerance");
    positiveFinite(o.absoluteMolesTolerance, "Invalid absolute mole tolerance");
    positiveFinite(o.absoluteAtomsTolerance, "Invalid absolute atom tolerance");
    if(o.relativeTolerance>=0.1 ||
       o.minimumStepSeconds>o.initialStepSeconds ||
       o.initialStepSeconds>o.maximumStepSeconds ||
       o.maximumAcceptedSteps==0 || o.maximumRejectedSteps==0)
        throw std::invalid_argument("Inconsistent adaptive integration configuration");
}

double totalAtoms(const RadiolyticHydrogenVesselState& s) {
    double sum = 0.0;
    for(const auto& entry : s.vessel.nuclideAtoms)
        sum += entry.second;
    if(!std::isfinite(sum))
        throw std::runtime_error("Nuclide inventory atom-count overflow");
    return sum;
}

double pressurePa(const RadiolyticHydrogenVesselState& state,
                  const RadiolyticHydrogenVesselModel& model) {
    const double p=state.vessel.headspaceGasMol*kGasConstant*
                   state.vessel.temperatureK/model.vessel.headspaceVolumeM3;
    if(!std::isfinite(p)||p<0.0)
        throw std::runtime_error("Invalid headspace pressure during integration");
    return p;
}

double normalizedDifference(double a, double b, double relative,
                            double absolute) {
    if(!std::isfinite(a) || !std::isfinite(b))
        throw std::runtime_error("Nonfinite state for adaptive error estimation");
    const double scale=absolute+relative*std::max(std::abs(a),std::abs(b));
    return std::abs(a-b)/scale;
}

double estimateError(const RadiolyticHydrogenVesselState& full,
                     const RadiolyticHydrogenVesselState& half,
                     const RadiolyticHydrogenVesselModel& model,
                     const AdaptiveVesselOptions& o) {
    double error=0;
    const auto measure=[&](double a,double b,double absTol) {
        error=std::max(error,normalizedDifference(a,b,o.relativeTolerance,absTol));
    };
    measure(full.vessel.temperatureK,half.vessel.temperatureK,
            o.absoluteTemperatureToleranceK);
    measure(full.vessel.headspaceGasMol,half.vessel.headspaceGasMol,
            o.absoluteMolesTolerance);
    measure(full.hydrogen.dissolvedMol,half.hydrogen.dissolvedMol,
            o.absoluteMolesTolerance);
    measure(full.hydrogen.headspaceMol,half.hydrogen.headspaceMol,
            o.absoluteMolesTolerance);
    measure(full.hydrogen.cumulativeConsumedMol,
            half.hydrogen.cumulativeConsumedMol,o.absoluteMolesTolerance);
    measure(pressurePa(full,model),pressurePa(half,model),
            o.absolutePressureTolerancePa);
    if(full.vessel.nuclideAtoms.size()!=half.vessel.nuclideAtoms.size())
        throw std::runtime_error("Nuclide inventory topology changed between trials");
    for(const auto& [name,count]:full.vessel.nuclideAtoms) {
        const auto it=half.vessel.nuclideAtoms.find(name);
        if(it==half.vessel.nuclideAtoms.end())
            throw std::runtime_error("Missing nuclide from accepted half-step");
        measure(count,it->second,o.absoluteAtomsTolerance);
    }
    return error;
}

double deposited(const RadiolyticHydrogenVesselStep& step) {
    double result=0.0;
    for(const auto& entry:step.vesselStep.decay.depositedEnergyJByPhase)
        result+=entry.second;
    if(!std::isfinite(result))
        throw std::runtime_error("Nonfinite deposited energy total");
    return result;
}

bool retryableFailure(const std::invalid_argument& exception) {
    const std::string message=exception.what();
    // Configuration and thermodynamic-domain failures are never retried.
    // Only time-step-specific numerical overflow/stiffness is retryable.
    return message.find("Decay step too stiff")!=std::string::npos ||
           message.find("Hydrogen kinetics time step/rates overflow")!=std::string::npos;
}

} // namespace

AdaptiveVesselTrajectory integrateRadiolyticHydrogenVessel(
    const RadiolyticHydrogenVesselState& initial,
    double endTimeSeconds,
    const DecayNetwork& network,
    const RadiolyticHydrogenVesselModel& model,
    const AdaptiveVesselOptions& options,
    const std::function<void(const AdaptiveVesselSample&)>& acceptedObserver) {
    checkOptions(options);
    if(!std::isfinite(initial.vessel.timeSeconds) ||
       !std::isfinite(endTimeSeconds) ||
       endTimeSeconds<=initial.vessel.timeSeconds)
        throw std::invalid_argument("End time must exceed finite initial time");

    // Zero-duration call checks preconditions without changing inventories.
    const auto validation=advanceRadiolyticHydrogenVessel(initial,0.0,network,model);
    (void)validation;
    AdaptiveVesselTrajectory trajectory;
    RadiolyticHydrogenVesselState accepted=initial;
    AdaptiveVesselSample start;
    start.state=initial;
    start.pressurePa=pressurePa(initial,model);
    if(acceptedObserver)acceptedObserver(start);
    trajectory.history.push_back(start);
    double nextDt=options.initialStepSeconds;

    while(accepted.vessel.timeSeconds < endTimeSeconds) {
        if(trajectory.acceptedSteps >= options.maximumAcceptedSteps)
            throw std::runtime_error("Adaptive integration exceeded accepted-step budget");
        const double remaining=endTimeSeconds-accepted.vessel.timeSeconds;
        if(!std::isfinite(remaining) || remaining<=0)
            throw std::runtime_error("Adaptive integration time grid became invalid");
        double dt=std::min(nextDt,remaining);
        if(accepted.vessel.timeSeconds+dt<=accepted.vessel.timeSeconds)
            throw std::runtime_error("Adaptive time grid cannot progress at floating-point precision");
        // The last fractional step may legitimately be smaller than minStep.
        const double minimumAllowed=std::min(options.minimumStepSeconds,remaining);
        double error=std::numeric_limits<double>::infinity();
        RadiolyticHydrogenVesselStep firstHalf,secondHalf;
        std::string rejectionReason;
        bool computed=false;
        try {
            const auto full=advanceRadiolyticHydrogenVessel(
                accepted,dt,network,model);
            firstHalf=advanceRadiolyticHydrogenVessel(
                accepted,dt/2.0,network,model);
            secondHalf=advanceRadiolyticHydrogenVessel(
                firstHalf.state,dt/2.0,network,model);
            error=estimateError(full.state,secondHalf.state,model,options);
            computed=std::isfinite(error);
            if(!computed)rejectionReason="Nonfinite adaptive error";
        } catch(const std::invalid_argument& e) {
            if(!retryableFailure(e))throw;
            rejectionReason=e.what();
        } catch(const std::runtime_error& e) {
            rejectionReason=e.what();
        }
        if(computed && error<=1.0) {
            const double atomsBefore=totalAtoms(accepted);
            const double atomsAfter=totalAtoms(secondHalf.state);
            const double nucleusClosure=atomsAfter-atomsBefore;
            if(std::abs(nucleusClosure)>
                1e-9*std::max(1.0,atomsBefore))
                throw std::runtime_error("Adaptive nuclide conservation audit failed");
            const double H2Before=accepted.hydrogen.dissolvedMol+
                                  accepted.hydrogen.headspaceMol;
            const double H2After=secondHalf.state.hydrogen.dissolvedMol+
                                 secondHalf.state.hydrogen.headspaceMol;
            const double produced=firstHalf.hydrogenStep.grossGeneratedMol+
                                  secondHalf.hydrogenStep.grossGeneratedMol;
            const double consumed=firstHalf.hydrogenStep.consumedDuringStepMol+
                                  secondHalf.hydrogenStep.consumedDuringStepMol;
            const double closure=H2Before+produced-H2After-consumed;
            if(!std::isfinite(closure)||
               std::abs(closure)>1e-9*std::max(1.0,H2Before+produced))
                throw std::runtime_error("Adaptive hydrogen conservation audit failed");
            const double expectedEnd=accepted.vessel.timeSeconds+dt;
            // Exact requested horizon when the last step ends there.
            auto candidate=secondHalf.state;
            candidate.vessel.timeSeconds=
                (dt>=remaining)?endTimeSeconds:expectedEnd;
            AdaptiveVesselSample sample;
            sample.state=candidate;
            sample.acceptedStepSeconds=dt;
            sample.pressurePa=pressurePa(candidate,model);
            sample.estimatedNormalizedError=error;
            sample.intervalEmittedEnergyJ=
                firstHalf.vesselStep.decay.emittedEnergyJ+
                secondHalf.vesselStep.decay.emittedEnergyJ;
            sample.intervalDepositedEnergyJ=
                deposited(firstHalf)+deposited(secondHalf);
            sample.intervalGrossH2Mol=produced;
            sample.intervalConsumedH2Mol=consumed;
            sample.intervalH2ClosureResidualMol=closure;
            sample.intervalNuclideAtomClosureResidual=nucleusClosure;
            // Reject invalid chemical projections before committing this
            // accepted candidate; errors abort without a partial result.
            if(acceptedObserver)acceptedObserver(sample);
            accepted=std::move(candidate);
            trajectory.history.push_back(sample);
            trajectory.attempts.push_back({
                trajectory.history[trajectory.history.size()-2].state.vessel.timeSeconds,
                dt,error,AdaptiveStepDisposition::Accepted,""});
            ++trajectory.acceptedSteps;
            // BE local error estimate is O(h^2), so step change scales as
            // sqrt(target / observed normalized error).
            const double factor=error<=1e-14?2.0:
                std::clamp(0.9/std::sqrt(error),0.5,2.0);
            nextDt=std::clamp(dt*factor,options.minimumStepSeconds,
                              options.maximumStepSeconds);
        } else {
            trajectory.attempts.push_back({
                accepted.vessel.timeSeconds,dt,error,
                AdaptiveStepDisposition::Rejected,
                rejectionReason.empty()?"Local error exceeds tolerance":rejectionReason});
            ++trajectory.rejectedSteps;
            if(trajectory.rejectedSteps>options.maximumRejectedSteps)
                throw std::runtime_error("Adaptive integration exceeded rejection budget");
            if(dt<=minimumAllowed*(1.0+1e-12))
                throw std::runtime_error("Adaptive integration could not meet minimum step");
            const double shrink=(computed && error>1.0)
                ?std::clamp(0.8/std::sqrt(error),0.1,0.5):0.5;
            nextDt=std::max(minimumAllowed,dt*shrink);
        }
    }
    trajectory.finalState=accepted;
    if(trajectory.acceptedSteps+1!=trajectory.history.size())
        throw std::runtime_error("Adaptive history is inconsistent");
    return trajectory;
}
} // namespace thermo
