#include "adaptive_vessel_integration.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace thermo;
static void near(double value,double expected,double tol,const char* label){
    if(!std::isfinite(value) ||
       std::abs(value-expected)>tol*std::max(1.0,std::abs(expected)))
        throw std::runtime_error(label);
}
static void ensure(bool condition,const char* label){
    if(!condition)throw std::runtime_error(label);
}
template<class F> static void expectsInvalid(F operation,const char* label){
    bool rejected=false;
    try{operation();}catch(const std::invalid_argument&){rejected=true;}
    ensure(rejected,label);
}
template<class F> static void expectsRuntime(F operation,const char* label){
    bool rejected=false;
    try{operation();}catch(const std::runtime_error&){rejected=true;}
    ensure(rejected,label);
}
static DecayNetwork stableNetwork(){
    // Zero radiological generation but a known deposition phase must be
    // represented so gross G(H2) can evaluate to exactly zero.
    return DecayNetwork({{"StableX","X",0,0,{},{{"water",1.0}}}});
}
static RadiolyticHydrogenVesselModel vessel(){
    RadiolyticHydrogenVesselModel m;
    m.liquidPhaseId="water";
    m.vessel.headspaceVolumeM3=0.01;
    m.vessel.sensibleHeatFraction=1.0;
    m.vessel.thermal={1.0,1000.0,0.0,298.15};
    m.vessel.radiolysisYields={{"water","H2",1.0}};
    m.kinetics.transferRatePerSecond=1.0;
    m.kinetics.dissolvedConsumptionRatePerSecond=0.0;
    m.kinetics.equilibriumDissolvedToGasMoleRatio=0.0;
    m.kinetics.substeps=1; // exercise controller's independent adaptive refinement
    return m;
}
static RadiolyticHydrogenVesselState state(){
    RadiolyticHydrogenVesselState s;
    s.vessel.timeSeconds=0.0;
    s.vessel.temperatureK=298.15;
    s.vessel.headspaceGasMol=0.1;
    s.vessel.nuclideAtoms={{"StableX",1000.0}};
    s.hydrogen.dissolvedMol=1.0;
    s.hydrogen.headspaceMol=0;
    return s;
}
static AdaptiveVesselOptions options(){
    AdaptiveVesselOptions o;
    o.initialStepSeconds=3.0;
    o.minimumStepSeconds=1e-7;
    o.maximumStepSeconds=3.0;
    o.relativeTolerance=5e-4;
    o.absoluteMolesTolerance=1e-11;
    o.absolutePressureTolerancePa=1e-4;
    o.maximumAcceptedSteps=10000;
    o.maximumRejectedSteps=10000;
    return o;
}
static void testAdaptiveReductionAndRollback(){
    const auto initial=state();
    const auto model=vessel();
    const auto result=integrateRadiolyticHydrogenVessel(
        initial,3.0,stableNetwork(),model,options());
    ensure(result.rejectedSteps>0,
           "Coarse initial time step should be rejected by local error");
    ensure(result.acceptedSteps>1,"Time grid should be refined");
    ensure(result.history.size()==result.acceptedSteps+1,"History count");
    ensure(result.attempts.size()==result.acceptedSteps+result.rejectedSteps,
           "Decision count");
    near(result.finalState.vessel.timeSeconds,3.0,1e-12,"Exact end horizon");
    near(initial.vessel.timeSeconds,0.0,1e-12,"No mutation of caller's state");
    near(initial.hydrogen.dissolvedMol,1.0,1e-12,"Immutable input inventory");
    const double target=std::exp(-3.0);
    near(result.finalState.hydrogen.dissolvedMol,target,0.002,
         "Adaptive solution approaches analytical transfer");
    near(result.finalState.hydrogen.headspaceMol,1.0-target,0.002,
         "Analytical headspace accumulation");
    near(result.finalState.vessel.headspaceGasMol,
         0.1+result.finalState.hydrogen.headspaceMol,1e-10,
         "Inert gas stays in headspace");
    near(result.finalState.hydrogen.cumulativeConsumedMol,0,1e-12,
         "No H2 consumption");
    for(std::size_t i=1;i<result.history.size();++i) {
        const auto& previous=result.history[i-1];
        const auto& current=result.history[i];
        ensure(current.state.vessel.timeSeconds>previous.state.vessel.timeSeconds,
               "Every accepted sample advances time");
        near(current.state.vessel.timeSeconds-previous.state.vessel.timeSeconds,
             current.acceptedStepSeconds,1e-11,"Sample time increment");
        ensure(current.estimatedNormalizedError<=1.0,
               "Accepted step violates local error criterion");
        near(current.intervalH2ClosureResidualMol,0,1e-10,"Hydrogen audit");
        near(current.intervalNuclideAtomClosureResidual,0,1e-9,"Nuclear audit");
    }
}
static void testIntegratedRadiolysisAndHeating(){
    const DecayNetwork network({
        {"Parent","X",10,2,{{"Daughter",1.0}},{{"water",0.5}}},
        {"Daughter","Y",0,0,{},{}}
    });
    auto model=vessel();
    model.vessel.sensibleHeatFraction=0.5;
    model.kinetics.transferRatePerSecond=0.2;
    RadiolyticHydrogenVesselState initial;
    initial.vessel.timeSeconds=0.0;
    initial.vessel.temperatureK=298.15;
    initial.vessel.headspaceGasMol=0.5;
    initial.vessel.nuclideAtoms={{"Parent",8.0}};
    auto o=options();
    o.initialStepSeconds=5.0;
    o.maximumStepSeconds=5.0;
    const auto result=integrateRadiolyticHydrogenVessel(
        initial,10.0,network,model,o);
    near(result.finalState.vessel.nuclideAtoms.at("Parent"),4.0,1e-9,
         "Parent half-life propagated over adaptive time grid");
    near(result.finalState.vessel.nuclideAtoms.at("Daughter"),4.0,1e-9,
         "Daughter produced over adaptive time grid");
    near(result.finalState.vessel.temperatureK,298.152,1e-9,
         "Integrated deposited sensible energy");
    double gross=0,emitted=0,absorbed=0;
    for(std::size_t i=1;i<result.history.size();++i) {
        gross+=result.history[i].intervalGrossH2Mol;
        emitted+=result.history[i].intervalEmittedEnergyJ;
        absorbed+=result.history[i].intervalDepositedEnergyJ;
    }
    near(emitted,8.0,1e-9,"Total decay emissions across accepted steps");
    near(absorbed,4.0,1e-9,"Total phase absorption across accepted steps");
    near(gross,4.0/(100*1.602176634e-19*6.02214076e23),1e-10,
         "Gross H2 production follows deposited energy");
    near(result.finalState.hydrogen.dissolvedMol+
         result.finalState.hydrogen.headspaceMol, gross,1e-10,
         "H2 balance across all accepted steps");
    near(result.finalState.vessel.headspaceGasMol,
         0.5+result.finalState.hydrogen.headspaceMol,1e-10,
         "Gas pressure inventory driven by H2 mass transfer");
    ensure(result.history.back().pressurePa>0,
           "Pressure must be positive");
}
static void testGuardsAndFailureDoesNotMutate(){
    const auto initial=state();
    const auto model=vessel();
    auto o=options();
    o.relativeTolerance=0;
    expectsInvalid([&]{(void)integrateRadiolyticHydrogenVessel(
        initial,3.0,stableNetwork(),model,o);},"Reject invalid tolerance");
    o=options();
    expectsInvalid([&]{(void)integrateRadiolyticHydrogenVessel(
        initial,0.0,stableNetwork(),model,o);},"Reject empty time horizon");
    o=options();
    o.maximumRejectedSteps=1;
    expectsRuntime([&]{(void)integrateRadiolyticHydrogenVessel(
        initial,3.0,stableNetwork(),model,o);},"Honor rejection budget");
    near(initial.vessel.timeSeconds,0,1e-12,"Failure leaves input state unmodified");
    near(initial.hydrogen.dissolvedMol,1,1e-12,"Failure preserves caller inventory");
    o=options();
    o.initialStepSeconds=0.01;
    o.maximumStepSeconds=0.01;
    o.maximumAcceptedSteps=1;
    expectsRuntime([&]{(void)integrateRadiolyticHydrogenVessel(
        initial,3.0,stableNetwork(),model,o);},"Honor accepted-step budget");
}
static void testStoppedAtNonzeroStart(){
    auto s=state();
    s.vessel.timeSeconds=4.0;
    auto o=options();
    o.initialStepSeconds=0.05;
    o.maximumStepSeconds=0.2;
    auto t=integrateRadiolyticHydrogenVessel(s,4.3,stableNetwork(),vessel(),o);
    near(t.history.front().state.vessel.timeSeconds,4,1e-12,
         "Preserve starting timestamp");
    near(t.history.back().state.vessel.timeSeconds,4.3,1e-12,
         "Respect absolute end horizon");
}
int main(){
    testAdaptiveReductionAndRollback();
    testIntegratedRadiolysisAndHeating();
    testGuardsAndFailureDoesNotMutate();
    testStoppedAtNonzeroStart();
    std::cout<<"Adaptive vessel step doubling and conservation checks passed\n";
}
