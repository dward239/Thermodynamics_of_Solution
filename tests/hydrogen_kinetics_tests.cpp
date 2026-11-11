#include "hydrogen_kinetics.h"
#include "radiolytic_hydrogen_vessel.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace thermo;
static void ensure(bool ok,const char* why){
    if(!ok)throw std::runtime_error(why);
}
static void near(double x,double target,double tol,const char* why){
    if(!std::isfinite(x)||
       std::abs(x-target)>tol*std::max(1.0,std::abs(target)))
        throw std::runtime_error(why);
}
template<class F> static void invalid(F f,const char* why){
    bool rejected=false;
    try{f();}catch(const std::invalid_argument&){rejected=true;}
    ensure(rejected,why);
}
static void testNoTransferAndGrossProduction(){
    HydrogenKineticsParameters p;
    auto result=advanceHydrogenKinetics({0.5,0.2,0},10.0,0.3,p);
    near(result.state.dissolvedMol,0.8,1e-11,"Gross H2 enters liquid");
    near(result.state.headspaceMol,0.2,1e-11,"No headspace transfer at zero k");
    near(result.consumedDuringStepMol,0,1e-11,"No consumption at zero k");
    near(result.closureResidualMol,0,1e-11,"No-transfer conservation");
}
static void testIrreversibleTransfer(){
    HydrogenKineticsParameters p;
    p.transferRatePerSecond=1.0;
    p.equilibriumDissolvedToGasMoleRatio=0.0;
    p.substeps=1000;
    auto result=advanceHydrogenKinetics({1.0,0.0,0},10.0,0.0,p);
    near(result.state.dissolvedMol,std::exp(-10.0),1e-4,
         "First-order gas transfer convergence");
    near(result.state.dissolvedMol+result.state.headspaceMol,1,1e-11,
         "Transfer conserves hydrogen");
    near(result.netTransferToHeadspaceMol,result.state.headspaceMol,1e-11,
         "Positive headspace transfer");
}
static void testReverseTransferAndPartitionRatio(){
    HydrogenKineticsParameters p;
    p.transferRatePerSecond=5.0;
    p.equilibriumDissolvedToGasMoleRatio=2.0;
    p.substeps=1000;
    auto result=advanceHydrogenKinetics({0,3.0,0},5,0,p);
    near(result.state.dissolvedMol,2.0,1e-8,
         "Liquid equilibrium inventory");
    near(result.state.headspaceMol,1.0,1e-8,
         "Gas equilibrium inventory");
    ensure(result.netTransferToHeadspaceMol<0,
           "Gas should dissolve into liquid if oversaturated");
}
static void testRapidTransferStability(){
    HydrogenKineticsParameters p;
    p.transferRatePerSecond=1e10;
    p.equilibriumDissolvedToGasMoleRatio=1.0;
    p.substeps=32;
    auto result=advanceHydrogenKinetics({0.8,0.2,0},1.0,0,p);
    near(result.state.dissolvedMol,0.5,1e-7,
         "Stiff liquid gas balance");
    near(result.state.headspaceMol,0.5,1e-7,
         "Stiff liquid gas redistribution");
    near(result.state.dissolvedMol+result.state.headspaceMol,1,1e-10,
         "Stiff transfer conservation");
}
static void testConsumptionConservation(){
    HydrogenKineticsParameters p;
    p.dissolvedConsumptionRatePerSecond=0.25;
    p.substeps=2000;
    auto result=advanceHydrogenKinetics({1,0,0},10,0,p);
    near(result.state.dissolvedMol,std::exp(-2.5),0.001,
         "First-order consumption");
    near(result.state.dissolvedMol+result.consumedDuringStepMol,1,1e-10,
         "Consumption tracked separately");
    near(result.state.cumulativeConsumedMol,result.consumedDuringStepMol,1e-11,
         "Cumulative pseudo-reactive sink");
    near(result.closureResidualMol,0,1e-10,"Hydrogen accounting closure");
}
static void testMixedProductionAndLoss(){
    HydrogenKineticsParameters p;
    p.transferRatePerSecond=1.0;
    p.dissolvedConsumptionRatePerSecond=0.4;
    p.equilibriumDissolvedToGasMoleRatio=0.5;
    p.substeps=100;
    auto result=advanceHydrogenKinetics({0.1,0.2,0.3},10,0.4,p);
    ensure(result.state.dissolvedMol>=0 && result.state.headspaceMol>=0,
           "Positive inventories");
    near(0.1+0.2+0.4,
         result.state.dissolvedMol+result.state.headspaceMol+
         result.consumedDuringStepMol,1e-10,"H2 generated-loss balance");
    near(result.state.cumulativeConsumedMol,0.3+result.consumedDuringStepMol,
         1e-10,"Cumulative consumption balance");
}
static void testVesselCoupling(){
    const DecayNetwork network({
        {"SyntheticX","X",10,2,{{"SyntheticY",1}},{{"water",0.5}}},
        {"SyntheticY","Y",0,0,{},{}}
    });
    RadiolyticHydrogenVesselModel model;
    model.vessel.headspaceVolumeM3=0.01;
    model.vessel.sensibleHeatFraction=0.5;
    model.vessel.thermal={2,500,0,298.15};
    model.vessel.radiolysisYields={{"water","H2",1.0}};
    model.liquidPhaseId="water";
    model.kinetics.transferRatePerSecond=1.0;
    model.kinetics.substeps=100;
    RadiolyticHydrogenVesselState state;
    state.vessel.nuclideAtoms={{"SyntheticX",8.0}};
    state.vessel.headspaceGasMol=1.0;
    state.vessel.temperatureK=298.15;
    const auto step=advanceRadiolyticHydrogenVessel(state,10,network,model);
    near(step.state.vessel.timeSeconds,10,1e-11,"Vessel time");
    near(step.state.vessel.temperatureK,298.152,1e-11,"Decay sensible heat");
    near(step.vesselStep.decay.emittedEnergyJ,8,1e-11,"Decay energy");
    near(step.vesselStep.decay.depositedEnergyJByPhase.at("water"),4,1e-11,
         "Phase absorbed energy");
    near(step.state.vessel.nuclideAtoms.at("SyntheticX"),4,1e-11,
         "Parent decay");
    ensure(step.hydrogenStep.grossGeneratedMol>0 &&
           step.state.hydrogen.headspaceMol>0,
           "Radiolytic H2 should reach gas headspace");
    near(step.state.vessel.headspaceGasMol,
         1.0+step.state.hydrogen.headspaceMol,1e-10,
         "Inert gas plus newly formed headspace H2");
    near(step.vesselStep.headspacePressurePa,
         step.state.vessel.headspaceGasMol*8.31446261815324*
         step.state.vessel.temperatureK/0.01,1e-10,
         "Updated pressure from gas and temperature");
    near(step.hydrogenStep.grossGeneratedMol,
         step.state.hydrogen.dissolvedMol+
         step.state.hydrogen.headspaceMol+
         step.hydrogenStep.consumedDuringStepMol,1e-10,
         "Gross production, net gas and sink bookkeeping");
    auto next=advanceRadiolyticHydrogenVessel(step.state,10,network,model);
    ensure(next.state.vessel.headspaceGasMol>step.state.vessel.headspaceGasMol,
           "Gas accumulates across time steps");
    near(next.state.vessel.nuclideAtoms.at("SyntheticX"),2,1e-10,
         "Two half-lives of decay");
}
static void testGuards(){
    HydrogenKineticsParameters p;
    invalid([&]{(void)advanceHydrogenKinetics({0,0,0},0,1,p);},
            "Cannot add H2 at zero duration");
    invalid([&]{(void)advanceHydrogenKinetics({0,0,0},1,-1,p);},
            "Gross H2 must not be negative");
    p.substeps=0;
    invalid([&]{(void)advanceHydrogenKinetics({0,0,0},1,0,p);},
            "Zero substeps must reject");
    p.substeps=32;
    p.transferRatePerSecond=-1;
    invalid([&]{(void)advanceHydrogenKinetics({0,0,0},1,0,p);},
            "Negative transfer rate must reject");
    const DecayNetwork network({
        {"X","X",10,2,{{"Y",1}},{{"water",0.5}}},
        {"Y","Y",0,0,{},{}}
    });
    RadiolyticHydrogenVesselModel model;
    model.vessel.headspaceVolumeM3=0.01;
    model.vessel.sensibleHeatFraction=1;
    model.vessel.thermal={1,1000,0,298.15};
    model.vessel.radiolysisYields={{"water","H2",1.0}};
    RadiolyticHydrogenVesselState initial;
    initial.vessel.nuclideAtoms={{"X",8.0}};
    initial.vessel.headspaceGasMol=0.1;
    initial.hydrogen.headspaceMol=0.2;
    invalid([&]{(void)advanceRadiolyticHydrogenVessel(
                initial,10,network,model);},
            "Subset H2 cannot exceed all headspace gas");
    initial.hydrogen.headspaceMol=0;
    model.vessel.radiolysisYields.push_back({"water","O2",0.2});
    invalid([&]{(void)advanceRadiolyticHydrogenVessel(
                initial,10,network,model);},
            "Other radiolytic products require additional kinetic models");
}
int main(){
    testNoTransferAndGrossProduction();
    testIrreversibleTransfer();
    testReverseTransferAndPartitionRatio();
    testRapidTransferStability();
    testConsumptionConservation();
    testMixedProductionAndLoss();
    testVesselCoupling();
    testGuards();
    std::cout<<"Radiolytic H2 kinetic vessel accounting passed\n";
}
