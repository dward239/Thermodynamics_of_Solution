#include "time_dependent_chemistry.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace thermo;
namespace {
constexpr double NA=6.02214076e23;
constexpr double heatingDecayEnergyJ=1e-20;
void near(double a,double b,double tol,const char* what){
    if(!std::isfinite(a)||std::abs(a-b)>tol*std::max(1.0,std::abs(b)))
        throw std::runtime_error(what);
}
void require(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
template<class F> void invalid(F f,const char* reason){
    bool yes=false;try{f();}catch(const std::invalid_argument&){yes=true;}
    require(yes,reason);
}
DecayNetwork network(double decayEnergyJ=1e-23){return DecayNetwork({
    {"X-parent","X",10,decayEnergyJ,{{"Y-daughter",1}},{{"water",0.5}}},
    {"Y-daughter","Y",0,0,{},{}}});}
RadiolyticHydrogenVesselModel vessel(){
    RadiolyticHydrogenVesselModel m;
    m.vessel.headspaceVolumeM3=.01;m.vessel.sensibleHeatFraction=.5;
    m.vessel.thermal={2,500,0,298.15};
    m.vessel.radiolysisYields={{"water","H2",0.5}};
    m.kinetics.transferRatePerSecond=.1;
    m.kinetics.equilibriumDissolvedToGasMoleRatio=.2;
    m.kinetics.substeps=16;
    return m;
}
RadiolyticHydrogenVesselState initial(){
    RadiolyticHydrogenVesselState s;
    s.vessel.nuclideAtoms={{"X-parent",NA}};
    s.vessel.headspaceGasMol=.1;
    return s;
}
AdaptiveVesselOptions options(){
    AdaptiveVesselOptions o;
    o.initialStepSeconds=2;o.minimumStepSeconds=1e-5;
    o.maximumStepSeconds=10;o.relativeTolerance=1e-4;
    o.maximumAcceptedSteps=20000;
    return o;
}
TimeDependentChemistryModel chem(){
    TimeDependentChemistryModel m;
    m.equilibriumTemplate.components={{"X",0,true},{"Y",0,true},{"B",0,true}};
    m.equilibriumTemplate.phases={
        {"water",1,1,IdealFluidPhaseKind::AqueousLiquid},
        {"organic",1,1,IdealFluidPhaseKind::OrganicLiquid}
    };
    m.equilibriumTemplate.species={
        {"X(aq)","water",{{"X",1}},0},
        {"X(org)","organic",{{"X",1}},-std::log(2.)},
        {"Y(aq)","water",{{"Y",1}},0},
        {"B(aq)","water",{{"B",1}},0}
    };
    m.backgroundComponentMol={{"B",.25}};
    return m;
}
ChemicalProjectionThermodynamics thermodynamicsFor(const TimeDependentChemistryModel& model) {
    ChemicalProjectionThermodynamics data;
    for(const auto& species:model.equilibriumTemplate.species)
        data.standardEnthalpyJPerMolById[species.id]=0.0;
    for(const auto& solid:model.candidateSolids)
        data.standardEnthalpyJPerMolById[solid.id]=0.0;
    return data;
}
RadiolyticHydrogenVesselModel heatedVessel() {
    auto model=vessel();
    model.vessel.thermal.heatCapacityJPerKgK=25.0;
    return model;
}
void phaseAndDaughter(){
    auto s=initial();auto n=network();auto v=vessel();auto c=chem();auto o=options();
    auto t=integrateVesselWithChemicalProjection(s,10,n,v,c,o);
    auto bare=integrateRadiolyticHydrogenVessel(s,10,n,v,o);
    require(t.samples.size()==bare.history.size(),"Chemistry sample per vessel time");
    for(std::size_t i=0;i<t.samples.size();++i) {
        near(t.samples[i].vessel.state.vessel.timeSeconds,
             bare.history[i].state.vessel.timeSeconds,1e-12,
             "Projection must not change accepted time grid");
        require(t.samples[i].maximumScaledComponentClosure<=1,
                "Component conservation audit");
    }
    auto first=t.samples.front(),last=t.samples.back();
    near(first.inputComponentMol.at("Y"),0,1e-12,"Absent daughter starts zero");
    near(first.fluid.speciesMoles.at("Y(aq)"),0,1e-12,
         "No artificial daughter seed concentration");
    near(first.fluid.speciesMoles.at("X(aq)"),1./3,1e-10,
         "Initial aqueous X");
    near(first.fluid.speciesMoles.at("X(org)"),2./3,1e-10,
         "Initial organic X");
    near(last.inputComponentMol.at("X"),.5,1e-10,
         "Parent element after half-life");
    near(last.inputComponentMol.at("Y"),.5,1e-10,
         "Daughter inventory after half-life");
    near(last.fluid.speciesMoles.at("X(aq)"),1./6,1e-10,
         "Re-equilibrated aqueous parent");
    near(last.fluid.speciesMoles.at("X(org)"),1./3,1e-10,
         "Re-equilibrated organic parent");
    near(last.fluid.speciesMoles.at("Y(aq)"),.5,1e-10,
         "Daughter species appears at later time");
    near(last.fluid.speciesMoles.at("B(aq)"),.25,1e-10,
         "Independent background component is conserved");
}
void dynamicSolidAppearance(){
    auto c=chem();
    c.equilibriumTemplate.species.erase(c.equilibriumTemplate.species.begin()+1);
    c.candidateSolids={{"XY(s)",{{"X",1},{"Y",1}},std::log(.01)}};
    auto t=integrateVesselWithChemicalProjection(
        initial(),10,network(),vessel(),c,options());
    near(t.samples.front().solidMoles.at("XY(s)"),0,1e-12,
         "No XY(s) when Y is absent");
    auto last=t.samples.back();
    near(last.fluid.speciesMoles.at("X(aq)"),.1,1e-8,
         "X saturation");
    near(last.fluid.speciesMoles.at("Y(aq)"),.1,1e-8,
         "Y saturation");
    near(last.solidMoles.at("XY(s)"),.4,1e-8,
         "Solid forms after daughter ingrowth");
    near(last.fluid.calculatedComponentMol.at("X"),.1,1e-8,
         "Reported fluid X amount excludes solid");
    near(last.fluid.calculatedComponentMol.at("Y"),.1,1e-8,
         "Reported fluid Y amount excludes solid");
    near(last.fluid.componentResidualMol.at("X"),0,1e-8,
         "Reported component closure includes solids");
    require(last.maximumScaledComponentClosure<=1,
            "Conservation across precipitated solids");
}
void temperatureDependentFluidPartition() {
    constexpr double R=8.31446261815324;
    constexpr double transferEnthalpy=18000.0;
    auto model=chem();
    model.equilibriumTemplate.phases.push_back({"gas",1,1,IdealFluidPhaseKind::Gas});
    model.equilibriumTemplate.species.push_back({"B(g)","gas",{{"B",1}},0});
    auto fixedStandardStates=integrateVesselWithChemicalProjection(
        initial(),10,network(heatingDecayEnergyJ),heatedVessel(),model,options());
    model.thermodynamics=thermodynamicsFor(model);
    model.thermodynamics->standardEnthalpyJPerMolById.at("X(org)")=transferEnthalpy;
    auto heated=integrateVesselWithChemicalProjection(
        initial(),10,network(heatingDecayEnergyJ),heatedVessel(),model,options());
    require(heated.samples.size()==fixedStandardStates.samples.size(),
            "Temperature-dependent chemistry preserves the vessel grid");
    require(heated.samples.back().vessel.state.vessel.temperatureK>310.0,
            "Decay heats the vessel sufficiently to test redistribution");
    const double reference=model.thermodynamics->referenceTemperatureK;
    for(std::size_t i=0;i<heated.samples.size();++i) {
        const auto& sample=heated.samples[i];
        const auto& legacy=fixedStandardStates.samples[i];
        const double temperature=sample.vessel.state.vessel.temperatureK;
        const auto& moles=sample.fluid.speciesMoles;
        near(sample.vessel.state.vessel.timeSeconds,legacy.vessel.state.vessel.timeSeconds,
             1e-12,"Chemistry preserves accepted time");
        near(moles.at("X(org)")/moles.at("X(aq)"),
             2.0*std::exp(transferEnthalpy/R*(1.0/reference-1.0/temperature)),
             1e-9,"Liquid transfer ratio follows van't Hoff");
        near(moles.at("B(g)")/moles.at("B(aq)"),reference/temperature,
             1e-9,"Ideal gas uses fixed reference pressure");
        require(sample.maximumScaledComponentClosure<=1,"Heated fluid component closure");
    }
    near(fixedStandardStates.samples.back().fluid.speciesMoles.at("B(g)"),.125,1e-9,
         "Legacy mode retains fixed gas standard concentration");
}
void temperatureDependentPrecipitation() {
    constexpr double R=8.31446261815324;
    constexpr double dissolutionEnthalpy=30000.0;
    auto model=chem();
    model.equilibriumTemplate.species.erase(model.equilibriumTemplate.species.begin()+1);
    model.candidateSolids={{"XY(s)",{{"X",1},{"Y",1}},std::log(.01)}};
    model.thermodynamics=thermodynamicsFor(model);
    model.thermodynamics->standardEnthalpyJPerMolById.at("XY(s)")=-dissolutionEnthalpy;
    auto trajectory=integrateVesselWithChemicalProjection(
        initial(),10,network(heatingDecayEnergyJ),heatedVessel(),model,options());
    near(trajectory.samples.front().solidMoles.at("XY(s)"),0,1e-12,
         "Temperature model does not seed daughter solids");
    const auto& last=trajectory.samples.back();
    const double temperature=last.vessel.state.vessel.temperatureK;
    const double reference=model.thermodynamics->referenceTemperatureK;
    const double saturationMol=.1*std::exp(dissolutionEnthalpy/(2.0*R)*
                                          (1.0/reference-1.0/temperature));
    near(last.fluid.speciesMoles.at("X(aq)"),saturationMol,1e-8,
         "Saturation amount follows temperature-dependent solubility product");
    near(last.solidMoles.at("XY(s)"),.5-saturationMol,1e-8,
         "Temperature-dependent precipitate conserves the daughter and parent");
    require(last.maximumScaledComponentClosure<=1,"Heated solid component closure");
}
void thermodynamicsValidation() {
    auto model=chem();
    model.thermodynamics=thermodynamicsFor(model);
    auto run=[&](double endTime) {
        return integrateVesselWithChemicalProjection(
            initial(),endTime,network(heatingDecayEnergyJ),heatedVessel(),model,options());
    };
    auto& data=*model.thermodynamics;
    auto& enthalpies=data.standardEnthalpyJPerMolById;
    enthalpies.erase("Y(aq)");
    invalid([&]{(void)run(.01);},"Missing enthalpy is rejected even for absent daughter");
    enthalpies["Y(aq)"]=std::numeric_limits<double>::infinity();
    invalid([&]{(void)run(.01);},"Nonfinite enthalpy is rejected");
    enthalpies["Y(aq)"]=0.0;
    enthalpies["unknown"]=0.0;
    invalid([&]{(void)run(.01);},"Unknown enthalpy entry is rejected");
    enthalpies.erase("unknown");
    data.maximumTemperatureK=290.0;
    invalid([&]{(void)run(.01);},"Reference temperature must lie within the bounds");
    data.maximumTemperatureK=323.15;
    data.minimumTemperatureK=data.referenceTemperatureK;
    near(run(.01).samples.front().vessel.state.vessel.temperatureK,data.minimumTemperatureK,
         1e-12,"An endpoint of the temperature range is inclusive");
    data.minimumTemperatureK=273.15;
    data.maximumTemperatureK=data.referenceTemperatureK;
    bool failed=false;
    try {
        (void)run(10);
    } catch(const std::domain_error&) {
        failed=true;
    }
    require(failed,"A heated candidate outside the data range aborts before commit");
}
void allZero(){
    auto s=initial();s.vessel.nuclideAtoms["X-parent"]=0;
    auto c=chem();c.backgroundComponentMol.clear();
    auto t=integrateVesselWithChemicalProjection(
        s,1,network(),vessel(),c,options());
    for(const auto& item:t.samples){
        near(item.fluid.speciesMoles.at("Y(aq)"),0,1e-12,"Zero daughter");
        near(item.fluid.speciesMoles.at("X(aq)"),0,1e-12,"Zero parent");
        near(item.fluid.speciesMoles.at("X(org)"),0,1e-12,"Zero organic");
        near(item.maximumScaledComponentClosure,0,1e-12,"Zero conservation");
    }
}
void abortWhenNewDaughterHasNoSpecies() {
    auto c=chem();
    c.equilibriumTemplate.species.erase(c.equilibriumTemplate.species.begin()+2);
    auto s=initial();
    bool failed=false;
    try {
        // Initial t0 projection succeeds; once daughter Y enters the system,
        // an accepted candidate must fail its chemistry projection.
        (void)integrateVesselWithChemicalProjection(
            s,10,network(),vessel(),c,options());
    } catch(const std::invalid_argument&) {
        failed=true;
    }
    require(failed,"Missing daughter phase representation must abort");
    near(s.vessel.timeSeconds,0,1e-12,"Failed run must not mutate caller time");
    near(s.vessel.nuclideAtoms.at("X-parent"),NA,1e-10,
         "Failed run must not mutate caller nuclides");
}
void invalidModels(){
    auto c=chem();
    c.backgroundComponentMol["B"]=-1;
    invalid([&]{(void)integrateVesselWithChemicalProjection(
        initial(),1,network(),vessel(),c,options());},
        "Reject negative background");
    c=chem();c.equilibriumTemplate.components[1].id="notY";
    invalid([&]{(void)integrateVesselWithChemicalProjection(
        initial(),1,network(),vessel(),c,options());},
        "Reject missing daughter element mapping");
    c=chem();c.equilibriumTemplate.species[0].charge=1;
    invalid([&]{(void)integrateVesselWithChemicalProjection(
        initial(),1,network(),vessel(),c,options());},
        "Reject unsupported charged chemistry");
    c=chem();c.equilibriumTemplate.components[0].totalMol=.1;
    invalid([&]{(void)integrateVesselWithChemicalProjection(
        initial(),1,network(),vessel(),c,options());},
        "Reject nonzero template placeholder");
    c=chem();c.equilibriumTemplate.phases[1].id="water";
    invalid([&]{(void)integrateVesselWithChemicalProjection(
        initial(),1,network(),vessel(),c,options());},
        "Reject duplicate phase names");
    c=chem();c.equilibriumTemplate.species[1].id="X(aq)";
    invalid([&]{(void)integrateVesselWithChemicalProjection(
        initial(),1,network(),vessel(),c,options());},
        "Reject duplicate species names");
    c=chem();c.candidateSolids={{"XY",{{"X",1},{"Y",1}},std::log(.01)},
                              {"XY",{{"X",1},{"Y",1}},std::log(.02)}};
    invalid([&]{(void)integrateVesselWithChemicalProjection(
        initial(),1,network(),vessel(),c,options());},
        "Reject duplicate solid names");
}
}
int main(){
    phaseAndDaughter();dynamicSolidAppearance();allZero();
    abortWhenNewDaughterHasNoSpecies();invalidModels();
    temperatureDependentFluidPartition();temperatureDependentPrecipitation();
    thermodynamicsValidation();
    std::cout<<"Time-dependent isotope-to-neutral-equilibrium projection passed\n";
}
