#include "time_dependent_chemistry.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace thermo;
namespace {
constexpr double NA=6.02214076e23;
void near(double a,double b,double tol,const char* what){
    if(!std::isfinite(a)||std::abs(a-b)>tol*std::max(1.0,std::abs(b)))
        throw std::runtime_error(what);
}
void require(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
template<class F> void invalid(F f,const char* reason){
    bool yes=false;try{f();}catch(const std::invalid_argument&){yes=true;}
    require(yes,reason);
}
DecayNetwork network(){return DecayNetwork({
    {"X-parent","X",10,1e-23,{{"Y-daughter",1}},{{"water",0.5}}},
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
    require(last.maximumScaledComponentClosure<=1,
            "Conservation across precipitated solids");
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
}
}
int main(){
    phaseAndDaughter();dynamicSolidAppearance();allZero();invalidModels();
    std::cout<<"Time-dependent isotope-to-neutral-equilibrium projection passed\n";
}
