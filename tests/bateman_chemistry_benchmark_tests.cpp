#include "time_dependent_chemistry.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace thermo;
namespace {
constexpr double NA=6.02214076e23;
void check(bool ok,const char* label){if(!ok)throw std::runtime_error(label);}
void near(double actual,double expected,double tol,const char* label){
    if(!std::isfinite(actual)||std::abs(actual-expected)>tol*std::max(1.0,std::abs(expected)))
        throw std::runtime_error(label);
}
DecayNetwork network(){return DecayNetwork({
    {"A-parent","A",10.0,1e-23,{{"B-middle",1.0}},{{"water",0.5}}},
    {"B-middle","B",20.0,1e-23,{{"C-stable",1.0}},{{"water",0.5}}},
    {"C-stable","C",0.0,0.0,{},{}}
});}
RadiolyticHydrogenVesselModel vessel(){
    RadiolyticHydrogenVesselModel m;
    m.vessel.headspaceVolumeM3=.01;
    m.vessel.sensibleHeatFraction=.5;
    m.vessel.thermal={2,500,0,298.15};
    m.vessel.radiolysisYields={{"water","H2",0.5}};
    m.kinetics.transferRatePerSecond=.1;
    m.kinetics.equilibriumDissolvedToGasMoleRatio=.2;
    m.kinetics.substeps=16;
    return m;
}
TimeDependentChemistryModel chemistry(){
    TimeDependentChemistryModel m;
    m.equilibriumTemplate.components={{"A",0,true},{"B",0,true},{"C",0,true},{"I",0,true}};
    m.equilibriumTemplate.phases={
        {"water",1,1,IdealFluidPhaseKind::AqueousLiquid},
        {"organic",2,1,IdealFluidPhaseKind::OrganicLiquid},
        {"gas",3,1,IdealFluidPhaseKind::Gas}
    };
    m.equilibriumTemplate.species={
        {"A(aq)","water",{{"A",1}},0},
        {"A(org)","organic",{{"A",1}},0},
        {"A(g)","gas",{{"A",1}},0},
        {"B(aq)","water",{{"B",1}},0},
        {"C(aq)","water",{{"C",1}},0},
        {"I(aq)","water",{{"I",1}},0}
    };
    m.backgroundComponentMol={{"I",.125}};
    return m;
}
void analyticalChainAndPartition(){
    RadiolyticHydrogenVesselState s;
    s.vessel.nuclideAtoms={{"A-parent",NA}};
    s.vessel.headspaceGasMol=.1;
    AdaptiveVesselOptions o;
    o.initialStepSeconds=1;
    o.minimumStepSeconds=1e-6;
    o.maximumStepSeconds=4;
    o.relativeTolerance=1e-5;
    o.maximumAcceptedSteps=20000;
    const auto run=integrateVesselWithChemicalProjection(s,20,network(),vessel(),chemistry(),o);
    check(run.samples.size()>2,"Expected accepted time samples");
    const double la=std::log(2.)/10.;
    const double lb=std::log(2.)/20.;
    bool sawDaughter=false;
    for(const auto& sample:run.samples){
        const double t=sample.vessel.state.vessel.timeSeconds;
        const double a=std::exp(-la*t);
        const double b=la/(lb-la)*(std::exp(-la*t)-std::exp(-lb*t));
        const double c=1-a-b;
        near(sample.inputComponentMol.at("A"),a,2e-9,"Bateman A");
        near(sample.inputComponentMol.at("B"),b,2e-9,"Bateman B");
        near(sample.inputComponentMol.at("C"),c,2e-9,"Bateman C");
        near(sample.fluid.speciesMoles.at("A(aq)"),a/6,2e-9,"A aqueous fraction");
        near(sample.fluid.speciesMoles.at("A(org)"),a/3,2e-9,"A organic fraction");
        near(sample.fluid.speciesMoles.at("A(g)"),a/2,2e-9,"A gas fraction");
        near(sample.fluid.speciesMoles.at("B(aq)"),b,2e-9,"B aqueous");
        near(sample.fluid.speciesMoles.at("C(aq)"),c,2e-9,"C aqueous");
        near(sample.fluid.speciesMoles.at("I(aq)"),.125,2e-9,"Background conservation");
        near(a+b+c,1.,2e-9,"Analytical nucleus closure");
        check(sample.maximumScaledComponentClosure<=1.,"Component closure tolerance");
        if(t>0 && c>1e-8) sawDaughter=true;
    }
    near(run.samples.front().fluid.speciesMoles.at("C(aq)"),0.,1e-12,"No artificial granddaughter seed");
    check(sawDaughter,"Granddaughter must appear during decay");
    near(run.samples.back().vessel.state.vessel.timeSeconds,20.,1e-12,"Final time");
}
}
int main(){analyticalChainAndPartition();std::cout<<"Analytical Bateman and three-phase benchmark passed\n";}
