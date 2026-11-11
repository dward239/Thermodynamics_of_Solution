#include "davies_multiphase_equilibrium.h"
#include "aqueous_engine.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace thermo;
static void ensure(bool okay,const char* message) {
    if(!okay)throw std::runtime_error(message);
}
static void near(double got,double expected,double tolerance,const char* message) {
    if(!std::isfinite(got)||
       std::abs(got-expected)>tolerance*std::max(1.0,std::abs(expected)))
        throw std::runtime_error(message);
}
template<class F> static void invalid(F f,const char* message) {
    bool rejected=false;
    try{f();}catch(const std::invalid_argument&){rejected=true;}
    ensure(rejected,message);
}
static IdealMultiphaseModel nacl(bool pair=false) {
    IdealMultiphaseModel model;
    model.components={{"Na",0.01,true},{"Cl",0.0,false}};
    model.phases={{"water",0.001,1000.0,IdealFluidPhaseKind::AqueousLiquid}};
    model.species={
        {"Na+","water",{{"Na",1}},0.0,+1},
        {"Cl-","water",{{"Cl",1}},0.0,-1}
    };
    if(pair)model.species.push_back(
        {"NaCl(aq)","water",{{"Na",1},{"Cl",1}},-std::log(100.0),0});
    model.electroneutralPhaseId="water";
    return model;
}
static void testSimpleSalt() {
    auto result=solveDaviesMultiphase(nacl());
    near(result.equilibrium.speciesMoles.at("Na+"),0.01,1e-10,
         "Na amount at equilibrium");
    near(result.equilibrium.speciesMoles.at("Cl-"),0.01,1e-10,
         "Electroneutral free Cl");
    near(result.ionicStrengthMolL,0.01,1e-10,"NaCl ionic strength");
    const double gamma=result.activityCoefficients.at("Na+");
    near(gamma,result.activityCoefficients.at("Cl-"),1e-12,"Symmetric ionic gamma");
    // Observed in the pinned NaCl PyPitzer comparison: the *molarity-based*
    // Davies estimate at equal numeric 0.01 is roughly 0.9018.
    ensure(gamma>0.899 && gamma<0.904,"Unexpected dilute Davies coefficient");
    near(result.activities.at("Na+"),gamma*0.01,1e-10,"Na activity");
    near(result.equilibrium.chargeResidualMol,0,1e-10,"NaCl electroneutrality");
    ensure(result.relativeIonicStrengthResidual<=1e-9,
           "Outer self-consistency residual");
}
static void testIonicAssociationAndLegacyComparison() {
    const auto model=nacl(true);
    const auto state=solveDaviesMultiphase(model);
    const auto ideal=solveIdealMultiphase(model);
    const double free=state.equilibrium.speciesMoles.at("Na+");
    ensure(free>ideal.speciesMoles.at("Na+"),
           "Davies ion-pair association should be weaker at low I");
    near(state.equilibrium.speciesMoles.at("Cl-"),free,1e-10,
         "Davies free ion equality");
    near(state.activities.at("NaCl(aq)")/
         (state.activities.at("Na+")*state.activities.at("Cl-")),
         100.0,1e-8,"Nonideal mass action");
    near(state.equilibrium.componentResidualMol.at("Na"),0,1e-10,
         "Nonideal component conservation");
    near(state.equilibrium.chargeResidualMol,0,1e-10,
         "Nonideal charge conservation");
    aqueous::Model old;
    old.temperatureK=298.15;
    old.referenceTemperatureK=298.15;
    old.activity=aqueous::ActivityModel::Davies;
    old.components={{"Na+",1,true,0.01},{"Cl-",-1,false,0.0}};
    old.complexes={{"NaCl(aq)",0,2.0,0.0,{{"Na+",1},{"Cl-",1}}}};
    const auto legacy=aqueous::solve(old);
    // This is cross-implementation numerical consistency, not experimental
    // validation; both use the same approximate Davies coefficient formula.
    near(free,legacy.concentrationMolL.at("Na+"),2e-7,
         "Legacy aqueous solver agreement");
    near(state.ionicStrengthMolL,legacy.ionicStrengthMolL,2e-7,
         "Legacy ionic strength agreement");
}
static void testNeutralOrganicAndGasPhases() {
    auto model=nacl();
    model.components.push_back({"X",0.004,true});
    model.phases.push_back({"organic",0.001,1000.0,
                            IdealFluidPhaseKind::OrganicLiquid});
    model.phases.push_back({"headspace",0.001,1000.0,
                            IdealFluidPhaseKind::Gas});
    model.species.push_back({"X(aq)","water",{{"X",1}},0.0,0});
    model.species.push_back({"X(org)","organic",{{"X",1}},
                             -std::log(2.0),0});
    model.species.push_back({"X(g)","headspace",{{"X",1}},
                             -std::log(1.0),0});
    const auto state=solveDaviesMultiphase(model);
    near(state.equilibrium.speciesMoles.at("X(aq)"),0.001,1e-10,
         "Neutral aqueous X");
    near(state.equilibrium.speciesMoles.at("X(org)"),0.002,1e-10,
         "Neutral organic X");
    near(state.equilibrium.speciesMoles.at("X(g)"),0.001,1e-10,
         "Neutral gas X");
    near(state.activityCoefficients.at("X(org)"),1.0,1e-12,
         "Neutral organic gamma unity");
}
static void testDomainAndOptionGuards() {
    auto tooStrong=nacl();
    tooStrong.components[0].totalMol=0.2;
    invalid([&]{(void)solveDaviesMultiphase(tooStrong);},
            "Davies model must reject ionic strength beyond 0.1 mol/L");
    auto options=DaviesMultiphaseOptions{};
    options.temperatureK=350;
    invalid([&]{(void)solveDaviesMultiphase(nacl(),options);},
            "Davies model must reject unsupported temperature");
    options=DaviesMultiphaseOptions{};
    options.relaxation=0.0;
    invalid([&]{(void)solveDaviesMultiphase(nacl(),options);},
            "Zero relaxation must reject");
    auto neutral=nacl();
    neutral.electroneutralPhaseId.clear();
    neutral.components[1].fixedTotal=true;
    neutral.components[1].totalMol=0.01;
    neutral.species[0].charge=0;
    neutral.species[1].charge=0;
    invalid([&]{(void)solveDaviesMultiphase(neutral);},
            "Nonideal electrolyte requires a charged phase");
}
int main() {
    testSimpleSalt();
    testIonicAssociationAndLegacyComparison();
    testNeutralOrganicAndGasPhases();
    testDomainAndOptionGuards();
    std::cout<<"Restricted Davies multi-fluid activity coupling checks passed\n";
}
