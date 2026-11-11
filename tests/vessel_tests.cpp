#include "vessel_simulation.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace thermo;
static void near(double a,double b,double tol,const char* msg){
    if(!std::isfinite(a)||std::abs(a-b)>tol*std::max(1.0,std::abs(b)))
        throw std::runtime_error(msg);
}
int main(){
    DecayNetwork network({{"X","Synthetic",10,2,{{"Y",1}},{{"water",0.5}}},
                          {"Y","Synthetic",0,0,{},{}}});
    VesselModel model;
    model.headspaceVolumeM3=0.01;
    model.sensibleHeatFraction=0.5;
    model.thermal={2,500,0,298.15};
    model.radiolysisYields={{"water","H2",1.0}};
    VesselState initial;
    initial.temperatureK=298.15;
    initial.nuclideAtoms={{"X",8}};
    initial.headspaceGasMol=1.0;
    auto result=advanceClosedVessel(initial,10,network,model);
    near(result.state.timeSeconds,10,1e-11,"Elapsed time");
    near(result.state.nuclideAtoms.at("X"),4,1e-11,"Parent inventory");
    near(result.state.nuclideAtoms.at("Y"),4,1e-11,"Daughter inventory");
    near(result.decay.emittedEnergyJ,8,1e-11,"Emitted energy");
    near(result.decay.depositedEnergyJByPhase.at("water"),4,1e-11,"Deposited energy");
    near(result.state.temperatureK,298.152,1e-11,"Lumped sensible heating");
    near(result.state.headspaceGasMol,1.0,1e-11,"Gross gas is not net gas");
    near(result.headspacePressurePa,8.31446261815324*298.152/0.01,1e-11,
         "Ideal gas pressure");
    near(result.grossRadiolysisMolByPhaseSpecies.at("water::H2"),
         4/(100.0*1.602176634e-19*6.02214076e23),1e-11,"Gross G conversion");
    auto added=advanceClosedVessel(initial,10,network,model,0.5);
    near(added.state.headspaceGasMol,1.5,1e-11,"Explicit net gas addition");
    bool rejected=false;
    try{advanceClosedVessel(initial,10,network,model,-2);}
    catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("Must reject gas overdraft");
    auto zero=advanceClosedVessel(initial,0,network,model);
    near(zero.state.temperatureK,initial.temperatureK,1e-11,"Zero duration temperature");
    std::cout<<"Closed-vessel coupling checks passed\n";
}
