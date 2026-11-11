#include "inventory_chemistry_bridge.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
using namespace thermo;
static void near(double a,double b,double tol,const char* msg) {
    if(!std::isfinite(a) || std::abs(a-b)>tol*std::max(1.0,std::abs(b)))
        throw std::runtime_error(msg);
}
int main() {
    constexpr double NA=6.02214076e23;
    DecayNetwork network({{"A-1","ElementA",10,0,{{"B-1",1}},{}},
                          {"A-2","ElementA",0,0,{},{}},
                          {"B-1","ElementB",0,0,{},{}}});
    const auto elemental=elementalMolesFromNuclides(network,
        {{"A-1",2*NA},{"A-2",3*NA},{"B-1",NA}});
    near(elemental.at("ElementA"),5,1e-12,"Isotopes sum to elemental moles");
    near(elemental.at("ElementB"),1,1e-12,"Daughter element moles");
    auto evolved=network.advance({{"A-1",2*NA},{"A-2",3*NA}},10.0);
    auto later=elementalMolesFromNuclides(network,evolved.finalAtoms);
    near(later.at("ElementA"),4,1e-11,"Elemental conversion after decay");
    near(later.at("ElementB"),1,1e-11,"Daughter element after decay");
    MultiphaseSystem system;
    system.phases={{"water",PhaseKind::AqueousLiquid,0.001,298.15}};
    system.species={{"A(aq)","water",0,{{"ElementA",1}}},
                    {"B(aq)","water",0,{{"ElementB",1}}}};
    validateElementCoverage(system,later);
    system.species.pop_back();
    bool rejected=false;
    try{validateElementCoverage(system,later);}
    catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("Missing daughter element species was accepted");
    rejected=false;
    try{elementalMolesFromNuclides(network,{{"Unknown",1}});}
    catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("Unknown nuclear identity was accepted");
    std::cout<<"Nuclear-to-element inventory interface checks passed\n";
}
