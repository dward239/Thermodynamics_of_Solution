#include "aqueous_solid.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace aqueous {
namespace {
double logSaturation(const Result& solution, const Solid& solid, double logKsp) {
    double logIap=0.0;
    for(const auto& [name,nu]:solid.dissolutionStoichiometry) {
        auto it=solution.activity.find(name);
        if(it==solution.activity.end() || !std::isfinite(it->second) || it->second<=0)
            throw std::runtime_error("Solid basis activity missing or invalid: "+name);
        logIap+=nu*std::log10(it->second);
    }
    if(!std::isfinite(logIap)) throw std::runtime_error("Nonfinite ion activity product");
    return logIap-logKsp;
}
}
Result solveSingleIdealSolid(const Model& model,const Options& options) {
    if(model.activity!=ActivityModel::Ideal)
        throw std::invalid_argument("Solid equilibrium currently supports ideal activities only");
    if(model.solids.size()!=1)
        throw std::invalid_argument("Exactly one solid phase is supported; competing solids unsupported");
    for(const auto& comp:model.components)
        if(!comp.fixedTotal)
            throw std::invalid_argument("Solid equilibrium does not yet support a charge-balance free component");

    const auto& solid=model.solids.front();
    if(solid.name.empty() || solid.dissolutionStoichiometry.empty() ||
       !std::isfinite(solid.logKspReference) || !std::isfinite(solid.deltaHkJMol))
        throw std::invalid_argument("Invalid solid definition");
    Complex dissolution;
    dissolution.logKReference=solid.logKspReference;
    dissolution.deltaHkJMol=solid.deltaHkJMol;
    const double logKsp=logKAtTemperature(dissolution,model.temperatureK,model.referenceTemperatureK);

    long long charge=0;
    double xmax=std::numeric_limits<double>::infinity();
    for(const auto& [name,nu]:solid.dissolutionStoichiometry) {
        if(nu<=0) throw std::invalid_argument("Solid dissolution stoichiometry must be positive");
        const auto it=std::find_if(model.components.begin(),model.components.end(),
                                  [&](const Component& c){return c.name==name;});
        if(it==model.components.end() || !it->fixedTotal)
            throw std::invalid_argument("Solid refers to unknown or free basis: "+name);
        if(!std::isfinite(it->totalMolL) || it->totalMolL<=0)
            throw std::invalid_argument("Solid basis total must be finite and positive");
        charge+=static_cast<long long>(nu)*it->charge;
        xmax=std::min(xmax,it->totalMolL/nu);
    }
    if(charge!=0) throw std::invalid_argument("Only electrically neutral precipitating solids are supported");
    for(const auto& c:model.components)
        if(c.name==solid.name)
            throw std::invalid_argument("Solid name collides with aqueous component");
    for(const auto& c:model.complexes)
        if(c.name==solid.name)
            throw std::invalid_argument("Solid name collides with aqueous complex");

    Model aqueous=model;
    aqueous.solids.clear(); // avoid recursion; never silently ignore a declared solid
    auto state=solve(aqueous,options);
    const double initialSI=logSaturation(state,solid,logKsp);
    // Require charge-neutral input inventories even when no precipitation occurs.
    double chargeScale=1e-9;
    for(const auto& c:model.components) chargeScale+=std::abs(c.charge)*c.totalMolL;
    if(std::abs(state.chargeResidualMolL)>1e-8*chargeScale)
        throw std::invalid_argument("Solid equilibrium requires charge-neutral initial solution");
    if(initialSI<=0) {
        state.solidPrecipitatedMolL[solid.name]=0.0;
        state.saturationIndex[solid.name]=initialSI;
        return state;
    }

    const auto evaluate=[&](double amount) -> Result {
        Model modified=aqueous;
        for(auto& c:modified.components) {
            const auto it=solid.dissolutionStoichiometry.find(c.name);
            if(it!=solid.dissolutionStoichiometry.end())
                c.totalMolL-=it->second*amount;
        }
        return solve(modified,options);
    };

    double lower=0.0;
    double upper=xmax*(1.0-1e-10);
    auto upperState=evaluate(upper);
    if(logSaturation(upperState,solid,logKsp)>0)
        throw std::runtime_error("Cannot bracket solid saturation before basis depletion");
    Result midpointState=state;
    double amount=0;
    double si=initialSI;
    for(int iteration=0;iteration<100;++iteration) {
        amount=0.5*(lower+upper);
        midpointState=evaluate(amount);
        si=logSaturation(midpointState,solid,logKsp);
        if(std::abs(si)<1e-9) break;
        if(si>0) lower=amount;
        else upper=amount;
    }
    if(std::abs(si)>1e-7)
        throw std::runtime_error("Solid saturation solve did not converge");
    midpointState.solidPrecipitatedMolL[solid.name]=amount;
    midpointState.saturationIndex[solid.name]=si;
    return midpointState;
}
