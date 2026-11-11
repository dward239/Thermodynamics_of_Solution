#include "aqueous_sweep.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>

namespace aqueous {
std::vector<SweepPoint> runSweep(
    const Model& base, const std::string& parameter,
    const std::string& componentName, double start, double stop,
    int steps, const Options& options) {
    if(parameter!="temperature" && parameter!="total")
        throw std::invalid_argument("Sweep parameter must be 'temperature' or 'total'");
    if(!std::isfinite(start)||!std::isfinite(stop)||steps<2||steps>1001 || start>=stop)
        throw std::invalid_argument("Sweep requires finite increasing endpoints and 2..1001 steps");
    if(parameter=="temperature" && (start<273.15 || stop>323.15))
        throw std::invalid_argument("Temperature sweep restricted to 273.15..323.15 K");
    if(parameter=="total") {
        if(start<=0)
            throw std::invalid_argument("Component total sweep requires positive concentrations");
        const auto it=std::find_if(base.components.begin(),base.components.end(),
                         [&](const Component& c){ return c.name==componentName; });
        if(it==base.components.end() || !it->fixedTotal)
            throw std::invalid_argument("Sweep total requires an existing fixed-total component");
    }
    std::vector<SweepPoint> series;
    series.reserve(static_cast<std::size_t>(steps));
    for(int i=0;i<steps;++i) {
        const double alpha=static_cast<double>(i)/(steps-1);
        const double value=start+(stop-start)*alpha;
        Model m=base;
        if(parameter=="temperature") m.temperatureK=value;
        else {
            for(auto& c:m.components)
                if(c.name==componentName) c.totalMolL=value;
        }
        try {
            series.push_back({value,solve(m,options)});
        } catch(const std::exception& err) {
            throw std::runtime_error("Sweep failed at step "+std::to_string(i+1)+
                                     ": "+err.what());
        }
    }
    return series;
}
namespace {
std::string quoteCsv(const std::string& input) {
    std::string out="\"";
    for(char c:input) { if(c=='"') out+='\"'; out+=c; }
    return out+"\"";
}
}
std::string sweepCsv(const std::vector<SweepPoint>& points,
                     const std::string& parameter) {
    if(points.empty()) throw std::invalid_argument("Sweep has no results");
    if(parameter!="temperature" && parameter!="total")
        throw std::invalid_argument("Invalid sweep label");
    std::set<std::string> species, solids;
    for(const auto& p:points) {
        for(const auto& [name,conc]:p.result.concentrationMolL) species.insert(name);
        for(const auto& [name,amount]:p.result.solidPrecipitatedMolL) solids.insert(name);
    }
    std::ostringstream out;out<<std::setprecision(15);
    out<<"step,parameter,value,temperature_K,normalized_residual,charge_residual_mol_L,ionic_strength_mol_L,pH";
    for(const auto& name:species) out<<","<<quoteCsv("concentration_"+name+"_mol_L");
    for(const auto& name:solids) {
        out<<","<<quoteCsv("solid_"+name+"_mol_L");
        out<<","<<quoteCsv("SI_"+name);
    }
    out<<"\n";
    for(std::size_t i=0;i<points.size();++i) {
        const auto& p=points[i];
        out<<(i+1)<<","<<parameter<<","<<p.parameterValue<<","
           <<p.result.temperatureK<<","<<p.result.normalizedResidual<<","
           <<p.result.chargeResidualMolL<<","<<p.result.ionicStrengthMolL<<",";
        if(p.result.hasPH) out<<p.result.pH;
        for(const auto& name:species) {
            const auto found=p.result.concentrationMolL.find(name);
            out<<",";
            if(found!=p.result.concentrationMolL.end()) out<<found->second;
        }
        for(const auto& name:solids) {
            auto s=p.result.solidPrecipitatedMolL.find(name);
            auto si=p.result.saturationIndex.find(name);
            out<<",";
            if(s!=p.result.solidPrecipitatedMolL.end()) out<<s->second;
            out<<",";
            if(si!=p.result.saturationIndex.end()) out<<si->second;
        }
        out<<"\n";
    }
    return out.str();
}
}
