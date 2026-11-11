#include "ideal_multiphase_equilibrium.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
namespace thermo {
namespace {
struct Prepared {
    std::vector<std::vector<double>> nu; // species-by-component stoichiometry
    std::vector<double> lnScale;          // ln(V * concentration standard) - mu0/RT
    std::vector<double> totals;
    std::vector<double> volumes;
};
struct Calculation {
    std::vector<double> moles;
    std::vector<double> residual;
    double norm=std::numeric_limits<double>::infinity();
};
Prepared validate(const IdealMultiphaseModel& model) {
    if(model.components.empty()||model.phases.empty()||model.species.empty())
        throw std::invalid_argument("Nonempty components, phases and species are required");
    std::map<std::string,std::size_t> componentIds;
    Prepared p;
    for(std::size_t i=0;i<model.components.size();++i) {
        const auto& c=model.components[i];
        if(c.id.empty()||!componentIds.emplace(c.id,i).second||
           !std::isfinite(c.totalMol)||c.totalMol<=0)
            throw std::invalid_argument("Invalid component name/positive total");
        p.totals.push_back(c.totalMol);
    }
    std::map<std::string,IdealPhase> phases;
    for(const auto& phase:model.phases) {
        if(phase.id.empty()||!phases.emplace(phase.id,phase).second||
           !std::isfinite(phase.volumeM3)||phase.volumeM3<=0||
           !std::isfinite(phase.standardConcentrationMolM3)||
           phase.standardConcentrationMolM3<=0)
            throw std::invalid_argument("Invalid phase, volume or concentration reference");
    }
    std::set<std::string> speciesIds;
    std::vector<bool> covered(model.components.size(),false);
    for(const auto& s:model.species) {
        if(s.id.empty()||!speciesIds.insert(s.id).second||
           !phases.count(s.phaseId)||!std::isfinite(s.standardMuOverRT)||
           s.componentStoichiometry.empty())
            throw std::invalid_argument("Invalid phase/species/chemical-potential definition");
        std::vector<double> n(model.components.size(),0.0);
        for(const auto& [name,count]:s.componentStoichiometry) {
            auto it=componentIds.find(name);
            if(it==componentIds.end()||count<=0)
                throw std::invalid_argument("Unknown/invalid component stoichiometry");
            n[it->second]=count;
            covered[it->second]=true;
        }
        const auto& phase=phases.at(s.phaseId);
        p.nu.push_back(std::move(n));
        p.volumes.push_back(phase.volumeM3);
        p.lnScale.push_back(std::log(phase.volumeM3)+
                            std::log(phase.standardConcentrationMolM3)-
                            s.standardMuOverRT);
    }
    for(bool found:covered)
        if(!found)throw std::invalid_argument("A component has no representing species");
    return p;
}
bool evaluate(const Prepared& p,const std::vector<double>& x,Calculation& result) {
    const std::size_t n=p.totals.size();
    result.moles.assign(p.nu.size(),0.0);
    std::vector<double> calculated(n,0.0);
    for(std::size_t s=0;s<p.nu.size();++s) {
        double ln=p.lnScale[s];
        for(std::size_t i=0;i<n;++i)ln+=p.nu[s][i]*x[i];
        if(!std::isfinite(ln)||ln>695.0)return false;
        const double amount=ln < -740.0 ? 0.0:std::exp(ln);
        if(!std::isfinite(amount))return false;
        result.moles[s]=amount;
        for(std::size_t i=0;i<n;++i)calculated[i]+=p.nu[s][i]*amount;
    }
    result.norm=0.0;
    result.residual.assign(n,0.0);
    for(std::size_t i=0;i<n;++i) {
        result.residual[i]=(calculated[i]-p.totals[i])/p.totals[i];
        if(!std::isfinite(result.residual[i]))return false;
        result.norm=std::max(result.norm,std::abs(result.residual[i]));
    }
    return true;
}
bool solveLinear(std::vector<std::vector<double>> a,
                 std::vector<double> b,std::vector<double>& x) {
    const std::size_t n=b.size();
    for(std::size_t k=0;k<n;++k) {
        std::size_t pivot=k;
        for(std::size_t i=k+1;i<n;++i)
            if(std::abs(a[i][k])>std::abs(a[pivot][k]))pivot=i;
        if(!std::isfinite(a[pivot][k])||std::abs(a[pivot][k])<1e-14)
            return false;
        std::swap(a[pivot],a[k]);
        std::swap(b[pivot],b[k]);
        for(std::size_t i=k+1;i<n;++i) {
            const double ratio=a[i][k]/a[k][k];
            for(std::size_t j=k+1;j<n;++j)a[i][j]-=ratio*a[k][j];
            b[i]-=ratio*b[k];
        }
    }
    x.assign(n,0.0);
    for(std::size_t i=n;i-- >0;) {
        double rhs=b[i];
        for(std::size_t j=i+1;j<n;++j)rhs-=a[i][j]*x[j];
        x[i]=rhs/a[i][i];
        if(!std::isfinite(x[i]))return false;
    }
    return true;
}
}
IdealEquilibriumResult solveIdealMultiphase(const IdealMultiphaseModel& model,
                                           int maxIterations,double tolerance) {
    if(maxIterations<=0||maxIterations>100000||!std::isfinite(tolerance)||
       tolerance<=0||tolerance>=0.01)
        throw std::invalid_argument("Invalid multiphase solver options");
    const auto p=validate(model);
    const std::size_t n=p.totals.size();
    std::vector<double> x(n,0.0);
    // A pure component monomer supplies a physically meaningful initial scale.
    for(std::size_t i=0;i<n;++i) {
        bool found=false;
        for(std::size_t s=0;s<p.nu.size();++s) {
            bool monomer=p.nu[s][i]==1;
            for(std::size_t j=0;j<n;++j)
                if(j!=i&&p.nu[s][j]!=0)monomer=false;
            if(monomer) {
                x[i]=std::log(p.totals[i])-p.lnScale[s];
                found=true;
                break;
            }
        }
        if(!found)x[i]=std::log(p.totals[i]);
    }
    Calculation current;
    if(!evaluate(p,x,current))throw std::runtime_error("Invalid initial ideal multiphase state");
    int iterations=0;
    while(current.norm>tolerance&&iterations<maxIterations) {
        std::vector<std::vector<double>> jac(n,std::vector<double>(n,0.0));
        for(std::size_t s=0;s<p.nu.size();++s)
            for(std::size_t i=0;i<n;++i)
                for(std::size_t j=0;j<n;++j)
                    jac[i][j]+=(p.nu[s][i]*p.nu[s][j]*current.moles[s])/p.totals[i];
        std::vector<double> rhs=current.residual;
        for(double& v:rhs)v=-v;
        std::vector<double> delta;
        if(!solveLinear(jac,rhs,delta))
            throw std::runtime_error("Ideal multiphase stoichiometric Jacobian singular");
        bool accepted=false;
        for(unsigned k=0;k<55;++k) {
            const double alpha=std::ldexp(1.0,-static_cast<int>(k));
            std::vector<double> trial=x;
            for(std::size_t i=0;i<n;++i)trial[i]+=alpha*delta[i];
            Calculation next;
            if(evaluate(p,trial,next)&&next.norm<current.norm) {
                x.swap(trial);
                current=std::move(next);
                accepted=true;
                break;
            }
        }
        if(!accepted)throw std::runtime_error("Ideal multiphase line search stalled");
        ++iterations;
    }
    if(current.norm>tolerance)
        throw std::runtime_error("Ideal multiphase equilibrium did not converge");
    IdealEquilibriumResult result;
    result.iterations=iterations;
    result.maxNormalizedResidual=current.norm;
    for(std::size_t s=0;s<p.nu.size();++s) {
        result.speciesMoles[model.species[s].id]=current.moles[s];
        result.concentrationsMolM3[model.species[s].id]=current.moles[s]/p.volumes[s];
    }
    for(std::size_t i=0;i<n;++i) {
        double calculated=0.0;
        for(std::size_t s=0;s<p.nu.size();++s)
            calculated+=p.nu[s][i]*current.moles[s];
        result.componentResidualMol[model.components[i].id]=calculated-p.totals[i];
    }
    return result;
}
} // namespace thermo
