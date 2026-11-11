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
    std::vector<double> charge;
    std::vector<bool> fixed;
    double chargeScale=1.0;
    int freeIndex=-1;
};
struct Calculation {
    std::vector<double> moles;
    std::vector<double> residual;
    std::vector<double> componentMoles;
    double chargeResidualMol=0.0;
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
           !std::isfinite(c.totalMol)||
           (c.fixedTotal && c.totalMol<=0.0) ||
           (!c.fixedTotal && c.totalMol!=0.0))
            throw std::invalid_argument("Invalid fixed total or free component: "+c.id);
        if(!c.fixedTotal) {
            if(p.freeIndex!=-1)
                throw std::invalid_argument("Exactly one free component is supported");
            p.freeIndex=static_cast<int>(i);
        } else p.chargeScale+=c.totalMol;
        p.totals.push_back(c.totalMol);
        p.fixed.push_back(c.fixedTotal);
    }
    std::map<std::string,IdealPhase> phases;
    for(const auto& phase:model.phases) {
        if(phase.id.empty()||!phases.emplace(phase.id,phase).second||
           !std::isfinite(phase.volumeM3)||phase.volumeM3<=0||
           !std::isfinite(phase.standardConcentrationMolM3)||
           phase.standardConcentrationMolM3<=0)
            throw std::invalid_argument("Invalid phase, volume or concentration reference");
    }
    if((p.freeIndex!=-1)!=(!model.electroneutralPhaseId.empty()))
        throw std::invalid_argument("Electroneutrality requires exactly one free component");
    if(!model.electroneutralPhaseId.empty() &&
       !phases.count(model.electroneutralPhaseId))
        throw std::invalid_argument("Unknown charge-balanced phase");
    if(p.freeIndex!=-1 && model.components.size()<2)
        throw std::invalid_argument("A charge-balanced model needs a fixed component");
    std::set<std::string> speciesIds;
    std::vector<bool> covered(model.components.size(),false);
    bool hasChargedSpecies=false;
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
        if(s.charge!=0) {
            hasChargedSpecies=true;
            if(s.phaseId!=model.electroneutralPhaseId)
                throw std::invalid_argument(
                    "Charged species outside designated charge-balanced phase: "+s.id);
        }
        const auto& phase=phases.at(s.phaseId);
        p.charge.push_back(static_cast<double>(s.charge));
        p.nu.push_back(std::move(n));
        p.volumes.push_back(phase.volumeM3);
        p.lnScale.push_back(std::log(phase.volumeM3)+
                            std::log(phase.standardConcentrationMolM3)-
                            s.standardMuOverRT);
    }
    for(bool found:covered)
        if(!found)throw std::invalid_argument("A component has no representing species");
    if(p.freeIndex!=-1 && !hasChargedSpecies)
        throw std::invalid_argument("Charge balance requires charged species");
    return p;
}
bool evaluate(const Prepared& p,const std::vector<double>& x,Calculation& result) {
    const std::size_t n=p.totals.size();
    result.moles.assign(p.nu.size(),0.0);
    std::vector<double> calculated(n,0.0);
    result.chargeResidualMol=0.0;
    for(std::size_t s=0;s<p.nu.size();++s) {
        double ln=p.lnScale[s];
        for(std::size_t i=0;i<n;++i)ln+=p.nu[s][i]*x[i];
        if(!std::isfinite(ln)||ln>695.0)return false;
        const double amount=ln < -740.0 ? 0.0:std::exp(ln);
        if(!std::isfinite(amount))return false;
        result.moles[s]=amount;
        result.chargeResidualMol+=p.charge[s]*amount;
        for(std::size_t i=0;i<n;++i)calculated[i]+=p.nu[s][i]*amount;
    }
    result.norm=0.0;
    result.residual.assign(n,0.0);
    result.componentMoles=calculated;
    for(std::size_t i=0;i<n;++i) {
        result.residual[i]=p.fixed[i]
            ?(calculated[i]-p.totals[i])/p.totals[i]
            :result.chargeResidualMol/p.chargeScale;
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
                for(std::size_t j=0;j<n;++j) {
                    const double ai=p.fixed[i]?p.nu[s][i]:p.charge[s];
                    const double scale=p.fixed[i]?p.totals[i]:p.chargeScale;
                    jac[i][j]+=ai*p.nu[s][j]*current.moles[s]/scale;
                }
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
    result.chargeResidualMol=current.chargeResidualMol;
    for(std::size_t i=0;i<n;++i)
        result.componentPotentialOverRT[model.components[i].id]=x[i];
    for(std::size_t s=0;s<p.nu.size();++s) {
        result.speciesMoles[model.species[s].id]=current.moles[s];
        result.concentrationsMolM3[model.species[s].id]=current.moles[s]/p.volumes[s];
    }
    for(std::size_t i=0;i<n;++i) {
        const auto& name=model.components[i].id;
        result.calculatedComponentMol[name]=current.componentMoles[i];
        if(p.fixed[i])
            result.componentResidualMol[name]=current.componentMoles[i]-p.totals[i];
    }
    return result;
}

IdealSolidEquilibriumResult solveIdealMultiphaseWithSolids(
    const IdealMultiphaseModel& fluidModel,
    const std::vector<IdealPureSolid>& candidates,
    int maxIterations,double tolerance) {
    if(maxIterations<1 || maxIterations>100000 ||
       !std::isfinite(tolerance) || tolerance<=0 || tolerance>=0.01)
        throw std::invalid_argument("Invalid solid equilibrium solver options");
    if(!fluidModel.electroneutralPhaseId.empty())
        throw std::invalid_argument(
            "Charge balance with solid phase enumeration is not yet supported");
    if(candidates.size()>8)
        throw std::invalid_argument("Solid candidate enumeration limited to eight phases");
    const auto prepared=validate(fluidModel);
    const auto base=solveIdealMultiphase(fluidModel,maxIterations,tolerance);
    IdealSolidEquilibriumResult output;
    output.fluidSolution=base;
    if(candidates.empty())return output;
    const std::size_t n=fluidModel.components.size();
    const std::size_t ns=candidates.size();
    std::map<std::string,std::size_t> indices;
    for(std::size_t i=0;i<n;++i)indices[fluidModel.components[i].id]=i;
    std::set<std::string> names;
    for(const auto& sp:fluidModel.species)names.insert(sp.id);
    std::vector<std::vector<double>> matrix(ns,std::vector<double>(n,0.0));
    for(std::size_t s=0;s<ns;++s) {
        const auto& solid=candidates[s];
        if(solid.id.empty() || !names.insert(solid.id).second ||
           solid.componentStoichiometry.empty() ||
           !std::isfinite(solid.standardMuOverRT))
            throw std::invalid_argument("Invalid/duplicate solid definition");
        for(const auto& [name,count]:solid.componentStoichiometry) {
            auto it=indices.find(name);
            if(it==indices.end()||count<=0)
                throw std::invalid_argument("Invalid solid stoichiometry: "+solid.id);
            matrix[s][it->second]=count;
        }
    }
    std::vector<double> start(n);
    for(std::size_t i=0;i<n;++i)
        start[i]=base.componentPotentialOverRT.at(fluidModel.components[i].id);
    const auto affinity=[&](const std::vector<double>& x,std::size_t s) {
        double value=-candidates[s].standardMuOverRT;
        for(std::size_t i=0;i<n;++i)value+=matrix[s][i]*x[i];
        return value;
    };
    double totalBound=0.0;
    for(double b:prepared.totals)totalBound+=b;
    const double qTolerance=std::max(1e-13,1e-8*totalBound);
    // Enumerate every subset of <=8 candidate solids, solving the equality
    // system and accepting only thermodynamically complementary solutions.
    for(std::size_t mask=0;mask<(std::size_t(1)<<ns);++mask) {
        std::vector<std::size_t> active;
        for(std::size_t s=0;s<ns;++s)
            if(mask&(std::size_t(1)<<s))active.push_back(s);
        if(active.size()>n)continue;
        const std::size_t k=active.size();
        std::vector<double> x=start,q(k,0.0);
        Calculation currentFluid;
        std::vector<double> residual;
        double norm=0.0;
        const auto assess=[&](const std::vector<double>& xx,
                              const std::vector<double>& qq,
                              Calculation& fluid,
                              std::vector<double>& residualOut,
                              double& normOut) {
            if(!evaluate(prepared,xx,fluid))return false;
            residualOut=fluid.residual;
            for(std::size_t j=0;j<k;++j)
                for(std::size_t i=0;i<n;++i)
                    residualOut[i]+=matrix[active[j]][i]*qq[j]/prepared.totals[i];
            for(std::size_t j=0;j<k;++j)
                residualOut.push_back(affinity(xx,active[j]));
            normOut=0.0;
            for(double v:residualOut) {
                if(!std::isfinite(v))return false;
                normOut=std::max(normOut,std::abs(v));
            }
            return true;
        };
        if(!assess(x,q,currentFluid,residual,norm))continue;
        int iterations=0;
        while(norm>tolerance && iterations<maxIterations) {
            std::vector<std::vector<double>> jac(n+k,std::vector<double>(n+k,0));
            for(std::size_t s=0;s<prepared.nu.size();++s)
                for(std::size_t i=0;i<n;++i)
                    for(std::size_t j=0;j<n;++j)
                        jac[i][j]+=(prepared.nu[s][i]*prepared.nu[s][j]*
                                    currentFluid.moles[s])/prepared.totals[i];
            for(std::size_t j=0;j<k;++j) {
                for(std::size_t i=0;i<n;++i) {
                    jac[i][n+j]=matrix[active[j]][i]/prepared.totals[i];
                    jac[n+j][i]=matrix[active[j]][i];
                }
            }
            auto rhs=residual;
            for(auto& v:rhs)v=-v;
            std::vector<double> delta;
            if(!solveLinear(jac,rhs,delta))break;
            bool accepted=false;
            for(unsigned trial=0;trial<55;++trial) {
                const double alpha=std::ldexp(1.0,-static_cast<int>(trial));
                auto xx=x,qq=q;
                for(std::size_t i=0;i<n;++i)xx[i]+=alpha*delta[i];
                for(std::size_t j=0;j<k;++j)qq[j]+=alpha*delta[n+j];
                Calculation ff;std::vector<double> rr;double nn=0.0;
                if(assess(xx,qq,ff,rr,nn)&&nn<norm) {
                    x.swap(xx);q.swap(qq);
                    currentFluid=std::move(ff);residual=std::move(rr);norm=nn;
                    accepted=true;
                    break;
                }
            }
            if(!accepted)break;
            ++iterations;
        }
        if(norm>tolerance)continue;
        bool complementary=true;
        for(double amount:q)
            if(!std::isfinite(amount) || amount < -qTolerance)
                complementary=false;
        for(std::size_t s=0;s<ns;++s) {
            const double a=affinity(x,s);
            const bool present=(mask&(std::size_t(1)<<s))!=0;
            if(!std::isfinite(a) ||
               (present && std::abs(a)>tolerance) ||
               (!present && a>tolerance))
                complementary=false;
        }
        if(!complementary)continue;
        output.fluidSolution=IdealEquilibriumResult{};
        output.fluidSolution.iterations=iterations;
        output.fluidSolution.maxNormalizedResidual=norm;
        for(std::size_t i=0;i<n;++i)
            output.fluidSolution.componentPotentialOverRT[fluidModel.components[i].id]=x[i];
        for(std::size_t s=0;s<fluidModel.species.size();++s) {
            output.fluidSolution.speciesMoles[fluidModel.species[s].id]=currentFluid.moles[s];
            output.fluidSolution.concentrationsMolM3[fluidModel.species[s].id]=
                currentFluid.moles[s]/prepared.volumes[s];
        }
        for(std::size_t i=0;i<n;++i)
            output.fluidSolution.componentResidualMol[fluidModel.components[i].id]=
                residual[i]*prepared.totals[i];
        for(std::size_t s=0;s<ns;++s) {
            output.solidMoles[candidates[s].id]=0.0;
            output.solidAffinityOverRT[candidates[s].id]=affinity(x,s);
        }
        for(std::size_t j=0;j<k;++j)
            output.solidMoles[candidates[active[j]].id]=std::max(0.0,q[j]);
        return output;
    }
    throw std::runtime_error(
        "No complementary ideal solid-phase equilibrium found; "
        "phase candidates may be inconsistent or the nonlinear solve may have failed");
}

} // namespace thermo
