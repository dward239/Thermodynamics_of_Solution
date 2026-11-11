#include "decay_network.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
namespace thermo {
namespace {
constexpr double kLn2=0.693147180559945309417232121458;
constexpr double kTolerance=1e-10;
constexpr double kMuPerChunk=12.0;
constexpr double kMaxChunks=10000.0;
void nonnegative(double v, const char* name) {
    if (!std::isfinite(v) || v<0.0) throw std::invalid_argument(std::string(name)+" must be finite and nonnegative");
}
}
DecayNetwork::DecayNetwork(std::vector<NuclideDefinition> definitions)
    : definitions_(std::move(definitions)) {
    if (definitions_.empty()) throw std::invalid_argument("Empty nuclide network");
    std::map<std::string,std::size_t> ids;
    for (std::size_t i=0;i<definitions_.size();++i) {
        const auto& n=definitions_[i];
        if(n.id.empty() || n.element.empty() || !ids.emplace(n.id,i).second)
            throw std::invalid_argument("Duplicate/empty nuclide or element id");
        nonnegative(n.halfLifeSeconds,"half life");
        nonnegative(n.emittedEnergyJPerDecay,"emitted energy per decay");
        if(n.halfLifeSeconds==0.0 && !n.daughters.empty())
            throw std::invalid_argument("Stable nuclide has daughter branches: "+n.id);
        if(n.halfLifeSeconds>0.0 && n.daughters.empty())
            throw std::invalid_argument("Unstable nuclide has no daughter branches: "+n.id);
        double sumDeposit=0.0;
        for(const auto& [phase,fraction]:n.depositionFractions) {
            if (phase.empty()) throw std::invalid_argument("Empty deposition phase");
            nonnegative(fraction,"deposition fraction");
            sumDeposit+=fraction;
        }
        if(!std::isfinite(sumDeposit) || sumDeposit>1.0+kTolerance)
            throw std::invalid_argument("Deposition fractions exceed unity: "+n.id);
    }
    for(const auto& n:definitions_) {
        double sum=0.0;
        std::set<std::string> seen;
        for(const auto& d:n.daughters) {
            if(!ids.count(d.daughter) || !seen.insert(d.daughter).second)
                throw std::invalid_argument("Unknown/duplicate daughter of "+n.id);
            nonnegative(d.fraction,"daughter branch fraction");
            if(d.fraction==0.0) throw std::invalid_argument("Zero branch fraction");
            sum+=d.fraction;
        }
        if(n.halfLifeSeconds>0.0 && (!std::isfinite(sum) || std::abs(sum-1.0)>kTolerance))
            throw std::invalid_argument("Branch probabilities do not sum to one: "+n.id);
    }
}
DecayStep DecayNetwork::advance(const std::map<std::string,double>& initialAtoms,
                                double durationSeconds) const {
    nonnegative(durationSeconds,"duration");
    const std::size_t n=definitions_.size();
    std::map<std::string,std::size_t> ids;
    for(std::size_t i=0;i<n;++i)ids[definitions_[i].id]=i;
    std::vector<double> state(2*n,0.0); // [atoms ..., cumulative decays ...]
    double totalInitially=0.0;
    for(const auto& [name,atoms]:initialAtoms) {
        nonnegative(atoms,"initial atoms");
        const auto it=ids.find(name);
        if(it==ids.end()) throw std::invalid_argument("Unknown inventory nuclide: "+name);
        state[it->second]=atoms;
        totalInitially+=atoms;
    }
    if(!std::isfinite(totalInitially)) throw std::invalid_argument("Inventory overflow");
    std::vector<double> lambda(n,0.0);
    double maxLambda=0.0;
    for(std::size_t i=0;i<n;++i) {
        if(definitions_[i].halfLifeSeconds>0) {
            lambda[i]=kLn2/definitions_[i].halfLifeSeconds;
            if(!std::isfinite(lambda[i]))
                throw std::invalid_argument("Half-life below representable decay rate");
            maxLambda=std::max(maxLambda,lambda[i]);
        }
    }
    if(durationSeconds>0.0 && maxLambda>0.0) {
        const double events=maxLambda*durationSeconds;
        if(!std::isfinite(events) || std::ceil(events/kMuPerChunk)>kMaxChunks)
            throw std::invalid_argument("Decay step too stiff for bounded uniformization; shorten dt");
        const auto chunks=static_cast<unsigned>(std::max(1.0,std::ceil(events/kMuPerChunk)));
        const double mu=events/chunks;
        for(unsigned chunk=0;chunk<chunks;++chunk) {
            std::vector<double> term=state;
            std::vector<double> sum(2*n,0.0);
            double weight=std::exp(-mu), cdf=0.0;
            bool converged=false;
            for(unsigned k=0;k<256;++k) {
                for(std::size_t i=0;i<2*n;++i)sum[i]+=weight*term[i];
                cdf+=weight;
                if(k>=static_cast<unsigned>(mu) && 1.0-cdf<=5e-15) {
                    converged=true;
                    break;
                }
                // P = I + generator/maxLambda: stochastic inventory transition,
                // augmented with cumulative event counts per nuclide.
                std::vector<double> next(2*n,0.0);
                for(std::size_t i=0;i<n;++i) {
                    const double p=lambda[i]/maxLambda;
                    next[i]+=(1.0-p)*term[i];
                    next[n+i]+=term[n+i]+p*term[i];
                    for(const auto& d:definitions_[i].daughters)
                        next[ids.at(d.daughter)]+=p*d.fraction*term[i];
                }
                term.swap(next);
                weight*=mu/(k+1.0);
            }
            if(!converged) throw std::runtime_error("Decay uniformization failed to converge");
            state.swap(sum);
        }
    }
    DecayStep out;
    double totalFinal=0.0, depositedTotal=0.0;
    for(std::size_t i=0;i<n;++i) {
        if(!std::isfinite(state[i]) || !std::isfinite(state[n+i]) ||
           state[i]<-1e-12 || state[n+i]<-1e-12)
            throw std::runtime_error("Decay produced invalid inventory/count");
        const double atoms=std::max(0.0,state[i]);
        const double decays=std::max(0.0,state[n+i]);
        totalFinal+=atoms;
        out.finalAtoms[definitions_[i].id]=atoms;
        out.decaysByNuclide[definitions_[i].id]=decays;
        const double energy=decays*definitions_[i].emittedEnergyJPerDecay;
        out.emittedEnergyJ+=energy;
        for(const auto& [phase,fraction]:definitions_[i].depositionFractions) {
            const double part=energy*fraction;
            out.depositedEnergyJByPhase[phase]+=part;
            depositedTotal+=part;
        }
    }
    if(!std::isfinite(out.emittedEnergyJ) || !std::isfinite(depositedTotal) ||
       depositedTotal>out.emittedEnergyJ+1e-10*std::max(1.0,out.emittedEnergyJ))
        throw std::runtime_error("Invalid radioactive energy bookkeeping");
    if(!std::isfinite(totalFinal) ||
       std::abs(totalFinal-totalInitially)>1e-10*std::max(1.0,totalInitially))
        throw std::runtime_error("Decay inventory conservation failed");
    return out;
}
} // namespace thermo
