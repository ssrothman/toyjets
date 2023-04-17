#include "gen.h"
#include <random>
#include <armadillo>
#include "Minuit2/FunctionMinimum.h"

static std::default_random_engine generator;
static std::uniform_real_distribution<float> prob(0.0f, 1.0f);
static std::normal_distribution<float> norm(0.0f, 1.0f);

bool check(float p){
    return prob(generator) < p;
}

void makegenpart(float recopt, float recoeta, float recophi,
                 float& genpt, float& geneta, float& genphi,
                 float ptsmear, float etasmear, float phismear){
    genpt = std::max(recopt * (norm(generator) * ptsmear + 1), 0.0001f);
    geneta = recoeta + norm(generator) * etasmear;
    genphi = recophi + norm(generator) * phismear;
}

arma::mat genJet(const jet& recojet, jet& jetout,
                                     float ptsmear,
                                     float etasmear,
                                     float phismear,
                                     float psplit,
                                     float pmerge,
                                     float ppu){
    std::geometric_distribution<int> geom(1-psplit);

    jetout.sumpt=0;

    jetout.eta.clear();
    jetout.phi.clear();
    jetout.pt.clear();

    jetout.eta.reserve(recojet.nPart);
    jetout.phi.reserve(recojet.nPart);
    jetout.pt.reserve(recojet.nPart);

    std::vector<unsigned> from(0);
    for(unsigned i=0; i<recojet.nPart; ++i){
        if(check(ppu)){ //PU particles have no gen counterpart
            continue;
        }

        int nsplit = geom(generator);

        std::vector<float> fracs(nsplit, 0);
        float sumfrac=0;
        for(unsigned isplit=0; isplit<nsplit; ++isplit){
            fracs[isplit] = prob(generator);
            sumfrac += fracs[isplit];
        }
        for(unsigned isplit=0; isplit<nsplit; ++isplit){
            fracs[isplit]/=sumfrac;

            float pt, eta, phi;
            makegenpart(recojet.pt[i]*fracs[isplit], recojet.eta[i], recojet.phi[i],
                        pt, eta, phi,
                        ptsmear, etasmear, phismear);
            jetout.pt.emplace_back(pt);
            jetout.eta.emplace_back(eta);
            jetout.phi.emplace_back(phi);
            from.emplace_back(i);
        }
    }


}
