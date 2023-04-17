#include "gen.h"
#include <random>

static std::default_random_engine generator;
static std::uniform_real_distribution<float> prob(0.0f, 1.0f);
static std::normal_distribution<float> norm(0.0f, 1.0f);

void genJet(const jet& recojet, jet& jetout,
                                float ptsmear,
                                float etasmear,
                                float phismear,
                                float psplit,
                                float pmerge,
                                float ppu){
    jetout.sumpt=0;
    jetout.eta.resize(recojet.nPart);
    jetout.phi.resize(recojet.nPart);
    jetout.pt.resize(recojet.nPart);

    for(unsigned i=0; i<recojet.nPart; ++i){
        float r = prob(generator);
        if(r < ppu){
            continue;
        }

        r = 
    }
}
