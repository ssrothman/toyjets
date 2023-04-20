#include "gaus.h"
#include "common.h"
#include <random>
#include <stdio.h>
#include <chrono>

static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());


void gausJet(unsigned nPart, jet& jetout,
                             float ptscale,
                             float angscale,
                             float ptsmear){
    jetout.sumpt = 0;
    jetout.nPart = nPart;
    jetout.eta.resize(nPart);
    jetout.phi.resize(nPart);
    jetout.pt.resize(nPart);

    std::normal_distribution<float> rng(0, angscale);
    std::normal_distribution<float> smear(1, ptsmear);

    for(unsigned i=0; i<nPart; ++i){
        jetout.eta[i] = rng(generator);
        jetout.phi[i] = rng(generator);
        jetout.pt[i] = std::max(normal_pdf(jetout.eta[i], 0, angscale) *
                                normal_pdf(jetout.phi[i], 0, angscale) * 
                                smear(generator), 
                                0.0001f);
        jetout.sumpt += jetout.pt[i];
    }

    float rescale = ptscale/jetout.sumpt;
    jetout.sumpt *= rescale;
    for (unsigned i=0; i<nPart; ++i){
        jetout.pt[i] *= rescale;
    }
}

