#include "gaus.h"
#include "common.h"
#include <random>
#include <stdio.h>
#include <math>

static std::default_random_engine generator;


bool gausJet(unsigned nPart, jet& jetout,
                             float ptscale,
                             float angscale,
                             float ptsmear){
    jetout.sumpt = 0;
    jetout.nPart = nPart;
    jetout.eta.resize(nPart);
    jetout.phi.resize(nPart);
    jetout.pt.resize(nPart);

    rng = std::normal_distribution(0, angscale);
    smear = std::normal_distribution(1, ptsmear);

    for(unsigned i=0; i<nPart; ++i){
        eta[i] = rng(generator);
        phi[i] = rng(generator);
        pt[i] = normal_pdf(eta[i]) * normal_pdf(phi[i]) * smear(generator);
        sumpt += pt[i];
    }

    float rescale = ptscale;
    for (unsigned i=0; i<nPart; ++i){
        pt[i] *= rescale;
    }
}

