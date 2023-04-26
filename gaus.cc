#include "gaus.h"
#include <random>
#include <stdio.h>
#include <chrono>

static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());


void gausJet(unsigned nPart, jet& jetout,
                             double ptscale,
                             double angscale,
                             double ptsmear){
    jetout.sumpt = 0;
    jetout.nPart = nPart;
    jetout.particles.resize(nPart);

    std::normal_distribution<double> rng(0, angscale);
    std::normal_distribution<double> smear(1, ptsmear);

    for(unsigned i=0; i<nPart; ++i){
        jetout.particles[i].eta = rng(generator);
        jetout.particles[i].phi = rng(generator);
        jetout.particles[i].pt = std::max(normal_pdf(jetout.particles[i].eta, 0, angscale) *
                                normal_pdf(jetout.particles[i].phi, 0, angscale) * 
                                smear(generator), 
                                0.0001);
        jetout.particles[i].dpt = ptsmear * jetout.particles[i].pt;
        jetout.particles[i].dphi = 0.1;
        jetout.particles[i].deta = 0.1;
        jetout.sumpt += jetout.particles[i].pt;
    }

    double rescale = ptscale/jetout.sumpt;
    jetout.sumpt *= rescale;
    for (unsigned i=0; i<nPart; ++i){
        jetout.particles[i].pt *= rescale;
    }
}

