#include "gen.h"
#include <random>
#include <armadillo>
#include <stdio.h>
#include <chrono>

static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
static std::uniform_real_distribution<float> prob(0.0f, 1.0f);
static std::normal_distribution<float> norm(0.0f, 1.0f);

bool check(float p){
    return prob(generator) < p;
}

float square(float x){
    return x*x;
}

void makegenpart(float recopt, float recoeta, float recophi,
                 float& genpt, float& geneta, float& genphi,
                 float ptsmear, float etasmear, float phismear){
    genpt = std::max(recopt * (norm(generator) * ptsmear + 1), 0.0001f);
    geneta = recoeta + norm(generator) * etasmear;
    genphi = recophi + norm(generator) * phismear;
}

arma::fmat genJet(const jet& recojet, jet& jetout,
                                     float ptsmear,
                                     float etasmear,
                                     float phismear,
                                     float psplit,
                                     float pmerge,
                                     float ppu,
                                     float pmiss,
                                     float mergethreshold){
    std::geometric_distribution<int> geom(1-psplit);
    std::geometric_distribution<int> geom2(1-pmiss);

    jetout.sumpt=0;

    jetout.eta.clear();
    jetout.phi.clear();
    jetout.pt.clear();

    jetout.eta.reserve(recojet.nPart);
    jetout.phi.reserve(recojet.nPart);
    jetout.pt.reserve(recojet.nPart);

    std::vector<std::vector<unsigned>> from(0);
    for(unsigned i=0; i<recojet.nPart; ++i){//for each particle in the reco jet
        if(check(ppu)){ //PU particles have no gen counterpart
            continue;
        }

        //add particles to genjet
        int nsplit = geom(generator)+1;
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
            jetout.nPart += 1;
            jetout.sumpt += pt;
            from.emplace_back(std::vector<unsigned>(1, i));
        }
    }//end for each particle in the reco jet
    
    //do merge
    std::vector<bool> todelete(jetout.nPart, false);
    for(unsigned i=0; i<jetout.nPart-1; ++i){
        if(todelete[i]){
            continue;
        }
        for(unsigned j=i+1; j<jetout.nPart; ++j){
            if(todelete[j]){
                continue;
            }
            float dR = std::sqrt(square(jetout.eta[i] - jetout.eta[j]) +
                                 square(jetout.phi[i] - jetout.phi[j]) );
            if(dR < mergethreshold && check(pmerge)){
                float pt = jetout.pt[i] + jetout.pt[j];
                float pt1 = jetout.pt[i]/pt;
                float pt2 = jetout.pt[j]/pt;
                jetout.pt[i] = pt;
                jetout.eta[i] = jetout.eta[i]*pt1 + jetout.eta[j]*pt2;
                jetout.phi[i] = jetout.phi[i]*pt1 + jetout.phi[j]*pt2;
                from[i].emplace_back(from[j][0]);

                todelete[j] = true;
            }
        }
    }

    //delete merged particles
    for(int i=jetout.nPart-1; i>=0; --i){
        if(todelete[i]){
            jetout.pt.erase(jetout.pt.begin()+i);
            jetout.eta.erase(jetout.eta.begin()+i);
            jetout.phi.erase(jetout.phi.begin()+i);
            from.erase(from.begin()+i);
            jetout.nPart-=1;
        }
    }

    //add non-reconstructed particles
    int nmiss = geom2(generator);
    for(unsigned i=0; i<nmiss; ++i){
        float nextPt = std::max(norm(generator)+1, 0.0f);
        jetout.pt.push_back(nextPt);
        jetout.eta.push_back(norm(generator));
        jetout.phi.push_back(norm(generator));
        from.push_back(std::vector<unsigned>(0));
        jetout.nPart+=1;
        jetout.sumpt += nextPt;
    }

    //make particle transfer matrix
    arma::fmat result(recojet.nPart, jetout.nPart, arma::fill::zeros);
    for(unsigned iGen=0; iGen<jetout.nPart; ++iGen){
        auto ip = std::unique(from[iGen].begin(), from[iGen].end());
        from[iGen].resize(std::distance(from[iGen].begin(), ip));
        for(unsigned iReco : from[iGen]){
            result(iReco, iGen) = 1;
        }
    }

    arma::fvec rowsum = arma::sum(result, 1);
    rowsum.replace(0.0f, 1.0f);
    arma::fvec rowfactor = arma::fvec(recojet.pt) / rowsum;
    result.each_col() %= rowfactor;

    arma::frowvec colfactor(jetout.pt);
    colfactor.replace(0.0f, 1.0f);
    result.each_row() /= colfactor;

    result *= jetout.sumpt / recojet.sumpt;

    return result;
}//end genJet()
