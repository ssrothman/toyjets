#include <random>
#include <armadillo>
#include <stdio.h>
#include <chrono>
#include "gen.h"

static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
static std::uniform_real_distribution<double> prob(0.0f, 1.0f);
static std::normal_distribution<double> norm(0.0f, 1.0f);

bool check(double p){
    return prob(generator) < p;
}

double square(double x){
    return x*x;
}

void makegenpart(double recopt, double recoeta, double recophi,
                 double& genpt, double& geneta, double& genphi,
                 double ptsmear, double etasmear, double phismear){
    genpt = std::max(recopt * (norm(generator) * ptsmear + 1), 0.0001);
    geneta = recoeta + norm(generator) * etasmear;
    genphi = recophi + norm(generator) * phismear;
}

arma::mat genJet(const jet& recojet, jet& jetout,
                                     double ptsmear,
                                     double etasmear,
                                     double phismear,
                                     double psplit,
                                     double pmerge,
                                     double ppu,
                                     double pmiss,
                                     double mergethreshold){
    std::geometric_distribution<int> geom(1-psplit);
    std::geometric_distribution<int> geom2(1-pmiss);

    jetout.sumpt=0;

    jetout.particles.clear();
    jetout.particles.reserve(recojet.nPart);

    std::vector<std::vector<unsigned>> from(0);
    for(unsigned i=0; i<recojet.nPart; ++i){//for each particle in the reco jet
        if(check(ppu)){ //PU particles have no gen counterpart
            continue;
        }

        //add particles to genjet
        unsigned nsplit = geom(generator)+1;
        std::vector<double> fracs(nsplit, 0);
        double sumfrac=0;
        for(unsigned isplit=0; isplit<nsplit; ++isplit){
            fracs[isplit] = prob(generator);
            sumfrac += fracs[isplit];
        }
        for(unsigned isplit=0; isplit<nsplit; ++isplit){
            fracs[isplit]/=sumfrac;

            double pt, eta, phi;
            makegenpart(recojet.particles[i].pt*fracs[isplit], recojet.particles[i].eta, recojet.particles[i].phi,
                        pt, eta, phi,
                        ptsmear, etasmear, phismear);
            jetout.particles.emplace_back(pt, eta, phi, 0, 0, 0, 0, 0);
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
            double dR = std::sqrt(square(jetout.particles[i].eta - jetout.particles[j].eta) +
                                 square(jetout.particles[i].phi - jetout.particles[j].phi));
            if(dR < mergethreshold && check(pmerge)){
                double pt = jetout.particles[i].pt + jetout.particles[j].pt;
                double pt1 = jetout.particles[i].pt/pt;
                double pt2 = jetout.particles[j].pt/pt;
                jetout.particles[i].pt = pt;
                jetout.particles[i].eta = jetout.particles[i].eta*pt1 + jetout.particles[j].eta*pt2;
                jetout.particles[i].phi = jetout.particles[i].phi*pt1 + jetout.particles[j].phi*pt2;
                from[i].emplace_back(from[j][0]);

                todelete[j] = true;
            }
        }
    }

    //delete merged particles
    for(int i=jetout.nPart-1; i>=0; --i){
        if(todelete[i]){
            jetout.particles.erase(jetout.particles.begin()+i);
            from.erase(from.begin()+i);
            jetout.nPart-=1;
        }
    }

    //add non-reconstructed particles
    unsigned nmiss = geom2(generator);
    for(unsigned i=0; i<nmiss; ++i){
        double nextPt = std::max(0.1, 0.0);
        jetout.particles.emplace_back(nextPt, norm(generator), norm(generator), 
                            0, 0, 0, 0, 0);
        from.push_back(std::vector<unsigned>(0));
        jetout.nPart+=1;
        jetout.sumpt += nextPt;
    }

    //make particle transfer matrix
    arma::mat result(recojet.nPart, jetout.nPart, arma::fill::zeros);
    for(unsigned iGen=0; iGen<jetout.nPart; ++iGen){
        auto ip = std::unique(from[iGen].begin(), from[iGen].end());
        from[iGen].resize(std::distance(from[iGen].begin(), ip));
        for(unsigned iReco : from[iGen]){
            result(iReco, iGen) = 1;
        }
    }

    arma::vec rowsum = arma::sum(result, 1);
    rowsum.replace(0.0f, 1.0f);
    arma::vec rowfactor(recojet.nPart, arma::fill::none); 
    for(unsigned i=0; i<recojet.nPart; ++i){
        rowfactor(i) = recojet.particles[i].pt;
    }
    rowfactor/=rowsum;
    result.each_col() %= rowfactor;

    arma::rowvec colfactor(jetout.nPart);
    for(unsigned i=0; i<jetout.nPart; ++i){
        colfactor(i) = jetout.particles[i].pt;
    }
    colfactor.replace(0.0f, 1.0f);
    result.each_row() /= colfactor;

    result *= jetout.sumpt / recojet.sumpt;

    return result;
}//end genJet()
