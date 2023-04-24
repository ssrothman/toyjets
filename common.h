#ifndef TOYJETS_COMMON_H
#define TOYJETS_COMMON_H

#include <vector>
#include <math.h>
#include <armadillo>

struct particle{
    double pt, eta, phi;
    double dpt, deta, dphi;
    unsigned pdgid; //absolute value
    int charge;

    particle(double pt, double eta, double phi,
             double dpt, double deta, double dphi,
             unsigned pdgid, int charge):
        pt(pt), eta(eta), phi(phi),
        dpt(dpt), deta(deta), dphi(dphi),
        pdgid(pdgid), charge(charge) {}

    particle() :
        pt(0), eta(0), phi(0),
        dpt(0), deta(0), dphi(0), 
        pdgid(0), charge(0) {}
};


struct jet{
    unsigned nPart;
    std::vector<particle> particles;
    double sumpt;

    arma::vec ptvec(){
        arma::vec ans(nPart, arma::fill::none);
        for(unsigned i=0; i<nPart; ++i){
            ans(i) = particles[i].pt;
        }
        return ans;
    }
};

constexpr double inv_sqrt_2pi = 1/sqrt(M_PI);
inline double normal_pdf(double x, double m, double s){
    double a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

#endif
