#ifndef TOYJETS_COMMON_H
#define TOYJETS_COMMON_H

#include <vector>
#include <math.h>
#include "SRothman/armadillo-12.2.0/include/armadillo"

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
    double pt, eta, phi;
    unsigned iJet;
};

arma::vec get_jet_pts(const jet& j);

constexpr double inv_sqrt_2pi = 0.3989422804014327;
inline double normal_pdf(double x, double m, double s){
    double a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

#endif
