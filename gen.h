#ifndef TOYJETS_GEN_H
#define TOYJETS_GEN_H

#include "common.h"
#include <armadillo>

arma::mat genJet(const jet& recojet, jet& jetout,
                                double ptsmear=0.05,  //multiplicative
                                double etasmear=0.05, //additive
                                double phismear=0.05, //additive
                                double psplit=0.10,   //probability
                                double pmerge=0.02,   //probability
                                double ppu=0.10,      //probability
                                double pmiss=0.10,    //probability
                                double mergethreshold=0.05); //dR

#endif 
