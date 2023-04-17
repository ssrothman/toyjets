#ifndef TOYJETS_GEN_H
#define TOYJETS_GEN_H

#include "common.h"

void genJet(const jet& recojet, jet& jetout,
                                float ptsmear=0.05,  //multiplicative
                                float etasmear=0.05, //additive
                                float phismear=0.05, //additive
                                float psplit=0.10,   //probability
                                float pmerge=0.02,   //probability
                                float ppu=0.10);     //probability

#endif 
