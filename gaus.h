#ifndef TOYJETS_GAUS_H
#define TOYJETS_GAUS_H

#include <vector>
#include "common.h"

/*
 * Create jet that looks like a bell curve
 * overwrites contents of eta, phi, pt vectors with jet constituents
 */
void gausJet(unsigned nPart, jet& jetout,
                             double ptscale=1,
                             double angscale=0.4,
                             double ptsmear=0.1);

#endif
