#ifndef TOYJETS_GAUS_H
#define TOYJETS_GAUS_H

#include <vector>
#include "common.h"

/*
 * Create jet that looks like a bell curve
 * overwrites contents of eta, phi, pt vectors with jet constituents
 */
void gausJet(int nPart, jet& jetout,
                        float ptscale=1,
                        float angscale=0.4,
                        float ptsmear=0.1);

#endif
