#include "common.h"

arma::vec get_jet_pts(const jet& j){
    arma::vec ans(j.nPart);
    for(unsigned i=0; i<j.nPart; ++i){
        ans(i) = j.particles[i].pt;
    }
    return ans;
}
