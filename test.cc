#include <armadillo>
#include <iostream>
#include <stdio.h>
#include <vector>
#include "gen.h"
#include "gaus.h"
#include "common.h"

int main(){
    jet result;
    unsigned npart = 20;
    gausJet(npart, result);
    printf("Made a jet with %d constituents\n", result.nPart);

    jet genjet;
    arma::fmat ans = genJet(result, genjet, 0.05, 0.05, 0.05, 0.2, 0.9, 0.05, 0.15);
    printf("Made corresponding genjet with %d consistuents\n", genjet.nPart);

    arma::fvec genpt(genjet.pt);
    arma::frowvec recopt(result.pt);
    std::cout << ans;

    printf("GEN PT\n");
    std::cout << arma::trans(genpt);
    printf("RECO PT\n");
    std::cout << recopt;
    printf("MAT * GEN\n");
    std::cout << arma::trans(ans * genpt);
    printf("ROWSUM\n");
    std::cout << arma::trans(arma::sum(ans, 1));
        
    return 0;
}
