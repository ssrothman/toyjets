#ifndef TOYJETS_COMMON_H
#define TOYJETS_COMMON_H

struct jet{
    unsigned nPart;
    std::vector<float> eta, phi, pt;
    float sumpt;
};

constexpr float inv_sqrt_2pi = 1/sqrt(M_PI);
float normal_pdf(float x, float m, float s){
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

#endif
