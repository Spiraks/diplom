#include <vector>
#include <cmath>
#include "FDTD.h"







float FDTD::GetSigma(float D, float dStep, float pmlG, float pmlR,int N) {
    return -(eps0 * c) / (2 * dStep) * log(pmlG) / (pow(pmlG, N) - 1) * log(pmlR) * pow(pmlG, D / dStep);
}

void FDTD::Pml(){
    PML_E();
    PML_H();
}


void FDTD::InitPML()
{
    for (size_t K = 0; K < plmLayerNumber; ++K) {
        pmlSigmaStarE[K] = this->GetSigma(K,dx,2.15,0.001,plmLayerNumber);
        pmlSigmaStarH[K] = (pmlSigmaStarE[K] *mu0)/eps0;
    }
}