#ifndef EXACTSFFTBYOFDM_2D_H_INCLUDED
#define EXACTSFFTBYOFDM_2D_H_INCLUDED

#include <iostream>


class ExactSfftByOfdm_2d
{
private:
    complex<float> * src;
    complex<float> * sparsefouries;
    int length;
    int width;
    int sparseK;
public:
    void ExactSfftByOfdm_2d();
    void ~ExactSfftByOfdm_2d();

    void Normalization(char * srchar);
    void FoldToBins(complex<float> * src,int Br,int Bc,int Tr,int Tc);
    void BasicExact2DSfft(complex<float> * src,int k);
    void BasicEstFreq(complex<float> * sam_srcX,complex<float> * sam_srcY,int T,bool isCol);

};

#endif // EXACTSFFTBYOFDM_2D_H_INCLUDED
