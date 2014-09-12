#ifndef EXACTSFFTBYOFDM_2D_H_INCLUDED
#define EXACTSFFTBYOFDM_2D_H_INCLUDED

#include <iostream>
#include <fftw3.h>


class ExactSfftByOfdm_2d
{
private:
    fftw_complex * src;
    fftw_complex * sparsefouries;
    int length;
    int width;
    int sparseK;
public:
    ExactSfftByOfdm_2d(Raw2D & raw);
    ~ExactSfftByOfdm_2d();

    void Normalization(char * srchar);
    void FoldToBins(fftw_complex * src,int Br,int Bc,int Tr,int Tc);
    void BasicExact2DSfft(fftw_complex * src,int k);
    void BasicEstFreq(fftw_complex * sam_srcX,fftw_complex * sam_srcY,int T,bool isCol);

};

#endif // EXACTSFFTBYOFDM_2D_H_INCLUDED
