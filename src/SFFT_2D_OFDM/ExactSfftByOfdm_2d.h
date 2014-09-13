#ifndef EXACTSFFTBYOFDM_2D_H_INCLUDED
#define EXACTSFFTBYOFDM_2D_H_INCLUDED

#include <iostream>
#include <fftw3.h>

struct TwoDFreElement
{
    long rowIndex;
    long colIndex;
    fftw_complex value;
};

class ExactSfftByOfdm_2d
{
private:
    double * src;
    fftw_complex * sparsefouries;
    int length;
    int width;
public:
    int sparseK;
    int sampleDims;
    ExactSfftByOfdm_2d(Raw2D & raw);
    ~ExactSfftByOfdm_2d();

    void setSparse(int kspase);
    void setSampleDims(int sampleT);
    void Normalization(char * srchar);
    fftw_complex * FoldToBins(double * src,int Br,int Bc,int Tr,int Tc,bool samDir);
    void BasicExact2DSfft(fftw_complex * src);
    TwoDFreElement* BasicEstFreq(fftw_complex ** sam_srcX,fftw_complex ** sam_srcY,int T,bool isCol);

};


#endif // EXACTSFFTBYOFDM_2D_H_INCLUDED
