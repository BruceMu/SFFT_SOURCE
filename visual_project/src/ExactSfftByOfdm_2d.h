#ifndef EXACTSFFTBYOFDM_2D_H_INCLUDED
#define EXACTSFFTBYOFDM_2D_H_INCLUDED

#include <iostream>
#include <fftw3.h>

/*
struct TwoDFreElement
{
    long rowIndex;
    long colIndex;
    fftw_complex value;
};
*/

class ExactSfftByOfdm_2d
{
private:
    double * src;
    fftw_complex * sparsefouries;
    int length;
    int width;
	int sparseK;
	int sampleDims;
public:
    ExactSfftByOfdm_2d(Raw2D & raw);
	ExactSfftByOfdm_2d(const double * src,int length,int width){  //测试用接口
		this->src = src;
		this->length = length;
		this->width = width;
	}

    ~ExactSfftByOfdm_2d();

    void setSparse(int kspase);
    void setSampleDims(int sampleT);
    void Normalization(char * srchar);//针对数字图像像素值做归一化;
    fftw_complex * FoldToBins(double * src,int Br,int Bc,int Tr,int Tc,bool samDirect);
    void BasicExact2DSfft(int C_LogN);
    TwoDFreElement* BasicEstFreq(const vector<fftw_complex *> sam_srcX,const vector<fftw_complex *> sam_srcY,int T,bool isCol);

};


#endif // EXACTSFFTBYOFDM_2D_H_INCLUDED
