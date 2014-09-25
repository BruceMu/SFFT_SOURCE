#ifndef EXACTSFFTBYOFDM_2D_H_INCLUDED
#define EXACTSFFTBYOFDM_2D_H_INCLUDED

#include <iostream>
#include <fftw3.h>
#include <Eigen/Dense>
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
    const unsigned char * src;
    fftw_complex * sparsefouries;
    int width;
	int height;
	int sparseK;
	int sampleDims;

	//使用C++矩阵运算库;
	MatrixXcd fourRowBase;
	MatrixXcd fourColBase;

public:
//    ExactSfftByOfdm_2d(Raw2D & raw);
	ExactSfftByOfdm_2d(char const * filename,int width,int height){  //测试用接口
		this->height = height;
		this->width = width;
		this->fourRowBase = MatrixXcd(this->sparseK,this->width);
		this->fourColBase = MatrixXcd(this->sparseK,this->height);
		this->src =(unsigned char *) fftw_malloc(sizeof(unsigned char)*width*height);
		ReadImage(src,filename,width*height,0);
		Inital(this->width,this->height);
	};

    ~ExactSfftByOfdm_2d();

    void setSparse(int kspase);
    void setSampleDims(int sampleT);
    void Normalization(char * srchar);//针对数字图像像素值做归一化;
    fftw_complex * FoldToBins(double * src,int Br,int Bc,int Tr,int Tc,bool samDirect);
    void BasicExact2DSfft(int C_LogN);
    void BasicEstFreq(const vector<fftw_complex *> sam_srcX,const vector<fftw_complex *> sam_srcY,bool isCol);

};


#endif // EXACTSFFTBYOFDM_2D_H_INCLUDED
