/************************************************************
*��ά��ɢ����Ҷ�任���˲���������˹��ͨ����ʵ��
*���ڼ���ͼ����ɢ���еĸ���ҶƵ��ͼ��������Ƶ��ͼ����
*�㷨������ʽ���ٸ���Ҷ�㷨�����Σ� ���ۻ������㷨����,ͼ����
*ʱ�临�Ӷȣ�һάO(NlogN)����άO(N^2)
*���ߣ�openCV��ȡͼ����չʾЧ��
*�汾������
*@auto:Bruce mu
************************************************************/
#include <stdio.h>
#include <iostream>
#include <complex>
#include <bitset>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <math.h>
#include "cv.h"
#include "highgui.h"
#define pi 3.1415926535
#define DELTA 8.0
using std::iostream;
using std::bitset;
using std::complex;
using namespace std;

/*******Ϊ�Ե����ϵĵ����������м���λ��************/
int rev(int k,int n)
{
     bitset<32> tmp(k);
     bitset<32> dist;
     for(int m = 0;m<n;m++)
     {
        if(tmp.test(m))
        {
            dist.set(n-1-m);
        }
     }
     int revk = (int)dist.to_ulong();
     return revk;
}
//������ʽ
complex<double>* FFT(const complex<double> *srcimg,int n)
{
    double flag = log10((double)n)/log10(2.0);
    int N = n;
    if(flag - (int)flag != 0){   //�ж��Ƿ�Ϊ2��ָ����
        cout <<"the length of srcimg is wrong"<< endl;
        /*����2��ָ����*/
        cout <<"need padding"<<endl;
        N = pow(2,(double)((int)flag+1));
        flag = log10((double)N)/log10(2.0);
    }
    /*�ı����ż˳��*/
    complex<double> *arr= new complex<double>[N];
    int sub;
    for(int k =0;k<N;k++)
    {
        sub =rev(k,(int)flag); 
        if(sub <= n-1){
            arr[k] = *(srcimg + sub);
        }else{
            complex<double> t = complex<double>(0,0);
            arr[k] = t;
        }
    }
     for(int s =1;s <= flag;s++)
    {
        int m = pow(2,(double)s);
        complex<double> wm = complex<double>(cos(2*pi/m),sin(2*pi/m));//wm1����W21��ʼ,���ڱ任
        for(int p = 0;p<N;p+=m)
        {
            complex<double> w(1,0);
            for(int j = 0;j<=m/2-1;j++)
            {
                complex<double> t = w * arr[p+j+m/2];
                complex<double> u = arr[p+j];
                arr[p+j] = u+t;
                arr[p+j+m/2] = u-t;
                w = w*wm;
            }
        }
    }
    return arr;
}
/***********һά���ٸ���Ҷ�任********************
*srcimg : ԭʼһά����                          *
*n     ��һά���еĳ���
*************************************************/
complex<double>* FFT(const double *srcimg,int n)
{
    double flag = log10((double)n)/log10(2.0);
    int N = n;
    if(flag - (int)flag != 0){   //�ж��Ƿ�Ϊ2��ָ����
//        cout <<"the length of srcimg is wrong"<< endl;
        /*����2��ָ����*/
//        cout <<"need padding"<<endl;
        N = pow(2,(double)((int)flag+1));
        flag = log10((double)N)/log10(2.0);
    }
    /*�ı����ż˳��*/
    complex<double> *arr= new complex<double>[N];
    int sub;
    for(int k =0;k<N;k++)
    {
        sub =rev(k,(int)flag); 
        if(sub <= n-1){
            arr[k] = complex<double>(*(srcimg + sub),0);
        }else{
            complex<double> t = complex<double>(0,0);
            arr[k] = t;
        }
    }
//    cout<<"------------after padding and retrival-----------------"<<endl;
//    for(size_t y=0;y<N;y++)
//    {
//        cout << arr[y].real()<<"  "<<arr[y].imag()<<endl;
//    }
    /*���ڵ����ĵ��ο��ٸ���Ҷ�任���Ե�����*/
     for(int s =1;s <= flag;s++)
    {
        int m = pow(2,(double)s);
        complex<double> wm = complex<double>(cos(2*pi/m),sin(2*pi/m));//wm1����W21��ʼ,���ڱ任
        for(int p = 0;p<N;p+=m)
        {
            complex<double> w(1,0);
            for(int j = 0;j<=m/2-1;j++)
            {
                complex<double> t = w * arr[p+j+m/2];
                complex<double> u = arr[p+j];
                arr[p+j] = u+t;
                arr[p+j+m/2] = u-t;
                w = w*wm;
            }
        }
    }
    return arr;
}

int countPadding(int n);
/*****************һά���ٸ���Ҷ��任********************/
/*fftimg:ԭʼһά����Ҷ����
  n     : ���г���
*******************************************************/
complex<double>* IFFT(const complex<double> *fftimg,int n)
{
    n = countPadding(n);
    double flag = log10((double)n)/log10(2.0);
    int N = n;
    if(flag - (int)flag != 0){   //�ж��Ƿ�Ϊ2��ָ����
        cout <<"the length of srcimg is wrong"<< endl;
        /*����2��ָ����*/
        cout <<"need padding"<<endl;
        N = pow(2,(double)((int)flag+1));
        flag = log10((double)N)/log10(2.0);
    }
    /*�ı����ż˳��*/
    complex<double> * spit = new complex<double>[N];
    int sub=0;
    for(int k =0;k<N;k++)
    {
        sub =rev(k,(int)flag); 
        if(sub < n){
            spit[k] = complex<double>(*(fftimg + sub));
        }else{
            spit[k] = complex<double>(0,0);
        }
    }

    for(int s =1;s <= flag;s++)
    {
        int m = pow(2,(double)s);
        complex<double> wm = complex<double>(cos(-2*pi/m),sin(-2*pi/m));//wm1����W2(-1)��ʼ
        for(int p = 0;p<N;p+=m)
        {
            complex<double> w(1,0);
            for(int j = 0;j<=m/2-1;j++)
            {
                complex<double> t = w * spit[p+j+m/2];
                complex<double> u = spit[p+j];
                spit[p+j] = u+t;
                spit[p+j+m/2] = u-t;
                w = w*wm;
            }
        }
    }

    for(size_t p =0;p<n;p++)
    {
        spit[p] = spit[p]/complex<double>(N,0);
    }
    return spit;
}
  
/*******ʹ�ù���Ŀ��ٸ���Ҷ��任**************/
complex<double>* IFFT2(const complex<double> *fftimg,int n)
{
    n = countPadding(n);
    complex<double> *gfftimg = new complex<double>[n];
    for(size_t i = 0;i<n;i++){
        gfftimg[i] = complex<double>(fftimg[i].real(),-fftimg[i].imag());
    }
    complex<double> *ifft = FFT(gfftimg,n); 
    for(size_t j = 0;j<n;j++)
    {
        ifft[j] = ifft[j]/complex<double>(n,0);
    }
    delete gfftimg;
    return ifft;
}

/*************��ά���ٸ���Ҷ�任**************************
*srcimg: ��һά��ʾ�Ķ�άԭʼ����
*width �����
*height���߶�
********************************************************/
complex<double>* twoDFFT(const double *srcimg,int width,int height)
{
    int w = countPadding(width);
    int h = countPadding(height);
    int pixes = w * h;
    complex<double> *hdirection = new complex<double>[w];
    complex<double> *vdirection = new complex<double>[h];
    complex<double> *fourier = new complex<double>[pixes];
    /*��άˮƽ����*/
    for(size_t i = 0;i<h;i++){
        for(size_t j = 0;j<w;j++){
            if(i>=height || j >=width){
                hdirection[j] = complex<double>(0,0);
            }else{
                hdirection[j] = complex<double>(srcimg[i*width + j],0);
            }
    //        cout << hdirection[j] << " ";
        }
    //    cout <<""<<endl;
        complex<double> *hfourier = FFT(hdirection,w);
        for(size_t m = 0;m<w;m++){
            fourier[i*w+m] = hfourier[m];
        }
        delete hfourier;
    }
    /*��ά��ֱ����*/
    for(size_t ii = 0;ii<w;ii++){
        for(size_t jj = 0;jj<h;jj++){
            vdirection[jj] = fourier[jj*w + ii];
        }
        complex<double> *vfourier = FFT(vdirection,h);
        for(size_t mm = 0;mm < h;mm++){
            fourier[mm*w +ii] = vfourier[mm];
        }
        delete vfourier;
    }
    delete hdirection;
    delete vdirection;
    return fourier;

}
/**************��ά���ٸ���Ҷ��任*************************
*fourier : һά��ʾ�Ķ�ά����Ҷ�任����                     *
*width   :��                                             *
*height  :��                                             *
***********************************************************/
complex<double>* twoDIFFT(const complex<double> *fourier,int width,int height)
{
    width = countPadding(width);
    height = countPadding(height);
    int fpoints = width * height;
    complex<double> *hdirection = new complex<double>[width];
    complex<double> *vdirection = new complex<double>[height];
    complex<double> *ifourier = new complex<double>[fpoints];

    for(size_t ii = 0;ii<height;ii++)
    {
        for(size_t jj = 0;jj<width;jj++){
            hdirection[jj] = fourier[ii*width+jj];
        }
        complex<double> *hifour = IFFT(hdirection,width);//��ʱ����
        for(size_t mm = 0;mm<width;mm++){
            ifourier[ii*width+mm] = hifour[mm];
        }
        delete hifour;
    }
    for(size_t i = 0;i<width;i++){
        for(size_t j = 0;j<height;j++){
            vdirection[j] = ifourier[j*width+i];
        }
        complex<double> *vifour = IFFT(vdirection,height);
        for(size_t m = 0;m<height;m++){
            ifourier[m*width+i] = vifour[m];
        }
        delete vifour;
    }    
    delete hdirection;
    delete vdirection;
    return ifourier;
}
/******************���������********************************************
*���θ���Ҷ�任�㷨ֻ����2��ָ���ε���ɢ���У����ڷ�2��ָ���ε����У�ʹ�������*
***********************************************************************/
inline int countPadding(int n)
{
    double lg = log10((double)n)/log10(2.0);
    if((lg - (int)lg) == 0){
        return n;
    }
    int N = pow(2.0,((int)lg+1));
    return N;
}

/*****��˹��ͨ�˲�����*************************
*src��   ����Ƶ��
*width�� ����Ƶ�׿��
*height������Ƶ�׸߶�
*D��     ��˹����������˲���ֵ
*ֻ������λ���к�ĸ���ҶƵ�׽��и�˹��ͨ�˲�
*/
void guass(complex<double> *src,int width,int height,double D)
{
    //find centre point
    int orgx = floor(width/2.0);
    int orgy = floor(height/2.0);
    double distence = 0;
    for(size_t i = 0;i<height;i++)
    {
        for(size_t j = 0;j<width;j++)
        {
            distence = sqrt(pow(abs((int)(i-orgy)),2.0)+pow(abs((int)(j-orgx)),2.0));
            src[i*width+j] = src[i*width+j] * exp(-distence/(2*pow(D,2.0)));
            
        }
        //cout<<i<<" "<<distence<<" "<<endl;
        
    }
//    cout<<orgx<<" "<<orgy<<endl;
}
/************��������ҶƵ������ת����256���Ҷ�����*****************/
void toGray(complex<double> *twodfourier,int pwid,int phei,char* pixval)
{
        double *vals = new double[pwid*phei];
        double max = 0;
        double min = 255;
        for(size_t p = 0;p<pwid*phei;p++){
            vals[p]=log(1+sqrt(pow(twodfourier[p].real(),2.0)+pow(twodfourier[p].imag(),2.0)));//�������ķ�����
            if(vals[p] > max){
                max = vals[p];
            }
            if(vals[p] < min){
                min = vals[p];
            }
        }
        cout<<min<< " "<<max<<endl;
        cout<<pwid<<" "<<phei<<endl;
        for(size_t s = 0;s<pwid*phei;s++){
            pixval[s] = (char)((vals[s]-min)/(max-min)*255);
        }
        delete vals;
}
/******************opencv ��ȡͼ����չʾЧ��***********************/
int StatisFreqFreqWithUdefine(int argc,char **argv)
{
    IplImage *img;
    if((img = cvLoadImage("E:\\code\\SFFT_SOURCE\\doc\\lena.jpg",0))!=0){
        int dim = img->imageSize;
        //��ͼ������Ƴ���ά���飬����Ƶ�׾��в�������Ӧ��д�ӿڣ����ص�ά���飻
        double * src = new double[dim];
        size_t j =0;
        for(int y =0;y<img->height;y++)
        {
            uchar * ptr = (uchar*)(img->imageData + y * img->widthStep);
            for(int x =0;x<img->width;x++){
                /***���뺯�����ԣ�-1���ģ�x+y���η���Ƶ�׽�����*/
                src[j] = (double)ptr[x]*pow(-1,(double)(x+y))/256;
                j++;
            }
        }
        int w = img->width;
        int h = img->height;
        int pwid = countPadding(w);
        int phei = countPadding(h);

        complex<double> *twodfourier = twoDFFT(src,w,h);
        char * pixval = new char[pwid*phei];
        CvMat freq;
        toGray(twodfourier,pwid,phei,pixval);//����Ƶ��ת256�Ҷ�
        cvInitMatHeader(&freq,pwid,phei,CV_8UC1,pixval);
        IplImage *ipl = cvCreateImage(cvGetSize(&freq),8,1);
        cvGetImage(&freq,ipl);//��ȡƵ��ͼ��

        guass(twodfourier,pwid,phei,DELTA);//����DELTA�ĸ�˹ƽ������˹��ͨ�˲���
        CvMat gaufre;
        char *pixvals = new char[pwid*phei];
        toGray(twodfourier,pwid,phei,pixvals);
        cvInitMatHeader(&gaufre,pwid,phei,CV_8UC1,pixvals);
        IplImage *gausf = cvCreateImage(cvGetSize(&gaufre),8,1);
        cvGetImage(&gaufre,gausf);//��˹ƽ�����Ƶ��ͼ��

        complex<double> *twodifourier = twoDIFFT(twodfourier,w,h);
        double ipix = 0;
        size_t po = 0;
        for(int y =0;y<img->height;y++)
        {
            uchar * ptr = (uchar*)(img->imageData + y * img->widthStep);
            for(int x =0;x<img->width;x++){
                ipix = twodifourier[po*pwid+x].real();
                ptr[x] = (uchar)(ipix * 256)*pow(-1,(double)(x+y));
            }
            po++;
        }
        cvNamedWindow("frequence_domain",CV_WINDOW_AUTOSIZE);
        cvShowImage("frequence_domain",ipl);
        cvNamedWindow("gauss_fre",CV_WINDOW_AUTOSIZE);
        cvShowImage("gauss_fre",gausf);
        cvNamedWindow("fourier_t",CV_WINDOW_AUTOSIZE);
        cvShowImage("fourier_t",img);
        cvWaitKey(0);
        cvReleaseImage(&ipl);
        cvReleaseImage(&gausf);
        cvReleaseImage(&img);
        cvDestroyWindow("fourier_t");
        cvDestroyWindow("frequence_domain");
        cvDestroyWindow("gauss_fre");
        delete pixval;
        delete pixvals;
        delete src;
        delete twodfourier;
        delete twodifourier;
        return 1;
    }
    return 0;
}