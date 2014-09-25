/************************************************************
*二维离散傅立叶变换和滤波函数（高斯低通）的实现
*用于计算图像离散序列的傅立叶频谱图，并用于频域图像处理
*算法：多项式快速傅立叶算法（蝶形） 理论基础：算法导论,图像处理
*时间复杂度：一维O(NlogN)，二维O(N^2)
*工具：openCV读取图像与展示效果
*版本：测试
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

/*******为自底向上的迭代重排序列计算位置************/
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
//重载形式
complex<double>* FFT(const complex<double> *srcimg,int n)
{
    double flag = log10((double)n)/log10(2.0);
    int N = n;
    if(flag - (int)flag != 0){   //判断是否为2的指数次
        cout <<"the length of srcimg is wrong"<< endl;
        /*填充成2的指数项*/
        cout <<"need padding"<<endl;
        N = pow(2,(double)((int)flag+1));
        flag = log10((double)N)/log10(2.0);
    }
    /*改变成奇偶顺序*/
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
        complex<double> wm = complex<double>(cos(2*pi/m),sin(2*pi/m));//wm1：从W21开始,周期变换
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
/***********一维快速傅立叶变换********************
*srcimg : 原始一维序列                          *
*n     ：一维序列的长度
*************************************************/
complex<double>* FFT(const double *srcimg,int n)
{
    double flag = log10((double)n)/log10(2.0);
    int N = n;
    if(flag - (int)flag != 0){   //判断是否为2的指数次
//        cout <<"the length of srcimg is wrong"<< endl;
        /*填充成2的指数项*/
//        cout <<"need padding"<<endl;
        N = pow(2,(double)((int)flag+1));
        flag = log10((double)N)/log10(2.0);
    }
    /*改变成奇偶顺序*/
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
    /*基于迭代的蝶形快速傅立叶变换，自底向上*/
     for(int s =1;s <= flag;s++)
    {
        int m = pow(2,(double)s);
        complex<double> wm = complex<double>(cos(2*pi/m),sin(2*pi/m));//wm1：从W21开始,周期变换
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
/*****************一维快速傅立叶逆变换********************/
/*fftimg:原始一维傅立叶序列
  n     : 序列长度
*******************************************************/
complex<double>* IFFT(const complex<double> *fftimg,int n)
{
    n = countPadding(n);
    double flag = log10((double)n)/log10(2.0);
    int N = n;
    if(flag - (int)flag != 0){   //判断是否为2的指数次
        cout <<"the length of srcimg is wrong"<< endl;
        /*填充成2的指数项*/
        cout <<"need padding"<<endl;
        N = pow(2,(double)((int)flag+1));
        flag = log10((double)N)/log10(2.0);
    }
    /*改变成奇偶顺序*/
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
        complex<double> wm = complex<double>(cos(-2*pi/m),sin(-2*pi/m));//wm1：从W2(-1)开始
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
  
/*******使用共轭的快速傅立叶逆变换**************/
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

/*************二维快速傅立叶变换**************************
*srcimg: 用一维表示的二维原始序列
*width ：宽度
*height：高度
********************************************************/
complex<double>* twoDFFT(const double *srcimg,int width,int height)
{
    int w = countPadding(width);
    int h = countPadding(height);
    int pixes = w * h;
    complex<double> *hdirection = new complex<double>[w];
    complex<double> *vdirection = new complex<double>[h];
    complex<double> *fourier = new complex<double>[pixes];
    /*二维水平方向*/
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
    /*二维垂直方向*/
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
/**************二维快速傅立叶逆变换*************************
*fourier : 一维表示的二维傅立叶变换序列                     *
*width   :宽                                             *
*height  :高                                             *
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
        complex<double> *hifour = IFFT(hdirection,width);//临时变量
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
/******************计算填充数********************************************
*蝶形傅立叶变换算法只计算2的指数次的离散序列，对于非2的指数次的序列，使用零填充*
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

/*****高斯低通滤波函数*************************
*src：   输入频谱
*width： 输入频谱宽度
*height：输入频谱高度
*D：     高斯函数方差，即滤波阈值
*只对于移位居中后的傅立叶频谱进行高斯低通滤波
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
/************复数傅立叶频谱数组转换成256级灰度数组*****************/
void toGray(complex<double> *twodfourier,int pwid,int phei,char* pixval)
{
        double *vals = new double[pwid*phei];
        double max = 0;
        double min = 255;
        for(size_t p = 0;p<pwid*phei;p++){
            vals[p]=log(1+sqrt(pow(twodfourier[p].real(),2.0)+pow(twodfourier[p].imag(),2.0)));//对数级的幅度铺
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
/******************opencv 读取图像与展示效果***********************/
int StatisFreqFreqWithUdefine(int argc,char **argv)
{
    IplImage *img;
    if((img = cvLoadImage("E:\\code\\SFFT_SOURCE\\doc\\lena.jpg",0))!=0){
        int dim = img->imageSize;
        //从图像矩阵复制出单维数组，并做频谱居中操作，对应改写接口，返回单维数组；
        double * src = new double[dim];
        size_t j =0;
        for(int y =0;y<img->height;y++)
        {
            uchar * ptr = (uchar*)(img->imageData + y * img->widthStep);
            for(int x =0;x<img->width;x++){
                /***输入函数乘以（-1）的（x+y）次方，频谱将居中*/
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
        toGray(twodfourier,pwid,phei,pixval);//复数频谱转256灰度
        cvInitMatHeader(&freq,pwid,phei,CV_8UC1,pixval);
        IplImage *ipl = cvCreateImage(cvGetSize(&freq),8,1);
        cvGetImage(&freq,ipl);//获取频谱图像

        guass(twodfourier,pwid,phei,DELTA);//方差DELTA的高斯平滑（高斯低通滤波）
        CvMat gaufre;
        char *pixvals = new char[pwid*phei];
        toGray(twodfourier,pwid,phei,pixvals);
        cvInitMatHeader(&gaufre,pwid,phei,CV_8UC1,pixvals);
        IplImage *gausf = cvCreateImage(cvGetSize(&gaufre),8,1);
        cvGetImage(&gaufre,gausf);//高斯平滑后的频谱图像

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