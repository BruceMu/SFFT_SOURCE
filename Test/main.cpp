#include <fftw3.h>
#include <iostream>
#include <complex>
#include <algorithm> 
#include "cv.h"
#include "highgui.h"
#define DELTA 8.0
using std::complex;
using std::cout;
using std::endl;
/*
int main()
{
	//使用OpenCV读取图像,抽取像素数组 赋值
	IplImage *img;
	if((img = cvLoadImage("E:\\code\\SFFT_SOURCE\\doc\\lena.jpg",0))!=0){
		int dim = img->imageSize;
		fftw_complex * src = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*img->width);

		size_t j =0;
		cout<<"10 origin pixs"<<endl;
		for(int y =100;y<101;y++)
		{
			uchar * ptr = (uchar*)(img->imageData + y * img->widthStep);
			for(int x =0;x<img->width;x++){
				src[j][0] = ((double)ptr[x]/255)*pow(-1,(double)(x+y));//归一化到[-1,1]区
				src[j][1] = 0;
				if(j < 10){
					cout<<src[j][0]<<endl;
				}
				j++;
			}
		}




		//计算单一序列傅立叶变换;
		fftw_complex * out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*img->width);
		fftw_plan p =  fftw_plan_dft_1d(img->width,src,out,FFTW_FORWARD,FFTW_ESTIMATE);
		fftw_execute(p);

		//抑制小幅值噪声频谱;
		int sumN = 0;
		double max=0,min = 0;
		double fremag = 0;
		for (int j = 0;j<img->width;j++)
		{
			fremag = sqrt(out[j][0]*out[j][0]+out[j][1]*out[j][1]);
			if (max < fremag)
			{
				max = fremag;
			}else if (min > fremag)
			{
				min = fremag;
			}
			
			if (sqrt(out[j][0]*out[j][0]+out[j][1]*out[j][1])<1)
			{
				out[j][0] = 0;
				out[j][1] = 0;
				sumN ++;
			}
			
		}

		cout<<"small magnitudes: "<<sumN<<endl;
		cout<<"max : "<<max<<endl;
		cout<<"min : "<<min<<endl;
		
		

		//测试抑制小噪声后对原信号的影响;
		fftw_complex * outNext = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*img->width);
		fftw_plan p2 = fftw_plan_dft_1d(img->width,out,outNext,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_execute(p2);


		cout<<"after recover pix"<<endl;
		for (int i = 0;i<10;i++)
		{
			cout<<(outNext[i][0])/img->width<<" "<<(outNext[i][1])/img->width<<" "<<
				sqrt(((outNext[i][0])/img->width)*((outNext[i][0])/img->width)+((outNext[i][1])/img->width)*((outNext[i][1])/img->width))<<endl;
		}

	}
	
}
*/
