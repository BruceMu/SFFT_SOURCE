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



//double转255灰度
void toGray(fftw_complex *twodfourier,int pwid,int phei,char* pixval)
{
	double *vals = new double[pwid*phei];
	double max = 0;
	double min = 255;
	for(size_t p = 0;p<pwid*phei;p++){
		vals[p]=log(1+sqrt(pow(twodfourier[p][0],2.0)+pow(twodfourier[p][1],2.0)));//对数级的幅度铺
		//vals[p]=sqrt(pow(twodfourier[p].real(),2.0)+pow(twodfourier[p].imag(),2.0));
		if(vals[p] > max){
			max = vals[p];
		}
		if(vals[p] < min){
			min = vals[p];
		}
	}
	cout<<min<< " "<<max<<endl;
	//	cout<<pwid<<" "<<phei<<endl;
	for(size_t s = 0;s<pwid*phei;s++){
		pixval[s] = (char)((vals[s]-min)/(max-min)*255);
	}
	delete vals;
}
//二维高斯
void guass(fftw_complex *src,int width,int height,double D)
{
	//find centre point
	//	int orgx = floor(width/2.0);
	int orgy = floor(height/2.0);
	double distence = 0;
	for(size_t i = 0;i<height;i++)
	{
		for(size_t j = 0;j<width;j++)
		{
			distence = sqrt(pow(abs((int)(i-orgy)),2.0)+pow(abs((int)(j-orgx)),2.0));
			//src[i*width+j] = src[i*width+j] * exp(-distence/(2*pow(D,2.0)));
			src[i*width+j][0] = src[i*width+j][0]*exp(-distence/(2*pow(D,2.0)));
			src[i*width+j][1] = src[i*width+j][1]*exp(-distence/(2*pow(D,2.0)));

		}
		//cout<<i<<" "<<distence<<" "<<endl;

	}
	//    cout<<orgx<<" "<<orgy<<endl;
}


//二维变换，在一维变换后在一维方向上加高斯，看最后的效果；
int main()
{
	//使用OpenCV读取图像,抽取像素数组 赋值
	IplImage *img;
	if((img = cvLoadImage("E:\\code\\SFFT_SOURCE\\doc\\lena.jpg",0))!=0){
		int dim = img->imageSize;
		fftw_complex * src = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*img->width);

		size_t j =0;
		for(int y =0;y<1;y++)
		{
			uchar * ptr = (uchar*)(img->imageData + y * img->widthStep);
			for(int x =0;x<img->width;x++){
				/***输入函数*/
				src[j][0] = (double)ptr[x]/255-0.5);//归一化到[-1,1]区
				src[j][1] = 0;
				cout<<src[j][0]<<endl;
				j++;
			}
		}



		//计算单一序列傅立叶变换
		fftw_complex * out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*img->width);
		fftw_plan p =  fftw_plan_dft_1d(img->width,src,out,FFTW_FORWARD,FFTW_ESTIMATE);
		fftw_execute(p);

		for (int j = 0;j<img->width;j++)
		{
			cout<<out[j][0]<<" "<<out[j][1]<<endl;
		}



/*
		//使用高斯函数过滤噪声；
		guass(out,img->width,img->height,DELTA);


		//统计每一行非零频谱的数量；
		int NonZore = 0;
		int matchline = 0;
		for(int i = 0;i<img->height;i++){
			NonZore = 0;
			for(int j = 0;j<img->width;j++){
				fftw_complex tmp = {out[i*img->width+j][0],out[i*img->width+j][1]};
				if (sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1])>25)
				{
					NonZore ++;
				}else{
					out[i*img->width+j][0]=0;//过滤小幅值噪声
					out[i*img->width+j][1]=0;
				}
			}
			if (NonZore <3)
			{
				cout<<"##############Match###############"<<endl;
				matchline ++;
			}
			else{
				cout<<i<<"  "<<NonZore<<endl;
			}
		}
		cout<<"matchline ;"<<matchline<<endl;
		cout<<"_________________________________________________"<<endl;
		cout<<"dim"<<img->width<<"*"<<img->height<<endl;

		//显示频谱图
		char * pixval = (char*)fftw_malloc(sizeof(char)*dim);

		toGray(out,img->width,img->height,pixval);
		//		fftw_free(values);

		CvMat freq;
		cvInitMatHeader(&freq,img->width,img->height,CV_8UC1,pixval);
		IplImage *ipl = cvCreateImage(cvGetSize(&freq),8,1);
		cvGetImage(&freq,ipl);//获取频谱图像




		cout<<"-----------------------------"<<endl;
		fftw_complex * timeRecove = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*dim);
		fftw_plan p2 = fftw_plan_dft_2d(img->width,img->height,out,timeRecove,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_execute(p2);

		double ipix = 0;
		size_t po = 0;
		for(int y =0;y<img->height;y++)
		{
			uchar * ptr = (uchar*)(img->imageData + y * img->widthStep);
			for(int x =0;x<img->width;x++){
				ipix = timeRecove[po*img->width+x][0];
				ptr[x] = (uchar)((ipix/dim+0.5)*pow(-1,(double)(x+y)) * 255);
			}
			po++;
		}
		cvNamedWindow("frequence_domain",CV_WINDOW_AUTOSIZE);
		cvShowImage("frequence_domain",ipl);
		cvNamedWindow("fourier_t",CV_WINDOW_AUTOSIZE);
		cvShowImage("fourier_t",img);
		cvWaitKey(0);
		cvReleaseImage(&img);
		cvDestroyWindow("fourier_t");
		cvReleaseImage(&ipl);
		cvDestroyWindow("frequence_domain");
		fftw_free(pixval);
		fftw_free(timeRecove);
		fftw_free(src);
		fftw_destroy_plan(p);
		*/
	}

}

