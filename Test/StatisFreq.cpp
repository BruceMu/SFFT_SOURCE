#include <fftw3.h>
#include <iostream>
#include <complex>
#include <algorithm> 
#include "cv.h"
#include "highgui.h"
#define DELTA 30
#define  PI 3.1415926535627
using std::complex;
using std::cout;
using std::endl;


bool compare(const complex<double> a,const complex<double> b){
	double model_a = sqrt(a.real()*a.real()+a.imag()*a.imag());
	double model_b = sqrt(b.real()*b.real()+b.imag()*b.imag());

	return model_a < model_b;
}
//自底向上归并排序频谱序列
void complex_sort(fftw_complex * queue,int dim){
	int length = 1;
}

void toGray(fftw_complex *twodfourier,int pwid,int phei,char* pixval)
{
	double *vals = new double[pwid*phei];
	double max = 0;
	double min = 255;
	for(size_t p = 0;p<pwid*phei;p++){
		vals[p]=sqrt(pow(twodfourier[p][0],2.0)+pow(twodfourier[p][1],2.0));//对数级的幅度铺
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

void guass(fftw_complex *src,int width,int height,double D)
{
	//find centre point
	int orgx = floor(width/2.0);
	int orgy = floor(height/2.0);
	double distence = 0;
	for(size_t i = 0;i<width;i++)
	{
		/*
		for(size_t j = 0;j<width;j++)
		{
			distence = (pow(abs((int)(i-orgy)),2.0)+pow(abs((int)(j-orgx)),2.0));
			//src[i*width+j] = src[i*width+j] * exp(-distence/(2*pow(D,2.0)));
			src[i*width+j][0] = src[i*width+j][0]*exp(-distence/(2*D*D));
			src[i*width+j][1] = src[i*width+j][1]*exp(-distence/(2*D*D));

		}*/
		distence = abs((int)(i-orgy))*abs((int)(i-orgy));
		for (size_t j = 0;j<height;j++)
		{
			src[i*width+j][0] *= exp(-distence/(2*D*D));
			src[i*width+j][1] *= exp(-distence/(2*D*D));
		}
		//cout<<i<<" "<<distence<<" "<<endl;

	}
	//    cout<<orgx<<" "<<orgy<<endl;
}
int main()
{
	//使用OpenCV读取图像,抽取像素数组 赋值
	IplImage *img;
	if((img = cvLoadImage("E:\\code\\SFFT_SOURCE\\doc\\lena.jpg",0))!=0){
		int dim = img->imageSize;
		fftw_complex * src = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*dim);

		size_t j =0;
		for(int y =0;y<img->height;y++)
		{
			uchar * ptr = (uchar*)(img->imageData + y * img->widthStep);
			for(int x =0;x<img->width;x++){
				/***输入函数*/
				src[j][0] = ((double)ptr[x]/255-0.5)*pow(-1,(double)(x+y));//归一化到[-1,1]区
				src[j][1] = 0;
				j++;
			}
		}


		//对图像二维矩阵序列做FFTW傅立叶变换
		fftw_complex * out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*dim);
		fftw_plan p =  fftw_plan_dft_2d(img->width,img->height,src,out,FFTW_FORWARD,FFTW_ESTIMATE);
		fftw_execute(p);

		//使用高斯函数过滤噪声；
		guass(out,img->width,img->height,DELTA);


		//统计每一行非零频谱的数量；
		int NonZore = 0;
		int smallNoise = 0;
		int matchline = 0;
		for(int i = 0;i<img->height;i++){
			NonZore = 0;
			smallNoise = 0;
			for(int j = 0;j<img->width;j++){
				fftw_complex tmp = {out[i*img->width+j][0],out[i*img->width+j][1]};
				if (sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1])>0)
				{
					NonZore ++;
					if (sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1])<20)
					{
						smallNoise++;
					}
				}/*else{
					if (sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1])>0)
					{
						smallNoise++;
						out[i*img->width+j][0]=0;//过滤小幅值噪声
						out[i*img->width+j][1]=0;
					}
				}*/
			}
			cout<<NonZore<<" NoneZeor freq and "<<smallNoise<<" small Noise inline"<<endl;
			/*
			if (NonZore <3)
			{
				cout<<"Match after filter "<<smallNoise<<" small noise"<<endl;
				matchline ++;
			}
			else{
				cout<<i<<"  "<<NonZore<<" after filter "<<smallNoise<<" small noise"<<endl;
			}*/
		}
		cout<<"matchline ;"<<matchline<<endl;
		cout<<"_________________________________________________"<<endl;
		cout<<"dim"<<img->width<<"*"<<img->height<<endl;

		//统计频谱数据:排序采样分位点
/*
		complex<double> * values = (complex<double>*)fftw_malloc(sizeof(complex<double>)*dim);

		for(int i = 0;i<dim;i++){
			values[i] = complex<double>(out[i][0],out[i][1]);
		}

		sort(values,values+dim,compare);

		int step = dim / 100;
		for(int i = 0;i<dim;i+=step){
			cout<<i<<".: "<<values[i].real()<<"  "<<values[i].imag()<<endl;
		}
*/
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

	}

}

