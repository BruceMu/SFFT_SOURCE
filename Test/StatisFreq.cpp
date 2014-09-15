#include <fftw3.h>
#include <iostream>
#include "cv.h"
#include "highgui.h"

int main()
{
	//ʹ��OpenCV��ȡͼ��,��ȡ��������
	IplImage *img;
	if((img = cvLoadImage("F:\\Fig0222(a)(face).tif",0))!=0){
		int dim = img->imageSize;
		double * src = (double *)fftw_malloc(sizeof(double)*dim);
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
		//��ͼ���ά������FFTW����Ҷ�任
		fftw_complex * out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*dim);
		fftw_plan p =  fftw_plan_dft_r2c_2d(img->width,img->height,src,out);
		fftw_execute(p);

		//ͳ��Ƶ������
		for(int i = 0;i<dim;i++){

		}

	}

}