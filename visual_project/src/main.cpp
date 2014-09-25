#include <iostream>
#include <fftw3.h>
using namespace std;

int main()
{
	int N = 64;
	fftw_complex * src = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N);
	fftw_complex * testSrc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N);
	fftw_plan p =  fftw_plan_dft_2d(8,8,src,testSrc,FFTW_BACKWARD,FFTW_ESTIMATE);

	//���쳤64ϡ��15��Ƶ��;
    ifstream fin("E:\\testDataSmall.txt",ios::binary);
	for(int i = 0;i<N;i++){
		fin>>src[i][0];
		src[i][1] = 0;
		cout<<src[i][0]<<endl;
	}


	fftw_execute(p);//��任��ʱ��,testSrcΪƵ��ϡ���ʱ���źţ�

	ExactSfftByOfdm_2d esfft(testSrc,8,8);
	esfft.setSampleDims(3);
	esfft.setSparse(15);
	esfft.BasicExact2DSfft(16);

	fftw_free(src);
	fftw_free(testSrc);
	fftw_destory_plan(p);
    return 0;
}
