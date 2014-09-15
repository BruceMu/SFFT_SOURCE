#include <iostream>
#include <fftw3.h>
using namespace std;

int main()
{
	int N = 64;
	fftw_complex * src = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N);
	fftw_complex * testSrc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N);
	fftw_plain p =  fftw_plan_dft_2d(8,8,src,testSrc,FFTW_BACKWARD,FFTW_ESTIMATE);


    ifstream fin("E:\\testDataSmall.txt",ios::binary);
	for(int i = 0;i<N;i++){
		fin>>src[i][0];
		src[i][1] = 0;
		cout<<src[i][0]<<endl;
	}


	fftw_execute(p);

	ExactSfftByOfdm_2d esfft = new ExactSfftByOfdm_2d(testSrc,8,8);
	esfft.setSampleDims(3);
	esfft.setSparse(15);
	esfft.BasicExact2DSfft(16);
    return 0;
}
