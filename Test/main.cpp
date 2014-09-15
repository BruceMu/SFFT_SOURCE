#include <iostream>
#include <fstream>
#include <fftw3.h>
#include <math.h>
#define PI 3.1415926535627
using namespace std;


int main()
{
    int N = 20;

    fftw_complex * src = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);

    ifstream fin("E:\\testOneSparse.txt",ios::binary);

    for(int i = 0;i<20;i++){
        fin>>src[i][0];
        fin >> src[i][1];
        cout<<src[i][0]<<" "<<src[i][1]<<endl;
    }

    cout<<"--------------------------------------------------------"<<endl;
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan p = fftw_plan_dft_1d(N,src,out,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p);


    for(int i = 0;i<20;i++){
        cout<<out[i][0]<<" "<<out[i][1]<<endl;
    }

    fftw_complex devide = {0};
    //a=out[1][0] b=out[1][1] c=out[0][0] d=out[0][1]
    //a+bi=out[1] c+di=out[0]  devide = out[1]/out[0]
    devide[0] = (out[0][0]*out[1][0]+out[0][1]*out[1][1])/(out[0][0]*out[0][0]+out[0][1]*out[0][1]);
    devide[1] = (out[1][1]*out[0][0]-out[1][0]*out[0][1])/(out[0][0]*out[0][0]+out[0][1]*out[0][1]);
    double arctan = atan(devide[1]/devide[0]);
    if(arctan < 0)
        arctan +=PI;
    double phase = arctan*N/(2*PI);
    cout<<"-------------------------------------------------"<<endl;
    cout<<devide[0]<<endl;
    cout<<devide[1]<<endl;
    cout << phase <<endl;

	cout<<"------------yanzheng------------------------------"<<endl;
	fftw_complex freValue = {out[0][0],out[0][1]};
	fftw_complex tmp = {0,0};
	double model = 0;
	for(int i = 0;i< 4;i++){
		fftw_complex wn = {cos(2*PI*i*phase/N),sin(2*PI*i*phase/N)};
		cout<<freValue[0]*wn[0]-freValue[1]*wn[1]<<"  "<<freValue[0]*wn[1]+freValue[1]*wn[0]<<endl;
		tmp[0] = out[i][0] - (freValue[0]*wn[0]-freValue[1]*wn[1]);
		tmp[1] = out[i][1] - (freValue[0]*wn[1]+freValue[1]*wn[0]);
		model += sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
		cout<<tmp[0]<<" "<<tmp[1]<<endl;
	}
	cout<<"-------------------------------------------"<<endl;
	cout<<model<<endl;
	if (model < 0.00000001)
	{
		cout<<"sccess!"<<endl;
	}
    fftw_destroy_plan(p);
    fftw_free(src);
    fftw_free(out);
}
