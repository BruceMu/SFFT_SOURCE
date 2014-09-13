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
        cout<<src[i][0]<<src[i][1]<<endl;
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
    //a+bi c+di
    devide[0] = (out[0][0]*out[1][0]+out[0][1]*out[1][1])/(out[0][0]*out[0][0]+out[0][1]*out[0][1]);
    devide[1] = (out[1][1]*out[0][0]-out[1][0]*out[0][1])/(out[0][0]*out[0][0]+out[0][1]*out[0][1]);
    double phase = atan(devide[1]/devide[0])*N/(2*PI);
    cout<<"-------------------------------------------------"<<endl;
    cout<<devide[0]<<endl;
    cout<<devide[1]<<endl;
    cout << phase <<endl;
    fftw_destroy_plan(p);
    fftw_free(src);
    fftw_free(out);
}
