#include <iostream>
#include <fstream>
#include <fftw3.h>
using namespace std;


int main()
{
    int N = 20;
    fftw_complex * out;
    double * in;
    fftw_plan p;
    double * src = new double[N];

    ifstream fin("E:\\testDataSmall.txt",ios::binary);

    for(int i = 0;i<20;i++){
        fin>>src[i];
        cout<<src[i]<<endl;
    }
    cout<<"-----------------------"<<endl;
//  in = (fftw_complex*)fftw_molloc(sizeof(fftw_complex) * N);
    in = (double *)fftw_malloc(sizeof(double)*N);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    double * out2 = (double *)fftw_malloc(sizeof(double)*N);
    fftw_plan p2 = fftw_plan_dft_c2r_1d(N,out,out2,FFTW_MEASURE);
    p = fftw_plan_dft_r2c_1d(N,in,out,FFTW_MEASURE);
    for(int i = 0;i< N;i++){
        in[i] = src[i];
    }
//  p = fftw_plan_dft_1d(N,in,out,FFTW_FORWORD,FFTW_MEASURE);
//    p = fftw_plan_dft_r2c_1d(N,in,out,FFTW_MEASURE);
    for(int i = 0;i<20;i++){
        cout<<in[i]<<endl;
    }
    cout<<"---------------------------------"<<endl;
    fftw_execute(p);
    fftw_complex inital = out[0];
    for(int i = 0;i<20;i++){
        if(out[i][0]> 0 or out[i][1]> 0 ){
            cout<<out[i][1]<<", "<<out[i][2]<<"i"<<out[i]-inital<<endl;

        }
    }

    fftw_execute(p2);
    cout << "------------------------"<<endl;
    for(int i = 0;i<20;i++){
        cout<<out2[i]<<endl;
    }
    delete [] src;
    fftw_destroy_plan(p);
    fftw_destroy_plan(p2);
    fftw_free(in);
    fftw_free(out);
    fftw_free(out2);
}
