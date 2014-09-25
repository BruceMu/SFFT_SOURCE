#include <iostream>
#include <fstream>
#include <fftw3.h>
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
    fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan p = fftw_plan_dft_1d(N,src,out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(p);


    for(int i = 0;i<20;i++){
        cout<<out[i][0]<<" "<<out[i][1]<<endl;
    }
    fftw_destroy_plan(p);
    fftw_free(src);
    fftw_free(out);
}
