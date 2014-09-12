#include <fftw3.h>
#include <iostream>
#include <fstream>
void main()
{
    int N = 0;
    fftw_complex * in, * out;
    fftw_plan p;
    char * src = new char[N];

    ifstream fin("d:\\test_data.txt",ios::nocreate|ios::binary);
    char c;
    int i = 0;
    while((c = fin.get()!=EOF){
        src[i] = c;
        i++;
    }
//  in = (fftw_complex*)fftw_molloc(sizeof(fftw_complex) * N);
    in = (double *)fftw_molloc(sizeof(double)*N);
    out = (fftw_complex*)fftw_molloc(sizeof(fftw_complex) * N);
    for(int i = 0;i< N;i++){
        in[i] = (double)src[i];
    }
//  p = fftw_plan_dft_1d(N,in,out,FFTW_FORWORD,FFTW_ESTIMATE);
    p = fftw_plan_dft_r2c_1d(N,in,out,FFTW_FORWORD,FFTW_ESTIMATE);
    fftw_execute(p);


    fftw_destory_plan(p);
    fftw_free(in);
    fftw_free(out);
}
