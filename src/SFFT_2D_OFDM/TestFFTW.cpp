#include <fftw3.h>
#include <iostream>
void main()
{
    fftw_complex * in, * out;
    fftw_plan p;

    in = (fftw_complex*)fftw_molloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*)fftw_molloc(sizeof(fftw_complex) * N);

    p = fftw_plan_dft_1d(N,in,out,FFTW_FORWORD,FFTW_ESTIMATE);
    fftw_execute(p);


    fftw_destory_plan(p);
    fftw_free(in);
    fftw_free(out);
}
