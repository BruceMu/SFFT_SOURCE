#include "ExactSfftByOfdm_2d.h"
#include <iostream>
#include <math.h>
#define PI 3.1415926535627


ExactSfftByOfdm_2d::ExactSfftByOfdm_2d(Raw2D & raw2)
{
    this.length = raw2.getXsize();
    this.width = raw.getYsize();
    src = (double *)fftw_malloc(sizeof(double)*length*width);
    for(long i = 0;i < length*width;i++){
        src[i] = raw2.getXY(i)/255;
    }
}

ExactSfftByOfdm_2d::~ExactSfftByOfdm_2d()
{
    fftw_free(this.src);
    fftw_free(this.sparsefouries)
}

/** \brief ����ϡ���
 * \param ϡ���
 */
ExactSfftByOfdm_2d::setSparse(int ksparse)
{
    this.sparseK = ksparse;
}
/** \brief ���ò������ޣ�
 * \param ������������
 */

ExactSfftByOfdm_2d::setSampleDims(int sampleT)
{
    this.sampleDims = sampleT;
}

fftw_complex * ExactSfftByOfdm_2d::FoldToBins(double * src,int Br,int Bc,int Tr,int Tc,bool samDir)
{
    long samlength = 0;
    double * sampleEles ;
    if(samDir){
        sampleEles = (double *)fftw_malloc(sizeof(double)*width);//���еĳ���Ϊ�����������������Ϊ�еĳ��ȡ�
        samlength = width;
    }else{
        sampleEles = (double *)fftw_malloc(sizeof(double)*length);//
        samlength = length;
    }
    //��ԭ�ź����в�����ֵ��
    long p = 0;
    for(int i =0;i<Br;i++){
        for(int j = 0;j<Bc;j++){
            sampleEles[p] = src[(i*(this.length/Br)+Tr)*this.length+(j*(this.width/Bc)+Tc)];
            p++;
        }
    }
    fftw_complex * dftrans_sam = fftw_malloc(sizeof(fftw_complex)*samlength);
    fftw_plan p = fftw_plan_dft_r2c_1d(samlength,sampleEles,dftrans_sam,FFTW_ESTIMATE);//FFTW������Ҷ�任;
    fftw_free(sampleEles);
    fftw_destroy_plan(p);
    return dftrans_sam;
}

TwoDFreElement* ExactSfftByOfdm_2d::BasicEstFreq(fftw_complex ** in_iter,fftw_complex ** in_update,int T,bool isCol)
{
    TwoDFreElement * spareEle = new TwoDFreElement();
    long iterlen  = length;
    long updatelen = width;
    if(isCol){
        iterlen = width;
        updatelen = length;
    }
    for(int i = 0;i<iterlen;i++){
        fftw_complex sample_sum = {0};
        for(int j = 0;j<T;j++){
            sample_sum[0] += in_iter[j][i][0];
            sample_sum[1] += in_iter[j][i][1];
        }
        if(sample_sum[0] ==0 && sample_sum[1] ==0){
          continue;
        }//�жϲ���Ƶ���Ƿ����0
        fftw_complex dived = {0};
        dived[0] = in_iter[1][i][0]*in_iter[0][i][0]+in_iter[1][i][1]*in_iter[0][i][1]/(in_iter[0][i][0]*in_iter[0][i][0]+in_iter[0][i][1]*in_iter[0][i][1]);
        dived[1] = in_iter[1][i][1]*in_iter[0][i][0]-in_iter[1][i][0]*in_iter[0][i][1]/(in_iter[0][i][0]*in_iter[0][i][0]+in_iter[0][i][1]*in_iter[0][i][1]);
        long phase = atan(dived[0]/dived[1])*updatelen/2*PI;
        fftw_complex freValue = {in_iter[0][i][0],in_iter[0][i][1]};
        fftw_complex tmp = {0};
        for(int j = 0;j<T;j++){
            tmp[0] = in_iter[j][i][0]-freValue[0];
            tmp[1] = in_iter[j][i][1]-freValue[1];
        }
    }
}

void ExactSfftByOfdm_2d::BasicExact2DSfft(double * src,int C_LogN)
{
    fftw_complex ** u_inCol;
    fftw_complex ** v_inRow;
    for(int t = 0;t<this.sampleDims;t++){
        u_inCol[t] = FoldToBins(src,this.length,1,0,t,true);//FoldToBinsӦ�з���ֵ�����ز����л��еĸ���Ҷ�任Ƶ��
        v_inRow[t] = FoldToBins(src,this.width,0,1,t,false);
    }
    this.sparsefouries = fftw_malloc(sizeof(fftw_complex)*width*length);

    for(int j = 0;j<C_LogN;j++){//C_LogN���õ�������
        TwoDFreElement * sparseFre = BasicEstFreq(u_inCol,v_inRow,true);
        long row = sparseFre->rowIndex;
        long col = sparseFre->colIndex;
        this.sparsefouries[row*this.length+col] = sparseFre->value;
        sparseFre = BasicEstFreq(v_inRow,u_inCol,false);
        row = sparseFre->rowIndex;
        col = sparseFre->colIndex;
        this.sparsefouries[row*this.length+col] = sparseFre->value;
    }
}
