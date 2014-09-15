#include "ExactSfftByOfdm_2d.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <bitset>
#define PI 3.1415926535627

using namespace std;


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
    fftw_free(this->src);
    fftw_free(this->sparsefouries);
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

fftw_complex * ExactSfftByOfdm_2d::FoldToBins(double * src,int Br,int Bc,int Tr,int Tc,bool samDirect)
{
    long samlength = 0;
    double * sampleEles ;
    if(samDirect){
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
            sampleEles[p] = src[(i*(this.length/Br)+Tr)*this->length+(j*(this->width/Bc)+Tc)];//�˴���������ڴ棬���ܻή�����ܣ������Ż���
            p++;
        }
    }
    fftw_complex * dftrans_sam = fftw_malloc(sizeof(fftw_complex)*samlength);
    fftw_plan p = fftw_plan_dft_r2c_1d(samlength,sampleEles,dftrans_sam,FFTW_ESTIMATE);//FFTW������Ҷ�任;
	fftw_execute(p);
    fftw_free(sampleEles);
    fftw_destroy_plan(p);
    return dftrans_sam;
}

void ExactSfftByOfdm_2d::BasicEstFreq(const vector<fftw_complex *> in_iter,const vector<fftw_complex*> in_update,int T,bool isCol)
{
//    TwoDFreElement * spareEle = new TwoDFreElement();
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
		//����  in_iter[1][i]/in_iter[1][0]
        dived[0] = in_iter[1][i][0]*in_iter[0][i][0]+in_iter[1][i][1]*in_iter[0][i][1]/(in_iter[0][i][0]*in_iter[0][i][0]+in_iter[0][i][1]*in_iter[0][i][1]);
        dived[1] = in_iter[1][i][1]*in_iter[0][i][0]-in_iter[1][i][0]*in_iter[0][i][1]/(in_iter[0][i][0]*in_iter[0][i][0]+in_iter[0][i][1]*in_iter[0][i][1]);
        double arctan = atan(devide[1]/devide[0]);
		if (arctan < 0)
		{
			arctan += PI;//����λ����ƽ��
		}
		long phase = arctan * updatelen / (2 * PI);  //ϡ����λ
        fftw_complex freValue = {in_iter[0][i][0],in_iter[0][i][1]}; //ϡ���ֵ
        fftw_complex tmp = {0};
		double Model = 0;
        for(int j = 0;j<T;j++){
			fftw_complex w_n_j = {cos(2*PI*phase*j/updatelen),sin(2*PI*phase*j/updatelen)};//������λ����;
			//���� freValue * w_n_j ,a=freValue[0] b=freValue[1],c=w_n_j[0],d=w_n_j[1]
			tmp[0] = in_iter[j][i][0]-(freValue[0]*w_n_j[0]-freValue[1]*w_n_j[1]);
			tmp[1] = in_iter[j][i][1]-(freValue[0]*w_n_j[1]+freValue[1]*w_n_j[0]);
			Model += sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
        }
		if (Model <=0.00000001) //Model��ֵ�ĵ����Ƿ�������õ���������;
		{
			if (isCol)
			{
				this->sparsefouries[i*length+phase] = freValue;
			}else{
				this->sparsefouries[phase*length+i] = freValue;
			}
			for (int j = 0;j<T;j++)
			{
				fftw_complex w_n_i = {cos(2*PI*i*j/iterlen,sin(2*PI*i*j/iterlen)};
				in_iter[j][i][0] = 0;
				in_iter[j][i][1] = 0;
				in_update[j][phase][0] -= (freValue[0]*w_n_i[0]-freValue[1]*w_n_i[1]);
				in_update[j][phase][1] -= (freValue[0]*w_n_i[1]+freValue[1]*w_n_i[0]);
			}
		}
    }
}


void ExactSfftByOfdm_2d::BasicExact2DSfft(int C_LogN)
{
/*	
	fftw_complex ** u_inCol = (fftw_complex *)fftw_malloc(sizeof(fftw_complex *)*sampleDims);//ָ������;
    fftw_complex ** v_inRow = (fftw_complex *)fftw_malloc(sizeof(fftw_complex *)*sampleDims);
*/
	vector<fftw_complex *> u_inCol(sampleDims);
	vector<fftw_complex *> v_inRow(sampleDims);
    for(int t = 0;t<this->sampleDims;t++){
        u_inCol[t] = FoldToBins(this->src,this->length,1,0,t,true);//FoldToBins����һ��Ƶ��ָ�룬ָ���л��в����õ���Ƶ������
        v_inRow[t] = FoldToBins(this->src,this->width,0,1,t,false);//FoldToBins����һ��Ƶ��ָ�룬ָ���л��в����õ���Ƶ������
    }


    this->sparsefouries = fftw_malloc(sizeof(fftw_complex)*width*length);//�洢����õ�ϡ��Ƶ��

    for(int j = 0;j<C_LogN;j++){//C_LogN���õ�������
//        TwoDFreElement * sparseFre = BasicEstFreq(u_inCol,v_inRow,true);
		BasicEstFreq(u_inCol,v_inRow,true);
//        long row = sparseFre->rowIndex;
//        long col = sparseFre->colIndex;
//        this.sparsefouries[row*this.length+col] = sparseFre->value;
//        sparseFre = BasicEstFreq(v_inRow,u_inCol,false);
		BasicEstFreq(v_inRow,u_inCol,false);
//        row = sparseFre->rowIndex;
//        col = sparseFre->colIndex;
 //       this.sparsefouries[row*this.length+col] = sparseFre->value;
    }
	for (int i = 0;i<this->sampleDims;i++)
	{
		fftw_free(u_inCol[i]);
		fftw_free(v_inRow[i]);
	}
}
