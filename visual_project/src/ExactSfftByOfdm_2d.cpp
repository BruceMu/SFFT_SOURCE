#include <iostream>
#include <math.h>
#include <vector>
#include <bitset>
#include <complex>
#include "ExactSfftByOfdm_2d.h"
#include <Eigen/Dense> 
#define PI 3.1415926535627
#define D 0.8;
#define NoiseThrold 25;

using namespace std;
using namespace Eigen;


void ExactSfftByOfdm_2d::ReadImage(unsigned char * buf,char const * file,int size,int start)
{
	FILE op = fopen(file,"rb");
	if (op == null)
	{
		printf("open file");
	}
	fseek(op,start,SEEK_SET);
	fread(buf,sizeof(unsigned char),size,op);
}
ExactSfftByOfdm_2d::~ExactSfftByOfdm_2d()
{
    fftw_free(this->src);
    fftw_free(this->sparsefouries);
}

/** \brief 设置稀疏度
 * \param 稀疏度
 */
void ExactSfftByOfdm_2d::setSparse(int ksparse)
{
    this->sparseK = ksparse;
}
/** \brief 设置采样上限；
 * \param 采样行列数。
 */

void ExactSfftByOfdm_2d::setSampleDims(int sampleT)
{
    this->sampleDims = sampleT;
}
/*
fftw_complex * ExactSfftByOfdm_2d::FoldToBins(const double * src,int Br,int Bc,int Tr,int Tc,bool samDirect)
{
    long samlength = 0;
    double * sampleEles ;
    if(samDirect){
        sampleEles = (double *)fftw_malloc(sizeof(double)*width);//以行的长度为间隔采样，样本长度为列的长度。
        samlength = width;
    }else{
        sampleEles = (double *)fftw_malloc(sizeof(double)*length);//
        samlength = length;
    }
    //对原信号序列采样赋值；
    long p = 0;
    for(int i =0;i<Br;i++){
        for(int j = 0;j<Bc;j++){
            sampleEles[p] = src[(i*(this->length/Br)+Tr)*this->length+(j*(this->width/Bc)+Tc)];//此处采样会跨内存，可能会降低性能，建议优化。
            p++;
        }
    }
    fftw_complex * dftrans_sam = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*samlength);
    fftw_plan p = fftw_plan_dft_r2c_1d(samlength,sampleEles,dftrans_sam,FFTW_ESTIMATE);//FFTW做傅立叶变换;
	fftw_execute(p);
    fftw_free(sampleEles);
    fftw_destroy_plan(p);
    return dftrans_sam;
}
*/
/*
按行采样二维序列，并对采样序列做傅立叶变换;
*/
void ExactSfftByOfdm_2d::SampleInRow(const char * src,int sampleNums){
	//sampleEles 用于保存采样序列
	fftw_complex * sampleEles = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*this->width);
	//sampleForiers用于保存采样序列的傅立叶变换序列;
	fftw_complex * sampleFouriers ;
	fftw_plan p =  fftw_plan_dft_1d(this->width,sampleEles,sampleFouriers,FFTW_FORWARD,FFTW_ESTIMATE);

	for (int i = 0;i<sampleNums;i++)
	{
		sampleFouriers = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*this->width);
		for (j = 0;j<this->width;j++
		{
			sampleEles[j][0] = ((double)src[i*this->width+j]/255-0.5)*pow(-1,(double)(i+j));
			sampleEles[j][1] = 0;
		}
		//FFTW快速傅立叶变换；
		fftw_execute(p);
		this->rowSampleFouries[i] = sampleFouriers;
	}

	//清理；
	fftw_free(sampleEles);
	fftw_destory_plan(p);
	//返回采样序列的傅立叶变换后序列的内存指针；

}

//按列采样二维序列，并对采用序列做快速傅立叶变换;
void ExactSfftByOfdm_2d::SampleInCol(const char * src,int sampleNums){
	fftw_complex * sampleEles = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*this->height);
	fftw_complex * sampleFouriers;
	fftw_plan p = fftw_plan_dft_1d(this->height,sampleEles,sampleFouriers,FFTW_FORWARD,FFTW_ESTIMATE);
	for (int i = 0;i<sampleNums;i++)
	{
		for (int j = 0;j<this->height;j++)
		{
			//此处跨内存访问，无法利用缓存局部性。建议优化;
			sampleEles[j][0] = ((double)src[j*this->width+i]/255-0.5)*pow(-1,(double)(i+j));
			sampleEles[j][1] = 0;
		}
		fftw_execute(p);

		this->colSampleFouriers[i] = sampleFouriers;
	}

	fftw_free(sampleEles);
	fftw_destory_plan(p);
}

//directFlag为True，高斯模糊采样列；去除部分噪声;
//changeLog:并入第一次迭代，以减少数组遍历次数;
/*
void ExactSfftByOfdm_2d::MakeGaussFilt(bool directFlag,int sampleNums;){
	if (directFlag)
	{
		for (int i = 0;i<this->height;i++)
		{
			for(int j = 0;j<sampleNums;j++)
			{
				this->colSampleFouriers[j][i][0] *= exp(-(((double)i-(double)this->height/2)*((double)i-(double)this->height/2)/(2*D*D))/sqrt(2*PI*D*D);
			}
		}
	}else{
		for (int i = 0;i<this->width;i++)
		{
			for(int j = 0;j<sampleNums;j++)
			{
				this->rowSampleFouriers[j][i][0] *= exp(-((double)i-(double)this->width/2)*((double)i-(double)this->width/2)/(2*D*D))/sqrt(2*PI*D*D);
			}
		}
	}
}
*/
//使用多线程方式计算，在第二次迭代前完成。逆矩阵的计算可人工完成;
void ExactSfftByOfdm_2d::Inital(const int width,const int height)
{
	//初始化行和列的稀疏标志位;
	static bitset<height> rowFlag;
	static bitset<width> colFlag;
	rowFlag.set();
	colFlag.set();
	//初始化傅立叶基矩阵;
	//初始化矩阵，采用逐个赋值的方式，这样的矩阵具有通用性，是否可以使用多线程创建;
	complex<double> wm(cos(-2*PI/this->width),sin(-2*PI/this->width));
	complex<double> wm1(conj(wm));
	for (int i = 0;i<this->sparseK;i++)
	{
		for (int j = 0;j<this->width;j++)
		{
			this->fourRowBase(i,j) = pow(wm1,i*j);
		}
	}

	complex<double> wn(cos(-2*PI/this->height,sin(-2*PI/this->height)));
	complex<double> wn1(conj(wn));

	for(int i =0;i<this->sparseK;i++)
	{
		for (int j = 0;j<this->height;j++)
		{
			this->fourColBase(i,j) = pow(wn1,i*j);
		}
	}

}
//计算高斯模糊后的恢复频谱是否为零；
double ExactSfftByOfdm_2d::SumMagnitude(int recoveIndex,bool direct,int testlen)
{
	//direct = True为使用列采样样本进行行恢复；
	double sum = 0;
	if (direct)
	{
		for(int i = 0;i<testlen;i++)
		{
			sum += sqrt(colSampleFouriers[i][recoveIndex][0]*colSampleFouriers[i][recoveIndex][0]
			+colSampleFouriers[i][recoveIndex][1]*colSampleFouriers[i][recoveIndex][1]);
		}
	}else{
		for (int i = 0;i<this->testlen;i++)
		{
			sum += sqrt(RowSampleFouriers[i][recoveIndex][0]*RowSampleFouriers[i][recoveIndex][0]
			+RowSampleFouriers[i][recoveIndex][1]*RowSampleFouriers[i][recoveIndex][1]);
		}
	}
	return sum;
}

void ExactSfftByOfdm_2d::ReocverByRow()
{
	for(int i = 0;i<this->height;i++)
	{
		if (rowFlag.test(i))
		{
			if (SumMagnitude(i,true,3) <= NoiseThrold)//根据一定的噪声阈值来控制模糊后的频谱稀疏;
			{
				rowFlag.reset(i);
				continue;
			}else{

			}
		}
	}
}

//One-sparse迭代
void ExactSfftByOfdm_2d::BasicEstFreq(const vector<fftw_complex *> in_iter,const vector<fftw_complex*> in_update,bool isCol)
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
        }//判断采样频谱是否大于0
        fftw_complex dived = {0};
		//计算  in_iter[1][i]/in_iter[1][0]
        dived[0] = in_iter[1][i][0]*in_iter[0][i][0]+in_iter[1][i][1]*in_iter[0][i][1]/(in_iter[0][i][0]*in_iter[0][i][0]+in_iter[0][i][1]*in_iter[0][i][1]);
        dived[1] = in_iter[1][i][1]*in_iter[0][i][0]-in_iter[1][i][0]*in_iter[0][i][1]/(in_iter[0][i][0]*in_iter[0][i][0]+in_iter[0][i][1]*in_iter[0][i][1]);
        double arctan = atan(devide[1]/devide[0]);
		if (arctan < 0)
		{
			arctan += PI;//负相位周期平移
		}
		long phase = arctan * updatelen / (2 * PI);  //稀疏相位
        fftw_complex freValue = {in_iter[0][i][0],in_iter[0][i][1]}; //稀疏幅值
        fftw_complex tmp = {0};
		double Model = 0;
        for(int j = 0;j<this->sampleDims;j++){
			fftw_complex w_n_j = {cos(2*PI*phase*j/updatelen),sin(2*PI*phase*j/updatelen)};//周期相位因子;
			//计算 freValue * w_n_j ,a=freValue[0] b=freValue[1],c=w_n_j[0],d=w_n_j[1]
			tmp[0] = in_iter[j][i][0]-(freValue[0]*w_n_j[0]-freValue[1]*w_n_j[1]);
			tmp[1] = in_iter[j][i][1]-(freValue[0]*w_n_j[1]+freValue[1]*w_n_j[0]);
			Model += sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
        }
		if (Model <=0.00000001) //Model阈值的调节是否可以适用到噪声情形;
		{
			if (isCol)
			{
				this->sparsefouries[i*length+phase] = freValue;
			}else{
				this->sparsefouries[phase*length+i] = freValue;
			}
			for (int j = 0;j<this->sampleDims;j++)
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
	fftw_complex ** u_inCol = (fftw_complex *)fftw_malloc(sizeof(fftw_complex *)*sampleDims);//指针数组;
    fftw_complex ** v_inRow = (fftw_complex *)fftw_malloc(sizeof(fftw_complex *)*sampleDims);
*/
	vector<fftw_complex *> u_inCol(sampleDims);
	vector<fftw_complex *> v_inRow(sampleDims);
    for(int t = 0;t<this->sampleDims;t++){
        u_inCol[t] = FoldToBins(this->src,this->length,1,0,t,true);//FoldToBins返回一个频谱指针，指向按行或列采样得到的频谱序列
        v_inRow[t] = FoldToBins(this->src,this->width,0,1,t,false);//FoldToBins返回一个频谱指针，指向按行或列采样得到的频谱序列
    }


    this->sparsefouries = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*width*length);//存储计算得到稀疏频谱

    for(int j = 0;j<C_LogN;j++){//C_LogN设置迭代次数
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
