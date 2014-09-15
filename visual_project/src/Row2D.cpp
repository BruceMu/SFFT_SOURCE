#include "Raw2D.h"
#include<math.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
//#include "Filter.h"
using namespace std;
#define M_EXP 2.7182818284590452353602874713527




static float lgtt=log10(2.0f);

Raw2D::Raw2D()
{
	xsize=0;
	ysize=0;
	data=NULL;
}

Raw2D::Raw2D(int xsize,int ysize,PIXTYPE *y)
{
	this->xsize=xsize;
	this->ysize=ysize;
	this->data=y;
}

Raw2D::Raw2D( int xsize,int ysize )
{
	this->xsize=xsize;
	this->ysize=ysize;
	this->data=new PIXTYPE[xsize*ysize];
}

Raw2D::Raw2D( Raw2D *r)
{
	this->xsize=r->xsize;
	this->ysize=r->ysize;
	this->data=new PIXTYPE[xsize*ysize];
	//for (int i=0;i<xsize;i++)
	//{
	// for (int j=0;j<ysize;j++)
	// {
	//	 data[i*ysize+j]=r->data[i*(r->ysize)+j];
	// }
	//}
	memcpy(this->data,r->data,sizeof(PIXTYPE)*xsize*ysize);
}

Raw2D::Raw2D(const Raw2D& r)
{
	//if (this->xsize != r.getXsize() && this->ysize != r.getYsize())
	//{
	this->xsize=r.xsize;
	this->ysize=r.ysize;
	//	//if (this->data != NULL)
	//	//	delete[] this->data;
	this->data=new PIXTYPE[xsize*ysize];
	//}
	memcpy(this->data, r.data, sizeof(PIXTYPE)*xsize*ysize);
}

Raw2D::~Raw2D(void)
{
	if(this->data!=NULL)
		delete [] this->data;
	data=NULL;
}
