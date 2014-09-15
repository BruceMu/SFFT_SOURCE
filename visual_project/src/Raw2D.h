#ifndef RAW2D_H_INCLUDED
#define RAW2D_H_INCLUDED
#include <iostream>
#include <string.h>
#define PIXTYPE float
class Raw2D
{
private:   			//-----------------DATA-----------------
	int xsize;		// # of pixels per scanline,
	int ysize;		// # of scanlines in this Raw2D.
	PIXTYPE *data;		// 1D array of PIXTYPE that are accessed as a 2D array.

	void putData(const PIXTYPE *res)
	{
		if (data != NULL)
			delete [] data;
		memcpy(data, res, xsize*ysize);
	}

public:
	const static int MAXPIXEL = 255;
	const static int  MINPIXEL = 0;

	Raw2D(int, int, PIXTYPE*);
	Raw2D(int, int);
	Raw2D(Raw2D* r);
	Raw2D(const Raw2D& r);
	Raw2D(void);		// constructor for 'empty' Raw2Ds
	~Raw2D(void);		// destructor; releases memory

	int size() const { return xsize*ysize; }


	void sizer(int ixsize, int iysize);	// get mem for rectangle of pixels
	void sizer(Raw2D* src);					// get same amt. of mem as 'src'
	int getXsize(void) const {return xsize;}		// get # pixels per scanline
	int getYsize(void) const {return ysize;}		// get # of scanlines.

	inline void put(int ix, int iy, PIXTYPE val)	// write 'val' at location ix,iy.
	{
#ifdef _DEBUG  //only check under debug mode
		if (iy + ysize*ix < size())
		{
#endif
			data[iy + ysize*ix] = val;
#ifdef _DEBUG
		}
		else
			cout<<"out of size put"<<endl;
#endif
	}

	inline PIXTYPE get(int ix, int iy) {	// read the value at ix,iy.
#ifdef _DEBUG
		if(iy + ysize*ix<=size())
		{
#endif
			return data[iy + ysize*ix];
#ifdef _DEBUG
		}
		else
		{
			cout<<"out of size get"<<endl;
		}
#endif
	}

	PIXTYPE getXY(int ixy)
	{		// read value at 1D address ixy
		if(ixy<xsize*ysize)
		{
			return data[ixy];
		}
		else
		{
			cout<<"out of size get "<<endl;
		}
	}

	void putXY(int ixy,PIXTYPE val){// write value at 1D address ixy
		if (ixy<xsize*ysize)
		{
			data[ixy] = val;
		}
		else cout<<"out of size putxy"<<endl;

	}
}

#endif // RAW2D_H_INCLUDED
