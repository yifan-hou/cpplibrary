#pragma once
#ifndef _YIFAN_ULTILITIES_H_
#define _YIFAN_ULTILITIES_H_

#include <math.h>
#include <conio.h>
#include <stdio.h>
#include <stdarg.h>


namespace YF
{
	/////////////////////////////////////////////////////////////////////////
	//                   types and static variables 
	/////////////////////////////////////////////////////////////////////////

	// typedef int (*printf_ptr) (const char *str, ...);
	// static printf_ptr _printf = printf;


	/////////////////////////////////////////////////////////////////////////
	//                          print
	/////////////////////////////////////////////////////////////////////////

	static void print_header(const char *info)
	{
		printf("**************************************************************************\n");
		printf("\t\t");
		printf(info);
		printf("**************************************************************************\n");
	}

	static void PressEnterToContinue()
	{
		int c;
		printf( "Press ENTER to continue... " );
		fflush( stdout );
		do c = getchar(); while ((c != '\n') && (c != EOF));
	}


	static int error (const char *fmt, ...)
	{
		printf("\nERROR INFO: \n");

		// The following four lines do the job of
		// printf(fmt,...);
		va_list args;
		va_start(args,fmt);
		int rt = vprintf(fmt,args);
		va_end(args);

		PressEnterToContinue();
		return rt;
	}

	 
	/////////////////////////////////////////////////////////////////////////
	//                          scalar
	/////////////////////////////////////////////////////////////////////////

	static void truncate(double *ele, const double _max, const double _min)
	{
		if ( (*ele) > _max)
			(*ele) = _max;
		else if ( (*ele) < _min)
			(*ele) = _min;
	}

	/////////////////////////////////////////////////////////////////////////
	//                          vector
	/////////////////////////////////////////////////////////////////////////


	static void buf_insert(const double ele, const int size, double * buf)
	{
		for (int i = 1; i < size; ++i)
		{
			buf[size - i] = buf[size - 1 - i];
		}
		buf[0] = ele;
	}

	static double vec_max(const double * vec, const int size)
	{
		double m = vec[0];
		double t1;
		for (int i = 0; i < size; ++i)
		{
			t1 = vec[i];
			if (t1 > m) m = t1;
		}
		return m;
	}

	static double vec_min(const double * vec, const int size)
	{
		double m = vec[0];
		double t1;
		for (int i = 0; i < size; ++i)
		{
			t1 = vec[i];
			if (t1 < m) m = t1;
		}
		return m;
	}

	static double vec_mean(const double * vec, const int size)
	{
		double sum = 0;
		for (int i = 0; i < size; ++i)
		{
			sum += vec[i];
		}
		return sum/double(size);
	}


	static double vec_slope(const double * x, const double * y,const int size)
    {
    	double avgX = vec_mean(x,size);
    	double avgY = vec_mean(y,size);

    	double numerator = 0.0;
    	double denominator = 0.0;

    	double xd = 0;
    	for(int i=0; i<size; ++i)
    	{
    		xd = x[i] - avgX;
    		numerator += (xd) * (y[i] - avgY);
    		denominator += xd * xd;
    	}
    	
    	return numerator / denominator;
    }

    // numerical differentiation with low pass filter
    // input x, calculate dx/dt
    // s/(as+1),
	static double diff_LPF(const double xdold, const double xnew, const double xold, const double T,const double a)
    {
		double As = exp(-T/a);
		return As*xdold + (1 - As)*((xnew - xold)/T);
    }



}
	
#endif // _YIFAN_ULTILITIES_H_
