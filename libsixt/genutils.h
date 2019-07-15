/***********************************************************************
   This file is part of SIXTE/SIRENA software.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

   Copyright 2014:  GENUTILS has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01, 
   ESP2013-48637-C2-1-P and ESP2014-53672-C3-1-P.

/***********************************************************************
*                      GENUTILS
*
*  File:       genutils.h
*  Developers: Beatriz Cobo
* 	       cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef GENUTILS_H_
#define GENUTILS_H_

// GSL

	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_blas.h>
	#include <gsl/gsl_multifit.h>
	#include <gsl/gsl_fft_complex.h>
	#include <gsl/gsl_complex_math.h>
	#include <gsl/gsl_sort.h>
	#include <gsl/gsl_sort_vector.h>
	#include <gsl/gsl_statistics.h>
	#include <gsl/gsl_interp.h>
	#include <gsl/gsl_spline.h>
	#include <fftw3.h>
	#include <gsl/gsl_eigen.h>
	#include <gsl/gsl_wavelet.h>
	#include <gsl/gsl_errno.h>

// General

	#include <fitsio.h>
	#include <math.h>
	#include <iostream>
	#include <boost/lexical_cast.hpp>
	#include <vector>
        #include <complex>
	#include <getopt.h> // For getopt module
	#include <stdarg.h>
	#include <string>
	#include <string.h>
	#include <stdio.h>
	#include <fstream>
	using std::ifstream;
	#include <sstream>
	#include <ctype.h>
	#include <sys/stat.h>
	#include "assert.h"

	#ifndef EPOK	//Event Processing OK
	#define EPOK                 (0)
	#endif

	#ifndef EPFAIL	//Event Processing Failure
	#define EPFAIL                 (1)
	#endif

	#define EP_EXIT_ERROR(msg,status) (exit_error(__func__, msg, status))

	#define EP_PRINT_ERROR(msg,status) (print_error(__func__, msg, status))

	using namespace std;

	const double pi = 4.0 * atan(1.0);
	
	//const int safetyMargin = 50; // In samples

	int polyFit (gsl_vector *x_fit, gsl_vector *y_fit, double *a, double *b, double *c);
	int polyFitLinear (gsl_vector *x_fit, gsl_vector *y_fit, double *a, double *b);
	int FFT(gsl_vector *invector,gsl_vector_complex *outvector,double STD);
	int FFTinverse(gsl_vector_complex *invector,gsl_vector *outvector,double STD);

	//GSL vectors
	void gsl_vector_sqrtIFCA(gsl_vector *cvnew,gsl_vector *cv);
	void gsl_vector_complex_absIFCA(gsl_vector *cvnew,gsl_vector_complex *cv);
	void gsl_vector_complex_scaleIFCA(gsl_vector_complex *cv,gsl_complex z);
	void gsl_vector_complex_argIFCA(gsl_vector *varg,gsl_vector_complex *vin);
	int gsl_vector_Sumsubvector(gsl_vector *invector, long offset, long n, double *sum);

	void print_error( const char* const func, string message, int status);
	void writeLog (FILE *fileRef, string type, int verbosity, string message);
	void exit_error(const char* const func, string msg,int status);
	bool fileExists (const std::string& name);
	
	int parabola3Pts (gsl_vector *x, gsl_vector *y, double *a, double *b, double *c);
        
        bool isNumber(string s);

#endif /*GENUTILS_H_*/
