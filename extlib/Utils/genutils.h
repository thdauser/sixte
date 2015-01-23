/***********************************************************************
   This file is part of SIXTE.

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

   Copyright 2014:  Trigger has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01 and
   ESP2013-48637-C2-1-P.

/***********************************************************************
*                      GENUTILS
*
*  File:      genutils.h
*  Developer: Beatriz Cobo Martí­n
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene González Pérez
*             José Ramón Rodón Ortiz
*                                                                     
*  Revision History:                                                  
*                                                                     
*  21/09/06	First version
*  20/10/06 Include convol function
*  14/12/06	Polyfit function include a new in/out parameters: chisq
*  09/01/07 Include convolfull
*  22/04/08	Including functions "getQuality" and "setQuality"
*  16/09/08	Adding function "polyFit" and "linear_fit"
* 			Adding function "findUpperK" and "findLowK"
* 			Adding function "baseline"
*  06/09/08	Adding new library "stdlib"
*  28/09/10	"fit_linear" function modified
*  07/01/11 Adding function "FFT"
*           Adding functions to handle vectors, complex vectors and complex matrix:
*           	gsl_vector_sqrtIFCA
*				gsl_vector_absIFCA
*				gsl_vector_powIFCA
*				gsl_vector_complex_addIFCA
*				gsl_vector_complex_mulIFCA
*				gsl_vector_complex_divIFCA
*				gsl_vector_complex_absIFCA
*				gsl_vector_complex_absRIFCA
*				gsl_vector_complex_conjugateIFCA
*				gsl_vector_complex_scaleIFCA
*				gsl_matrix_complex_absIFCA
*				gsl_matrix_complex_conjugateIFCA
*			Adding functions "Hw" and "absHw"
*  29/03/11 Updated .h
*  06/03/12 Adding "gsl_vector_complex_argIFCA" and "FFTinverse"
*           New libraries "vector" and "complex"
*    Dec/14 Deleted some non used functions
************************************************************************/
#ifndef GENUTILS_H_
#define GENUTILS_H_

// GSL

	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_multifit.h>
	#include <gsl/gsl_fft_complex.h>
	#include <gsl/gsl_complex_math.h>
	#include <gsl/gsl_sort.h>
	#include <gsl/gsl_sort_vector.h>
	#include <gsl/gsl_statistics.h>
	#include <gsl/gsl_interp.h>
	#include <gsl/gsl_spline.h>

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

	typedef struct
	{
		int errCODE;
		string errNAME;
		string errMESS;
	} EPerr;

	int polyFit (gsl_vector *x_fit, gsl_vector *y_fit, double *a, double *b, double *c);
	int FFT(gsl_vector *invector,gsl_vector_complex *outvector,double STD);
	int FFTinverse(gsl_vector_complex *invector,gsl_vector *outvector,double STD);
	int area0 (gsl_vector **input);

	//GSL vectors and matrix

	gsl_vector gsl_vector_sqrtIFCA(gsl_vector *cvnew,gsl_vector *cv);
	gsl_vector gsl_vector_absIFCA(gsl_vector *vnew,gsl_vector *v);
	gsl_vector gsl_vector_powIFCA(gsl_vector *cvnew,gsl_vector *cv, double exp);

	gsl_vector_complex gsl_vector_complex_addIFCA(gsl_vector_complex *cvnew,gsl_vector_complex *cv);
	gsl_vector_complex gsl_vector_complex_mulIFCA(gsl_vector_complex *cvnew,gsl_vector_complex *cv);
	gsl_vector_complex gsl_vector_complex_divIFCA(gsl_vector_complex *cvnew,gsl_vector_complex *cv);
	gsl_vector_complex gsl_vector_complex_absIFCA(gsl_vector_complex *cvnew,gsl_vector_complex *cv);
	gsl_vector gsl_vector_complex_absRIFCA(gsl_vector *cvnew,gsl_vector_complex *cv);
	gsl_vector_complex gsl_vector_complex_conjugateIFCA(gsl_vector_complex *cvnew,gsl_vector_complex *cv);
	gsl_vector_complex gsl_vector_complex_scaleIFCA(gsl_vector_complex *cv,gsl_complex z);
	gsl_matrix_complex gsl_matrix_complex_absIFCA(gsl_matrix_complex *cmnew,gsl_matrix_complex *cm);
	gsl_matrix_complex gsl_matrix_complex_conjugateIFCA(gsl_matrix_complex *cmnew,gsl_matrix_complex *cm);
	gsl_vector gsl_vector_complex_argIFCA(gsl_vector *varg,gsl_vector_complex *vin);
	int gsl_vector_Sumsubvector(gsl_vector *invector, long offset, long n, double *sum);

	EPerr readErrCodes(char *errorfile, int errNum);
	void EPexit(int instatus=0, string fname="");
	void print_error( const char* const func, string message, int status);
	void writeLog (FILE *fileRef, string type, int verbosity, string message);
	void exit_error(const char* const func, string msg,int status);
	bool fileExists (const std::string& name);
	
//using namespace std;

#endif /*GENUTILS_H_*/
