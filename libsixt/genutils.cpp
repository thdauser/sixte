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

   Copyright 2014:  Trigger has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01 and
   ESP2013-48637-C2-1-P.

/***********************************************************************
*                      GENUTILS
*
*  File:      genutils.cpp
*  Developer: Beatriz Cobo Martín
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene González Pérez
*             José Ramón Rodón Ortiz
*                                                                     
*  Revision History:                                                  
*                                                                     
*	version	1.0.0	21/09/06    First version               		  	
*	version	1.0.1	20/10/06    Include convol function		
* 	version	1.0.2	26/10/06    Modify critical points function
*	version 1.0.3	15/11/06	Adding new in/out put parameter in convol function:	
*								freqresult is the result in frequency
* 								Adding new function: Shiftleft	
*	version	1.0.4  	14/12/06	Polyfit function include a new in/out parameters: chisq
* 	version	1.0.5  	09/12/06    Modifying convol function adding padding to pulse vector
* 								Adding new function: convolFull: Is the same than convol adding edges.
*	version 1.0.6	31/05/07	Free allocate of memory
* 	version 1.1.0	22/04/08	Adding functions "getQuality" and "setQuality"
* 	version 1.4.0	16/09/08	Adding function "polyFit".
* 								Adding function "findUpperK" and "findLowK"
* 								Adding function "baseline" 
*   version 1.5.0	02/10/08    Adding function "fit_linear" 
* 	version 1.6.0	06/09/08	Adding new library "stdlib"
*  
* 	version 1.7.0	28/09/10	"fit_linear" function modified
*   07/01/11    Adding function "FFT"
*               Adding functions to handle vectors, complex vectors and complex matrix:
*           		gsl_vector_sqrtIFCA
*					gsl_vector_absIFCA
*					gsl_vector_powIFCA
*					gsl_vector_complex_addIFCA
*					gsl_vector_complex_mulIFCA
*					gsl_vector_complex_divIFCA
*					gsl_vector_complex_absIFCA
*					gsl_vector_complex_absRIFCA
*					gsl_vector_complex_conjugateIFCA
*					gsl_vector_complex_scaleIFCA
*					gsl_matrix_complex_absIFCA
*					gsl_matrix_complex_conjugateIFCA
*				Adding functions "Hw" and "absHw"
*   06/03/12    Adding "gsl_vector_complex_argIFCA" and "FFTinverse"
*     Dec/14    Deleted some non used functions
***********************************************************************/

/******************************************************************************
DESCRIPTION:
 
The objective of this package is to support general utilities.

MAP OF SECTIONS IN THIS FILE::

 - 1. polyFit
 - 2. FFT
 - 3. gsl_vector_sqrtIFCA
 - 4. gsl_vector_absIFCA
 - 5. gsl_vector_powIFCA
 - 6. gsl_vector_complex_addIFCA
 - 7. gsl_vector_complex_mulIFCA
 - 8. gsl_vector_complex_divIFCA
 - 9. gsl_vector_complex_absIFCA
 - 10. gsl_vector_complex_absRIFCA
 - 11. gsl_vector_complex_conjugateIFCA
 - 12. gsl_vector_complex_scaleIFCA
 - 13. gsl_matrix_complex_absIFCA
 - 14. gsl_matrix_complex_conjugateIFCA
 - 15. gsl_vector_complex_argIFCA
 - 16. FFTinverse
 - 17. exit_error
 - 18. print_error
 - 19. writeLog
 - 20. area0

*******************************************************************************/


/***** SECTION 1 ************************************
*       INCLUDE's
****************************************************/
#include "genutils.h"


/***** SECTION 1 ************************************************************
* polyFit: This function make a polynomial fit: ax² + bx + c using the regression quadratic analysis
****************************************************************************/
int polyFit (gsl_vector *x_fit, gsl_vector *y_fit, double *a, double *b, double *c)
{
	int status=EPOK;
	double sxx=0, sxy=0, sxx2=0, sx2y=0, sx2x2=0; 
	double x=0, x2=0, x3=0, x4=0; 
	double y=0, xy=0, x2y=0;
	int n = x_fit->size;
	
 	for (int i=0; i<n;i++)
 	{
		x = x + gsl_vector_get(x_fit,i);
		x2 = x2 + pow(gsl_vector_get(x_fit,i),2);
		x3 = x3 + pow(gsl_vector_get(x_fit,i),3);
		x4 = x4 + pow(gsl_vector_get(x_fit,i),4);
		
		y = y + gsl_vector_get(y_fit,i);
		
		xy = xy + ( gsl_vector_get(x_fit,i)* gsl_vector_get(y_fit,i));
		x2y = x2y + (pow(gsl_vector_get(x_fit,i),2) * gsl_vector_get(y_fit,i));				
	}
 	
 	sxx 	= x2 - (pow(x,2)/n);
 	sxy 	= xy - (x*y/n);
 	sxx2 	= x3 - (x*x2/n);
 	sx2y 	= x2y - (x2*y/n);
 	sx2x2	= x4 - (pow(x2,2)/n);

	*a = ((sx2y * sxx)  - (sxy  * sxx2))/ ((sxx * sx2x2)- pow(fabs(sxx2),2));
	*b = ((sxy  * sx2x2)- (sx2y * sxx2))/ ((sxx * sx2x2) -pow(fabs(sxx2),2));
	*c = (y/n) - ((*b)*x/n) - ((*a)*x2/n);
	 	
	return EPOK;
}
/*xxxx end of SECTION 1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 2 ************************************************************
* FFT: This function calculates the FFT of the elements of a vector
*
* GSL library (overview of FFTs):
*    For physical applications it is important to remember that the index appearing in the DFT does not correspond
*    directly to a physical frequency. If the time-step of the DFT is \Delta then the frequency-domain includes both
*    positive and negative frequencies, ranging from -1/(2\Delta) through 0 to +1/(2\Delta). The positive frequencies
*    are stored from the beginning of the array up to the middle, and the negative frequencies are stored backwards
*    from the end of the array.
*
*    Here is a table which shows the layout of the array data, and the correspondence between the time-domain data z,
*    and the frequency-domain data x.
*
*     index    z               x = FFT(z)
*
*     0        z(t = 0)        x(f = 0)
*     1        z(t = 1)        x(f = 1/(n Delta))
*     2        z(t = 2)        x(f = 2/(n Delta))
*     .        ........        ..................
*     n/2      z(t = n/2)      x(f = +1/(2 Delta),
*                                    -1/(2 Delta))
*     .        ........        ..................
*     n-3      z(t = n-3)      x(f = -3/(n Delta))
*     n-2      z(t = n-2)      x(f = -2/(n Delta))
*     n-1      z(t = n-1)      x(f = -1/(n Delta))
*
* Therefore, when frequency axis is built f = i/SelectedTimeDuration = i/(size/samprate)
* (n=size and Delta=1/samprate seg/sample) with i varying from 0 to size/2-1.
*****************************************************************************/
int FFT(gsl_vector *invector,gsl_vector_complex *outvector,double STD)
{
	int status=EPOK;

	//Declare variables
 	gsl_fft_complex_workspace * work;
 	gsl_fft_complex_wavetable * wavetable;
 	#define REAL(z,i) ((z)[2*(i)])
 	#define IMAG(z,i) ((z)[2*(i)+1])
 	double data[2*invector->size];

 	//FFT calculus
 	work = gsl_fft_complex_workspace_alloc(invector->size);
 	wavetable = gsl_fft_complex_wavetable_alloc(invector->size);
 		//In order to work with complex numbers
 	for (int i=0; i< invector->size; i++) {
 		REAL(data,i) = gsl_vector_get(invector,i);
 		IMAG(data,i) = 0.0;
 	}
 	gsl_fft_complex_forward(data,1,invector->size,wavetable,work);
 	for (int i=0; i< invector->size; i++) {
 		gsl_vector_complex_set(outvector,i,gsl_complex_rect(REAL(data,i),IMAG(data,i)));
 	}
 	gsl_vector_complex_scaleIFCA(outvector,gsl_complex_rect(1.0/sqrt(invector->size),0.0)); //Factor 1/sqrt(N), in the FFT expression
 	gsl_vector_complex_scaleIFCA(outvector,gsl_complex_rect(sqrt(2*STD),0.0));

 	gsl_fft_complex_wavetable_free(wavetable);
 	gsl_fft_complex_workspace_free(work);

 	return EPOK;
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 ************************************************************
* gsl_vector_sqrtIFCA: This function calculates the sqrt of the elements of a vector
*****************************************************************************/
gsl_vector gsl_vector_sqrtIFCA(gsl_vector *cvnew,gsl_vector *cv)
{
	size_t i;

	for (i = 0; i < cv->size; i++) {
	    gsl_vector_set(cvnew, i, sqrt(gsl_vector_get(cv,i)));
	}
}
/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
* gsl_vector_absIFCA: This function calculates the absolute value of the elements of a vector
*****************************************************************************/
gsl_vector gsl_vector_absIFCA(gsl_vector *vnew,gsl_vector *v)
{
	size_t i;

	for (i = 0; i < v->size; i++)
	{
		if (gsl_vector_get(v,i)<0) {gsl_vector_set(vnew,i,gsl_vector_get(v,i)*(-1));}
		else {gsl_vector_set(vnew,i,gsl_vector_get(v,i));}
	}
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* gsl_vector_powIFCA: This function calculates the power of the elements of a vector
*****************************************************************************/
gsl_vector gsl_vector_powIFCA(gsl_vector *cvnew,gsl_vector *cv, double exp)
{
	size_t i;

	for (i = 0; i < cv->size; i++)
	{
	    gsl_vector_set(cvnew, i, pow(gsl_vector_get(cv,i),exp));
	}
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
* gsl_vector_complex_addIFCA: This function adds the complex elements of two vectors
*****************************************************************************/
gsl_vector_complex gsl_vector_complex_addIFCA(gsl_vector_complex *cvnew,gsl_vector_complex *cv)
{
	gsl_complex z,z1;
	size_t i;

	for (i = 0; i < cv->size; i++)
	{
	    z = gsl_vector_complex_get(cv, i);
	    z1 = gsl_vector_complex_get(cvnew, i);
	    gsl_vector_complex_set(cvnew, i, gsl_complex_add(z,z1));
	}
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* gsl_vector_complex_mulIFCA: This function multiplies the complex elements of two vectors
*****************************************************************************/
gsl_vector_complex gsl_vector_complex_mulIFCA(gsl_vector_complex *cvnew,gsl_vector_complex *cv)
{
	gsl_complex z,z1;
	size_t i;

	for (i = 0; i < cv->size; i++)
	{
	    z = gsl_vector_complex_get(cv, i);
	    z1 = gsl_vector_complex_get(cvnew, i);
	    gsl_vector_complex_set(cvnew, i, gsl_complex_mul(z,z1));
	}
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
* gsl_vector_complex_divIFCA: This function divides the complex elements of two vectors
*****************************************************************************/
gsl_vector_complex gsl_vector_complex_divIFCA(gsl_vector_complex *cvnew,gsl_vector_complex *cv)
{
	gsl_complex z,z1;
	size_t i;

	for (i = 0; i < cv->size; i++)
	{
	    z = gsl_vector_complex_get(cvnew, i);
	    z1 = gsl_vector_complex_get(cv, i);
	    gsl_vector_complex_set(cvnew, i, gsl_complex_div(z,z1));
	}
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 9 ************************************************************
* gsl_vector_complex_absIFCA: This function calculates the magnitude of the complex elements of a vector.
*                             It returns a complex number whose real part is the abs and whose imaginary part is zero.
*****************************************************************************/
gsl_vector_complex gsl_vector_complex_absIFCA(gsl_vector_complex *cvnew,gsl_vector_complex *cv)
{
	gsl_complex z;
	size_t i;

	for (i = 0; i < cv->size; i++)
	{
	    z = gsl_vector_complex_get(cv, i);
	    gsl_vector_complex_set(cvnew, i, gsl_complex_rect(gsl_complex_abs(z),0.0));
	}
}
/*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 10 ************************************************************
* gsl_vector_complex_absRIFCA: This function calculates the magnitude of the complex elements of a vector.
*****************************************************************************/
gsl_vector gsl_vector_complex_absRIFCA(gsl_vector *cvnew,gsl_vector_complex *cv)
{
	gsl_complex z;
	size_t i;

	for (i = 0; i < cv->size; i++)
	{
	    z = gsl_vector_complex_get(cv, i);
	    gsl_vector_set(cvnew, i, gsl_complex_abs(z));
	}
}
/*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 11 ************************************************************
* gsl_vector_complex_conjugateIFCA: This function conjugates the complex elements of a vector
*****************************************************************************/
gsl_vector_complex gsl_vector_complex_conjugateIFCA(gsl_vector_complex *cvnew,gsl_vector_complex *cv)
{
	gsl_complex z;
	size_t i;

	for (i = 0; i < cv->size; i++)
	{
	    z = gsl_vector_complex_get(cv, i);
	    gsl_vector_complex_set(cvnew, i, gsl_complex_conjugate(z));
	}
}
/*xxxx end of SECTION 11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 12 ************************************************************
* gsl_vector_complex_scaleIFCA: This function multiplies the complex elements of a vector by a complex number
*****************************************************************************/
gsl_vector_complex gsl_vector_complex_scaleIFCA(gsl_vector_complex *cv,gsl_complex z)
{
	gsl_complex z_aux;
	size_t i;

	for (i = 0; i < cv->size; i++)
	{
		z_aux = gsl_complex_mul(gsl_vector_complex_get(cv, i),z);
		gsl_vector_complex_set(cv, i, z_aux);
	}
}
/*xxxx end of SECTION 12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 13 ************************************************************
* gsl_matrix_complex_absIFCA: This function calculates the magnitude of the complex elements of a matrix.
*                             It returns a complex number whose real part is the abs and whose imaginary part is zero.
*****************************************************************************/
gsl_matrix_complex gsl_matrix_complex_absIFCA(gsl_matrix_complex *cmnew,gsl_matrix_complex *cm)
{
	gsl_complex z;
	size_t i, j;

	for (i = 0; i < cm->size1; i++)
	{
		for (j = 0; j < cm->size2; j++)
		{
			z = gsl_matrix_complex_get(cm, i, j);
			gsl_matrix_complex_set(cmnew, i, j, gsl_complex_rect(gsl_complex_abs(z),0.0));
		}
	}
}
/*xxxx end of SECTION 13 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 14 ************************************************************
* gsl_matrix_complex_conjugateIFCA: This function conjugates the complex elements of a matrix
*****************************************************************************/
gsl_matrix_complex gsl_matrix_complex_conjugateIFCA(gsl_matrix_complex *cmnew,gsl_matrix_complex *cm)
{
	gsl_complex z;
	size_t i,j;

	for (i = 0; i < cm->size1; i++)
	{
		for (j = 0; j < cm->size2; j++)
		{
			z = gsl_matrix_complex_get(cm, i, j);
			gsl_matrix_complex_set(cmnew, i, j, gsl_complex_conjugate(z));
		}
	}
}
/*xxxx end of SECTION 14 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 15 ************************************************************
* gsl_vector_complex_argIFCA: This function calculates the arguments of the complex elements of a vector
*****************************************************************************/
gsl_vector gsl_vector_complex_argIFCA(gsl_vector *varg,gsl_vector_complex *vin)
{
	gsl_complex z;
	size_t i;

	for (i = 0; i < vin->size; i++)
	{
		z = gsl_vector_complex_get(vin, i);
		gsl_vector_set(varg, i, gsl_complex_arg(z));
	}
}
/*xxxx end of SECTION 15 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 16 ************************************************************
* gsl_vector_Sumsubvector: This function sums some of the elements of the input vector and it writes the result in 'sum'.
*
* The starting element of the sum is 'offset' from the start of the input vector.
* It is going to be summed 'n' elements.
*
* Offset can take values from 0 to invector->size.
*
* Arguments:
*      - invector: Input vector
*      - offset: It is the first element to sum
*      - n: Number of elements to sum
*      - sum: Calculated output value (sum of the corresponding elements)
*      - status: Auxiliary variable which has the function status: Error, Warning or OK
*
****************************************/
int gsl_vector_Sumsubvector(gsl_vector *invector, long offset, long n, double *sum)
{
    int status=EPOK;
	*sum = 0.0;

    for (int i=offset; i<offset+n; i++)
    {
    	*sum = *sum + gsl_vector_get(invector,i);
    }

    return(EPOK);
}
/*xxxx end of SECTION 16 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 17 ************************************************************
* FFTinverse: This function calculates the inverse FFT of the elements of a vector
*****************************************************************************/
int FFTinverse(gsl_vector_complex *invector,gsl_vector *outvector,double STD)
{
	//Declare variables
	int status=EPOK;
	gsl_fft_complex_workspace * work (gsl_fft_complex_workspace_alloc(invector->size));
	gsl_fft_complex_wavetable * wavetable (gsl_fft_complex_wavetable_alloc(invector->size));
	std::vector<std::complex<double> >cd (invector->size,std::complex<double> (0, 0));
	gsl_complex_packed_array cpa (reinterpret_cast<double *>(&cd[0]));

	for (int i=0;i<invector->size;i++) cd[i]=std::complex<double> (GSL_REAL(gsl_vector_complex_get(invector,i)),GSL_IMAG(gsl_vector_complex_get(invector,i)));

	//Inverse FFT calculus
	gsl_fft_complex_inverse(cpa,1,invector->size,wavetable,work);

	for (int i=0;i<invector->size;i++)
	{
		gsl_vector_set(outvector,i,cd[i].real());
	}
	gsl_vector_scale(outvector,sqrt(invector->size)); //Factor 1/sqrt(N), in the FFT expression
	gsl_vector_scale(outvector,1/sqrt(2*STD));

	gsl_fft_complex_wavetable_free(wavetable);
	gsl_fft_complex_workspace_free(work);

 	return EPOK;
}
/*xxxx end of SECTION 17 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 18 ************************************************************
* exit_error function: This function prints out error messages and exits program
*****************************************************************************/
void exit_error(const char* const func, string msg, int status)
{
    cout<<"Error in " <<string(func) << ": " << msg<<endl;
    if (status > 10)  // CFITSIO Error
    {
	fits_report_error(stderr, status); /* print error report */
       //exit( status );    /* terminate the program, returning error status */
    }
    exit(EPFAIL);
}
/*xxxx end of SECTION 18 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 19 ************************************************************
* print_error function: This function prints out error messages
*****************************************************************************/
void print_error( const char* const func, string message, int status)
{
    if (status > 0 && status <= 10)
    {
	cout<<"Error in "<<string(func)<<": "<<message<<endl;
    }
    else if (status > 10)  // CFITSIO Error
    {
    	fits_report_error(stderr, status); /* print error report */
    }
    else if (status == -999)
    {
        cout<<"Warning in "<<string(func)<<": "<<message<<endl;
    }
    return;
}
/*xxxx end of SECTION 19 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 20 ************************************************************
* writeLog function: This function includes the processing of the each level of message in the log file
*                    and the output screen
*
* Verbosity = 0 --> The log file and the output screen include Errors
* Verbosity = 1 --> The log file and the output screen include Errors and Warnings
* Verbosity = 2 --> The log file and the output screen include Errors, Warnings and Alerts
* Verbosity = 3 --> The log file and the output screen include Errors, Warnings, Alerts and Log messages
*
* Parameters:
* - fileRef (in): File reference to log file
* - type (in): String to indicate error type "Error", "Warning", "Alert","Log" or "OK"
* - verbosity (in): Integer value for verbosity
* -	message (in): String message to print
*****************************************************************************/
void writeLog (FILE *fileRef, string type, int verbosity, string message)
{
	if  ( type == "Error" || type == "OK"      ||
	    ((type == "Warning") && (verbosity>0)) ||
	    ((type == "Alert") && (verbosity>1))   ||
	    ((type == "Log") && (verbosity==3)) )
	{
		cout << type << ": " << message << "\n";
	    fprintf(fileRef,"%s: %s\n",type.c_str(),message.c_str());
	}

	return;
}
/*xxxx end of SECTION 20 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 21 ************************************
* area0 function: This function shifts the filter template vertically to set its area to zero (with equal areas below
*                 and above the x-axis). Therefore, in the frequency domain it means to remove the bin 0, i.e., the
*                 baseline, essentially in order to be baseline insensitive, to calculate properly the pulse energy
*                 in the next tasks.
****************************************************/
int area0 (gsl_vector **input)
{
	int status = EPOK;
	string message = "";

	double sum=0;
	int size_input = (*input)->size;
	gsl_vector *sumgsl=gsl_vector_alloc(size_input);
	gsl_vector *inputsize=gsl_vector_alloc(size_input);

	for (int i=0; i<size_input; i++)
	{
		sum=sum+gsl_vector_get(*input,i);
	}
	gsl_vector_set_all(sumgsl,sum);

	gsl_vector_set_all(inputsize,size_input);
	gsl_vector_div(sumgsl,inputsize);
	gsl_vector_sub(*input,sumgsl);

	// Free allocate of GSL vectors
	gsl_vector_free(sumgsl);
	gsl_vector_free(inputsize);

	return EPOK;
}
/*xxxx end of SECTION 21  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/* Check for file existence
 * Filename is a char
 */
bool fileExists(const std::string& name) {
  struct stat buffer;   
  return (stat(name.c_str(), &buffer) == 0); 
}
