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
*  File:       genutils.cpp
*  Developers: Beatriz Cobo
* 	       cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

/******************************************************************************
DESCRIPTION:
 
The objective of this package is to support general utilities.

MAP OF SECTIONS IN THIS FILE::

 - 1. polyFit
 - 2. polyFitLinear
 - 3. FFT
 - 4. FFTinverse
 - 5. gsl_vector_sqrtIFCA
 - 6. gsl_vector_complex_absIFCA
 - 7. gsl_vector_complex_scaleIFCA
 - 8. gsl_vector_complex_argIFCA
 - 9. gsl_vector_Sumsubvector
 - 10. exit_error
 - 11. print_error
 - 12. writeLog
 - 13. fileExists
 - 14. parabola3Pts

*******************************************************************************/

#include "genutils.h"

/***** SECTION 1 ************************************************************
* polyFit: This function makes a polynomial fitting: ax² + bx + c using the regression quadratic analysis
*          To measure how well the model agrees with the data, the chi-square merit function is used
* 
* Parameters:
* - x_fit: Input GSL x vector
* - y_fit: Input GSL y vector
* - a: Fit coefficient of the quadratic term
* - b: Fit coefficient of the linear term
* - c: Fit coefficient (independent term)
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
* polyFitLinear: This function makes a linear fitting: ax + b using the regression linear analysis
*                To measure how well the model agrees with the data, the chi-square merit function is used
* 
* Parameters:
* - x_fit: Input GSL x vector
* - y_fit: Input GSL y vector
* - a: Fit coefficient of the quadratic term
* - b: Fit coefficient of the linear term
****************************************************************************/
int polyFitLinear (gsl_vector *x_fit, gsl_vector *y_fit, double *a, double *b)
{
	int status=EPOK;
	double sxx=0, sxy=0, sxx2=0, sx2y=0, sx2x2=0; 
	double x=0, x2=0; 
	double y=0, xy=0;
	int n = x_fit->size;
	
 	for (int i=0; i<n;i++)
 	{
		x = x + gsl_vector_get(x_fit,i);
		x2 = x2 + pow(gsl_vector_get(x_fit,i),2);
		
		y = y + gsl_vector_get(y_fit,i);
		
		xy = xy + ( gsl_vector_get(x_fit,i)* gsl_vector_get(y_fit,i));
	}
	
	*a = (n*xy-x*y) / (n*x2-pow(x,2));
        *b = y/n - (*a)*x/n;
        	 	
	return EPOK;
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 ************************************************************
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
* The frequency axis will be as built f = i/STD = i/(size/samprate) with i varying from 0 to size/2-1 (n=size and Delta=1/samprate sec/sample).
* 
* Parameters:
* - invector: Input GSL vector
* - outvector: Output GSL complex vector with the FFT of invector
* - STD: SelectedTimeDuration=(Size of invector)/samprate
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

 	gsl_fft_complex_wavetable_free(wavetable); wavetable = 0;
 	gsl_fft_complex_workspace_free(work); work = 0;

 	return EPOK;
}
/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
* FFTinverse: This function calculates the inverse FFT of the elements of a vector
* 
* Parameters:
* - invector: Input GSL complex vector
* - outvector: Output GSL vector with the inverse FFT of invector
* - STD: SelectedTimeDuration=(Size of invector)/samprate
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

	gsl_fft_complex_wavetable_free(wavetable); wavetable = 0;
	gsl_fft_complex_workspace_free(work); work = 0;

 	return EPOK;
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* gsl_vector_sqrtIFCA: This function calculates the sqrt of the elements of a vector
* 
* Parameters:
* - cv: Input GSL vector
* - cvnew: Output GSL vector with the square root values of the elements of cv
*****************************************************************************/
void gsl_vector_sqrtIFCA(gsl_vector *cvnew,gsl_vector *cv)
{
	size_t i;

	for (i = 0; i < cv->size; i++) {
	    gsl_vector_set(cvnew, i, sqrt(gsl_vector_get(cv,i)));
	}
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
* gsl_vector_complex_absIFCA: This function calculates the magnitude of the complex elements of a vector (real part)
* 
* Parameters:
* - cv: Input GSL complex vector
* - cvnew: Output GSL vector with the absolute values of the elements of cv
*****************************************************************************/
void gsl_vector_complex_absIFCA(gsl_vector *cvnew,gsl_vector_complex *cv)
{
	gsl_complex z;
	size_t i;

	for (i = 0; i < cv->size; i++)
	{
	    z = gsl_vector_complex_get(cv, i);
	    gsl_vector_set(cvnew, i, gsl_complex_abs(z));
	}
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* gsl_vector_complex_scaleIFCA: This function multiplies the complex elements of a vector by a complex number
* 
* Parameters:
* - cv: Input/Output (scaled) GSL complex vector
* - z: Input GSL complex number 
*****************************************************************************/
void gsl_vector_complex_scaleIFCA(gsl_vector_complex *cv,gsl_complex z)
{
	gsl_complex z_aux;
	size_t i;

	for (i = 0; i < cv->size; i++)
	{
		z_aux = gsl_complex_mul(gsl_vector_complex_get(cv, i),z);
		gsl_vector_complex_set(cv, i, z_aux);
	}
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
* gsl_vector_complex_argIFCA: This function calculates the arguments of the complex elements of a vector
* 
* Parameters: 
* - vin: Input GSL complex vector
* - varg: Output GSL vector with the arguments of the elements of vin
*****************************************************************************/
void gsl_vector_complex_argIFCA(gsl_vector *varg,gsl_vector_complex *vin)
{
	gsl_complex z;
	size_t i;

	for (i = 0; i < vin->size; i++)
	{
		z = gsl_vector_complex_get(vin, i);
		gsl_vector_set(varg, i, gsl_complex_arg(z));
	}
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 9 ************************************************************
* gsl_vector_Sumsubvector: This function returns the sum of some elements of the input vector.
*
* The starting element of the sum is 'offset' from the start of the input vector.
* It will sum up n elements.
*
* Offset can take values from 0 to invector->size.
* 
* Parameters:
* - invector: Input GSL vector
* - offset: It is the first element to be summed
* - n: Number of elements in the sum
* - sum: Calculated output value (sum of the corresponding elements)
****************************************/
int gsl_vector_Sumsubvector(gsl_vector *invector, long offset, long n, double *sum)
{
	int status=EPOK;
	string message = "";
	char valERROR[256];
	
	*sum = 0.0;

	if ((offset < 0) || (offset+n > invector->size))
	{
		sprintf(valERROR,"%d",__LINE__+7);
		string str(valERROR);
		message = "Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	for (int i=offset; i<offset+n; i++)
	{
		*sum = *sum + gsl_vector_get(invector,i);
	}

	return(EPOK);
}
/*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 10 ************************************************************
* exit_error function: This function prints out error messages and exits program
* 
* Parameters:
* - func: Function name whose error is printed
* - msg: Error message to be printed
* - status: Status
*****************************************************************************/
void exit_error(const char* const func, string msg, int status)
{
	cout<<"Error in " <<string(func) << ": " << msg<<endl;
	if (status > 10)  // CFITSIO Error
	{
		fits_report_error(stderr, status); // Print error report 
		//exit( status );    		   // Terminate the program, returning error status 
	}
	exit(EPFAIL);
}
/*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 11 ************************************************************
* print_error function: This function prints out error messages
* 
* Parameters:
* - func: Function name whose error is printed
* - message: Error message to be printed
* - status: Status
*****************************************************************************/
void print_error(const char* const func, string message, int status)
{
	if (status > 0 && status <= 10)
	{
		cout<<"Error in "<<string(func)<<": "<<message<<endl;
	}
	else if (status > 10)  // CFITSIO Error
	{
		fits_report_error(stderr, status); // print error report 
	}
	else if (status == -999)
	{
		cout<<"Warning in "<<string(func)<<": "<<message<<endl;
	}
	return;
}
/*xxxx end of SECTION 11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 12 ************************************************************
* writeLog function: This function includes the processing of the each level of message in the log file
*                    and the output screen
*
* Verbosity = 0 --> The log file and the output screen include Errors
* Verbosity = 1 --> The log file and the output screen include Errors and Warnings
* Verbosity = 2 --> The log file and the output screen include Errors, Warnings and Alerts
* Verbosity = 3 --> The log file and the output screen include Errors, Warnings, Alerts and Log messages
*
* Parameters:
* - fileRef: File reference to log file
* - type: String to indicate error type "Error", "Warning", "Alert","Log" or "OK"
* - verbosity: Integer value for verbosity
* - message: String message to print
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
/*xxxx end of SECTION 12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 13 ************************************************************
* fileExists function: This function checks for file existence returning a boolean value
*                       Filename is a char
* Parameters:
* - name: Filename
*****************************************************************************/
bool fileExists(const std::string& name) 
{
	struct stat buffer;   
	return (stat(name.c_str(), &buffer) == 0); 
}
/*xxxx end of SECTION 13 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 14 ************************************************************
* parabola3Pts: This function calculates the equation of a parabola given 3 points
* 
* Parameters:
* - x: Input GSL x vector
* - y: Input GSL y vector
* - a: Fit coefficient of the quadratic term
* - b: Fit coefficient of the linear term
* - c: Fit coefficient (independent term)
****************************************************************************/
int parabola3Pts (gsl_vector *x, gsl_vector *y, double *a, double *b, double *c)
{
	int status=EPOK;
	double x1,x2,x3,y1,y2,y3;
	
	x1 = gsl_vector_get(x,0);
	x2 = gsl_vector_get(x,1);
	x3 = gsl_vector_get(x,2);
	y1 = gsl_vector_get(y,0);
	y2 = gsl_vector_get(y,1);
	y3 = gsl_vector_get(y,2);	
	
	*a = (((y1-y2)/(x1-x2))-((y2-y3)/(x2-x3)))/(((pow(x1,2)-pow(x2,2))/(x1-x2))-((pow(x2,2)-pow(x3,2))/(x2-x3)));
	*b = (y1-y2-(*a)*(pow(x1,2)-pow(x2,2)))/(x1-x2);
	*c = y1-(*a)*pow(x1,2)-(*b)*x1;
	 	
	return EPOK;
}
/*xxxx end of SECTION 14 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
