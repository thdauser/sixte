
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

   Copyright 2014:  PULSEPROCESS has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish (
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01, 
   ESP2013-48637-C2-1-P, ESP2014-53672-C3-1-P and RTI2018-096686-B-C21.

***********************************************************************
*                      PULSEPROCESS
*
*  File:       pulseprocess.cpp
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

/******************************************************************************
DESCRIPTION:

The purpose of this package is to support utilities related to find pulses.

MAP OF SECTIONS IN THIS FILE:

 - 1. lpf_boxcar
 - 2. differentiate
 - 3. findMeanSigma
 - 4. medianKappaClipping
 - 5. getB
 - 6. getPulseHeight
 - 7. RS_filter
 - 8. find_model_energies
 - 9. find_model_maxDERs
 - 10. interpolate_model
 - 11. findPulsesCAL
 - 12. findTstartCAL
 - 13. InitialTriggering
 - 14. FindSecondaries
 - 15. find_model_samp1DERsNoReSCLD
 - 16. smoothDerivative
 - 17. noDetect

*******************************************************************************/

#include "pulseprocess.h"
#include <limits>

/***** SECTION 1 ************************************************************
* lpf_boxcar function: This function implements a low pass filtering as a box-car function in time
*
* The box-car function is a temporal average window:
*
* x_(i-1) = (sum(i=0,n-1)(Ii))/n	x_(i) = (sum(i=1,n)(Ii))/n
*
* If the cut frequency of the filter is fc, the box-car length (n) is (1/fc)*samprate
*
* According to Jan:
* 	sinc(f1)=0.6 where f1=1/(2pi*scaleFactor)
*	sinc(0.5)~0.6 => f1~0.5
*   	sinc(fc)=0, sinc(1)=0 => fc=1
*   	fc=kf1 => fc~2f1
*
* - Declare variables
* - Define the LPF (frequency domain) and the box-car function (time domain)
* - It is going to work with a longer vector to not have fake results for the last boxLength windows
* - Apply the box-car window by shifting it along the (lengthened) input vector
* - Free allocated GSL vectors
*
*  The function returns:
*    1: Function cannot run
*    3: Cut-off frequency too high => Equivalent to not filter
*    4: Cut-off frequency too low
*
* Parameters:
* - invector: Input/Output vector (non-filtered input vector/filtered input vector)
* - szVct: Size of the invector
* - scaleFactor: Scale factor
* - sampleRate: Sampling frequency (samples per second)
******************************************************************************/
int lpf_boxcar (gsl_vector **invector, int szVct, double scaleFactor, int sampleRate)
{
	string message = "";
	char valERROR[256];
	
	// Declare variables
	gsl_vector *invectorAux;
	gsl_vector *invectorAux1;
	double cutFreq;	//Frequency domain
	int boxLength;	//Time domain
	double boxSum = 0.0;

	// Define the LPF (frequency domain) and the box-car function (time domain)
	cutFreq = 2 * (1/(2*pi*scaleFactor));	//According to Jan, sinc(f1)=0.6 where f1=1/(2pi*scaleFactor)
						//sinc(0.5)~0.6 => f1~0.5
						//sinc(fc)=0, sinc(1)=0 => fc=1
						//fc=kf1 => fc~2f1
	boxLength =(int) ((1/cutFreq) * sampleRate);

	if (boxLength < 1)		boxLength = 1;
	if (boxLength >= szVct)	        return(4);

	// It is going to work with a longer vector to not have fake results for the last boxLength windows
	// (due to not having n samples to take into account)
	// It is not necessary to check the allocation because 'boxLength' is 1 at least
	invectorAux = gsl_vector_alloc(szVct+boxLength);
	invectorAux1 = gsl_vector_alloc(szVct+boxLength);
	for (int i=0;i<szVct;i++)
	{
		gsl_vector_set(invectorAux,i,gsl_vector_get(*invector,i));	// 'i' already correct in the for
	}
	double value = gsl_vector_get(*invector,szVct-1);
	for (int i=szVct;i<szVct+boxLength;i++)
	{
		if (value > 0) 		value = value-0.01*value;
		else if (value< 0)	value = value+0.01*value;
		gsl_vector_set(invectorAux,i,value);				// 'i' already correct in the for
	}

	// Apply the box-car window by shifting it along the (lengthened) input vector 
	for (int i=0;i<boxLength;i++)
	{
		boxSum = boxSum + gsl_vector_get(invectorAux,i);
	}
	gsl_vector_set(invectorAux1,0,boxSum/boxLength);			// 'i' already correct in the for

	if (invectorAux->size-boxLength-1+1 > invectorAux1->size-1)
	{
		sprintf(valERROR,"%d",__LINE__+8);
		string str(valERROR);
		message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
        str.clear();
		EP_PRINT_ERROR(message,EPFAIL);
	}
	for (int i=0;i<invectorAux->size-boxLength;i++)
	{
		boxSum = boxSum - gsl_vector_get(invectorAux,i) + gsl_vector_get(invectorAux,i+boxLength);
		gsl_vector_set(invectorAux1,i+1,boxSum/boxLength);		
	}

	if (invectorAux->size-boxLength+boxLength-2+1 > invectorAux1->size-1)
	{
		sprintf(valERROR,"%d",__LINE__+8);
		string str(valERROR);
		message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
        str.clear();
		EP_PRINT_ERROR(message,EPFAIL);
	}
	for (int i=0;i<boxLength-1;i++)
	{
		boxSum = boxSum - gsl_vector_get(invectorAux,invectorAux->size-boxLength+i);
		gsl_vector_set(invectorAux1,invectorAux->size-boxLength+i+1,boxSum/(boxLength-i-1));
	}

	for (int i=0;i<szVct;i++)
	{
		gsl_vector_set(*invector,i,gsl_vector_get(invectorAux1,i));	// 'i' already correct in the for
	}

	// Free allocated GSL vectors
	gsl_vector_free(invectorAux); invectorAux = 0;
	gsl_vector_free(invectorAux1); invectorAux1 = 0;
	
	if (boxLength == 1)	return(3);

    message.clear();
    
	return (EPOK);
}
/*xxxx end of SECTION 1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 2 ************************************************************
* differentiate function: This function applies the derivative method (x_i-x_(i-1)) to the input vector
*
* The derivative method provides more sensitivity to handle with pilot pulses.
* Moreover, little variations of the baseline will not affect.
*
* Parameters:
* - invector: Input/Ouput GSL vector (non-differentiated input vector/differentiated input vector)
* - szVct: Size of invector
******************************************************************************/
int differentiate (gsl_vector **invector,int szVct)
{
	for (int i=0; i<szVct-1; i++)
	{
		gsl_vector_set(*invector,i,gsl_vector_get(*invector,i+1)-gsl_vector_get(*invector,i));	
	}
	gsl_vector_set(*invector,szVct-1,gsl_vector_get(*invector,szVct-2));
	//gsl_vector_set(*invector,szVct-1,-999);

	return (EPOK);
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 ************************************************************
* findMeanSigma function: This function calculates the mean and the standard deviation of the input vector
* 
* Parameters:
* - invector: Input GSL vector
* - mean: Mean of the elements of 'invector'
* - sigma: Standard deviation of the elements of 'invector'
******************************************************************************/
int findMeanSigma (gsl_vector *invector, double *mean, double *sigma)
{
	// Declare variables
	int size = invector->size; // Size of the input vector
	// To calculate the mean
	bool veryBig = false;
    double suma = 0.0;

	// Mean
	for (int i=0;i<size;i++)
	{
        if (gsl_vector_get(invector,i)>1e10)
		{
			veryBig =true;
			break;
		}
	}
	
	if (veryBig == false)
	{
        *mean = 0.0;
        for (int i=0;i<size;i++)
        {
            *mean = *mean + gsl_vector_get(invector,i);
        }
        *mean = *mean/size;
		// Standard deviation
        for (int i=0;i<size;i++)
        {
            suma = suma + pow(gsl_vector_get(invector,i)-*mean,2.0);
        }
        *sigma = sqrt(suma/(size-1));
	}
	else	// To avoid an inf in IO or a NAN in JUPITER
	{
		*mean = 1e10;
		*sigma = 1e10;
	}

	return EPOK;
}
/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
* medianKappaClipping function: This function calculates a threshold in the first derivative of the record by using
*                               a Kappa-clipping method (replacing points beyond mean+-Kappa*sigma with the median).
*
* First, mean and sigma are calculated and invector values out of (mean+Kappa*sigma,mean-Kappa*sigma) are replaced
* with the median (it is trying to look for the baseline). And this process is iteratively repeated until there are
* no points beyond mean+-Kappa *sigma. Finally, the threshold is calculated as mean+nSigmas*sigma ('+' is used because
* if there are pulses in the input invector they are always positive).
*
* - Declare variables
* - Calculate the median
* - Iterate until there are no points out of the maximum excursion (kappa*sigma)
* - Establish the threshold as mean+nSigmas*sigma
*
* Parameters:
* - invector: First derivative of the (filtered) record
* - Kappa: To establish the range around of the mean
* - stopCriteria: It is given in %
* - nSigmas: Times sigma to calculate threshold (mean+nSigmas*sigma)
* - boxLPF: Length of the low-pass filtering box-car
* - threshold: Calculated threshold
******************************************************************************/
int medianKappaClipping (gsl_vector *invector, double kappa, double stopCriteria, double nSigmas, int boxLPF, double *threshold)
{
	string message = "";
	char valERROR[256];

	// Declare variables
	int size = invector->size; // Size of the input vector
	double mean1, sg1;
	double mean2, sg2;
	gsl_vector_view temp;
	// Variables to remove input vector elements higher than the maximum excursion (kappa*sg)
	int i;							// To go through the elements of a vector
	int cnt;						// Number of points inside the excursion (mean+-excursion)
	// It is not necessary to check the allocation because 'invector' size must already be > 0
	gsl_vector *invectorNew = gsl_vector_alloc(size);	// Auxiliary vector
        
	// To calculate the median
	//double data[size];					// Intermediate value to use 'gsl_stats_median_from_sorted_data'
	double median;

	// Median
    gsl_vector_memcpy(invectorNew,invector);
    gsl_sort_vector(invectorNew);
    if (size%2 == 0)	//Even
    {
        median = (gsl_vector_get(invectorNew,(int) (size/2)-1)+gsl_vector_get(invectorNew,(int) (size/2)))/2;
    }
    else                    //Odd
    {
        median = gsl_vector_get(invectorNew,(int) (size/2));
    }

	gsl_vector_memcpy(invectorNew,invector);

	// Iterate until no points out of the maximum excursion (kappa*sigma)
	do
	{
		if ((size-boxLPF-1 < 1) || (size-boxLPF-1 >invectorNew->size))
		{
			sprintf(valERROR,"%d",__LINE__+5);
			string str(valERROR);
			message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
            str.clear();
			EP_PRINT_ERROR(message,EPFAIL);
		}
		temp = gsl_vector_subvector(invectorNew,0,size-boxLPF-1);
		if (findMeanSigma (&temp.vector, &mean1, &sg1))
		{
			message = "Cannot run findMeanSigma routine for kappa-sigma iteration";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		i = 0;
		cnt = 0;

		while (i<invectorNew->size)
		{
			if ((gsl_vector_get(invectorNew,i) >= mean1 + kappa*sg1) || (gsl_vector_get(invectorNew,i) <= mean1 - kappa*sg1))
			// HARDPOINT!!!!!!!!!!!!!!!!!!! (kappa)
			{
				if ((i < 0) || (i >(invectorNew)->size-1))
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "Setting with <= 0 size in line " + str + " (" + __FILE__ + ")";
                    str.clear();
					EP_PRINT_ERROR(message,EPFAIL);
				}
				gsl_vector_set(invectorNew,i,median);
				cnt++;
			}
			i++;
		}

		if (cnt != 0)
		// Some points of the invector have been replaced with the median
		{
			temp = gsl_vector_subvector(invectorNew,0,size-boxLPF-1);
			if (findMeanSigma (&temp.vector, &mean2, &sg2))
			{
				message = "Cannot run findMeanSigma routine for kappa-sigma iteration after replacement with the median";
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
		}
		else
		// No points of the invector have been replaced with the median
		{
			mean2 =mean1;
			sg2 = sg1;
		}

	} while (fabs((mean2-mean1)/mean1)>(stopCriteria/100.0));	// HARDPOINT!!!!!!!!!!!!!!!!!!! (stopCriteria)

	// Establish the threshold as mean+nSigmas*sigma
	*threshold = mean2+nSigmas*sg2;	// HARDPOINT!!!!!!!!!!!!!!!!!!! (nSigmas)
	
	gsl_vector_free(invectorNew); invectorNew= 0;
    
    message.clear();

	return EPOK;
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* getB function: This function calculates the sum, B, of lb digitized data samples of a pulse-free interval immediately
*       before the each pulse. If the pulse-free interval before the current pulse is lower than lb, B is calculated with the available
*       number of samples. If there is not a pulse-free interval before the pulse, it is looked for it after the current pulse.
*       The number of samples of the pulse-free interval used to calculate B is stored in the lb vector.
*
* First of all, the auxiliary variable 'Baux' is initialized to -999 and all the elements of the lb vector are equal to the 'LbT' input variable in samples.
* Then, the code is divided into 2 *if* statements:
* 
*  - When the current pulse is the first pulse into the record:
*  	- tstart>=lb => Sum lb samples
*  	- 0<tstart<lb => Sum the available number of samples (although the available number of samples was lower than lb)
*  	- tstart=0 => If there is not a pulse-free interval before the pulse, it is looked for it after the current pulse
*  - When the current pulse is not the first pulse into the record:
*  	- tstart_(i)-tend_(i-1)>=lb => Sum lb samples
*  	- 0<tstart_(i)-tend_(i-1)<lb => Sum the available number of samples (although the available number of samples was lower than lb)
*  	- If there is not a pulse-free interval before the pulse, it is looked for it after the current pulse
*
* If 'Baux' is still -999, a pulse-free interval can not be found to apply the running sum filter. This has to be taken into account,
* out of the function, to try to get a usable B.
* 
* Parameters:
* - vectorin: Input record
* - tstart: Starting time of the pulses into the record
* - nPulses: Number of pulses into the record
* - lb: Vector containing the baseline averaging length used for each pulse
* - sizePulse:  Size of the pulse in samples
* - B: In general, sum of the Lb digitized data samples of a pulse-free interval immediately before the current pulse
* - rmsB: In general, rms of the baseline related to a pulse-free interval immediately before the current pulse
****************************************/
int getB(gsl_vector *vectorin, gsl_vector *tstart, int nPulses, gsl_vector **lb, int sizepulse, gsl_vector **B, gsl_vector **rmsB)
{
	string message = "";
	char valERROR[256];

	// Declare variables
	// It is not necessary to check the allocation because 'tstart' size must already be > 0
	*B = gsl_vector_alloc(tstart->size);
	gsl_vector_set_all(*B,-999);
    *rmsB = gsl_vector_alloc(tstart->size);
	gsl_vector_set_all(*rmsB,-999);
	double Baux = -999;
	double tendprev;

	// Auxiliary variables
	gsl_vector *input;
	gsl_vector_view temp;		// In order to handle with gsl_vector_view (subvectors)
	
	double mean, sigma;

    if (nPulses != 0)
    {
        if ((nPulses-1 >(tstart)->size-1))
        {
            sprintf(valERROR,"%d",__LINE__+8);
            string str(valERROR);
            message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
            str.clear();
            EP_PRINT_ERROR(message,EPFAIL);
        }
        
        for (int i=0;i<nPulses;i++)
        {
            gsl_vector_set(tstart,i,(int)(gsl_vector_get(tstart,i)));	
        }
        
        for (int i=0;i<nPulses;i++)
        {
            if (i == 0)     // First pulse into a record		
                //  Current pulse
                //      //\\          /\     /\       /\
                // (1) //  \\  (2)   /  \(3)/  \ (4) /  \    (5)
                // ----      --------    ---    -----    ----------
            {
                if (gsl_vector_get(tstart,0)>=gsl_vector_get(*lb,0))
                    // tstart>=lb => Sum lb samples
                    // length_(1)>=lb
                {
                    if ((input = gsl_vector_alloc(gsl_vector_get(*lb,0))) == 0)
                    {
                        sprintf(valERROR,"%d",__LINE__-2);
                        string str(valERROR);
                        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                    }
                    
                    if ((gsl_vector_get(tstart,0)-gsl_vector_get(*lb,0) < 0) || (gsl_vector_get(tstart,0)-gsl_vector_get(*lb,0) > vectorin->size-2)
                        || (gsl_vector_get(*lb,0) < 1) || (gsl_vector_get(*lb,0) >vectorin->size-(gsl_vector_get(tstart,0)-gsl_vector_get(*lb,0))))
                    {
                        sprintf(valERROR,"%d",__LINE__+6);
                        string str(valERROR);
                        message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL);
                    }
                    temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,0)-gsl_vector_get(*lb,0),gsl_vector_get(*lb,0));
                    if (gsl_vector_memcpy(input, &temp.vector) != 0)
                    {
                        sprintf(valERROR,"%d",__LINE__-2);
                        string str(valERROR);
                        message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL);
                    }
                    
                    // Sum all the elements of 'input'
                    if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                    {
                        message = "Cannot run gsl_vector_Sumsubvector routine when tstart>=lb";
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                    }
                    //cout<<"Baux: "<<Baux<<endl;
                    if (findMeanSigma (input, &mean, &sigma))
                    {
                        message = "Cannot run findMeanSigma in getB";
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                    }
                    gsl_vector_free(input); input = 0;
                }
                else if ((gsl_vector_get(tstart,0)<gsl_vector_get(*lb,0)) && (gsl_vector_get(tstart,0)>1))
                    // 0<tstart<lb => Sum the available number of samples (although the available number of samples was lower than lb)
                    // 0<length_(1)<lb
                {
                    if ((input = gsl_vector_alloc(gsl_vector_get(tstart,0))) == 0)
                    {
                        sprintf(valERROR,"%d",__LINE__-2);
                        string str(valERROR);
                        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                    }
                    
                    if ((gsl_vector_get(tstart,0) < 1) || (gsl_vector_get(tstart,0) >vectorin->size))
                    {
                        sprintf(valERROR,"%d",__LINE__+5);
                        string str(valERROR);
                        message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL);
                    }
                    temp = gsl_vector_subvector(vectorin,0,gsl_vector_get(tstart,0));
                    if (gsl_vector_memcpy(input, &temp.vector) != 0)
                    {
                        sprintf(valERROR,"%d",__LINE__-2);
                        string str(valERROR);
                        message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL);
                    }
                    
                    // Sum all the elements of 'input'
                    if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                    {
                        message = "Cannot run gsl_vector_Sumsubvector routine when tstart<lb";
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                    }
                    gsl_vector_set(*lb,0, gsl_vector_get(tstart,0));
                    if (findMeanSigma (input, &mean, &sigma))
                    {
                        message = "Cannot run findMeanSigma in getB";
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                    }
                    gsl_vector_free(input); input = 0;
                }
                else	// If there is not a pulse-free interval before the pulse, it is looked for it after the current pulse
                {
                    for (int j=1;j<nPulses;j++)
                        // (2),(3),(4) and (5) are analyzed
                        // When 0<length_(j)<lb => 'break' => Out of the 'for' loop
                    {
                        tendprev = gsl_vector_get(tstart,j-1)+sizepulse-1;
                        if (tendprev >= vectorin->size)
                        {
                            tendprev = vectorin->size-1;
                        }
                        
                        if (gsl_vector_get(tstart,j)-tendprev > 0)
                        {
                            if (gsl_vector_get(tstart,j)-tendprev >= gsl_vector_get(*lb,0))
                                // length_(j)>=lb (j/=nPulses)
                            {
                                if ((input = gsl_vector_alloc(gsl_vector_get(*lb,0))) == 0)
                                {
                                    sprintf(valERROR,"%d",__LINE__-2);
                                    string str(valERROR);
                                    message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                }
                                
                                if ((gsl_vector_get(tstart,j)-gsl_vector_get(*lb,0) < 0) || (gsl_vector_get(tstart,j)-gsl_vector_get(*lb,0) > vectorin->size-2)
                                    || (gsl_vector_get(*lb,0) < 1) || (gsl_vector_get(*lb,0) >vectorin->size-(gsl_vector_get(tstart,j)-gsl_vector_get(*lb,0))))
                                {
                                    sprintf(valERROR,"%d",__LINE__+5);
                                    string str(valERROR);
                                    message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,j)-gsl_vector_get(*lb,0),gsl_vector_get(*lb,0));
                                if (gsl_vector_memcpy(input, &temp.vector) != 0)
                                {
                                    sprintf(valERROR,"%d",__LINE__-2);
                                    string str(valERROR);
                                    message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                
                                // Sum all the elements of input
                                if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                                {
                                    message = "Cannot run gsl_vector_Sumsubvector routine when no pulse free interval before the pulse";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                if (findMeanSigma (input, &mean, &sigma))
                                {
                                    message = "Cannot run findMeanSigma in getB";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                gsl_vector_free(input); input = 0;
                                
                                break;
                            }
                            else if ((gsl_vector_get(tstart,j)-tendprev < gsl_vector_get(*lb,0)) && (gsl_vector_get(tstart,j)-tendprev > 1))
                                // 0<length_(j)<lb (j/=nPulses)
                            {
                                if ((input = gsl_vector_alloc(gsl_vector_get(tstart,j)-tendprev-1)) == 0)
                                {
                                    sprintf(valERROR,"%d",__LINE__-2);
                                    string str(valERROR);
                                    message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                }
                                
                                if ((tendprev+1 < 0) || (tendprev+1 > vectorin->size-2)
                                    || (gsl_vector_get(tstart,j)-tendprev-1 < 1) || (gsl_vector_get(tstart,j)-tendprev-1 >vectorin->size-(tendprev+1)))
                                {
                                    sprintf(valERROR,"%d",__LINE__+5);
                                    string str(valERROR);
                                    message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,j)-tendprev-1);
                                if (gsl_vector_memcpy(input, &temp.vector) != 0)
                                {
                                    sprintf(valERROR,"%d",__LINE__-2);
                                    string str(valERROR);
                                    message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                gsl_vector_set(*lb,0,gsl_vector_get(tstart,j)-tendprev-1);
                                
                                // Sum all the elements of input
                                if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                                {
                                    message = "Cannot run gsl_vector_Sumsubvector routine when no pulse free interval before the pulse";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                if (findMeanSigma (input, &mean, &sigma))
                                {
                                    message = "Cannot run findMeanSigma in getB";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                gsl_vector_free(input); input = 0;
                                
                                break;
                            }
                            else
                            {
                                for (int j=0;j<nPulses;j++)	// From the current pulse
                                    // (2),(3),(4) and (5) are analyzed
                                    // When 0<length_(j)<lb => 'break' => Out of the 'for' loop
                                {
                                    tendprev = gsl_vector_get(tstart,j)+sizepulse-1;
                                    if (tendprev >= vectorin->size)
                                    {
                                        tendprev = vectorin->size-1;
                                    }
                                    if ((j < nPulses-1) && (gsl_vector_get(tstart,j+1)-tendprev > 0)) // Not last pulse into a row=event
                                    {
                                        if (gsl_vector_get(tstart,j+1)-tendprev >= gsl_vector_get(*lb,0))
                                        {
                                            if ((input = gsl_vector_alloc(gsl_vector_get(*lb,0))) == 0)
                                            {
                                                sprintf(valERROR,"%d",__LINE__-2);
                                                string str(valERROR);
                                                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                            }
                                            
                                            if ((gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,0) < 0) || (gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,0) > vectorin->size-2)
                                                || (gsl_vector_get(*lb,0) < 1) || (gsl_vector_get(*lb,0) > vectorin->size-(gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,0))))
                                            {
                                                sprintf(valERROR,"%d",__LINE__+5);
                                                string str(valERROR);
                                                message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL);
                                            }
                                            temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,0),gsl_vector_get(*lb,0));
                                            if (gsl_vector_memcpy(input, &temp.vector) != 0)
                                            {
                                                sprintf(valERROR,"%d",__LINE__-2);
                                                string str(valERROR);
                                                message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL);
                                            }
                                            
                                            // Sum all the elements of input
                                            if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                                            {
                                                message = "Cannot run gsl_vector_Sumsubvector routine when no first pulse in row & tstart-tendprev >= lb";
                                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                            }
                                            if (findMeanSigma (input, &mean, &sigma))
                                            {
                                                message = "Cannot run findMeanSigma in getB";
                                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                            }
                                            gsl_vector_free(input); input = 0;
                                            
                                            break;
                                        }
                                        else if ((gsl_vector_get(tstart,j+1)-tendprev < gsl_vector_get(*lb,0)) && ((gsl_vector_get(tstart,j+1)-tendprev>1)))
                                        {
                                            if ((input = gsl_vector_alloc(gsl_vector_get(tstart,j+1)-tendprev-1)) == 0)
                                            {
                                                sprintf(valERROR,"%d",__LINE__-2);
                                                string str(valERROR);
                                                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                            }
                                            
                                            if ((tendprev+1 < 0) || (tendprev+1 > vectorin->size-2)
                                                || (gsl_vector_get(tstart,j+1)-tendprev-1 < 1) || (gsl_vector_get(tstart,j+1)-tendprev-1 > vectorin->size-(tendprev+1)))
                                            {
                                                sprintf(valERROR,"%d",__LINE__+5);
                                                string str(valERROR);
                                                message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL);
                                            }
                                            temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,j+1)-tendprev-1);
                                            if (gsl_vector_memcpy(input, &temp.vector) != 0)
                                            {
                                                sprintf(valERROR,"%d",__LINE__-2);
                                                string str(valERROR);
                                                message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL);
                                            }
                                            
                                            // Sum all the elements of input
                                            if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                                            {
                                                message = "Cannot run gsl_vector_Sumsubvector routine when no first pulse in row & tstart-tendprev < lb";
                                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                            }
                                            if (/*(0 < 0) || */(0 >(*lb)->size-1))
                                            {
                                                sprintf(valERROR,"%d",__LINE__+5);
                                                string str(valERROR);
                                                message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL);
                                            }
                                            gsl_vector_set(*lb,0,gsl_vector_get(tstart,j+1)-tendprev-1);
                                            if (findMeanSigma (input, &mean, &sigma))
                                            {
                                                message = "Cannot run findMeanSigma in getB";
                                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                            }
                                            gsl_vector_free(input); input = 0;
                                            
                                            break;
                                        }
                                    }
                                    
                                    if (j == nPulses-1)	// Last pulse into the record
                                        // (5) is analyzed
                                    {
                                        if (vectorin->size-tendprev >= gsl_vector_get(*lb,0))
                                        {
                                            if ((input = gsl_vector_alloc(gsl_vector_get(*lb,0))) == 0)
                                            {
                                                sprintf(valERROR,"%d",__LINE__-2);
                                                string str(valERROR);
                                                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                            }
                                            
                                            if ((tendprev < 0) || (tendprev > vectorin->size-2)
                                                || (gsl_vector_get(*lb,0) < 1) || (gsl_vector_get(*lb,0) > vectorin->size-tendprev))
                                            {
                                                sprintf(valERROR,"%d",__LINE__+5);
                                                string str(valERROR);
                                                message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL);
                                            }
                                            temp = gsl_vector_subvector(vectorin,tendprev,gsl_vector_get(*lb,0));
                                            if (gsl_vector_memcpy(input, &temp.vector) != 0)
                                            {
                                                sprintf(valERROR,"%d",__LINE__-2);
                                                string str(valERROR);
                                                message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL);
                                            }
                                            
                                            // Sum all the elements of input
                                            if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                                            {
                                                message = "Cannot run gsl_vector_Sumsubvector routine when no pulse-free interval before pulse & last pulse in row";
                                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                            }
                                            if (findMeanSigma (input, &mean, &sigma))
                                            {
                                                message = "Cannot run findMeanSigma in getB";
                                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                            }
                                            gsl_vector_free(input); input = 0;
                                            
                                            break;
                                        }
                                        else if ((vectorin->size-tendprev < gsl_vector_get(*lb,0)) && (vectorin->size-tendprev > 1))
                                        {
                                            if ((input = gsl_vector_alloc(vectorin->size-tendprev-1)) == 0)
                                            {
                                                sprintf(valERROR,"%d",__LINE__-2);
                                                string str(valERROR);
                                                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                            }
                                            
                                            if ((tendprev+1 < 0) || (tendprev+1 > vectorin->size-2)
                                                || (vectorin->size-tendprev-1 < 1))
                                            {
                                                sprintf(valERROR,"%d",__LINE__+5);
                                                string str(valERROR);
                                                message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL);
                                            }
                                            temp = gsl_vector_subvector(vectorin,tendprev+1,vectorin->size-tendprev-1);
                                            if (gsl_vector_memcpy(input, &temp.vector) != 0)
                                            {
                                                sprintf(valERROR,"%d",__LINE__-2);
                                                string str(valERROR);
                                                message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL);
                                            }
                                            if (/*(i < 0) || */(0 >(*lb)->size-1))
                                            {
                                                sprintf(valERROR,"%d",__LINE__+5);
                                                string str(valERROR);
                                                message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                                                str.clear();
                                                EP_PRINT_ERROR(message,EPFAIL);
                                            }
                                            gsl_vector_set(*lb,0,vectorin->size-tendprev-1);
                                            
                                            // Sum all the elements of input
                                            if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                                            {
                                                message = "Cannot run gsl_vector_Sumsubvector routine when no pulse-free interval before pulse & last pulse in row";
                                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                            }
                                            if (findMeanSigma (input, &mean, &sigma))
                                            {
                                                message = "Cannot run findMeanSigma in getB";
                                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                            }
                                            gsl_vector_free(input); input = 0;
                                            
                                            break;
                                        }
                                    }
                                }//for 0 to N
                            }
                        }
                    }
                }
            }// end of the first pulse scope
            else    // Not first pulse into a record
                //             Current pulse
                //      /\          //\\     /\       /\
                // (1) /  \  (2)   //  \\(3)/  \ (4) /  \    (5)
                // ----    --------      ---    -----    ----------
            {	
                
                tendprev = gsl_vector_get(tstart,i-1)+sizepulse-1;
                if (gsl_vector_get(tstart,i)-tendprev >= gsl_vector_get(*lb,i))
                    // tstart_(i)-tend_(i-1)>=lb => Sum lb samples
                    // length_(2)>=lb
                {
                    if ((input = gsl_vector_alloc(gsl_vector_get(*lb,i))) == 0)
                    {
                        sprintf(valERROR,"%d",__LINE__-2);
                        string str(valERROR);
                        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                    }
                    
                    if ((gsl_vector_get(tstart,i)-gsl_vector_get(*lb,i) < 0) || (gsl_vector_get(tstart,i)-gsl_vector_get(*lb,i) > vectorin->size-2)
                        || (gsl_vector_get(*lb,i) < 1) || (gsl_vector_get(*lb,i) > vectorin->size-(gsl_vector_get(tstart,i)-gsl_vector_get(*lb,i))))
                    {
                        sprintf(valERROR,"%d",__LINE__+5);
                        string str(valERROR);
                        message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL);
                    }
                    temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,i)-gsl_vector_get(*lb,i),gsl_vector_get(*lb,i));
                    if (gsl_vector_memcpy(input, &temp.vector) != 0)
                    {
                        sprintf(valERROR,"%d",__LINE__-2);
                        string str(valERROR);
                        message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL);
                    }
                    
                    // Sum all the elements of input
                    if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                    {
                        message = "Cannot run gsl_vector_Sumsubvector routine length_(2)>=lb";
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                    }
                    if (findMeanSigma (input, &mean, &sigma))
                    {
                        message = "Cannot run findMeanSigma in getB";
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                    }
                    gsl_vector_free(input); input = 0;
                }
                else if ((gsl_vector_get(tstart,i)-tendprev<gsl_vector_get(*lb,i)) && (gsl_vector_get(tstart,i)-tendprev>1))
                    // 0<tstart_(i)-tend_(i-1)<lb => Sum the available number of samples (although the available number of samples was lower than lb)
                    // 0<length_(2)<lb
                {
                    if ((input = gsl_vector_alloc(gsl_vector_get(tstart,i)-tendprev-1)) == 0)
                    {
                        sprintf(valERROR,"%d",__LINE__-2);
                        string str(valERROR);
                        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                    }
                    
                    if ((tendprev+1 < 0) || (tendprev+1 > vectorin->size-2)
                        || (gsl_vector_get(tstart,i)-tendprev-1 < 1) || (gsl_vector_get(tstart,i)-tendprev-1 > vectorin->size-(tendprev+1)))
                    {
                        sprintf(valERROR,"%d",__LINE__+5);
                        string str(valERROR);
                        message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL);
                    }
                    temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,i)-tendprev-1);
                    if (gsl_vector_memcpy(input, &temp.vector) != 0)
                    {
                        sprintf(valERROR,"%d",__LINE__-2);
                        string str(valERROR);
                        message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL);
                    }
                    
                    // Sum all the elements of input
                    if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                    {
                        message = "Cannot run gsl_vector_Sumsubvector routine when 0<length_(2)<lb";
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                    }
                    if ((i < 0) || (i >(*lb)->size-1))
                    {
                        sprintf(valERROR,"%d",__LINE__+5);
                        string str(valERROR);
                        message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL);
                    }
                    gsl_vector_set(*lb,i,gsl_vector_get(tstart,i)-tendprev-1);
                    if (findMeanSigma (input, &mean, &sigma))
                    {
                        message = "Cannot run findMeanSigma in getB";
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                    }
                    gsl_vector_free(input); input = 0;
                }
                else	// If there is not a pulse-free interval before the pulse, it is looked for it after the current pulse
                {
                    for (int j=i;j<nPulses;j++)	// From the current pulse
                        // (3),(4) and (5) are analyzed
                        // When 0<length_(j)<lb => 'break' => Out of the 'for' loop
                    {
                        tendprev = gsl_vector_get(tstart,j)+sizepulse-1;
                        if (tendprev >= vectorin->size)
                        {
                            tendprev = vectorin->size-1;
                        }
                        if ((j < nPulses-1) && (gsl_vector_get(tstart,j+1)-tendprev > 0)) // Not last pulse into a row=event
                        {
                            if (gsl_vector_get(tstart,j+1)-tendprev >= gsl_vector_get(*lb,i))
                            {
                                if ((input = gsl_vector_alloc(gsl_vector_get(*lb,i))) == 0)
                                {
                                    sprintf(valERROR,"%d",__LINE__-2);
                                    string str(valERROR);
                                    message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                }
                                
                                if ((gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i) < 0) || (gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i) > vectorin->size-2)
                                    || (gsl_vector_get(*lb,i) < 1) || (gsl_vector_get(*lb,i) > vectorin->size-(gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i))))
                                {
                                    sprintf(valERROR,"%d",__LINE__+5);
                                    string str(valERROR);
                                    message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i),gsl_vector_get(*lb,i));
                                if (gsl_vector_memcpy(input, &temp.vector) != 0)
                                {
                                    sprintf(valERROR,"%d",__LINE__-2);
                                    string str(valERROR);
                                    message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                
                                // Sum all the elements of input
                                if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                                {
                                    message = "Cannot run gsl_vector_Sumsubvector routine when no first pulse in row & tstart-tendprev >= lb";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                if (findMeanSigma (input, &mean, &sigma))
                                {
                                    message = "Cannot run findMeanSigma in getB";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                gsl_vector_free(input); input = 0;
                                
                                break;
                            }
                            else if ((gsl_vector_get(tstart,j+1)-tendprev < gsl_vector_get(*lb,i)) && ((gsl_vector_get(tstart,j+1)-tendprev>1)))
                            {
                                if ((input = gsl_vector_alloc(gsl_vector_get(tstart,j+1)-tendprev-1)) == 0)
                                {
                                    sprintf(valERROR,"%d",__LINE__-2);
                                    string str(valERROR);
                                    message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                }
                                
                                if ((tendprev+1 < 0) || (tendprev+1 > vectorin->size-2)
                                    || (gsl_vector_get(tstart,j+1)-tendprev-1 < 1) || (gsl_vector_get(tstart,j+1)-tendprev-1 > vectorin->size-(tendprev+1)))
                                {
                                    sprintf(valERROR,"%d",__LINE__+5);
                                    string str(valERROR);
                                    message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,j+1)-tendprev-1);
                                if (gsl_vector_memcpy(input, &temp.vector) != 0)
                                {
                                    sprintf(valERROR,"%d",__LINE__-2);
                                    string str(valERROR);
                                    message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                
                                // Sum all the elements of input
                                if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                                {
                                    message = "Cannot run gsl_vector_Sumsubvector routine when no first pulse in row & tstart-tendprev < lb";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                if ((i < 0) || (i >(*lb)->size-1))
                                {
                                    sprintf(valERROR,"%d",__LINE__+5);
                                    string str(valERROR);
                                    message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                gsl_vector_set(*lb,i,gsl_vector_get(tstart,j+1)-tendprev-1);
                                if (findMeanSigma (input, &mean, &sigma))
                                {
                                    message = "Cannot run findMeanSigma in getB";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                gsl_vector_free(input); input = 0;
                                
                                break;
                            }
                        }
                        
                        if (j == nPulses-1)	// Last pulse into the record
                            // (5) is analyzed
                        {
                            if (vectorin->size-tendprev >= gsl_vector_get(*lb,i))
                            {
                                if ((input = gsl_vector_alloc(gsl_vector_get(*lb,i))) == 0)
                                {
                                    sprintf(valERROR,"%d",__LINE__-2);
                                    string str(valERROR);
                                    message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                }
                                
                                //if ((tendprev < 0) || (tendprev > vectorin->size-2)                                                                   || (gsl_vector_get(*lb,i) < 1) || (gsl_vector_get(*lb,i) > vectorin->size-tendprev))
                                if ((tendprev < 0) || (gsl_vector_get(*lb,i) < 1) || (gsl_vector_get(*lb,i) > vectorin->size-tendprev))
                                {
                                    sprintf(valERROR,"%d",__LINE__+5);
                                    string str(valERROR);
                                    message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                temp = gsl_vector_subvector(vectorin,tendprev,gsl_vector_get(*lb,i));
                                if (gsl_vector_memcpy(input, &temp.vector) != 0)
                                {
                                    sprintf(valERROR,"%d",__LINE__-2);
                                    string str(valERROR);
                                    message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                
                                // Sum all the elements of input
                                if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                                {
                                    message = "Cannot run gsl_vector_Sumsubvector routine when no pulse-free interval before pulse & last pulse in row";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                if (findMeanSigma (input, &mean, &sigma))
                                {
                                    message = "Cannot run findMeanSigma in getB";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                gsl_vector_free(input); input = 0;
                                
                                break;
                            }
                            else if ((vectorin->size-tendprev < gsl_vector_get(*lb,i)) && (vectorin->size-tendprev > 1))
                            {
                                if ((input = gsl_vector_alloc(vectorin->size-tendprev-1)) == 0)
                                {
                                    sprintf(valERROR,"%d",__LINE__-2);
                                    string str(valERROR);
                                    message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                }
                                
                                //if ((tendprev+1 < 0) || (tendprev+1 > vectorin->size-2)                                                                    || (vectorin->size-tendprev-1 < 1))
                                if ((tendprev+1 < 0) || (vectorin->size-tendprev-1 < 1))
                                {
                                    sprintf(valERROR,"%d",__LINE__+5);
                                    string str(valERROR);
                                    message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                temp = gsl_vector_subvector(vectorin,tendprev+1,vectorin->size-tendprev-1);
                                if (gsl_vector_memcpy(input, &temp.vector) != 0)
                                {
                                    sprintf(valERROR,"%d",__LINE__-2);
                                    string str(valERROR);
                                    message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                if ((i < 0) || (i >(*lb)->size-1))
                                {
                                    sprintf(valERROR,"%d",__LINE__+5);
                                    string str(valERROR);
                                    message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                gsl_vector_set(*lb,i,vectorin->size-tendprev-1);
                                
                                // Sum all the elements of input
                                if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
                                {
                                    message = "Cannot run gsl_vector_Sumsubvector routine when no pulse-free interval before pulse & last pulse in row";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                if (findMeanSigma (input, &mean, &sigma))
                                {
                                    message = "Cannot run findMeanSigma in getB";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                gsl_vector_free(input); input = 0;
                                
                                break;
                            }
                            else if ((vectorin->size-tendprev < gsl_vector_get(*lb,i)) && (vectorin->size-tendprev <= 1))
                            {
                                if ((i < 0) || (i >(*lb)->size-1))
                                {
                                    sprintf(valERROR,"%d",__LINE__+5);
                                    string str(valERROR);
                                    message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                gsl_vector_set(*lb,i,gsl_vector_get(*lb,i-1));
                                
                                break;
                            }
                        }
                    }
                }
            }
            
            if ((i < 0) || (i >(*B)->size-1))
            {
                sprintf(valERROR,"%d",__LINE__+5);
                string str(valERROR);
                message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                str.clear();
                EP_PRINT_ERROR(message,EPFAIL);
            }
            gsl_vector_set(*B,i,Baux);
            gsl_vector_set(*rmsB,i,sigma);
        }//for
    }
    
    message.clear();
        
	return(EPOK);
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
* getPulseHeight function: This function estimates the pulse height of a pulse by using a running sum filter.
*                 It extracts from the record, vectorin, the pulse whose pulse height is going to be estimated
*                 by using RS_filter.
*
* - Declare variables
* - Extracting from the record the pulse whose pulse height is going to be estimated
* - Apply the running sum filter
*
* Parameters:
* - vectorin: Not filtered record
* - tstart: Starting time of the pulse whose pulse height is going to be estimated
* - tstartnext: Starting time of the next pulse whose pulse height is going to be estimated
* - lastPulse: It is 1 if the pulse is the last one into the record or the only one
* - lrs: Running sum length (equal to the 'Lrs' global_variable)
* - lb: Baseline averaging length used for the pulse whose pulse height is going to be estimated
* - B: In general, sum of the Lb digitized data samples of a pulse-free interval immediately before
*      the current pulse
* - sizepulse: Size of the pulse in samples ('PulseLength')
* - pulseheight: Estimated pulse height of the pulse
****************************************/
int getPulseHeight(gsl_vector *vectorin, double tstart, double tstartnext, int lastPulse, double lrs, double lb, double B, int sizepulse, double *pulseheight)
{
	string message = "";
	char valERROR[256];

	// Declare variables
	long tend;	// Ending time of the pulse
	double ph;	// Pulse height
	
	// Auxiliary variables
	gsl_vector *input;	// Segment of 'vectorin' where is the pulse
	gsl_vector_view temp;	// In order to handle with gsl_vector_view (subvectors)

	// Extracting from the record the pulse whose pulse height is going to be estimated
	if (lastPulse == 1)
	{
		tend = tstart+sizepulse-1;
		if (tend >= vectorin->size)	// Truncated pulses at the end of the record
		{
			tend = vectorin->size-1;
		}
	}
	else
	{
		tend = tstartnext-1;
	}

	if (vectorin->size-tstart-(vectorin->size-tend) <= 0)
	{
		*pulseheight = -999;
	}
	else
	{
		// It is not necessary to check allocation bacuse 'vectorin->size-tstart-(vectorin->size-tend)' has been checked previously
		input = gsl_vector_alloc(vectorin->size-tstart-(vectorin->size-tend));	// Only the pulse
		temp = gsl_vector_subvector(vectorin,tstart,input->size);
		if (gsl_vector_memcpy(input, &temp.vector) != 0)
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);
			message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
            str.clear();
			EP_PRINT_ERROR(message,EPFAIL);
		}

		// Apply the running sum filter
		if (RS_filter (input, lrs, lb, B, &ph))
		{
			message = "Cannot run RS_filter routine when size-tstart>size-tend";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		*pulseheight = ph;

		gsl_vector_free(input); input = 0;
	}

	message.clear();
    
	return(EPOK);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* RS_filter function: This function uses the running sum filter to find the pulse height.
*                     It always works in time domain.
*
* A running sum filter, RS, is the sum of lrs digitized data samples. It is continuously updated upon the arrival of
* new data point. Simultaneously a baseline filter, B, is the sum of lb digitized data samples without pulses.
*
* The algorithm looks for the time when RS/lrs reaches its maximum. At that time RS is stored, RS_max, and the baseline
* is scaled with lrs, Bp (Bp=B*lrs/lb). Then, the pulse height related to the pulse pseudoenergy is given by:
*
*                     RS_max-Bp
*    Pulse height = -------------
*                        lrs
*
* Parameters:
* - vector: Not filtered pulse (extracted from the record in 'getPulseHeight')
* - lrs: Running sum length (samples)
* - lb: Baseline averaging length (samples)
* - B: In general, sum of the lb digitized data samples of a pulse-free interval immediately before
*      the current pulse
* - pulseheight: Pulse height of the pulse
*****************************************/
int RS_filter (gsl_vector *vector, double lrs, double lb, double B, double *pulseheight)
{
	string message = "";

	// Declare variables
	double Rs;		// Sum of lrs digitized data samples
	double Rs_max = -1e20;	// Rs is continuously updated upon the arrival of a new data point if the new value is
	                        // higher than the old one
	double Bp;

	if (vector->size<lrs)
	{
		if (gsl_vector_Sumsubvector(vector,0,vector->size,&Rs))
		{
			message = "Cannot run gsl_vector_Sumsubvector routine when size<lrs";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		if (Rs > Rs_max)	Rs_max = Rs;
	}
	else
	{
		for (int i=0;i<vector->size;i++)
		{
			if (i+lrs > vector->size)	break;
			if (gsl_vector_Sumsubvector(vector,i,lrs,&Rs))
			{
				message = "Cannot run gsl_vector_Sumsubvector routine when size>lrs in iteration i=" + i;
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			if (Rs > Rs_max)	{Rs_max = Rs;}
		}
	}

	Bp = B*lrs/lb;
	*pulseheight = (Rs_max-Bp)/lrs;
    
    message.clear();

	return EPOK;
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
* find_model_energies function: This function uses 'energy' in order to choose the proper pulse template (pulse_templates_B0)
*                               of the calibration library.
*
* In general, it finds the two energies which straddle 'energy' in the calibration library and interpolates ('interpolate_model')
*   - If 'energy' is lower than the lowest energy in the library => The model with the lowest energy in the library is chosen
*   - If 'energy' is higher than the highest energy in the library => The model with the highest energy in the library is chosen
*
* Parameters:
* - energy: Energy of the pulse whose pulse template is being sought
* - reconstruct_init: Member of ReconstructInitSIRENA structure to initialize the reconstruction parameters (pointer and values). 
*                     In particular, this function uses the energies of the models ('energies') and their templates ('pulse_templates'), the number
*                     of templates in the library ('ntemplates'), the template duration ('template_duration') and the 'pulse_templates_B0'.
* - modelFound: Found template of the pulse whose energy is 'energy'
****************************************/
int find_model_energies(double energy, ReconstructInitSIRENA *reconstruct_init,gsl_vector **modelFound)
{
	string message = "";
	char valERROR[256];
	
	gsl_vector *modelFound_aux;
    if ((*modelFound)->size <= reconstruct_init->library_collection->pulse_templates[0].template_duration)
    {
        if ((modelFound_aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration)) == 0)
        {
            sprintf(valERROR,"%d",__LINE__-2);
            string str(valERROR);
            message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
            str.clear();
            EP_PRINT_ERROR(message,EPFAIL);
        }
    }
    else if ((*modelFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration)
    {
        if ((modelFound_aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration)) == 0)
        {
            sprintf(valERROR,"%d",__LINE__-2);
            string str(valERROR);
            message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
            str.clear();
            EP_PRINT_ERROR(message,EPFAIL);
        }
    } 

	long nummodels = reconstruct_init->library_collection->ntemplates;

	if (energy < gsl_vector_get(reconstruct_init->library_collection->energies,0))
	{
        if ((*modelFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration)
        {
            gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].ptemplate);
            gsl_vector_add_constant(modelFound_aux,-1.0*reconstruct_init->library_collection->baseline);
            gsl_vector_scale(modelFound_aux,energy/gsl_vector_get(reconstruct_init->library_collection->energies,0));
        }
        else
        {
            gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_B0[0].ptemplate);
            gsl_vector_scale(modelFound_aux,energy/gsl_vector_get(reconstruct_init->library_collection->energies,0));
        }
	}
	else if (energy > gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1))
	{
        if ((*modelFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration)
        {
            gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[nummodels-1].ptemplate);
            gsl_vector_add_constant(modelFound_aux,-1.0*reconstruct_init->library_collection->baseline);
            gsl_vector_scale(modelFound_aux,energy/gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1));
        }
        else
        {
            gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_B0[nummodels-1].ptemplate);
            gsl_vector_scale(modelFound_aux,energy/gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1));
        }
	}
	else
	{
        gsl_matrix *pulse_templatesMaxLengthFixedFilter_B0;
        
        if ((*modelFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration)
        {
            pulse_templatesMaxLengthFixedFilter_B0 = gsl_matrix_alloc(reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration, nummodels);
            for (int i=0;i<nummodels-1;i++)
            {
                gsl_matrix_set_row(pulse_templatesMaxLengthFixedFilter_B0,i,reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[i].ptemplate);
            }
            gsl_matrix_add_constant(pulse_templatesMaxLengthFixedFilter_B0,-1.0*reconstruct_init->library_collection->baseline);
        }
        
		for (int i=0;i<nummodels-1;i++)
		{
			if (fabs(energy-gsl_vector_get(reconstruct_init->library_collection->energies,i))<1e-6)
			{
                if ((*modelFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration)
                {
                    gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[i].ptemplate);
                    gsl_vector_add_constant(modelFound_aux,-1.0*reconstruct_init->library_collection->baseline);
                }
                else
                {
                    gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_B0[i].ptemplate);
                }
				gsl_vector_scale(modelFound_aux,energy/gsl_vector_get(reconstruct_init->library_collection->energies,i));

				break;
			}
			else if ((energy > gsl_vector_get(reconstruct_init->library_collection->energies,i)) && (energy < gsl_vector_get(reconstruct_init->library_collection->energies,i+1)))
			{
				// Interpolate between the two corresponding rows in "models"
				// It is not necessary to check the allocation because the same kind of allocation has been checked previously
				gsl_vector *modelAux;
				if ((*modelFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration)
                {
                    modelAux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration);
                }
                else
                {
                    modelAux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);
                    
                }
				gsl_vector_set_zero(modelAux);

                if ((*modelFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration)
                {
                    gsl_vector *pulse_templatesMaxLengthFixedFilter_B0row_i = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration);
                    gsl_vector *pulse_templatesMaxLengthFixedFilter_B0row_iplus1 = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration); 
                    gsl_matrix_get_row(pulse_templatesMaxLengthFixedFilter_B0row_i, pulse_templatesMaxLengthFixedFilter_B0,i);
                    gsl_matrix_get_row(pulse_templatesMaxLengthFixedFilter_B0row_iplus1, pulse_templatesMaxLengthFixedFilter_B0,i+1);
                    
                    if (interpolate_model(&modelAux,energy,pulse_templatesMaxLengthFixedFilter_B0row_i,gsl_vector_get(reconstruct_init->library_collection->energies,i),
					pulse_templatesMaxLengthFixedFilter_B0row_iplus1,gsl_vector_get(reconstruct_init->library_collection->energies,i+1)))
                    {
                        message = "Cannot run interpolate_model with two rows in models";
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                    }
                    gsl_vector_free(pulse_templatesMaxLengthFixedFilter_B0row_i); pulse_templatesMaxLengthFixedFilter_B0row_i = 0;
                    gsl_vector_free(pulse_templatesMaxLengthFixedFilter_B0row_iplus1); pulse_templatesMaxLengthFixedFilter_B0row_iplus1 = 0;
                }
                else
                {
                    if (interpolate_model(&modelAux,energy,reconstruct_init->library_collection->pulse_templates_B0[i].ptemplate,gsl_vector_get(reconstruct_init->library_collection->energies,i),
					reconstruct_init->library_collection->pulse_templates_B0[i+1].ptemplate,gsl_vector_get(reconstruct_init->library_collection->energies,i+1)))
                    {
                        message = "Cannot run interpolate_model with two rows in models";
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                    }
                }

				gsl_vector_memcpy(modelFound_aux,modelAux);
                //cout<<"find_model_energiesD"<<endl;
				gsl_vector_free(modelAux); modelAux = 0;

				break;
			}
		}

        if ((*modelFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration)
        {
            gsl_matrix_free(pulse_templatesMaxLengthFixedFilter_B0); pulse_templatesMaxLengthFixedFilter_B0 = 0;
        }

	}

	gsl_vector_view temp;
	temp = gsl_vector_subvector(modelFound_aux,0,(*modelFound)->size);
	gsl_vector_memcpy(*modelFound,&temp.vector);

	gsl_vector_free(modelFound_aux); modelFound_aux = 0;
    
    message.clear();

    return(EPOK);
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION xx ************************************************************
* find_model_maxDERS function: This function uses the maximum of the derivative of the filtered pulse (maxDER) in order to choose the proper
*                              pulse template (pulse_templates_filder) of the calibration library.
*
* In general, it finds the two maxDER which straddle 'maxDER' in the calibration library and interpolates ('interpolate_model'):
*   - If maxDER is lower than the lowest maxDER in the library => The model with the lowest maxDER in the library is chosen
*   - If maxDER is higher than the highest maxDER in the library => The model with the highest maxDER in the library is chosen
*
* Parameters:
* - maxDER: Maximum of the derivative of the filtered pulse whose pulse template is being sought
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values). 
*                     In particular, this function uses the number of templates in the library ('ntemplates'), the template duration ('template_duration'), 
*                     the filtered and differentiated templates ('pulse_templates_filder') and the 'maxDERs' of the templates
* - modelFound: Found template of the pulse whose maximum of the derivative of the filtered version is 'maxDER'
****************************************/
int find_model_maxDERs(double maxDER, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound)
{
	string message = "";
	char valERROR[256];

	long nummodels = reconstruct_init->library_collection->ntemplates;
        
    gsl_vector_view temp;

	if (maxDER < gsl_vector_get(reconstruct_init->library_collection->maxDERs,0))
	{
                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[0].ptemplate,0,(*modelFound)->size);
                gsl_vector_memcpy(*modelFound,&temp.vector);
		gsl_vector_scale(*modelFound,maxDER/gsl_vector_max(*modelFound));
	}
	else if (maxDER > gsl_vector_get(reconstruct_init->library_collection->maxDERs,nummodels-1))
	{
                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[nummodels-1].ptemplate,0,(*modelFound)->size);
                gsl_vector_memcpy(*modelFound,&temp.vector);
		gsl_vector_scale(*modelFound,maxDER/gsl_vector_max(*modelFound));
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
			if (fabs(maxDER-gsl_vector_get(reconstruct_init->library_collection->maxDERs,i))<1e-6)
			{
                                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[i].ptemplate,0,(*modelFound)->size);
                                gsl_vector_memcpy(*modelFound,&temp.vector);
                                gsl_vector_scale(*modelFound,maxDER/gsl_vector_max(*modelFound));

				break;
			}
			else if ((maxDER > gsl_vector_get(reconstruct_init->library_collection->maxDERs,i)) && (maxDER < gsl_vector_get(reconstruct_init->library_collection->maxDERs,i+1)))
			{
				// Interpolate between the two corresponding rows in "models"
                                
                                gsl_vector *modelA = gsl_vector_alloc((*modelFound)->size);
                                gsl_vector *modelB = gsl_vector_alloc((*modelFound)->size);
                                gsl_vector_set_zero(*modelFound);
                                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[i].ptemplate,0,(*modelFound)->size);
                                gsl_vector_memcpy(modelA,&temp.vector);
                                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[i+1].ptemplate,0,(*modelFound)->size);
                                gsl_vector_memcpy(modelB,&temp.vector);

                                if (interpolate_model(modelFound,maxDER,modelA,gsl_vector_get(reconstruct_init->library_collection->maxDERs,i),modelB,gsl_vector_get(reconstruct_init->library_collection->maxDERs,i+1)))
				{
					message = "Cannot run interpolate_model with two rows in models";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				
				gsl_vector_free(modelA); modelA = 0;
                gsl_vector_free(modelB); modelB = 0;

				break;
			}
		}
	}
	
	message.clear();

    return(EPOK);
}
/*xxxx end of SECTION xx xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 9 ************************************************************
* find_model_samp1DERS function: This function uses the 1st sample of the derivative of the filtered pulse (samp1DER) in order to choose the proper
*                              pulse template (pulse_templates_filder) of the calibration library.
*
* It finds the two samp1DER closer in the calibration library and interpolates ('interpolate_model')
*   - If samp1DER is lower than the lower samp1DER in the pulse models library => The model with
*     a lower samp1DER in the pulse models library is chosen
*   - If samp1DER is higher than the higher samp1DER in the pulse models library => The model with
*     a higher samp1DER in the pulse models library is chosen
*
* Parameters:
* - samp1DER: 1st sample of the derivative of the filtered pulse whose pulse template is looking for
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values). 
*                     In particular, this function uses the 1st samples of the derivative of the models ('samp1DERs'), the filter and differentiated 
*                     templates ('pulse_templates_filder'), the number of templates in the library ('ntemplates') and the template duration ('template_duration').
* - modelFound: Found template of the pulse whose 1st sample of the derivative of the filtered pulse is 'samp1DER'
****************************************/
int find_model_samp1DERs(double samp1DER, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound)
{
	string message = "";
	char valERROR[256];
	
	long nummodels = reconstruct_init->library_collection->ntemplates;
        
        gsl_vector_view temp;

	if (samp1DER < gsl_vector_get(reconstruct_init->library_collection->samp1DERs,0))
	{
                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[0].ptemplate,0,(*modelFound)->size);
                gsl_vector_memcpy(*modelFound,&temp.vector);
                gsl_vector_scale(*modelFound,samp1DER/gsl_vector_get(*modelFound,0));
	}
	else if (samp1DER > gsl_vector_get(reconstruct_init->library_collection->samp1DERs,nummodels-1))
	{
                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[nummodels-1].ptemplate,0,(*modelFound)->size);
                gsl_vector_memcpy(*modelFound,&temp.vector);
                gsl_vector_scale(*modelFound,samp1DER/gsl_vector_get(*modelFound,0));  ///Creo que había un error antes escalando con 'gsl_vector_get(modelFound_aux,nummodels-1)'
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
			if ((samp1DER >= gsl_vector_get(reconstruct_init->library_collection->samp1DERs,i)) && (samp1DER < gsl_vector_get(reconstruct_init->library_collection->samp1DERs,i+1)))
			{
				// Interpolate between the two corresponding rows in "models"
                                gsl_vector *modelA = gsl_vector_alloc((*modelFound)->size);
                                gsl_vector *modelB = gsl_vector_alloc((*modelFound)->size);
                                gsl_vector_set_zero(*modelFound);
                                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[i].ptemplate,0,(*modelFound)->size);
                                gsl_vector_memcpy(modelA,&temp.vector);
                                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[i+1].ptemplate,0,(*modelFound)->size);
                                gsl_vector_memcpy(modelB,&temp.vector);

                                if (interpolate_model(modelFound,samp1DER,modelA,gsl_vector_get(reconstruct_init->library_collection->samp1DERs,i),modelB,gsl_vector_get(reconstruct_init->library_collection->samp1DERs,i+1)))
				{
					message = "Cannot run interpolate_model with two rows in models";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				gsl_vector_free(modelA); modelA = 0;
                gsl_vector_free(modelB); modelB = 0;


				break;
			}
		}
	}
	
	message.clear();
    
    return(EPOK);
}
/*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 10 ************************************************************
* interpolate_model: This function interpolates the pulse model, p(t,E), between two models of the pulse models library,
*                    p(t,E1) and p(t,E2), being E1<E<E2.
*
* According to the interpolation method:
*
*                p(t,E1)+p(t,E2)
* 	1. p(t,E)= -----------------
* 	                    2
* 	           E2-E           E-E1
* 	2. p(t,E)=-------p(t,E1)+-------p(t,E2)
*                  E2-E1          E2-E1
*
* The more intelligent averaging (2) is used instead the simplest method (1).
*
* Parameters:
* - modelFound: Found model of the pulse whose energy or maxDER is p_model
* - p_model: Parameter (energy or maxDER) of the pulse whose model is looking for
* - modelIn1: Model of the pulse whose parameter (energy or maxDER) is immediately lower than p_model in the library FITS file
* - p_modelIn1: Parameter (energy or maxDER) immediately lower than p_model in the library FITS file
* - modelIn2: Model of the pulse whose parameter (energy or maxDER) is immediately greater than p_model in the library FITS file
* - p_modelIn2: Parameter (energy or maxDER) immediately greater than p_model in the library FITS file
****************************************/
int interpolate_model(gsl_vector **modelFound, double p_model, gsl_vector *modelIn1, double p_modelIn1, gsl_vector *modelIn2, double p_modelIn2)
{
	// Declare variables
	double factor1, factor2;
	
	// It is not necessary to check the allocation because 'modelIn1' and 'modelIn2' sizes must already be > 0
	gsl_vector *modelIn1Aux = gsl_vector_alloc(modelIn1->size);
	gsl_vector *modelIn2Aux = gsl_vector_alloc(modelIn2->size);
	gsl_vector_memcpy(modelIn1Aux,modelIn1);
	gsl_vector_memcpy(modelIn2Aux,modelIn2);

	// Method 1: The simplest method
	/*gsl_vector_add(*modelFound,modelIn1);
	gsl_vector_add(*modelFound,modelIn2);
	gsl_vector_scale(*modelFound,0.5);*/

	// Method 2: A bit more intelligent averaging
	factor1 = (p_modelIn2-p_model)/(p_modelIn2-p_modelIn1);
	factor2 = (p_model-p_modelIn1)/(p_modelIn2-p_modelIn1);
	gsl_vector_scale(modelIn1Aux,factor1);
	gsl_vector_scale(modelIn2Aux,factor2);
        gsl_vector_memcpy(*modelFound,modelIn1Aux);
	gsl_vector_add(*modelFound,modelIn2Aux);
	
	gsl_vector_free(modelIn1Aux); modelIn1Aux = 0;
	gsl_vector_free(modelIn2Aux); modelIn2Aux = 0;

    return(EPOK);
}
/*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 11 ************************************************************
* findPulsesCAL: This function is going to find the pulses in a record (in the CALibration mode) by using the function findTstartCAL
*
* - Declare variables
* - Establish the threshold (call medianKappaClipping)
* - Find pulses (call findTstartCAL)
* - If at least a pulse is found
* 	- Get the 'pulseheight' of each found pulse (in order to be used to build the pulse templates library)
* - Free allocated GSL vectors
*
* Parameters:
* - vectorin: Not filtered record
* - vectorinDER: Derivative of the low-pass filtered 'vectorin'
* - tstart: Starting time of the found pulses into the record (in samples)
* - quality: Quality of the found pulses into the record
* - pulseheight: Pulse height of the found pulses into the record
* - maxDERgsl: Maximum of the derivative of the found low-pass filtered pulses into the record
* - nPulses: Number of found pulses
* - threshold: Threshold used to find the pulses (output parameter because it is necessary out of the function)
* - scalefactor: Scale factor to calculate the LPF box-car length
* - samplingRate: Sampling rate
* - samplesup: Number of consecutive samples over the threshold to locate a pulse ('samplesUp')
* - nsgms: Number of Sigmas to establish the threshold
* - lb: Vector containing the baseline averaging length used for each pulse
* - lrs: Running sum length (equal to the 'Lrs' global_variable)
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'maxPulsesPerRecord' and 'pulse_length' (and 'findTstartCAL' uses 'tstartPulse1', 'tstartPulse2' and 'tstartPulse3')
* - stopCriteriamkc: Used in medianKappaClipping (%)
* - kappamkc: Used in medianKappaClipping
****************************************/
int findPulsesCAL
(
	gsl_vector *vectorin,
	gsl_vector *vectorinDER,
	gsl_vector **tstart,
	gsl_vector **quality,
	gsl_vector **pulseheight,
	gsl_vector **maxDERgsl,

	int *nPulses,
	double *threshold,

	double scalefactor,
	double samplingRate,

	int samplesUp,
        double nsgms,

	double lb,
	double lrs,

	ReconstructInitSIRENA *reconstruct_init,

	double stopcriteriamkc,
	double kappamkc)
{
	string message = "";

	const double pi = 4.0 * atan(1.0);

	int pulseFound;
	double thresholdmediankappa;		// Threshold to look for pulses in the first derivative

	gsl_vector_set_zero(*quality);
	gsl_vector_set_zero(*pulseheight);	// Pulse height of the single pulses

	// Establish the threshold
	if (medianKappaClipping (vectorinDER, kappamkc, stopcriteriamkc, nsgms, (int)(pi*samplingRate*scalefactor), &thresholdmediankappa))
	{
		message = "Cannot run medianKappaClipping looking for single pulses";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	*threshold = thresholdmediankappa;
    //cout<<*threshold<<endl;
        
    // Find pulses
	if (findTstartCAL (reconstruct_init->maxPulsesPerRecord, vectorinDER, thresholdmediankappa, samplesUp, reconstruct_init, nPulses, tstart, quality, maxDERgsl))
	{
		message = "Cannot run findTstartCAL";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	if (*nPulses != 0)
	{
		// It is not necessary to check the allocation because 'reconstruct_init->maxPulsesPerRecord'='EventListSize'(input parameter) must already be > 0
		gsl_vector *Lbgsl = gsl_vector_alloc(reconstruct_init->maxPulsesPerRecord);	// If there is no free-pulses segments longer than Lb=>
		gsl_vector_set_all(Lbgsl,lb);                                  			// segments shorter than Lb will be used and its length (< Lb)
		                                                               			// must be used instead Lb in RS_filter
		gsl_vector *Bgsl;
                gsl_vector *sigmagsl;

		if (getB(vectorin, *tstart, *nPulses, &Lbgsl, reconstruct_init->pulse_length, &Bgsl, &sigmagsl))
		{
			message = "Cannot run getB";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		double energyaux = gsl_vector_get(*pulseheight,0);
		for (int i=0;i<*nPulses;i++)
		{
			if (i != *nPulses-1)	// Not last pulse in the record
			{
				if (getPulseHeight(vectorin, gsl_vector_get(*tstart,i), gsl_vector_get(*tstart,i+1), 0, lrs, gsl_vector_get(Lbgsl,i), gsl_vector_get(Bgsl,i), reconstruct_init->pulse_length, &energyaux))
				{
					message = "Cannot run getPulseHeight routine when pulse i=" + boost::lexical_cast<std::string>(i) + " is not the last pulse";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
			}
			else
			{
				// 'gsl_vector_get(*tstart,i+1)' is always 0 ('tstart' has 'maxPulsesPerRecord' size but only the fisrt 'nPulses' elements are filled with correct values)
				// In spite of 'gsl_vector_get(*tstart,i+1)' = 0, 'getPulseHeight' knows how to work with the last pulse
				if (getPulseHeight(vectorin, gsl_vector_get(*tstart,i), gsl_vector_get(*tstart,i+1), 1, lrs, gsl_vector_get(Lbgsl,i), gsl_vector_get(Bgsl,i), reconstruct_init->pulse_length, &energyaux))
				{
					message = "Cannot run getPulseHeight routine when pulse i=" + boost::lexical_cast<std::string>(i) + " is the last pulse";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
			}
			gsl_vector_set(*pulseheight,i,energyaux);
		}

		gsl_vector_free(Lbgsl); Lbgsl = 0;
		gsl_vector_free(Bgsl); Bgsl = 0;
        gsl_vector_free(sigmagsl); sigmagsl = 0;
	}
	
	message.clear();
	
	return(EPOK);
}
/*xxxx end of SECTION 11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 12 ************************************************************
* findTstartCAL function: This function finds the pulses tstarts in the input vector (first derivative of the filtered record)
*
* This function scans all values the derivative of the (low-pass filtered) record until it finds nSamplesUp consecutive values (due to the noise more than 1 value is
* required) over the threshold. To look for more pulses, it finds nSamplesUp consecutive values
* (due to the noise) under the threshold and then, it starts to scan again.
*
* - Declare variables
* - Allocate GSL vectors
* - It is possible to find the tstarts...
* 	- Obtain tstart of each pulse in the derivative:
* 		- If der_i>threshold and foundPulse=false, it looks for nSamplesUp consecutive samples over the threshold
*   			- If not, it looks again for a pulse crossing over the threshold
*	        	- If yes, a pulse is found (truncated if it is at the beginning)
*	   	- If der_i>threshold and foundPulse=true, it looks for a sample under the threshold
*   	  		- If not, it looks again for a sample under the threshold
*       		- If yes, it looks for nSamplesUp consecutive samples under the threshold and again it starts to look for a pulse
* - ...or to use the tstart provided as input parameters
* 	- Obtain the maxDERs of the pulses whose tstarts have been provided
*
* Parameters:
* - maxPulsesPerRecord: Expected maximum number of pulses per record in order to not allocate the GSL variables with the size of the input vector
* - der: First derivative of the (low-pass filtered) record
* - adaptativethreshold: Threshold
* - nSamplesUp: Number of consecutive samples over the threshold to 'find' a pulse
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'pulse_length', 'tstartPulse1', 'tstartPulse2' and 'tstartPulse3'
* - numberPulses: Number of found pulses
* - tstartgsl: Pulses tstart (in samples)
* - flagTruncated: Flag indicating if the pulse is truncated (inside this function only initial truncated pulses are classified)
* - maxDERgsl: Maximum of the first derivative of the (low-pass filtered) record inside each found pulse
******************************************************************************/
int findTstartCAL
(
	int maxPulsesPerRecord,

	gsl_vector *der,
	double adaptativethreshold,
	int nSamplesUp,

	ReconstructInitSIRENA *reconstruct_init,

	int *numberPulses,
	
	gsl_vector **tstartgsl,
	gsl_vector **flagTruncated,
	gsl_vector **maxDERgsl)
{
	string message="";
	char valERROR[256];

	// Declare variables
	int szRw = der->size;	 // Size of segment to process
	*numberPulses = 0;
	bool foundPulse = false;
	int i = 0;		// To go through the elements of a vector
	
	gsl_vector_view temp;	// In order to handle with gsl_vector_view (subvectors)
	
	// Allocate GSL vectors
	// It is not necessary to check the allocation because 'maxPulsesPerRecord'='EventListSize'(input parameter) must already be > 0
	if (*tstartgsl)		{gsl_vector_free(*tstartgsl); *tstartgsl = 0;}
	*tstartgsl = gsl_vector_alloc(maxPulsesPerRecord);	
	if (*flagTruncated)	{gsl_vector_free(*flagTruncated); *flagTruncated = 0;}
	*flagTruncated = gsl_vector_alloc(maxPulsesPerRecord);
	gsl_vector_set_zero(*flagTruncated);
	if (*maxDERgsl)	{gsl_vector_free(*maxDERgsl); *maxDERgsl = 0;}
	*maxDERgsl = gsl_vector_alloc(maxPulsesPerRecord);	// Maximum of the first derivative
	gsl_vector_set_all(*maxDERgsl,-1E3);

	int cntUp = 0;
	int cntDown = 0;
	double possibleTstart;
	double possiblemaxDER;
        
	// To provide the tstarts (or not)
	bool findTstarts = true;
        if ((isNumber(reconstruct_init->tstartPulse1)) && (atoi(reconstruct_init->tstartPulse1) != 0)) findTstarts = false;
        
    /*cout<<"adaptativethreshold: "<<adaptativethreshold<<endl;
    for (int i=0;i<der->size;i++)
        cout<<i<<" "<<gsl_vector_get(der,i)<<endl;*/
	
	if (findTstarts == true)
	{
		// It looks for a pulse
		// If a pulse is found (foundPulse==true) => It looks for another pulse
		do
		{
			foundPulse = false;
			
			// It looks for a pulse since the beginning (or the previous pulse) to the end of the record
			while (i < szRw-1)
			{
				if (foundPulse == false)
				{
					// The first condition to detect a pulse is that the adjustedDerivative was over the threshold
					if (gsl_vector_get(der,i) > adaptativethreshold)
					{
                        cntUp++;
						if (cntUp == 1)
						{
							possibleTstart = i;
							possiblemaxDER = gsl_vector_get(der,i);
						}
						else if (cntUp == nSamplesUp)
						{
							if (*numberPulses == maxPulsesPerRecord)
							{
								sprintf(valERROR,"%d",__LINE__+5);
								string str(valERROR);
								message = "Found pulses in record>'EventListSize'(input parameter) => Change EventListSize or check if the threshold is too low => Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                                str.clear();
								EP_PRINT_ERROR(message,EPFAIL);
							}
							gsl_vector_set(*maxDERgsl,*numberPulses,possiblemaxDER);
							gsl_vector_set(*tstartgsl,*numberPulses,possibleTstart);
							if (possibleTstart == 0)	gsl_vector_set(*flagTruncated,*numberPulses,1);
							*numberPulses = *numberPulses +1;
							foundPulse = true;
							cntUp = 0;
						}
						cntDown = 0;
					}
					else cntUp = 0;
				}
				else
				{
					if (gsl_vector_get(der,i) > gsl_vector_get(*maxDERgsl,*numberPulses-1))
					{
						gsl_vector_set(*maxDERgsl,*numberPulses-1,gsl_vector_get(der,i));
					}
					else if (gsl_vector_get(der,i) < adaptativethreshold)
					{
						cntUp = 0;
						cntDown++;
						if (cntDown == nSamplesUp) // nSamplesUp samples under the threshold in order to look for another pulse
						{
							foundPulse = false; 
							cntDown = 0;
						}
					}
					else if (gsl_vector_get(der,i) > adaptativethreshold)
					{
						cntDown = 0;
					}
				}
				i++;
			}
		} while (foundPulse == true);
	}
	else
	{
		gsl_vector_view temp;
		// It is not necessary to check the allocation because 'reconstruct_init->pulse_length'='PulseLength'(input parameter) has been checked previously
		gsl_vector *model = gsl_vector_alloc(reconstruct_init->pulse_length);

		gsl_vector *tstartPulsei = gsl_vector_alloc(3);
//		gsl_vector_set(tstartPulsei,0,reconstruct_init->tstartPulse1);
		gsl_vector_set(tstartPulsei,1,reconstruct_init->tstartPulse2);
		gsl_vector_set(tstartPulsei,2,reconstruct_init->tstartPulse3);

		if (reconstruct_init->tstartPulse2 == 0) 	*numberPulses = 1;
		else if (reconstruct_init->tstartPulse3 == 0) 	*numberPulses = 2;
		else						*numberPulses = 3;

		double maxPulse;
		int cnt;
		for (int i=0;i<*numberPulses;i++)
		{
			gsl_vector_set(*tstartgsl,i,gsl_vector_get(tstartPulsei,i));

			if (i != *numberPulses-1)	
			{
				if ((gsl_vector_get(tstartPulsei,i) < 0) || (gsl_vector_get(tstartPulsei,i) > der->size-2)
					|| (gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i) < 1) || (gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i) > der->size-gsl_vector_get(tstartPulsei,i)))
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                    str.clear();
					EP_PRINT_ERROR(message,EPFAIL);
				}
				temp = gsl_vector_subvector(der,gsl_vector_get(tstartPulsei,i),gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i));
			}
			else
			{
				if ((gsl_vector_get(tstartPulsei,i) < 0) || (gsl_vector_get(tstartPulsei,i) > der->size-2)
					|| (szRw-gsl_vector_get(tstartPulsei,i) < 1) || (szRw-gsl_vector_get(tstartPulsei,i) > der->size-gsl_vector_get(tstartPulsei,i)))
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                    str.clear();
					EP_PRINT_ERROR(message,EPFAIL);
				}
				temp = gsl_vector_subvector(der,gsl_vector_get(tstartPulsei,i),szRw-gsl_vector_get(tstartPulsei,i));
			}
			
			gsl_vector_set(*maxDERgsl,i,gsl_vector_max(&temp.vector));
		}

		gsl_vector_free(tstartPulsei); tstartPulsei = 0;
		gsl_vector_free(model); model = 0;
	}
	
	message.clear();

	return (EPOK);
}
/*xxxx end of SECTION 12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 13 ************************************************************
* InitialTriggering function: This function finds the first pulse in the input vector, first derivative of the (low-pass filtered) record
*
* - Declare variables
* - Stablish the threshold
* - It is necessary to find the tstart of the first pulse...
*   Obtain tstart of the first pulse in the derivative if der_i>threshold
* - ...or to use the tstart provided as input parameter
*
* Parameters:
* - derivative: First derivative of the (low-pass filtered) record
* - nSgms: Number of Sigmas to establish the threshold 
* - scalefactor: Scale factor to calculate the LPF box-car length ('scaleFactor')
* - samplingRate: Sampling rate
* - stopCriteriamkc: Used in medianKappaClipping (%)
* - kappamkc: Used in medianKappaClipping
* - triggerCondition: true -> the algorithm has found the first event
*                     false -> the algorithm has not found any event
* - tstart: First event tstart (in samples)
* - flagTruncated: Flag indicating if the event is truncated (inside this function only initial truncated events are classified)
* - threshold: Calculated threshold  (output parameter because it is necessary out of the function)
****************************************/
int InitialTriggering
(
	gsl_vector *derivative,

	double nSgms,
	double scalefactor,
	double samplingRate,
	double stopcriteriamkc,
	double kappamkc,

	bool *triggerCondition,
	int *tstart,
	int *flagTruncated,

	double *threshold)
{
	string message = "";

	// Declare variables
	const double pi = 4.0 * atan(1.0);
	int sizeRecord = derivative->size;	// Size of segment to process
	int i = 0;

	*triggerCondition = false;

	// Stablish the threshold
	if (medianKappaClipping (derivative, kappamkc, stopcriteriamkc, nSgms, (int)(pi*samplingRate*scalefactor), threshold))
	{
		message = "Cannot run medianKappaClipping doing the initial triggering";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	//*threshold = 30; // samprate = 156.250 kHz
	//*threshold = 55; // samprate/2 = 78.125 kHz
	
	while (i < sizeRecord-1)
        {
            if (gsl_vector_get(derivative,i) > *threshold)
            {
                *triggerCondition = true;
                *tstart = i;
                if (i == 0)	*flagTruncated = 1;
                else		*flagTruncated = 0;
                
                break;
            }
            
            i++;
        }
        
    message.clear();

	return(EPOK);
}
/*xxxx end of SECTION 13 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 14 ************************************************************
* FindSecondaries function: This function runs after InitialTriggering to find all the events (except the first one) in the record first derivative of the (low-pass filtered) record by using the Adjusted Derivative detection method
*
* - Declare variables
* - Establishing the criteria of the slope of the derivative depending on the sampling rate
* - It looks for an event and if a pulse is found (foundPulse==true), it looks for another event
* 
* 	- It looks for an event since the beginning (or the previous event) to the end of the record
*         The first condition to detect an event is that the adjustedDerivative was over the threshold
* 
*           - Select the model of the found pulse from the libary by using the 1st sample of the derivative (samp1DER)
*           - Dot product between the detected pulse and the pulse template in 3 different lags
*            
*               - If maximum of the dot product found => Stop calculating dot products in more lags
*               - If maximum of the dot product not found => Calculate dot products in more lags (number of lags is limited to 5)
*            
*           - If maximum of the dot product not found => tstart is the first sample crossing above the threshold (without jitter)
*              
*               - Average of the first 4 samples of the derivative
*               - Find model in order to subtract
*
*           - If maximum of the dot product found => Parabola analytically defined => Locate the maximum => New tstart (with jitter)
*            
*               - Iterative process in order to extract the best template from the library:
*                   - samp1DER correction
*                   - Find the model from the libary by using the corrected samp1DER
*                   - Dot product in 3 lags
*                   - Locate the maximum of the parabola
*               - samp1DER correction
*               - Find model in order to subtract
*               - Template correction
*               - Average of the first 4 samples of the derivative
*                
*           - The second condition to detect an event is meeting the criteria of the slope of the derivative
* 
* 	- Subtract the model from the adjusted derivative
* 
* 		- Select the model of the found event from the libary by using the first sample of the derivative
* 		- Subtract
*  
* - Free allocated GSL vectors
*
* Parameters:
* - maxPulsesPerRecord: Expected maximum number of events per record in order to not allocate the GSL variables with the size of the input vector ('EventListSize')
*
* - adjustedDerivative: First derivative of the (low-pass filtered) record 
* - adaptativethreshold: Threshold
* - samprate: Sampling rate
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses some parameters ('pulse_length', 'tstartPulsex', 'EnergyMethod') and the templates
* - tstartFirstEvent: Tstart of the first event of the record (in samples) found by 'InitialTriggering'
* - numberPulses: Number of found events
* - tstartgsl: Staring time of the found events (in samples)
* - flagTruncated: Flag indicating if the event is truncated (inside this function only initial truncated pulses are classified)
* - maxDERgsl: Maximum of the derivative of the event
*              It is going to be used by 'find_matchfilter', 'find_matchfilterDAB', 'find_optimalfilter', 'find_optimalfilterDAB' and 'find_Esboundary'
* - samp1DERgsl: Average of the first 4 samples of the derivative of the event
* - lagsgsl: Number of necessary lags to establish the tstart (currently limited to 5)
****************************************/
int FindSecondaries
(
	int maxPulsesPerRecord,

	gsl_vector *adjustedDerivative,
	double adaptativethreshold,
        double samprate,

	ReconstructInitSIRENA *reconstruct_init,

	int tstartFirstEvent,

	int *numberPulses,
	gsl_vector **tstartgsl,
	gsl_vector **flagTruncated,
	gsl_vector **maxDERgsl,
	gsl_vector **samp1DERgsl,
    gsl_vector **lagsgsl)
{
	string message = "";
	char valERROR[256];

	// Declare variables
	bool foundPulse = false;
	int sizeRecord = adjustedDerivative->size;		// Size of segment to process
	*numberPulses = 0;
	gsl_vector_set_all(*maxDERgsl,-1E3);
	gsl_vector_set_all(*tstartgsl,-1E3);
    gsl_vector_set_all(*lagsgsl,-1E3);
	int i = tstartFirstEvent;
	// It is not necessary to check the allocation because 'reconstruct_init->pulse_length'='PulseLength'(input parameter) has been checked previously
    int pulse_length_ToConvolve = 25; // Instead of using pulse_length
    int pulse_length_ToSubtract = 100; // Instead of using pulse_length, because the results are equal
    gsl_vector *modelToConvolve = gsl_vector_alloc(pulse_length_ToConvolve);
    gsl_vector *modelToSubtract = gsl_vector_alloc(pulse_length_ToSubtract);
    double tstartJITTER;
    int indexMinNew,indexMaxNew;
    gsl_vector *indexMin;
    gsl_vector *indexMax;
    int i_indexMin, i_indexMax;
    int indexM = 0;
    double samp1DER_Aux;
    double m;
    double prev_samp1DER_Aux, next_samp1DER_Aux;
	
	// To provide the tstarts (or not)
	//bool findTstarts = true;
	//if (reconstruct_init->tstartPulse1 != 0)	findTstarts = false;
	
	// It is not necessary to check the allocation because 'maxPulsesPerRecord'='EventListSize'(input parameter) must already be > 0
	gsl_vector *index_maxDERgsl = gsl_vector_alloc(maxPulsesPerRecord);	// Index where the maximum of the first derivative of the (low-pass filtered) event is
        int ider;
	
	int numlags;	// Or 1 if NO lags 
	gsl_vector *lags_vector;
	gsl_vector *convolutionLags;
	int indexmax;
	
	//gsl_vector_view temp;
	double a,b,c;
	double xmax;
        bool exitLags = false;
        int newLag;
        double newconvolutionLags;
        int indexLags;
        gsl_vector *sublags_vector = gsl_vector_alloc(3);
        gsl_vector *subconvolutionLags_vector = gsl_vector_alloc(3);
        
        gsl_vector *ThreePoints_x = gsl_vector_alloc(3);
        gsl_vector_set(ThreePoints_x,0,1);
        gsl_vector_set(ThreePoints_x,1,2);
        gsl_vector_set(ThreePoints_x,2,3);
        gsl_vector *ThreePoints_y = gsl_vector_alloc(3);
        double a3Points,b3Points;
        double angleStart0 = -999.0, angleStart1 = -999.0;
        
        int previouslyFalsePulse = -1;
        
        // Establishing the criteria of the slope of the derivative depending on the sampling rate
        double criteriaDER_value;
        if (samprate == 156250)
            //criteriaDER_value = 86.7; // In degrees,  samprate = 156250 Hz
            criteriaDER_value = 88.3; // In degrees,  samprate = 156250 Hz
        else if (samprate == 156250/2)
            //criteriaDER_value = 87.5; // In degrees, samprate/2 = 78125 Hz
            criteriaDER_value = 89; // In degrees, samprate/2 = 78125 Hz      
        
        double sum_samp1DER;
        int limitMin, limitMax;
      
        gsl_vector_view temp;
        
        numlags = 3;
        lags_vector = gsl_vector_alloc(numlags);
        convolutionLags = gsl_vector_alloc(numlags);
        
        int tstartWITHOUTLags;
        
        // It looks for a pulse
        // If a pulse is found (foundPulse==true) => It looks for another pulse
        do
        {	
            foundPulse = false;
            
            // It looks for a pulse since the beginning (or the previous pulse) to the end of the record
            while (i < sizeRecord-1)
            {
                if (foundPulse == false)
                {	
                    // The first condition to detect a pulse is that the adjustedDerivative was over the threshold
                    if (gsl_vector_get(adjustedDerivative,i) > adaptativethreshold)
                        //if ((gsl_vector_get(adjustedDerivative,i) > adaptativethreshold) && (i>0))
                    {
                        //if (i == 0) cout<<"gsl_vector_get(adjustedDerivative,0): "<<gsl_vector_get(adjustedDerivative,0)<<" adaptativethreshold: "<<adaptativethreshold<<endl;
                        if (*numberPulses == (*maxDERgsl)->size)
                        {
                            sprintf(valERROR,"%d",__LINE__+5);
                            string str(valERROR);
                            message = "Found pulses in record>'EventListSize'(input parameter) => Change EventListSize or check if the threshold is too low => Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                            str.clear();
                            EP_PRINT_ERROR(message,EPFAIL);
                        }
                        gsl_vector_set(*maxDERgsl,*numberPulses,gsl_vector_get(adjustedDerivative,i));
                        gsl_vector_set(index_maxDERgsl,*numberPulses,i);
                        gsl_vector_set(*samp1DERgsl,*numberPulses,gsl_vector_get(adjustedDerivative,i));
                        gsl_vector_set(*tstartgsl,*numberPulses,i);
                        tstartWITHOUTLags = i;
                        if (i == 0)	gsl_vector_set(*flagTruncated,*numberPulses,1);
                        *numberPulses = *numberPulses +1;
                        foundPulse = true;
                        //cout<<"Supera el umbral en "<<i<<" numberPulses="<<*numberPulses<<endl;
                        //cout<<"Supera el umbral en "<<i<<" "<<gsl_vector_get(adjustedDerivative,i-1)<<" "<<gsl_vector_get(adjustedDerivative,i)<<" "<<adaptativethreshold<<endl;
                        if (i == 0) 
                            EP_PRINT_ERROR("FIRST sample is being above the threshold",-999);
                        
                        indexM = 0;
                        indexMin = gsl_vector_alloc(10);
                        indexMax = gsl_vector_alloc(10);
                        gsl_vector_set_all(indexMin,9999);
                        gsl_vector_set_all(indexMax,9999);
                    }
                    i++;
                }
                else
                {
                    gsl_vector_set_zero(convolutionLags);
                    exitLags = false;
                    if ((((strcmp(reconstruct_init->EnergyMethod,"I2RALL") == 0) && (gsl_vector_get(index_maxDERgsl,*numberPulses-1)-gsl_vector_get(*tstartgsl,*numberPulses-1) >= 0)) || ((strcmp(reconstruct_init->EnergyMethod,"I2RALL") != 0) && 
                        (gsl_vector_get(index_maxDERgsl,*numberPulses-1)-gsl_vector_get(*tstartgsl,*numberPulses-1) >= 0))) && (i < sizeRecord-1))
                    {
                        ider = gsl_vector_get(index_maxDERgsl,*numberPulses-1)+1;
                        while (gsl_vector_get(adjustedDerivative,ider) > gsl_vector_get(*maxDERgsl,*numberPulses-1))
                        {
                            gsl_vector_set(*maxDERgsl,*numberPulses-1,gsl_vector_get(adjustedDerivative,ider));
                            ider++;
                            if (ider == adjustedDerivative->size)       break;
                            
                        }
                        
                        // Select the model of the found pulse from the libary by using the 1st sample of the derivative (samp1DER)
                        if (find_model_samp1DERsNoReSCLD(gsl_vector_get(*samp1DERgsl,*numberPulses-1), reconstruct_init, &modelToConvolve, &indexMinNew, &indexMaxNew))
                        {
                            message = "Cannot run find_model_samp1DERs routine";
                            EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                        }
                        
                        if (numlags != 1)
                        {
                            for (int i=0;i<numlags;i++)	gsl_vector_set(lags_vector,i,-numlags/2+i);
                            
                            indexLags = 0;
                            
                            // Dot product between the detected pulse and the pulse template in 3 different lags
                            for (int j=0;j<numlags;j++)
                            {
                                for (int k=0;k<pulse_length_ToConvolve;k++)
                                {
                                    if (gsl_vector_get(*tstartgsl,*numberPulses-1)+gsl_vector_get(lags_vector,j) < 0) 
                                        gsl_vector_set(convolutionLags,j,-999.0);
                                    else if (gsl_vector_get(*tstartgsl,*numberPulses-1)+gsl_vector_get(lags_vector,j) > sizeRecord) 
                                        gsl_vector_set(convolutionLags,j,-999.0);
                                    else 	
                                    {
                                        if (gsl_vector_get(*tstartgsl,*numberPulses-1)+gsl_vector_get(lags_vector,j)+k < sizeRecord)
                                        {                                                                                    
                                            gsl_vector_set(convolutionLags,j,gsl_vector_get(convolutionLags,j)+gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)+gsl_vector_get(lags_vector,j)+k)*gsl_vector_get(modelToConvolve,k));
                                        }
                                        else
                                            break;
                                    }
                                }
                                
                                //cout<<"lag="<<gsl_vector_get(lags_vector,j)<<", convolution="<<gsl_vector_get(convolutionLags,j)<<endl;
                            }
                            
                            indexmax = gsl_vector_max_index(convolutionLags);
                            
                            if (indexmax == 1)
                                // If maximum of the dot product found => Stop calculating dot products in more lags
                            {
                                gsl_vector_set(*lagsgsl,*numberPulses-1,numlags);
                                exitLags = true;
                            }
                            else
                                // If maximum of the dot product not found => Calculate dot products in more lags (number of lags is limited to 5)
                            {
                                do 
                                {   
                                    indexLags = indexLags + 1;
                                    if (indexmax == 0)  
                                    {       
                                        newLag = gsl_vector_get(lags_vector,0)-indexLags;
                                        gsl_vector_set(convolutionLags,2,gsl_vector_get(convolutionLags,1));
                                        gsl_vector_set(convolutionLags,1,gsl_vector_get(convolutionLags,0));
                                    }
                                    else    
                                    {
                                        newLag = gsl_vector_get(lags_vector,2)+indexLags;
                                        gsl_vector_set(convolutionLags,0,gsl_vector_get(convolutionLags,1));
                                        gsl_vector_set(convolutionLags,1,gsl_vector_get(convolutionLags,2));
                                    }
                                    
                                    if (gsl_vector_get(*tstartgsl,*numberPulses-1)+newLag < 0)  break;
                                    
                                    newconvolutionLags = 0.0;
                                    for (int k=0;k<pulse_length_ToConvolve;k++)
                                    {
                                        if (gsl_vector_get(*tstartgsl,*numberPulses-1)+newLag < 0) 
                                            newconvolutionLags = -999.0;
                                        else if (gsl_vector_get(*tstartgsl,*numberPulses-1)+newLag > sizeRecord) 
                                            newconvolutionLags = -999.0;
                                        else
                                        {
                                            if (gsl_vector_get(*tstartgsl,*numberPulses-1)+newLag+k < sizeRecord)
                                            {
                                                newconvolutionLags = newconvolutionLags +                                                                                       +gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)+newLag+k)*gsl_vector_get(modelToConvolve,k);
                                                //cout<<k<<" "<<gsl_vector_get(*tstartgsl,*numberPulses-1)<<" "<<newLag<<" "<<gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)+newLag+k)<<" "<<gsl_vector_get(model,k)<<endl;
                                            }
                                            else
                                                break;
                                        }
                                    }
                                    //cout<<"lag="<<newLag<<", convolution="<<newconvolutionLags<<endl;
                                    
                                    if (indexmax == 0) 
                                    {
                                        if (newconvolutionLags < gsl_vector_get(convolutionLags,1))
                                        {
                                            exitLags = true;
                                            gsl_vector_set(lags_vector,0,newLag);
                                            gsl_vector_set(lags_vector,1,newLag+1);
                                            gsl_vector_set(lags_vector,2,newLag+2);
                                        }
                                        gsl_vector_set(convolutionLags,0,newconvolutionLags);
                                        gsl_vector_set(*lagsgsl,*numberPulses-1,fabs(newLag)+2);
                                    }
                                    else    
                                    {       
                                        if (newconvolutionLags < gsl_vector_get(convolutionLags,1))
                                        {
                                            exitLags = true;
                                            gsl_vector_set(lags_vector,2,newLag);
                                            gsl_vector_set(lags_vector,1,newLag-1);
                                            gsl_vector_set(lags_vector,0,newLag-2);
                                        }
                                        gsl_vector_set(convolutionLags,2,newconvolutionLags);
                                        gsl_vector_set(*lagsgsl,*numberPulses-1,fabs(newLag)+2);
                                    }
                                } while ((exitLags == false) && (gsl_vector_get(*lagsgsl,*numberPulses-1) < 5));
                            }
                            //cout<<"lags: "<<gsl_vector_get(*lagsgsl,*numberPulses-1)<<endl;
                            
                            if (exitLags == false)
                                // If maximum of the dot product not found
                            {
                                // tstart is the first sample crossing above the threshold (without jitter)
                                tstartJITTER = gsl_vector_get(*tstartgsl,*numberPulses-1);
                                //cout<<"tstartJITTER: "<<tstartJITTER<<endl;
                                
                                // Average of the first 4 samples of the derivative
                                sum_samp1DER = 0.0;
                                limitMin = 0;
                                limitMax = 3;
                                for (int index_samp1DER=limitMin;index_samp1DER<=limitMax;index_samp1DER++)
                                {
                                    if (gsl_vector_get(*tstartgsl,*numberPulses-1)+index_samp1DER > sizeRecord-1)
                                    {
                                        limitMax = index_samp1DER-1;
                                        limitMin = limitMax-3;
                                    }
                                }
                                
                                for (int index_samp1DER=limitMin;index_samp1DER<=limitMax;index_samp1DER++)
                                {
                                    sum_samp1DER = sum_samp1DER + gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)+index_samp1DER);
                                }
                                gsl_vector_set(*samp1DERgsl,*numberPulses-1,sum_samp1DER/4.0);
                                
                                // Find model in order to subtract
                                if (find_model_samp1DERsNoReSCLD(gsl_vector_get(*samp1DERgsl,*numberPulses-1), reconstruct_init, &modelToSubtract,&indexMinNew,&indexMaxNew))
                                {
                                    message = "Cannot run find_model_samp1DERs routine";
                                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                            }
                            else
                                // If maximum of the dot product found
                            {
                                // Parabola analytically defined
                                if (parabola3Pts (lags_vector, convolutionLags, &a, &b, &c))
                                {
                                    message = "Cannot run routine parabola3Pts";
                                    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                }
                                
                                // Maximum of the parabola
                                xmax = -b/(2*a);
                                
                                // New tstart (with jitter)
                                tstartJITTER = gsl_vector_get(*tstartgsl,*numberPulses-1)+xmax;
                                if (xmax >= 0)
                                {
                                    gsl_vector_set(*tstartgsl,*numberPulses-1,gsl_vector_get(*tstartgsl,*numberPulses-1)+floor(xmax));
                                    xmax = xmax - floor(xmax);
                                }
                                else
                                {
                                    gsl_vector_set(*tstartgsl,*numberPulses-1,gsl_vector_get(*tstartgsl,*numberPulses-1)+ceil(xmax));
                                    xmax = xmax - ceil(xmax);
                                }
                                //cout<<"tstart0: "<<gsl_vector_get(*tstartgsl,*numberPulses-1)<<endl;
                                //cout<<"xmax0: "<<xmax<<endl;
                                //cout<<"tstartJITTER0: "<<tstartJITTER<<endl;
                                
                                // Iterative process in order to extract the best template from the library
                                do 
                                {
                                    // samp1DER correction
                                    samp1DER_Aux = gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1));
                                    if (xmax < 0)
                                    {
                                        prev_samp1DER_Aux = gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)-1);
                                        m = (samp1DER_Aux-prev_samp1DER_Aux);   // m=(y1-y0)/(x1-x0)    x1-x0=1
                                        gsl_vector_set(*samp1DERgsl,*numberPulses-1,m*(1+xmax)+prev_samp1DER_Aux);
                                    }
                                    else if (xmax > 0)
                                    {
                                        next_samp1DER_Aux = gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)+1);
                                        m = (next_samp1DER_Aux-samp1DER_Aux);   // m=(y1-y0)/(x1-x0)    x1-x0=1
                                        gsl_vector_set(*samp1DERgsl,*numberPulses-1,m*xmax+samp1DER_Aux);
                                    }
                                    else
                                    {
                                        gsl_vector_set(*samp1DERgsl,*numberPulses-1,gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-
                                        1)));
                                    }
                                    
                                    gsl_vector_set(indexMin,indexM,indexMinNew);
                                    gsl_vector_set(indexMax,indexM,indexMaxNew);
                                    
                                    // Find the model from the libary by using the corrected samp1DER
                                    if (find_model_samp1DERsNoReSCLD(gsl_vector_get(*samp1DERgsl,*numberPulses-1), reconstruct_init, &modelToConvolve,&indexMinNew,&indexMaxNew))
                                    {
                                        message = "Cannot run find_model_samp1DERs routine";
                                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                    }
                                    
                                    // Dot product in 3 lags
                                    for (int i=0;i<numlags;i++)	gsl_vector_set(lags_vector,i,-numlags/2+i);
                                    gsl_vector_set_zero(convolutionLags);
                                    for (int j=0;j<numlags;j++)
                                    {
                                        for (int k=0;k<pulse_length_ToConvolve;k++)
                                        {
                                            if (gsl_vector_get(*tstartgsl,*numberPulses-1)+gsl_vector_get(lags_vector,j) < 0) 
                                                gsl_vector_set(convolutionLags,j,-999.0);
                                            else if (gsl_vector_get(*tstartgsl,*numberPulses-1)+gsl_vector_get(lags_vector,j) > sizeRecord) 
                                                gsl_vector_set(convolutionLags,j,-999.0);
                                            else 	
                                            {
                                                if (gsl_vector_get(*tstartgsl,*numberPulses-1)+gsl_vector_get(lags_vector,j)+k < sizeRecord)
                                                {                                                                                    
                                                    gsl_vector_set(convolutionLags,j,gsl_vector_get(convolutionLags,j)+gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)+gsl_vector_get(lags_vector,j)+k)*gsl_vector_get(modelToConvolve,k));
                                                }
                                                else
                                                    break;
                                            }
                                        }
                                        
                                        //cout<<"lagCONV2="<<gsl_vector_get(lags_vector,j)<<", convolution="<<gsl_vector_get(convolutionLags,j)<<endl;
                                    }
                                    
                                    // Locate the maximum of the parabola
                                    if (parabola3Pts (lags_vector, convolutionLags, &a, &b, &c))
                                    {
                                        message = "Cannot run routine parabola3Pts";
                                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                    }
                                    
                                    xmax = -b/(2*a);
                                    //cout<<"xmax: "<<xmax<<endl;
                                    //gsl_vector_set(*tstartgsl,*numberPulses-1,gsl_vector_get(*tstartgsl,*numberPulses-1)+round(xmax));  
                                    if (xmax >= 0)
                                    {
                                        gsl_vector_set(*tstartgsl,*numberPulses-1,gsl_vector_get(*tstartgsl,*numberPulses-1)+floor(xmax));
                                        xmax = xmax - floor(xmax);
                                    }
                                    else
                                    {
                                        gsl_vector_set(*tstartgsl,*numberPulses-1,gsl_vector_get(*tstartgsl,*numberPulses-1)+ceil(xmax));
                                        xmax = xmax - ceil(xmax);
                                    }
                                    tstartJITTER = gsl_vector_get(*tstartgsl,*numberPulses-1)+xmax;
                                    //cout<<"tstartgsl: "<<gsl_vector_get(*tstartgsl,*numberPulses-1)<<endl;
                                    //cout<<"tstartJITTER: "<<tstartJITTER<<endl;
                                    if ((tstartJITTER < 0) || (tstartJITTER >= sizeRecord))   break;
                                    samp1DER_Aux = gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1));
                                    
                                    i_indexMin = 999;
                                    for (int k=0;k<=indexM;k++)
                                    {
                                        if (gsl_vector_get(indexMin,k) == indexMinNew)
                                        {
                                            i_indexMin = k;
                                            break;    
                                        }
                                    }
                                    i_indexMax = -999;
                                    for (int k=0;k<=indexM;k++)
                                    {
                                        if (gsl_vector_get(indexMax,k) == indexMaxNew)
                                        {
                                            i_indexMax = k;
                                            break;    
                                        }
                                    }
                                    
                                    indexM = indexM + 1;
                                } while ((i_indexMin != i_indexMax) && (gsl_vector_get(*tstartgsl,*numberPulses-1)+1 < sizeRecord));
                                
                                
                                //cout<<"tstartJITTER: "<<tstartJITTER<<endl;
                                
                                if (indexMin != NULL) {gsl_vector_free(indexMin); indexMin = 0;}
                                if (indexMax != NULL) {gsl_vector_free(indexMax); indexMax = 0;}
                                
                                if ((tstartJITTER >= 0) && (tstartJITTER+1 < sizeRecord))
                                {
                                    // samp1der correction
                                    samp1DER_Aux = gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1));
                                    //cout<<"Muestra SIN CORREGIR con la que buscar el modelo: "<<samp1DER_Aux<<endl;
                                    if (xmax < 0)
                                    {
                                        prev_samp1DER_Aux = gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)-1);
                                        m = (samp1DER_Aux-prev_samp1DER_Aux);   // m=(y1-y0)/(x1-x0)    x1-x0=1
                                        gsl_vector_set(*samp1DERgsl,*numberPulses-1,m*(1+xmax)+prev_samp1DER_Aux);
                                    }
                                    else if (xmax > 0)
                                    {
                                        next_samp1DER_Aux = gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)+1);
                                        m = (next_samp1DER_Aux-samp1DER_Aux);   // m=(y1-y0)/(x1-x0)    x1-x0=1
                                        gsl_vector_set(*samp1DERgsl,*numberPulses-1,m*xmax+samp1DER_Aux);
                                    }
                                    else
                                    {
                                        gsl_vector_set(*samp1DERgsl,*numberPulses-1,gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)));
                                    }
                                    
                                    if (reconstruct_init->detectSP == 1)
                                    {
                                        // Find model in order to subtract
                                        gsl_vector *modelToSubtract_Aux = gsl_vector_alloc(pulse_length_ToSubtract);
                                        //cout<<"Muestra con la que buscar el modelo: "<<gsl_vector_get(*samp1DERgsl,*numberPulses-1)<<endl;
                                        if (find_model_samp1DERs(gsl_vector_get(*samp1DERgsl,*numberPulses-1), reconstruct_init, &modelToSubtract_Aux))
                                        {
                                            message = "Cannot run find_model_samp1DERs routine";
                                            EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                        }
                                        
                                        // Template correction
                                        for (int j=0;j<pulse_length_ToSubtract;j++)
                                        {
                                            if (xmax < 0)
                                            {
                                                if (j != pulse_length_ToSubtract-1)
                                                    gsl_vector_set(modelToSubtract,j,(gsl_vector_get(modelToSubtract_Aux,j+1)-      gsl_vector_get(modelToSubtract_Aux,j))*(-xmax)+gsl_vector_get(modelToSubtract_Aux,j));
                                                
                                                else 
                                                    gsl_vector_set(modelToSubtract,j,gsl_vector_get(modelToSubtract_Aux,j)); //?????????????????????
                                            }
                                            else if (xmax > 0)
                                            {
                                                if (j == 0)
                                                {
                                                    gsl_vector_set(modelToSubtract,j,(gsl_vector_get(modelToSubtract_Aux,j)-0)*(1-xmax)+0);
                                                }
                                                else 
                                                {
                                                    gsl_vector_set(modelToSubtract,j,(gsl_vector_get(modelToSubtract_Aux,j)-gsl_vector_get(modelToSubtract_Aux,j-1))*(1-xmax)+gsl_vector_get(modelToSubtract_Aux,j-1));
                                                }
                                            }
                                            else
                                            {
                                                gsl_vector_memcpy(modelToSubtract,modelToSubtract_Aux);
                                            }
                                        }
                                        
                                        //cout<<"modelToSubtract"<<endl;
                                        /*for (int j=0;j<10;j++)
                                         *                                       cout<<j<<" "<<gsl_vector_get(modelToSubtract_Aux,j)<<" "<<gsl_vector_get(modelToSubtract,j)<<endl;*/
                                        
                                        gsl_vector_free(modelToSubtract_Aux); modelToSubtract_Aux = 0;
                                    }
                                    else
                                        gsl_vector_set_all(modelToSubtract,1e6);
                                    
                                    // Average of the first 4 samples of the derivative
                                    sum_samp1DER = 0.0;
                                    limitMin = 0;
                                    limitMax = 3;
                                    for (int index_samp1DER=limitMin;index_samp1DER<=limitMax;index_samp1DER++)
                                    {
                                        if (gsl_vector_get(*tstartgsl,*numberPulses-1)+index_samp1DER > sizeRecord-1)
                                        {
                                            limitMax = index_samp1DER-1;
                                            limitMin = limitMax-3;
                                        }
                                    }
                                    
                                    for (int index_samp1DER=limitMin;index_samp1DER<=limitMax;index_samp1DER++)
                                    {
                                        sum_samp1DER = sum_samp1DER + gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)+index_samp1DER);
                                    }
                                    gsl_vector_set(*samp1DERgsl,*numberPulses-1,sum_samp1DER/4.0);
                                }
                            }
                        }
                        
                        // Slope of the straight line define by the tstart (digitized sample), and the preceding sample and the following one 
                        if ((tstartJITTER >= 0) && (tstartJITTER+2 < sizeRecord))
                        {
                            
                            if (gsl_vector_get(*tstartgsl,*numberPulses-1) == 0)
                            {
                                gsl_vector_set(ThreePoints_y,0,gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)));
                                gsl_vector_set(ThreePoints_y,1,gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)+1));
                                gsl_vector_set(ThreePoints_y,2,gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)+2));
                            }
                            else if (gsl_vector_get(*tstartgsl,*numberPulses-1) == sizeRecord-1)
                            {
                                gsl_vector_set(ThreePoints_y,0,gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)-2));
                                gsl_vector_set(ThreePoints_y,1,gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)-1));
                                gsl_vector_set(ThreePoints_y,2,gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)));
                            }
                            else
                            {
                                gsl_vector_set(ThreePoints_y,0,gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)-1));
                                gsl_vector_set(ThreePoints_y,1,gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)));
                                gsl_vector_set(ThreePoints_y,2,gsl_vector_get(adjustedDerivative,gsl_vector_get(*tstartgsl,*numberPulses-1)+1));
                            }
                            
                            //ax+b
                            if (polyFitLinear(ThreePoints_x,ThreePoints_y,&a3Points,&b3Points))
                            {
                                message = "Cannot run polyFitLinear routine";
                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                            }
                            angleStart1 = atan(a3Points)*180/pi; // In degrees
                        }
                        else
                            angleStart1 = -999;
                        
                        //cout<<*numberPulses<<" angleS: "<<angleStart1<<endl;
                        
                        if ((*numberPulses > 1) && ((tstartJITTER<= gsl_vector_get(*tstartgsl,*numberPulses-2)) || (fabs(tstartJITTER*(1/samprate)-gsl_vector_get(*tstartgsl,*numberPulses-2)*(1/samprate)) < 10e-6)) || (tstartJITTER < 0) || (tstartJITTER >= sizeRecord) ||
                            (angleStart1<criteriaDER_value))
                            //(angleStart1<criteriaDER_value) || (gsl_vector_get(convolutionLags,1) < 0))    
                        {
                            //cout<<"angleStart1<criteriaDER_value"<<endl;
                            previouslyFalsePulse = gsl_vector_get(*tstartgsl,*numberPulses-1);
                            *numberPulses = *numberPulses-1;
                            gsl_vector_set(*flagTruncated,*numberPulses,0);
                            foundPulse = false;
                        }
                        else
                        {
                            //cout<<*numberPulses<<" angleS: "<<angleStart1<<endl;
                            if ((*numberPulses == 1) || ((*numberPulses > 1 ) && (tstartJITTER > gsl_vector_get(*tstartgsl,*numberPulses-2)))) 
                            {
                                // Subtract the pulse template
                                if ((gsl_vector_get(*tstartgsl,*numberPulses-1) < 0) || (min(gsl_vector_get(*tstartgsl,*numberPulses-1)+pulse_length_ToSubtract,(double) sizeRecord)-1 > adjustedDerivative->size-1))
                                {
                                    sprintf(valERROR,"%d",__LINE__+9);
                                    string str(valERROR);
                                    message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                                    str.clear();
                                    EP_PRINT_ERROR(message,EPFAIL);
                                }
                                
                                for (int j=gsl_vector_get(*tstartgsl,*numberPulses-1);j<min(gsl_vector_get(*tstartgsl,*numberPulses-1)+pulse_length_ToSubtract,(double) sizeRecord);j++)
                                {
                                    gsl_vector_set(adjustedDerivative,j,gsl_vector_get(adjustedDerivative,j)-gsl_vector_get(modelToSubtract,j-gsl_vector_get(*tstartgsl,*numberPulses-1)));
                                }
                                /*cout<<"afterSubtracting"<<endl;
                                 *               for (int j=0;j<10;j++)
                                 *                   cout<<j<<" "<<gsl_vector_get(adjustedDerivative,j)<<" "<<gsl_vector_get(adjustedDerivative,j)<<endl;*/
                                gsl_vector_set(*tstartgsl,*numberPulses-1,tstartJITTER);  // This should be the new tstart
                                
                                if (gsl_vector_get(*flagTruncated,*numberPulses-1) == 1)	i = 0;
                                else
                                {
                                    if (tstartWITHOUTLags >= floor(gsl_vector_get(*tstartgsl,*numberPulses-1))) 
                                        i = tstartWITHOUTLags + 1;
                                    else
                                        i = floor(gsl_vector_get(*tstartgsl,*numberPulses-1))+1;
                                }
                                foundPulse = false; 
                            }
                        }
                    }
                    else 
                    {
                        if (strcmp(reconstruct_init->EnergyMethod,"I2RALL") != 0)
                        {
                            previouslyFalsePulse = gsl_vector_get(*tstartgsl,*numberPulses-1);
                            *numberPulses = *numberPulses-1;
                            gsl_vector_set(*flagTruncated,*numberPulses,0);
                            foundPulse = false;
                        }
                    }
                }
            }
        } while (foundPulse == true);
        
        gsl_vector_free(lags_vector); lags_vector = 0;
        gsl_vector_free(convolutionLags); convolutionLags = 0;

	// Free allocated GSL vectors
	//gsl_vector_free(model); model = 0;
	gsl_vector_free(index_maxDERgsl); index_maxDERgsl = 0;
	
	gsl_vector_free(modelToSubtract);modelToSubtract = 0;
    gsl_vector_free(modelToConvolve);modelToConvolve = 0;
        	
    gsl_vector_free(sublags_vector); sublags_vector = 0;
    gsl_vector_free(subconvolutionLags_vector); subconvolutionLags_vector = 0;
    gsl_vector_free(sublags_vector); sublags_vector = 0;
    gsl_vector_free(subconvolutionLags_vector); subconvolutionLags_vector = 0;
        
    gsl_vector_free(ThreePoints_x); ThreePoints_x = 0;
    gsl_vector_free(ThreePoints_y); ThreePoints_y = 0;
        
    message.clear();

	return(EPOK);
}
/*xxxx end of SECTION 14 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 15 ************************************************************
* find_model_samp1DERsNoReSCLD function: This function 
* 
******************************************************************************/
int find_model_samp1DERsNoReSCLD(double samp1DER, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound, int *indexMin, int *indexMax)
{
	string message = "";
	char valERROR[256];

	long nummodels = reconstruct_init->library_collection->ntemplates;
        
        gsl_vector_view temp;

	if (samp1DER < gsl_vector_get(reconstruct_init->library_collection->samp1DERs,0))
	{
                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[0].ptemplate,0,(*modelFound)->size);
                gsl_vector_memcpy(*modelFound,&temp.vector);
                
                *indexMin = -999;
                *indexMax = 0;
	}
	else if (samp1DER > gsl_vector_get(reconstruct_init->library_collection->samp1DERs,nummodels-1))
	{
                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[nummodels-1].ptemplate,0,(*modelFound)->size);
                gsl_vector_memcpy(*modelFound,&temp.vector);
                
                *indexMin = nummodels-1;
                *indexMax = -999;
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
			if ((samp1DER >= gsl_vector_get(reconstruct_init->library_collection->samp1DERs,i)) && (samp1DER < gsl_vector_get(reconstruct_init->library_collection->samp1DERs,i+1)))
			{
				// Interpolate between the two corresponding rows in "models"
	                        *indexMin = i;
                                *indexMax = i+1;
                                
                                gsl_vector *modelA = gsl_vector_alloc((*modelFound)->size);
                                gsl_vector *modelB = gsl_vector_alloc((*modelFound)->size);
                                gsl_vector_set_zero(*modelFound);
                                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[i].ptemplate,0,(*modelFound)->size);
                                gsl_vector_memcpy(modelA,&temp.vector);
                                temp = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_filder[i+1].ptemplate,0,(*modelFound)->size);
                                gsl_vector_memcpy(modelB,&temp.vector);

                                if (interpolate_model(modelFound,samp1DER,modelA,gsl_vector_get(reconstruct_init->library_collection->samp1DERs,i),modelB,gsl_vector_get(reconstruct_init->library_collection->samp1DERs,i+1)))
				{
					message = "Cannot run interpolate_model with two rows in models";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				
				gsl_vector_free(modelA); modelA = 0;
                gsl_vector_free(modelB); modelB = 0;

				break;
			}
		}
	}
	
	message.clear();

    return(EPOK);
}
/*xxxx end of SECTION 15 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 16 ************************************************************
* smoothDerivative function: This function applies the smooth derivative to the input vector
*
* Parameters:
* - invector: Input/Ouput GSL vector (input vector/smooth differentiated input vector)
* - N: box-car length (it must be an even number)
******************************************************************************/
int smoothDerivative (gsl_vector **invector, int N)
{
        if (N%2 != 0)   // Odd number
        {
                string message = "";
                message = "In the smoothDerivative function, N must be an even number)";
                EP_PRINT_ERROR(message,EPFAIL);
                message.clear();
        }
        else            // Even number
        {
                //gsl_vector *invectorAux = gsl_vector_alloc((*invector)->size+N);
                int szVct = (*invector)->size+N;
                //int szVct = invectorAux->size;
                //for (int i=0;i<N;i++)       gsl_vector_set(invectorAux,i,0.0);
                //for (int i=N;i<szVct;i++)   gsl_vector_set(invectorAux,i,gsl_vector_get(*invector,i-N));

                //cout<<"smoothDerivative0:"<<endl;
                //for (int i=0;i<szVct-N;i++)   cout<<i<<" "<<gsl_vector_get(*invector,i)<<endl;
                
                gsl_vector *window = gsl_vector_alloc(N);
                gsl_vector_set_all(window,1.0);
                for (int i=0;i<N/2;i++)     gsl_vector_set(window,i,-1.0);
                
                gsl_vector *conv = gsl_vector_alloc(szVct);
                gsl_vector_set_zero(conv);
                
            
                
                int index;
                for (int i=0;i<szVct;i++)
                {
                        //cout<<"i: "<<i<<endl;
                        index = N-1;
                        for (int j=i;j>=i-N+1;j--)
                        {
                                //cout<<"j: "<<j<<endl;
                                //if ((j>=0) && (j<szVct))
                                if ((j>=0) && (j<(*invector)->size))
                                {
                                        //cout<<i<<" "<<j<<" "<<index<<endl;
                                        gsl_vector_set(conv,i,gsl_vector_get(conv,i)+gsl_vector_get(*invector,j)*gsl_vector_get(window,index));
                                        //gsl_vector_set(conv,i,gsl_vector_get(conv,i)+gsl_vector_get(invectorAux,j)*gsl_vector_get(window,index));
                                        index = index - 1;
                                        //cout<<"if "<<i<<" "<<j<<" "<<index<<" conv="<<gsl_vector_get(conv,i)<<endl;
                                }
                        }
                }
                
                gsl_vector_free(window); window = 0;
                //cout<<"smoothDerivative1:"<<endl;
                //for (int i=0;i<szVct;i++)   cout<<i<<" "<<gsl_vector_get(conv,i)<<endl;
                
                //cout<<"(*invector)->size: "<<(*invector)->size<<endl;
                //cout<<"conv->size: "<<conv->size<<endl;
                gsl_vector_view temp;
                temp = gsl_vector_subvector(conv,0,szVct-N);
		gsl_vector_memcpy(*invector,&temp.vector);
                //for (int i=0;i<N-1;i++)     
                //        gsl_vector_set(*invector,i,0.0);
                //cout<<"sm0: "<<gsl_vector_get(*invector,(*invector)->size-3)<<" "<<gsl_vector_get(*invector,(*invector)->size-2)<<" "<<gsl_vector_get(*invector,(*invector)->size-1)<<endl;
                
                //gsl_vector_set(*invector,(*invector)->size-1,gsl_vector_get(*invector,(*invector)->size-2));
                
                //cout<<"sm1: "<<gsl_vector_get(*invector,(*invector)->size-3)<<" "<<gsl_vector_get(*invector,(*invector)->size-2)<<" "<<gsl_vector_get(*invector,(*invector)->size-1)<<endl;
                
                //cout<<"smoothDerivative2:"<<endl;
                //for (int i=0;i<(*invector)->size;i++)   cout<<i<<" "<<gsl_vector_get(*invector,i)<<endl;
                
                //gsl_vector_memcpy(invectorAux,conv);
                
                /*gsl_vector_view temp;
                temp = gsl_vector_subvector(conv,N,szVct-N);
		gsl_vector_memcpy(*invector,&temp.vector);*/
                
                gsl_vector_free(conv); conv = 0;
        }

	return (EPOK);
}
/*xxxx end of SECTION 16 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 17 ************************************************************
* FindSecondariesSTC function: This function runs after InitialTriggering to find all the events (except the first one) in the record first derivative 
*                           of the (low-pass filtered) record by using the Single Threshold Crossing detection method
*
* This function scans all values the derivative of the (low-pass filtered) record until it finds nSamplesUp consecutive values 
* (due to the noise more than 1 value is required) over the threshold. To look for more pulses, it finds nSamplesUp consecutive values
* (due to the noise) under the threshold and then, it starts to scan again.
*
* - Declare variables
* - It looks for an event and if a pulse is found, it looks for another event
*    
*       - It looks for an event since the beginning (or the previous event) to the end of the record. 
*         The condition to detect an event is that the adjustedDerivative was over the threshold at least nSamplesUp consecutive samples
*
* Parameters:
* - maxPulsesPerRecord: Expected maximum number of pulses per record in order to not allocate the GSL variables with the size of the input vector
* - der: First derivative of the (low-pass filtered) record
* - adaptativethreshold: Threshold
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'pulse_length', 'tstartPulse1', 'tstartPulse2' and 'tstartPulse3'
* - tstartFirstEvent: Tstart of the first event of the record (in samples) found by 'InitialTriggering'
* - numberPulses: Number of found pulses
* - tstartgsl: Pulses tstart (in samples)
* - flagTruncated: Flag indicating if the pulse is truncated (inside this function only initial truncated pulses are classified)
* - maxDERgsl: Maximum of the derivative of the pulse
* - samp1DERgsl: Average of the first 4 samples of the derivative of the event
******************************************************************************/
int FindSecondariesSTC
(
	int maxPulsesPerRecord,

	gsl_vector *der,
	double adaptativethreshold,

	ReconstructInitSIRENA *reconstruct_init,
 
        int tstartFirstEvent,

	int *numberPulses,
	
	gsl_vector **tstartgsl,
	gsl_vector **flagTruncated,
	gsl_vector **maxDERgsl,
        gsl_vector **samp1DERgsl
)
{
	string message="";
	char valERROR[256];

	// Declare variables
	int szRw = der->size;	 // Size of segment to process
	*numberPulses = 0;
	bool foundPulse = false;
	//int i = 0;		// To go through the elements of a vector
        int i = tstartFirstEvent;
	
	int nSamplesUp = reconstruct_init->samplesUp;
	int nSamplesDown = reconstruct_init->samplesDown;
	
	gsl_vector_view temp;	// In order to handle with gsl_vector_view (subvectors)
	
	// Allocate GSL vectors
	// It is not necessary to check the allocation because 'maxPulsesPerRecord'='EventListSize'(input parameter) must already be > 0
	if (*tstartgsl)		{gsl_vector_free(*tstartgsl);	*tstartgsl = 0;}
	*tstartgsl = gsl_vector_alloc(maxPulsesPerRecord);	
	if (*flagTruncated)	{gsl_vector_free(*flagTruncated); *flagTruncated = 0;}
	*flagTruncated = gsl_vector_alloc(maxPulsesPerRecord);
	gsl_vector_set_zero(*flagTruncated);
	if (*maxDERgsl)	{gsl_vector_free(*maxDERgsl); *maxDERgsl = 0;}
	*maxDERgsl = gsl_vector_alloc(maxPulsesPerRecord);	// Maximum of the first derivative
	gsl_vector_set_all(*maxDERgsl,-1E3);

	int cntUp = 0;
	int cntDown = 0;
	double possibleTstart;
	double possiblemaxDER;
        double possiblesamp1DER;
        
        double sum_samp1DER;
        int limitMin, limitMax;
        
        int nodetectSecondaries = 1;
        	
        // It looks for &tstartgsl,&qualitygsl, &maxDERgsl,&samp1DERgsla pulse
        // If a pulse is found (foundPulse==true) => It looks for another pulse
        do
        {
            foundPulse = false;
            
            // It looks for a pulse since the beginning (or the previous pulse) to the end of the record
            while ((i < szRw-1) && (nodetectSecondaries == 1))
            {
                if (foundPulse == false)
                {
                    // The condition to detect a pulse is that the adjustedDerivative was over the threshold at least nSamplesUp consecutive samples
                    if (gsl_vector_get(der,i) > adaptativethreshold)
                    {
                        cntUp++;
                        if (cntUp == 1)
                        {
                            possibleTstart = i;
                            //cout<<"cntUp: "<<cntUp<<" "<<i<<endl;
                            possiblemaxDER = gsl_vector_get(der,i);
                            possiblesamp1DER = gsl_vector_get(der,i);
                        }
                        //else if (cntUp == nSamplesUp)
                        if (cntUp == nSamplesUp)
                        {
                            if (*numberPulses == maxPulsesPerRecord)
                            {
                                sprintf(valERROR,"%d",__LINE__+5);
                                string str(valERROR);
                                message = "Found pulses in record>'EventListSize'(input parameter) => Change EventListSize or check if the threshold is too low => Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                                str.clear();
                                EP_PRINT_ERROR(message,EPFAIL);
                            }
                            gsl_vector_set(*maxDERgsl,*numberPulses,possiblemaxDER);
                            gsl_vector_set(*tstartgsl,*numberPulses,possibleTstart);
                            //cout<<"tstart0: "<<gsl_vector_get(*tstartgsl,*numberPulses)<<endl;
                            gsl_vector_set(*samp1DERgsl,*numberPulses,possiblesamp1DER);
                            
                            // Average of the first 4 samples of the derivative
                            sum_samp1DER = 0.0;
                            limitMin = 0;
                            limitMax = 3;
                            for (int index_samp1DER=limitMin;index_samp1DER<=limitMax;index_samp1DER++)
                            {
                                if (gsl_vector_get(*tstartgsl,*numberPulses)+index_samp1DER > szRw-1)
                                {
                                    limitMax = index_samp1DER-1;
                                    limitMin = limitMax-3;
                                }
                            }
                            
                            for (int index_samp1DER=limitMin;index_samp1DER<=limitMax;index_samp1DER++)
                            {
                                sum_samp1DER = sum_samp1DER + gsl_vector_get(der,gsl_vector_get(*tstartgsl,*numberPulses)+index_samp1DER);
                            }
                            gsl_vector_set(*samp1DERgsl,*numberPulses,sum_samp1DER/4.0);
                            
                            if (possibleTstart == 0)	gsl_vector_set(*flagTruncated,*numberPulses,1);
                            *numberPulses = *numberPulses +1;
                            foundPulse = true;
                            cntUp = 0;
                            nodetectSecondaries = reconstruct_init->detectSP;
                        }
                        cntDown = 0;
                    }
                    else cntUp = 0;
                }
                else
                {
                    if (gsl_vector_get(der,i) > gsl_vector_get(*maxDERgsl,*numberPulses-1))
                    {
                        gsl_vector_set(*maxDERgsl,*numberPulses-1,gsl_vector_get(der,i));
                    }
                    else if (gsl_vector_get(der,i) < adaptativethreshold)
                    {
                        cntUp = 0;
                        cntDown++;
                        if (cntDown == nSamplesDown) // nSamplesUp sample under the threshold in order to look for another pulse
                        {
                            foundPulse = false; 
                            cntDown = 0;
                        }
                    }
                    else if (gsl_vector_get(der,i) > adaptativethreshold)
                    {
                        cntDown = 0;
                    }
                }
                i++;
            }
        } while (foundPulse == true);

    message.clear();
    
	return (EPOK);
}
/*xxxx end of SECTION 17 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 18 ************************************************************
* noDetect function: This function runs if the starting time of the pulses are agiven as input parameters (tstartPulse1 != 0). 
*                    It looks for the maximum of the derivative of the pulse and the average of the first 4 samples of the derivative of the pulse.
*
* Parameters:
* - der: First derivative of the (low-pass filtered) record
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'maxDERs', 'tstartPulse1', 'tstartPulse2' and 'tstartPulse3'
* - numberPulses: Number of pulses
* - tstartgsl: Pulses tstart (in samples)
* - flagTruncated: Flag indicating if the pulse is truncated (inside this function only initial truncated pulses are classified)
* - maxDERgsl: Maximum of the derivative of the pulse
* - samp1DERgsl: Average of the first 4 samples of the derivative of the event
* - num_previousDetectedPulses: Number of previous detected pulses (to know the index to get the proper element from tstartPulse1_i in case tstartPulse1=nameFile)
* - samprate: Sampling rate
******************************************************************************/
int noDetect(gsl_vector *der, ReconstructInitSIRENA *reconstruct_init, int *numberPulses, gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl, gsl_vector **samp1DERgsl, long num_previousDetectedPulses, double samprate, double tstartRecord)
//int noDetect(gsl_vector *der, ReconstructInitSIRENA *reconstruct_init, int *numberPulses, gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl, gsl_vector **samp1DERgsl)
{
	string message="";
	char valERROR[256];
                
        gsl_vector_view temp;
        
        int sizeRecord = der->size;
        double sum_samp1DER;
        
        int pulse_length_ToSubtract = 100; // Instead of using pulse_length, because the results are equal
        gsl_vector *modelToSubtract = gsl_vector_alloc(pulse_length_ToSubtract);

	gsl_vector *tstartPulsei = gsl_vector_alloc(3);
	//gsl_vector_set(tstartPulsei,0,reconstruct_init->tstartPulse1);
        double tstartPulse1_seconds = gsl_vector_get(reconstruct_init->tstartPulse1_i,num_previousDetectedPulses);
        //std::cout << "tstartPulse1_seconds:" << std::setprecision(15) << tstartPulse1_seconds << '\n';
        //std::cout << "tstartRecord:" <<std::setprecision(15) << tstartRecord<< '\n';
        int tstartPulse1_samples;
        if ((strcmp(reconstruct_init->EnergyMethod,"I2R") == 0) || (strcmp(reconstruct_init->EnergyMethod,"I2RALL") == 0) || (strcmp(reconstruct_init->EnergyMethod,"I2RNOL") == 0) || (strcmp(reconstruct_init->EnergyMethod,"I2RFITTED") == 0))
        {
            //tstartPulse1_samples = ceil((tstartPulse1_seconds-tstartRecord)*samprate)-2;
            tstartPulse1_samples = ceil((tstartPulse1_seconds-tstartRecord)*samprate)-1;
        }
        else
        {
            //tstartPulse1_samples = ceil((tstartPulse1_seconds-tstartRecord)*samprate)-1;
            tstartPulse1_samples = ceil((tstartPulse1_seconds-tstartRecord)*samprate);
        }
        //std::cout << "(tstartPulse1_seconds-tstartRecord)*samprate: " << std::setprecision(15) << (tstartPulse1_seconds-tstartRecord)*samprate<< '\n';
        //std::cout << "tstartPulse1_samples: " << std::setprecision(15) << tstartPulse1_samples<< '\n';
        gsl_vector_set(tstartPulsei,0,tstartPulse1_samples);
	gsl_vector_set(tstartPulsei,1,reconstruct_init->tstartPulse2);
	gsl_vector_set(tstartPulsei,2,reconstruct_init->tstartPulse3);
        
        gsl_vector_set_zero(*flagTruncated);

	if (reconstruct_init->tstartPulse2 == 0) 	*numberPulses = 1;
	else if (reconstruct_init->tstartPulse3 == 0) 	*numberPulses = 2;
	else						*numberPulses = 3;

        for (int i=0;i<*numberPulses;i++)
        {
            gsl_vector_set(*tstartgsl,i,gsl_vector_get(tstartPulsei,i));
            
            if (i != *numberPulses-1) 	
            {
                if ((gsl_vector_get(tstartPulsei,i) < 0) || (gsl_vector_get(tstartPulsei,i) > sizeRecord-2)
                    || (gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i) < 1) || (gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i) > sizeRecord-gsl_vector_get(tstartPulsei,i)))
                {
                    sprintf(valERROR,"%d",__LINE__+5);
                    string str(valERROR);
                    message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                    str.clear();
                    EP_PRINT_ERROR(message,EPFAIL);
                }
                temp = gsl_vector_subvector(der,gsl_vector_get(tstartPulsei,i),gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i));
            }
            else	
            {
                if ((gsl_vector_get(tstartPulsei,i) < 0) || (gsl_vector_get(tstartPulsei,i) > sizeRecord-2)
                    || (sizeRecord-gsl_vector_get(tstartPulsei,i) < 1) || (sizeRecord-gsl_vector_get(tstartPulsei,i) > sizeRecord-gsl_vector_get(tstartPulsei,i)))
                {
                    sprintf(valERROR,"%d",__LINE__+5);
                    string str(valERROR);
                    message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                    str.clear();
                    EP_PRINT_ERROR(message,EPFAIL);
                }
                temp = gsl_vector_subvector(der,gsl_vector_get(tstartPulsei,i),sizeRecord-gsl_vector_get(tstartPulsei,i));
            }
            
            if (i == 0)
            {	
                
                gsl_vector_set(*maxDERgsl,i,gsl_vector_max(&temp.vector));
                gsl_vector_set(*samp1DERgsl,i,gsl_vector_get(der,gsl_vector_get(tstartPulsei,i)));
                
                // Average of the first 4 samples of the derivative
                sum_samp1DER = 0.0;
                for (int index_samp1DER=0;index_samp1DER<4;index_samp1DER++)
                {
                    sum_samp1DER = sum_samp1DER + gsl_vector_get(der,gsl_vector_get(tstartPulsei,i)+index_samp1DER);
                    
                }                                                                                   
                gsl_vector_set(*samp1DERgsl,i,sum_samp1DER/4.0);
                //cout<<"samp1DERgsl/4: "<<gsl_vector_get(*samp1DERgsl,i)<<endl;
            }
            else
            {
                if (find_model_maxDERs(gsl_vector_get(*maxDERgsl,i-1), reconstruct_init, &modelToSubtract))
                {
                    message = "Cannot run find_model routine for pulse i=" + boost::lexical_cast<std::string>(i) + " when newPulses = 1";
                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                }
                
                if ((gsl_vector_get(*tstartgsl,i-1) < 0) || (min(gsl_vector_get(*tstartgsl,i-1)+pulse_length_ToSubtract,(double) sizeRecord)-1 > der->size-1))
                {
                    sprintf(valERROR,"%d",__LINE__+7);
                    string str(valERROR);
                    message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                    str.clear();
                    EP_PRINT_ERROR(message,EPFAIL);
                }
                
                for (int j=gsl_vector_get(*tstartgsl,i-1);j<min(gsl_vector_get(*tstartgsl,i-1)+pulse_length_ToSubtract,(double) sizeRecord);j++)
                {
                    gsl_vector_set(der,j,gsl_vector_get(der,j)-gsl_vector_get(modelToSubtract,j-gsl_vector_get(*tstartgsl,i-1)));
                }
                
                if (i != *numberPulses-1)	
                {
                    if ((gsl_vector_get(tstartPulsei,i) < 0) || (gsl_vector_get(tstartPulsei,i) > der->size-2)
                        || (gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i) < 1) || (gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i) > der->size-gsl_vector_get(tstartPulsei,i)))
                    {
                        sprintf(valERROR,"%d",__LINE__+5);
                        string str(valERROR);
                        message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL);
                    }
                    temp = gsl_vector_subvector(der,gsl_vector_get(tstartPulsei,i),gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i));
                }
                else
                {
                    if ((gsl_vector_get(tstartPulsei,i) < 0) || (gsl_vector_get(tstartPulsei,i) > der->size-2)
                        || (sizeRecord-gsl_vector_get(tstartPulsei,i) < 1) || (sizeRecord-gsl_vector_get(tstartPulsei,i) > der->size-gsl_vector_get(tstartPulsei,i)))
                    {
                        sprintf(valERROR,"%d",__LINE__+5);
                        string str(valERROR);
                        message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                        str.clear();
                        EP_PRINT_ERROR(message,EPFAIL);
                    }
                    temp = gsl_vector_subvector(der,gsl_vector_get(tstartPulsei,i),sizeRecord-gsl_vector_get(tstartPulsei,i));
                }
                
                gsl_vector_set(*maxDERgsl,i,gsl_vector_max(&temp.vector));
                // Average of the first 4 samples of the derivative
                sum_samp1DER = 0.0;
                for (int index_samp1DER=0;index_samp1DER<4;index_samp1DER++)
                {
                    sum_samp1DER = sum_samp1DER + gsl_vector_get(&temp.vector,gsl_vector_get(tstartPulsei,i)+index_samp1DER);
                    
                }                                                                                   
                gsl_vector_set(*samp1DERgsl,i,sum_samp1DER/4.0);
            }
        }
        
        message.clear();
        
        return (EPOK);
}
/*xxxx end of SECTION 18 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
