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
   ESP2013-48637-C2-1-P and ESP2014-53672-C3-1-P.

/***********************************************************************
*                      PULSEPROCESS
*
*  File:       pulseprocess.cpp
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

*******************************************************************************/

#include "pulseprocess.h"

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
* 	sinc(f1)=0.6 where f1=1/(2pi*tau_fall)
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
* - tau_fall: Value related to the fall times of the pulses (seconds) (It's really tauFALL*sacaleFactor)
* - sampleRate: Sampling frequency (samples per second)
******************************************************************************/
int lpf_boxcar (gsl_vector **invector, int szVct, double tau_fall, int sampleRate)
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
	cutFreq = 2 * (1/(2*pi*tau_fall));	//According to Jan, sinc(f1)=0.6 where f1=1/(2pi*tau_fall)
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
	gsl_vector_free(invectorAux);
	gsl_vector_free(invectorAux1);
	
	if (boxLength == 1)	return(3);

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
	double data[size];	// Intermediate value to use 'gsl_stats_median_from_sorted_data'
	bool veryBig = false;

	// Mean
	for (int i=0;i<size;i++)
	{
		data[i] = gsl_vector_get(invector,i);
		if (data[i]>1e10)
		{
			veryBig =true;
			break;
		}
	}
	if (veryBig == false)
	{
		*mean = gsl_stats_mean(data, 1, size);
		// Standard deviation
		*sigma= gsl_stats_sd_m(data, 1, size, *mean);
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
	double data[size];					// Intermediate value to use 'gsl_stats_median_from_sorted_data'
	double median;

	// Median
	for (int i=0;i<size;i++)
	{
		data[i] = gsl_vector_get(invector,i);
	}
	gsl_sort(data,1,size);
	median = gsl_stats_median_from_sorted_data (data,1,size);

	gsl_vector_memcpy(invectorNew,invector);

	// Iterate until no points out of the maximum excursion (kappa*sigma)
	do
	{
		if ((size-boxLPF-1 < 1) || (size-boxLPF-1 >invectorNew->size))
		{
			sprintf(valERROR,"%d",__LINE__+5);
			string str(valERROR);
			message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
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

	gsl_vector_free(invectorNew);

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
****************************************/
int getB(gsl_vector *vectorin, gsl_vector *tstart, int nPulses, gsl_vector **lb, int sizepulse, gsl_vector **B)
{
	string message = "";
	char valERROR[256];

	// Declare variables
	// It is not necessary to check the allocation because 'tstart' size must already be > 0
	*B = gsl_vector_alloc(tstart->size);
	gsl_vector_set_all(*B,-999);
	double Baux = -999;
	double tendprev;

	// Auxiliary variables
	gsl_vector *input;
	gsl_vector_view temp;		// In order to handle with gsl_vector_view (subvectors)

	if (nPulses-1 >(tstart)->size-1)
	{
		sprintf(valERROR,"%d",__LINE__+7);
		string str(valERROR);
		message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL);
	}
	for (int i=0;i<nPulses;i++)
	{
		gsl_vector_set(tstart,i,(int)(gsl_vector_get(tstart,i)));	
	}

	for (int i=0;i<nPulses;i++)
	{
		if (i == 0)		// First pulse into a record
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
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				
				if ((gsl_vector_get(tstart,i)-gsl_vector_get(*lb,0) < 0) || (gsl_vector_get(tstart,i)-gsl_vector_get(*lb,0) > vectorin->size-2)
					|| (gsl_vector_get(*lb,0) < 1) || (gsl_vector_get(*lb,0) >vectorin->size-(gsl_vector_get(tstart,i)-gsl_vector_get(*lb,0))))
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL);
				}
				temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,i)-gsl_vector_get(*lb,0),gsl_vector_get(*lb,0));
				if (gsl_vector_memcpy(input, &temp.vector) != 0)
				{
					sprintf(valERROR,"%d",__LINE__-2);
					string str(valERROR);
					message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL);
				}

				// Sum all the elements of 'input'
				if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
				{
					message = "Cannot run gsl_vector_Sumsubvector routine when tstart>=lb";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				gsl_vector_free(input);
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
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
			  
				if ((gsl_vector_get(tstart,0) < 1) || (gsl_vector_get(tstart,0) >vectorin->size))
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL);
				}
				temp = gsl_vector_subvector(vectorin,0,gsl_vector_get(tstart,0));
				if (gsl_vector_memcpy(input, &temp.vector) != 0)
				{
					sprintf(valERROR,"%d",__LINE__-2);
					string str(valERROR);
					message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL);
				}

				// Sum all the elements of 'input'
				if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
				{
					message = "Cannot run gsl_vector_Sumsubvector routine when tstart<lb";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}				
				gsl_vector_set(*lb,0, gsl_vector_get(tstart,0));
				gsl_vector_free(input);
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
								EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
							}
				
							if ((gsl_vector_get(tstart,j)-gsl_vector_get(*lb,0) < 0) || (gsl_vector_get(tstart,j)-gsl_vector_get(*lb,0) > vectorin->size-2)
								|| (gsl_vector_get(*lb,0) < 1) || (gsl_vector_get(*lb,0) >vectorin->size-(gsl_vector_get(tstart,j)-gsl_vector_get(*lb,0))))
							{
								sprintf(valERROR,"%d",__LINE__+5);
								string str(valERROR);
								message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,j)-gsl_vector_get(*lb,0),gsl_vector_get(*lb,0));
							if (gsl_vector_memcpy(input, &temp.vector) != 0)
							{
								sprintf(valERROR,"%d",__LINE__-2);
								string str(valERROR);
								message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							
							// Sum all the elements of input
							if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
							{
								message = "Cannot run gsl_vector_Sumsubvector routine when no pulse free interval before the pulse";
								EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
							}
							gsl_vector_free(input);

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
								EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
							}
						  
							if ((tendprev+1 < 0) || (tendprev+1 > vectorin->size-2)
								|| (gsl_vector_get(tstart,j)-tendprev-1 < 1) || (gsl_vector_get(tstart,j)-tendprev-1 >vectorin->size-(tendprev+1)))
							{
								sprintf(valERROR,"%d",__LINE__+5);
								string str(valERROR);
								message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,j)-tendprev-1);
							if (gsl_vector_memcpy(input, &temp.vector) != 0)
							{
								sprintf(valERROR,"%d",__LINE__-2);
								string str(valERROR);
								message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							gsl_vector_set(*lb,0,gsl_vector_get(tstart,j)-tendprev-1);

							// Sum all the elements of input
							if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
							{
								message = "Cannot run gsl_vector_Sumsubvector routine when no pulse free interval before the pulse";
								EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
							}
							gsl_vector_free(input);

							break;
						}
						else
						{
							for (int j=i;j<nPulses;j++)	// From the current pulse
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
									if (gsl_vector_get(tstart,j+1)-tendprev >= gsl_vector_get(*lb,i))
									{
										if ((input = gsl_vector_alloc(gsl_vector_get(*lb,i))) == 0)
										{
											sprintf(valERROR,"%d",__LINE__-2);
											string str(valERROR);
											message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
											EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
										}
				
										if ((gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i) < 0) || (gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i) > vectorin->size-2)
											|| (gsl_vector_get(*lb,i) < 1) || (gsl_vector_get(*lb,i) > vectorin->size-(gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i))))
										{
											sprintf(valERROR,"%d",__LINE__+5);
											string str(valERROR);
											message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
											EP_PRINT_ERROR(message,EPFAIL);
										}
										temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i),gsl_vector_get(*lb,i));
										if (gsl_vector_memcpy(input, &temp.vector) != 0)
										{
											sprintf(valERROR,"%d",__LINE__-2);
											string str(valERROR);
											message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
											EP_PRINT_ERROR(message,EPFAIL);
										}

										// Sum all the elements of input
										if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
										{
											message = "Cannot run gsl_vector_Sumsubvector routine when no first pulse in row & tstart-tendprev >= lb";
											EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
										}
										gsl_vector_free(input);

										break;
									}
									else if ((gsl_vector_get(tstart,j+1)-tendprev < gsl_vector_get(*lb,i)) && ((gsl_vector_get(tstart,j+1)-tendprev>1)))
									{
										if ((input = gsl_vector_alloc(gsl_vector_get(tstart,j+1)-tendprev-1)) == 0)
										{
											sprintf(valERROR,"%d",__LINE__-2);
											string str(valERROR);
											message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
											EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
										}
										
										if ((tendprev+1 < 0) || (tendprev+1 > vectorin->size-2)
											|| (gsl_vector_get(tstart,j+1)-tendprev-1 < 1) || (gsl_vector_get(tstart,j+1)-tendprev-1 > vectorin->size-(tendprev+1)))
										{
											sprintf(valERROR,"%d",__LINE__+5);
											string str(valERROR);
											message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
											EP_PRINT_ERROR(message,EPFAIL);
										}
										temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,j+1)-tendprev-1);
										if (gsl_vector_memcpy(input, &temp.vector) != 0)
										{
											sprintf(valERROR,"%d",__LINE__-2);
											string str(valERROR);
											message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
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
											EP_PRINT_ERROR(message,EPFAIL);
										}
										gsl_vector_set(*lb,i,gsl_vector_get(tstart,j+1)-tendprev-1);
										gsl_vector_free(input);

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
											EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
										}
										
										if ((tendprev < 0) || (tendprev > vectorin->size-2)
											|| (gsl_vector_get(*lb,i) < 1) || (gsl_vector_get(*lb,i) > vectorin->size-tendprev))
										{
											sprintf(valERROR,"%d",__LINE__+5);
											string str(valERROR);
											message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
											EP_PRINT_ERROR(message,EPFAIL);
										}
										temp = gsl_vector_subvector(vectorin,tendprev,gsl_vector_get(*lb,i));
										if (gsl_vector_memcpy(input, &temp.vector) != 0)
										{
											sprintf(valERROR,"%d",__LINE__-2);
											string str(valERROR);
											message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
											EP_PRINT_ERROR(message,EPFAIL);
										}

										// Sum all the elements of input
										if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
										{
											message = "Cannot run gsl_vector_Sumsubvector routine when no pulse-free interval before pulse & last pulse in row";
											EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
										}
										gsl_vector_free(input);

										break;
									}
									else if ((vectorin->size-tendprev < gsl_vector_get(*lb,i)) && (vectorin->size-tendprev > 1))
									{
										if ((input = gsl_vector_alloc(vectorin->size-tendprev-1)) == 0)
										{
											sprintf(valERROR,"%d",__LINE__-2);
											string str(valERROR);
											message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
											EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
										}
										
										if ((tendprev+1 < 0) || (tendprev+1 > vectorin->size-2)
											|| (vectorin->size-tendprev-1 < 1))
										{
											sprintf(valERROR,"%d",__LINE__+5);
											string str(valERROR);
											message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
											EP_PRINT_ERROR(message,EPFAIL);
										}
										temp = gsl_vector_subvector(vectorin,tendprev+1,vectorin->size-tendprev-1);
										if (gsl_vector_memcpy(input, &temp.vector) != 0)
										{
											sprintf(valERROR,"%d",__LINE__-2);
											string str(valERROR);
											message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
											EP_PRINT_ERROR(message,EPFAIL);
										}
										if ((i < 0) || (i >(*lb)->size-1))
										{
											sprintf(valERROR,"%d",__LINE__+5);
											string str(valERROR);
											message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
											EP_PRINT_ERROR(message,EPFAIL);
										}
										gsl_vector_set(*lb,i,vectorin->size-tendprev-1);

										// Sum all the elements of input
										if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
										{
											message = "Cannot run gsl_vector_Sumsubvector routine when no pulse-free interval before pulse & last pulse in row";
											EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
										}
										gsl_vector_free(input);

										break;
									}
								}
							}
						}
					}
				}
			}
		}
		else	// Not first pulse into a record
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
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
			  
				if ((gsl_vector_get(tstart,i)-gsl_vector_get(*lb,i) < 0) || (gsl_vector_get(tstart,i)-gsl_vector_get(*lb,i) > vectorin->size-2)
					|| (gsl_vector_get(*lb,i) < 1) || (gsl_vector_get(*lb,i) > vectorin->size-(gsl_vector_get(tstart,i)-gsl_vector_get(*lb,i))))
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL);
				}
				temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,i)-gsl_vector_get(*lb,i),gsl_vector_get(*lb,i));
				if (gsl_vector_memcpy(input, &temp.vector) != 0)
				{
					sprintf(valERROR,"%d",__LINE__-2);
					string str(valERROR);
					message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL);
				}

				// Sum all the elements of input
				if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
				{
					message = "Cannot run gsl_vector_Sumsubvector routine length_(2)>=lb";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				gsl_vector_free(input);
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
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				
				if ((tendprev+1 < 0) || (tendprev+1 > vectorin->size-2)
					|| (gsl_vector_get(tstart,i)-tendprev-1 < 1) || (gsl_vector_get(tstart,i)-tendprev-1 > vectorin->size-(tendprev+1)))
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL);
				}
				temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,i)-tendprev-1);
				if (gsl_vector_memcpy(input, &temp.vector) != 0)
				{
					sprintf(valERROR,"%d",__LINE__-2);
					string str(valERROR);
					message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
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
					EP_PRINT_ERROR(message,EPFAIL);
				}
				gsl_vector_set(*lb,i,gsl_vector_get(tstart,i)-tendprev-1);
				gsl_vector_free(input);
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
								EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
							}
						  
							if ((gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i) < 0) || (gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i) > vectorin->size-2)
								|| (gsl_vector_get(*lb,i) < 1) || (gsl_vector_get(*lb,i) > vectorin->size-(gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i))))
							{
								sprintf(valERROR,"%d",__LINE__+5);
								string str(valERROR);
								message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i),gsl_vector_get(*lb,i));
							if (gsl_vector_memcpy(input, &temp.vector) != 0)
							{
								sprintf(valERROR,"%d",__LINE__-2);
								string str(valERROR);
								message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}

							// Sum all the elements of input
							if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
							{
								message = "Cannot run gsl_vector_Sumsubvector routine when no first pulse in row & tstart-tendprev >= lb";
								EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
							}
							gsl_vector_free(input);

							break;
						}
						else if ((gsl_vector_get(tstart,j+1)-tendprev < gsl_vector_get(*lb,i)) && ((gsl_vector_get(tstart,j+1)-tendprev>1)))
						{
							if ((input = gsl_vector_alloc(gsl_vector_get(tstart,j+1)-tendprev-1)) == 0)
							{
								sprintf(valERROR,"%d",__LINE__-2);
								string str(valERROR);
								message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
							}
							
							if ((tendprev+1 < 0) || (tendprev+1 > vectorin->size-2)
								|| (gsl_vector_get(tstart,j+1)-tendprev-1 < 1) || (gsl_vector_get(tstart,j+1)-tendprev-1 > vectorin->size-(tendprev+1)))
							{
								sprintf(valERROR,"%d",__LINE__+5);
								string str(valERROR);
								message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,j+1)-tendprev-1);
							if (gsl_vector_memcpy(input, &temp.vector) != 0)
							{
								sprintf(valERROR,"%d",__LINE__-2);
								string str(valERROR);
								message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
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
								EP_PRINT_ERROR(message,EPFAIL);
							}
							gsl_vector_set(*lb,i,gsl_vector_get(tstart,j+1)-tendprev-1);
							gsl_vector_free(input);

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
								EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
							}
							
							if ((tendprev < 0) || (tendprev > vectorin->size-2)
								|| (gsl_vector_get(*lb,i) < 1) || (gsl_vector_get(*lb,i) > vectorin->size-tendprev))
							{
								sprintf(valERROR,"%d",__LINE__+5);
								string str(valERROR);
								message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							temp = gsl_vector_subvector(vectorin,tendprev,gsl_vector_get(*lb,i));
							if (gsl_vector_memcpy(input, &temp.vector) != 0)
							{
								sprintf(valERROR,"%d",__LINE__-2);
								string str(valERROR);
								message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}

							// Sum all the elements of input
							if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
							{
								message = "Cannot run gsl_vector_Sumsubvector routine when no pulse-free interval before pulse & last pulse in row";
								EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
							}
							gsl_vector_free(input);

							break;
						}
						else if ((vectorin->size-tendprev < gsl_vector_get(*lb,i)) && (vectorin->size-tendprev > 1))
						{
							if ((input = gsl_vector_alloc(vectorin->size-tendprev-1)) == 0)
							{
								sprintf(valERROR,"%d",__LINE__-2);
								string str(valERROR);
								message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
							}
							
							if ((tendprev+1 < 0) || (tendprev+1 > vectorin->size-2)
								|| (vectorin->size-tendprev-1 < 1))
							{
								sprintf(valERROR,"%d",__LINE__+5);
								string str(valERROR);
								message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							temp = gsl_vector_subvector(vectorin,tendprev+1,vectorin->size-tendprev-1);
							if (gsl_vector_memcpy(input, &temp.vector) != 0)
							{
								sprintf(valERROR,"%d",__LINE__-2);
								string str(valERROR);
								message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							if ((i < 0) || (i >(*lb)->size-1))
							{
								sprintf(valERROR,"%d",__LINE__+5);
								string str(valERROR);
								message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							gsl_vector_set(*lb,i,vectorin->size-tendprev-1);

							// Sum all the elements of input
							if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
							{
								message = "Cannot run gsl_vector_Sumsubvector routine when no pulse-free interval before pulse & last pulse in row";
								EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
							}
							gsl_vector_free(input);

							break;
						}
						else if ((vectorin->size-tendprev < gsl_vector_get(*lb,i)) && (vectorin->size-tendprev <= 1))
						{
							if ((i < 0) || (i >(*lb)->size-1))
							{
								sprintf(valERROR,"%d",__LINE__+5);
								string str(valERROR);
								message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
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
			EP_PRINT_ERROR(message,EPFAIL);
		}
		gsl_vector_set(*B,i,Baux);
	}

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
* - sizepulse: Size of the pulse in samples, ntaus * tauFALL in samples (equal to 'sizePulse_b' global variable)
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
			EP_PRINT_ERROR(message,EPFAIL);
		}

		// Apply the running sum filter
		if (RS_filter (input, lrs, lb, B, &ph))
		{
			message = "Cannot run RS_filter routine when size-tstart>size-tend";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		*pulseheight = ph;

		gsl_vector_free(input);
	}

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
	if ((modelFound_aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration)) == 0)
	{
		sprintf(valERROR,"%d",__LINE__-2);
		string str(valERROR);
	        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL);
	}

	long nummodels = reconstruct_init->library_collection->ntemplates;

	if (energy < gsl_vector_get(reconstruct_init->library_collection->energies,0))
	{
		gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_B0[0].ptemplate);
		gsl_vector_scale(modelFound_aux,energy/gsl_vector_get(reconstruct_init->library_collection->energies,0));
	}
	else if (energy > gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1))
	{
		gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_B0[nummodels-1].ptemplate);
		gsl_vector_scale(modelFound_aux,energy/gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1));
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
			if (fabs(energy-gsl_vector_get(reconstruct_init->library_collection->energies,i))<1e-6)
			{
				gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_B0[i].ptemplate);
				gsl_vector_scale(modelFound_aux,energy/gsl_vector_get(reconstruct_init->library_collection->energies,i));

				break;
			}
			else if ((energy > gsl_vector_get(reconstruct_init->library_collection->energies,i)) && (energy < gsl_vector_get(reconstruct_init->library_collection->energies,i+1)))
			{
				// Interpolate between the two corresponding rows in "models"
				// It is not necessary to check the allocation because the same kind of allocation has been checked previously
				gsl_vector *modelAux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);
				gsl_vector_set_zero(modelAux);

				if (interpolate_model(&modelAux,energy,reconstruct_init->library_collection->pulse_templates_B0[i].ptemplate,gsl_vector_get(reconstruct_init->library_collection->energies,i),
					reconstruct_init->library_collection->pulse_templates_B0[i+1].ptemplate,gsl_vector_get(reconstruct_init->library_collection->energies,i+1)))
				{
					message = "Cannot run interpolate_model with two rows in models";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				gsl_vector_memcpy(modelFound_aux,modelAux);
				gsl_vector_free(modelAux);

				break;
			}
		}
	}

	gsl_vector_view temp;
	temp = gsl_vector_subvector(modelFound_aux,0,(*modelFound)->size);
	gsl_vector_memcpy(*modelFound,&temp.vector);

	gsl_vector_free(modelFound_aux);

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

	gsl_vector *modelFound_aux;
	if ((modelFound_aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration)) == 0)
	{
		sprintf(valERROR,"%d",__LINE__-2);
		string str(valERROR);
	        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL);
	}

	long nummodels = reconstruct_init->library_collection->ntemplates;

	if (maxDER < gsl_vector_get(reconstruct_init->library_collection->maxDERs,0))
	{
		gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_filder[0].ptemplate);
		gsl_vector_scale(modelFound_aux,maxDER/gsl_vector_max(modelFound_aux));
	}
	else if (maxDER > gsl_vector_get(reconstruct_init->library_collection->maxDERs,nummodels-1))
	{
		gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_filder[nummodels-1].ptemplate);
		gsl_vector_scale(modelFound_aux,maxDER/gsl_vector_max(modelFound_aux));
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
			if (fabs(maxDER-gsl_vector_get(reconstruct_init->library_collection->maxDERs,i))<1e-6)
			{
				gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_filder[i].ptemplate);
				gsl_vector_scale(modelFound_aux,maxDER/gsl_vector_max(modelFound_aux));

				break;
			}
			else if ((maxDER > gsl_vector_get(reconstruct_init->library_collection->maxDERs,i)) && (maxDER < gsl_vector_get(reconstruct_init->library_collection->maxDERs,i+1)))
			{
				// Interpolate between the two corresponding rows in "models"
				gsl_vector_set_zero(modelFound_aux);

				if (interpolate_model(&modelFound_aux,maxDER,reconstruct_init->library_collection->pulse_templates_filder[i].ptemplate,gsl_vector_get(reconstruct_init->library_collection->maxDERs,i),
						reconstruct_init->library_collection->pulse_templates_filder[i+1].ptemplate,gsl_vector_get(reconstruct_init->library_collection->maxDERs,i+1)))
				{
					message = "Cannot run interpolate_model with two rows in models";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}

				break;
			}
		}
	}

	gsl_vector_view temp;
	temp = gsl_vector_subvector(modelFound_aux,0,(*modelFound)->size);
	gsl_vector_memcpy(*modelFound,&temp.vector);

	gsl_vector_free(modelFound_aux);

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
	
	gsl_vector *modelFound_aux;
	if ((modelFound_aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration)) == 0)
	{
		sprintf(valERROR,"%d",__LINE__-2);
		string str(valERROR);
	        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL);
	}

	long nummodels = reconstruct_init->library_collection->ntemplates;

	if (samp1DER < gsl_vector_get(reconstruct_init->library_collection->samp1DERs,0))
	{
		gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_filder[0].ptemplate);
		gsl_vector_scale(modelFound_aux,samp1DER/gsl_vector_get(modelFound_aux,0));
	}
	else if (samp1DER > gsl_vector_get(reconstruct_init->library_collection->samp1DERs,nummodels-1))
	{
		gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_filder[nummodels-1].ptemplate);
		gsl_vector_scale(modelFound_aux,samp1DER/gsl_vector_get(modelFound_aux,nummodels-1));
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
			if ((samp1DER >= gsl_vector_get(reconstruct_init->library_collection->samp1DERs,i)) && (samp1DER < gsl_vector_get(reconstruct_init->library_collection->samp1DERs,i+1)))
			{
				// Interpolate between the two corresponding rows in "models"
				gsl_vector_set_zero(modelFound_aux);

				if (interpolate_model(&modelFound_aux,samp1DER,reconstruct_init->library_collection->pulse_templates_filder[i].ptemplate,gsl_vector_get(reconstruct_init->library_collection->samp1DERs,i),
					reconstruct_init->library_collection->pulse_templates_filder[i+1].ptemplate,gsl_vector_get(reconstruct_init->library_collection->samp1DERs,i+1)))
				{
					message = "Cannot run interpolate_model with two rows in models";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}

				break;
			}
		}
	}

	gsl_vector_view temp;
	temp = gsl_vector_subvector(modelFound_aux,0,(*modelFound)->size);
	gsl_vector_memcpy(*modelFound,&temp.vector);

	gsl_vector_free(modelFound_aux);

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
*              E2-E1          E2-E1
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
	gsl_vector_set_zero(*modelFound);
	
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
	gsl_vector_add(*modelFound,modelIn1Aux);
	gsl_vector_add(*modelFound,modelIn2Aux);
	
	gsl_vector_free(modelIn1Aux);
	gsl_vector_free(modelIn2Aux);

    return(EPOK);
}
/*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 11 ************************************************************
* findPulsesCAL: This function is going to find the pulses in a record (in the CALibration mode) by using the function findTstartCAL
*
* - Declare variables
* - First step to look for single pulses (call 'medianKappaClipping')
*   	- Call 'findTstartCAL'
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
* - taufall: Fall time of the pulses (seconds)
* - scalefactor: Scale factor to apply to the fall time of the pulses in order to calculate the LPF box-car length
* - samplingRate: Sampling rate
* - samplesup: Number of consecutive samples over the threshold to locate a pulse
* - nsgms: Number of Sigmas to establish the threshold
* - lb: Vector containing the baseline averaging length used for each pulse
* - lrs: Running sum length (equal to the 'Lrs' global_variable)
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'maxPulsesPerRecord' and 'pulse_length'
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

	double taufall,
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

	// First step to look for single pulses
	if (medianKappaClipping (vectorinDER, kappamkc, stopcriteriamkc, nsgms, (int)(pi*samplingRate*taufall*scalefactor), &thresholdmediankappa))
	{
		message = "Cannot run medianKappaClipping looking for single pulses";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	*threshold = thresholdmediankappa;

	if (findTstartCAL (reconstruct_init->maxPulsesPerRecord, vectorinDER, thresholdmediankappa, samplesUp, reconstruct_init, nPulses, tstart, quality, maxDERgsl))
	{
		message = "Cannot run findTstartCAL";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	if (*nPulses != 0)
	{
		// It is not necessary to check the allocation because 'reconstruct_init->maxPulsesPerRecord'='EventListSize'(input parameter) must already be > 0
		gsl_vector *Lbgsl = gsl_vector_alloc(reconstruct_init->maxPulsesPerRecord);	// If there is no free-pulses segments longer than Lb=>
		gsl_vector_set_all(Lbgsl,lb);                                  			// segments shorter than Lb will be useed and its length (< Lb)
		                                                               			// must be used instead Lb in RS_filter
		gsl_vector *Bgsl;

		if (getB(vectorin, *tstart, *nPulses, &Lbgsl, reconstruct_init->pulse_length, &Bgsl))
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

		gsl_vector_free(Lbgsl);
		gsl_vector_free(Bgsl);
	}
	return(EPOK);
}
/*xxxx end of SECTION 11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 12 ************************************************************
* findTstartCAL function: This function finds the pulses tstarts in the input vector (first derivative of the filtered record)
*
* This function scans all values the derivative of the (low-pass filtered) record until it finds nSamplesUp consecutive values (due to the noise more than 1 value is
* required) over the threshold. To look for more pulses, it finds nSamplesUp consecutive values
* (due to the noise) below the threshold and then, it starts to scan again.
*
* In general, the tstarts found in the derivative do not coincide with the tstart of the not filtered record
* (by the spreading caused by the low-pass filtered).
*
* - Declare variables
* - Allocate GSL vectors
* - It is possible to find the tstarts...
* 	- Obtain tstart of each pulse in the derivative:
* 		- If der_i>threshold and prevPulse=false, it looks for nSamplesUp consecutive samples over the threshold
*   			- If not, it looks again for a pulse crossing over the threshold
*	        	- If yes, a pulse is found (truncated if it is at the beginning)
*	   	- If der_i>threshold and prevPulse=true, it looks for a sample below the threshold
*   	  		- If not, it looks again for a sample below the threshold
*       		- If yes, it looks for nSamplesUp consecutive samples below the threshold and again it starts to look for a pulse
*	- Indentify if there is a truncated pulse at the end (it has not found its tend)
*	- Obtain tstart of each pulse in the not filtered signal??? CURRENTLY NOT USED
* - ...or to use the tstart provided as input parameters
* 	- Obtain the maxDERs of the pulses whose tstarts have been provided
*
* Parameters:
* - maxPulsesPerRecord: Expected maximum number of pulses per record in order to not allocate the GSL variables with the size of the input vector
* - der: First derivative of the (low-pass filtered) record
* - adaptativethreshold: Threshold
* - nSamplesUp: Number of consecutive samples over the threshold to 'find' a pulse ('samplesUp')
* - allPulsesMode: 0-> If it finds a pulse the function returns
*                  1-> It finds all pulses of the record
* - sampling: Sampling rate => Not neccessary => To be deleted
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
* - numberPulses: Number of found pulses
* - thereIsPulse: 0 -> The algorithm has not found any pulse
*                 1 -> The algorithm has found unless one pulse
* - tstartgsl: Pulses tstart (in samples)
* - flagTruncated: Flag indicating if the pulse is truncated (inside this function only initial truncated pulses are classified)
* - maxDERgsl: Maximum of the first derivative of the low-pass filtered record inside each found pulse
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
	*tstartgsl = gsl_vector_alloc(maxPulsesPerRecord);		// tstarts referred to the not filtered event
	*flagTruncated = gsl_vector_alloc(maxPulsesPerRecord);
	gsl_vector_set_zero(*flagTruncated);
	*maxDERgsl = gsl_vector_alloc(maxPulsesPerRecord);		// Maximum of the first derivative
	gsl_vector_set_all(*maxDERgsl,-1E3);

	// To provide the tstarts (or not)
	bool findTstarts = true;
	if (reconstruct_init->tstartPulse1 != 0)	findTstarts = false;
	
	//for (int i=990; i<(990+30);i++) cout<<i<<" "<<gsl_vector_get(der,i)<<endl;
	//for (int i=20995; i<(20995+10);i++) cout<<i<<" "<<gsl_vector_get(der,i)<<endl;
	//cout<<"threshold: "<<adaptativethreshold<<endl;

	int cntUp = 0;
	int cntDown = 0;
	double possibleTstart;
	double possiblemaxDER;
	
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
								EP_PRINT_ERROR(message,EPFAIL);
							}
							gsl_vector_set(*maxDERgsl,*numberPulses,possiblemaxDER);
							gsl_vector_set(*tstartgsl,*numberPulses,possibleTstart);
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
					else if  (gsl_vector_get(der,i) < adaptativethreshold)
					{
						cntUp = 0;
						cntDown++;
						if (cntDown == nSamplesUp) // nSamplesUp samples under the threshold in order to look for another pulse
						{
							foundPulse = false; 
							cntDown = 0;
						}
					}
					else if  (gsl_vector_get(der,i) > adaptativethreshold)
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
		gsl_vector_set(tstartPulsei,0,reconstruct_init->tstartPulse1);
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
					EP_PRINT_ERROR(message,EPFAIL);
				}
				temp = gsl_vector_subvector(der,gsl_vector_get(tstartPulsei,i),szRw-gsl_vector_get(tstartPulsei,i));
			}
			
			gsl_vector_set(*maxDERgsl,i,gsl_vector_max(&temp.vector));
		}

		gsl_vector_free(tstartPulsei);
		gsl_vector_free(model);
	}

	return (EPOK);
}
/*xxxx end of SECTION 12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 13 ************************************************************
* InitialTriggering function: This function is going to find the first pulse in the input vector (first derivative of the filtered record)
*
* - Declare variables
* - Stablish the threshold
* - It is possible to find the tstart of the first pulse...
*   Obtain tstart of the first pulse in the derivative. If der_i>threshold, it looks for nSamplesUp consecutive samples over the threshold
* 	- If not, it looks again for a pulse crossing over the threshold
*       - If yes, a pulse is found (truncated if it is at the beginning)
* - ...or to use the tstart provided as input parameter
*
* Parameters:
* - derivative: First derivative of the (low-pass filtered) record
* - nSamplesUp: Number of consecutive samples over the threshold to 'find' an event ('samplesUp')
* - nSgms: Number of Sigmas to establish the threshold 
* - taufall: Fall time of the pulses (seconds) ('tauFall')
* - scalefactor: Scale factor to apply to the fall time of the pulses in order to calculate the LPF box-car length ('scaleFactor')
* - samplingRate: Sampling rate
* - stopCriteriamkc: Used in medianKappaClipping (%)
* - kappamkc: Used in medianKappaClipping
* - triggerCondition: true -> the algorithm has not found any event
*                     false -> the algorithm has found the first event
* - tstart: First event tstart (in samples)
* - flagTruncated: Flag indicating if the event is truncated (inside this function only initial truncated events are classified)
* - threshold: Calculated threshold  (output parameter because it is necessary out of the function)
* - tstartProvided: Tstart of the first pulse provided as input parameter
****************************************/
int InitialTriggering
(
	gsl_vector *derivative,

	int nSamplesUp,
	double nSgms,

	double taufall,
	double scalefactor,
	double samplingRate,

	double stopcriteriamkc,
	double kappamkc,

	bool *triggerCondition,
	int *tstart,
	int *flagTruncated,

	double *threshold,

	int tstartProvided)
{
	string message = "";

	// Declare variables
	const double pi = 4.0 * atan(1.0);
	int sizeRecord = derivative->size;	// Size of segment to process
	int cntUp = 0;
	int i = 0;
	int possibleTstart;

	*triggerCondition = false;

	// Stablish the threshold
	if (medianKappaClipping (derivative, kappamkc, stopcriteriamkc, nSgms, (int)(pi*samplingRate*taufall*scalefactor), threshold))
	{
		message = "Cannot run medianKappaClipping doing the initial triggering";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	if (tstartProvided == 0)
	{
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
	}
	else
	{
		*triggerCondition = true;
		*tstart = tstartProvided;
		if (*tstart == 0)	*flagTruncated = 1;
		else			*flagTruncated = 0;
	}

	return(EPOK);
}
/*xxxx end of SECTION 13 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 14 ************************************************************
* FindSecondaries function: This function runs after InitialTriggering to find all the events (except the first one) in the record first derivative 
*                           of the (filtered) record
*
* - Declare variables
* - It is possible to find the tstarts...
* - It looks for an event and if a pulse is found (foundPulse==true), it looks for another event
* 	- It looks for an event since the beginning (or the previous event) to the end of the record
*         The first condition to detect an event is that the adjustedDerivative was over the threshold
* 	  In general, nSamplesUp over the threshold to 'find' an event
* 	- Subtract the model from the adjusted derivative
* 		- Select the model of the found event from the libary by using the maximum of the derivative (maxDER)
* 		- Subtract
* - ...or to use the tstart provided as input parameters
* 	- Obtain the maxDERs of the events whose tstarts have been provided
* - Free allocated GSL vectors
*
* Parameters:
* - maxPulsesPerRecord: Expected maximum number of events per record in order to not allocate the GSL variables with the size of the input vector ('EventListSize')
*
* - adjustedDerivative: First derivative of the (low-pass filtered) record 
* - adaptativethreshold: Threshold
* - nSamplesUp: Number of consecutive samples over the threshold to 'find' an event or pulse ('samplesUp')
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses some parameters ('pulse_length', 'tstartPulsex',...), the templates...
* - tstartFirstEvent: Tstart of the first event of the record (in samples) found by 'InitialTriggering'
* - numberPulses: Number of found events
* - tstartgsl: Events tstart (in samples)
* - flagTruncated: Flag indicating if the event is truncated (inside this function only initial truncated pulses are classified)
* - maxDERgsl: Maximum of the first derivative of the (low-pass filtered) record inside each found event
****************************************/
int FindSecondaries
(
	int maxPulsesPerRecord,

	gsl_vector *adjustedDerivative,
	double adaptativethreshold,

	int nSamplesUp,

	ReconstructInitSIRENA *reconstruct_init,

	int tstartFirstEvent,

	int *numberPulses,
	gsl_vector **tstartgsl,
	gsl_vector **flagTruncated,
	gsl_vector **maxDERgsl,
	gsl_vector **samp1DERgsl)
{
	string message = "";
	char valERROR[256];

	// Declare variables
	bool foundPulse = false;
	int sizeRecord = adjustedDerivative->size;		// Size of segment to process
	*numberPulses = 0;
	gsl_vector_set_all(*maxDERgsl,-1E3);
	gsl_vector_set_all(*tstartgsl,-1E3);
	int i = tstartFirstEvent;
	// It is not necessary to check the allocation because 'reconstruct_init->pulse_length'='PulseLength'(input parameter) has been checked previously
	gsl_vector *model = gsl_vector_alloc(reconstruct_init->pulse_length);
	
	// To provide the tstarts (or not)
	bool findTstarts = true;
	if (reconstruct_init->tstartPulse1 != 0)	findTstarts = false;
	
	//for (int i=990; i<(990+25);i++) cout<<i<<" "<<gsl_vector_get(adjustedDerivative,i)<<endl;
	//for (int i=20990; i<(20990+25);i++) cout<<i<<" "<<gsl_vector_get(adjustedDerivative,i)<<endl;
	//adaptativethreshold = 40.0;
	//cout<<"Threshold: "<<adaptativethreshold<<endl;
	
	// It is not necessary to check the allocation because 'maxPulsesPerRecord'='EventListSize'(input parameter) must already be > 0
	gsl_vector *index_maxDERgsl = gsl_vector_alloc(maxPulsesPerRecord);
	
	if (findTstarts == true)
	{
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
					{
						if (*numberPulses == (*maxDERgsl)->size)
						{
							sprintf(valERROR,"%d",__LINE__+5);
							string str(valERROR);
							message = "Found pulses in record>'EventListSize'(input parameter) => Change EventListSize or check if the threshold is too low => Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
							EP_PRINT_ERROR(message,EPFAIL);
						}
						gsl_vector_set(*maxDERgsl,*numberPulses,gsl_vector_get(adjustedDerivative,i));
						gsl_vector_set(index_maxDERgsl,*numberPulses,i);
						gsl_vector_set(*samp1DERgsl,*numberPulses,gsl_vector_get(adjustedDerivative,i));
						gsl_vector_set(*tstartgsl,*numberPulses,i);
						if (i == 0)	gsl_vector_set(*flagTruncated,*numberPulses,1);
						*numberPulses = *numberPulses +1;
						foundPulse = true;
						//cout<<"Supera el umbral en "<<i<<" "<<gsl_vector_get(adjustedDerivative,i)<<endl;
					}
					i++;
				}
				else
				{
					if (gsl_vector_get(adjustedDerivative,i) > gsl_vector_get(*maxDERgsl,*numberPulses-1))
					{
						gsl_vector_set(*maxDERgsl,*numberPulses-1,gsl_vector_get(adjustedDerivative,i));
						gsl_vector_set(index_maxDERgsl,*numberPulses-1,i);
						//cout<<"Nuevo maxDER "<<gsl_vector_get(adjustedDerivative,i)<<endl;
						
						i++;
					}
					else
					{
						if (((strcmp(reconstruct_init->EnergyMethod,"I2RBISALL") == 0) && (gsl_vector_get(index_maxDERgsl,*numberPulses-1)-gsl_vector_get(*tstartgsl,*numberPulses-1) >= 0)) ||
							((strcmp(reconstruct_init->EnergyMethod,"I2RBISALL") != 0) && (gsl_vector_get(index_maxDERgsl,*numberPulses-1)-gsl_vector_get(*tstartgsl,*numberPulses-1) > 0)))
						{
							// Select the model of the found pulse from the libary by using the 1st sample of the derivative (samp1DER)
							if (find_model_samp1DERs(gsl_vector_get(*samp1DERgsl,*numberPulses-1), reconstruct_init, &model))
							{
								message = "Cannot run find_model_samp1DERs routine";
								EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
							}

							if ((gsl_vector_get(*tstartgsl,*numberPulses-1) < 0) || (min(gsl_vector_get(*tstartgsl,*numberPulses-1)+reconstruct_init->pulse_length,(double) sizeRecord)-1 > adjustedDerivative->size-1))
							{
								sprintf(valERROR,"%d",__LINE__+9);
								string str(valERROR);
								message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							for (int j=gsl_vector_get(*tstartgsl,*numberPulses-1);j<min(gsl_vector_get(*tstartgsl,*numberPulses-1)+reconstruct_init->pulse_length,(double) sizeRecord);j++)
							{
								//if (j < gsl_vector_get(*tstartgsl,*numberPulses-1)+10)
								//cout<<j<<" "<<gsl_vector_get(adjustedDerivative,j)<<" "<<gsl_vector_get(model,j-gsl_vector_get(*tstartgsl,*numberPulses-1))<<endl;
								gsl_vector_set(adjustedDerivative,j,gsl_vector_get(adjustedDerivative,j)-gsl_vector_get(model,j-gsl_vector_get(*tstartgsl,*numberPulses-1)));
							}
							
							if (gsl_vector_get(*flagTruncated,*numberPulses-1) == 1)	i = 0;
							else								i = gsl_vector_get(*tstartgsl,*numberPulses-1) + 1;
							foundPulse = false; 
						}
						else 
						{
							if (strcmp(reconstruct_init->EnergyMethod,"I2RBISALL") != 0)
							{
								*numberPulses = *numberPulses-1;
								foundPulse = false;
							}
							i++;
						}
					}
				}
			}
		} while (foundPulse == true);
	}
	else // Use the tstartPulsei provided as input parameters
	{
		gsl_vector_view temp;

		gsl_vector *tstartPulsei = gsl_vector_alloc(3);
		gsl_vector_set(tstartPulsei,0,reconstruct_init->tstartPulse1);
		gsl_vector_set(tstartPulsei,1,reconstruct_init->tstartPulse2);
		gsl_vector_set(tstartPulsei,2,reconstruct_init->tstartPulse3);

		if (reconstruct_init->tstartPulse2 == 0) 	*numberPulses = 1;
		else if (reconstruct_init->tstartPulse3 == 0) 	*numberPulses = 2;
		else						*numberPulses = 3;

		for (int i=0;i<*numberPulses;i++)
		{
			gsl_vector_set(*tstartgsl,i,gsl_vector_get(tstartPulsei,i));

			if (i != *numberPulses-1) 	
			{
				if ((gsl_vector_get(tstartPulsei,i) < 0) || (gsl_vector_get(tstartPulsei,i) > adjustedDerivative->size-2)
					|| (gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i) < 1) || (gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i) > adjustedDerivative->size-gsl_vector_get(tstartPulsei,i)))
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL);
				}
				temp = gsl_vector_subvector(adjustedDerivative,gsl_vector_get(tstartPulsei,i),gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i));
			}
			else	
			{
				if ((gsl_vector_get(tstartPulsei,i) < 0) || (gsl_vector_get(tstartPulsei,i) > adjustedDerivative->size-2)
					|| (sizeRecord-gsl_vector_get(tstartPulsei,i) < 1) || (sizeRecord-gsl_vector_get(tstartPulsei,i) > adjustedDerivative->size-gsl_vector_get(tstartPulsei,i)))
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL);
				}
				temp = gsl_vector_subvector(adjustedDerivative,gsl_vector_get(tstartPulsei,i),sizeRecord-gsl_vector_get(tstartPulsei,i));
			}

			if (i == 0)	gsl_vector_set(*maxDERgsl,i,gsl_vector_max(&temp.vector));
			else
			{
				if (find_model_maxDERs(gsl_vector_get(*maxDERgsl,i-1), reconstruct_init, &model))
				{
					message = "Cannot run find_model routine for pulse i=" + boost::lexical_cast<std::string>(i) + " when newPulses = 1";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}

				if ((gsl_vector_get(*tstartgsl,i-1) < 0) || (min(gsl_vector_get(*tstartgsl,i-1)+reconstruct_init->pulse_length,(double) sizeRecord)-1 > adjustedDerivative->size-1))
				{
					sprintf(valERROR,"%d",__LINE__+7);
					string str(valERROR);
					message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL);
				}
				for (int j=gsl_vector_get(*tstartgsl,i-1);j<min(gsl_vector_get(*tstartgsl,i-1)+reconstruct_init->pulse_length,(double) sizeRecord);j++)
				{
					gsl_vector_set(adjustedDerivative,j,gsl_vector_get(adjustedDerivative,j)-gsl_vector_get(model,j-gsl_vector_get(*tstartgsl,i-1)));
				}

				if (i != *numberPulses-1)	
				{
					if ((gsl_vector_get(tstartPulsei,i) < 0) || (gsl_vector_get(tstartPulsei,i) > adjustedDerivative->size-2)
						|| (gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i) < 1) || (gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i) > adjustedDerivative->size-gsl_vector_get(tstartPulsei,i)))
					{
						sprintf(valERROR,"%d",__LINE__+5);
						string str(valERROR);
						message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
						EP_PRINT_ERROR(message,EPFAIL);
					}
					temp = gsl_vector_subvector(adjustedDerivative,gsl_vector_get(tstartPulsei,i),gsl_vector_get(tstartPulsei,i+1)-gsl_vector_get(tstartPulsei,i));
				}
				else
				{
					if ((gsl_vector_get(tstartPulsei,i) < 0) || (gsl_vector_get(tstartPulsei,i) > adjustedDerivative->size-2)
						|| (sizeRecord-gsl_vector_get(tstartPulsei,i) < 1) || (sizeRecord-gsl_vector_get(tstartPulsei,i) > adjustedDerivative->size-gsl_vector_get(tstartPulsei,i)))
					{
						sprintf(valERROR,"%d",__LINE__+5);
						string str(valERROR);
						message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
						EP_PRINT_ERROR(message,EPFAIL);
					}
					temp = gsl_vector_subvector(adjustedDerivative,gsl_vector_get(tstartPulsei,i),sizeRecord-gsl_vector_get(tstartPulsei,i));
				}
				
				gsl_vector_set(*maxDERgsl,i,gsl_vector_max(&temp.vector));
			}
		}

		gsl_vector_free(tstartPulsei);
	}

	// Free allocated GSL vectors
	gsl_vector_free(model);
	gsl_vector_free(index_maxDERgsl);

	return(EPOK);
}
/*xxxx end of SECTION 14 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/