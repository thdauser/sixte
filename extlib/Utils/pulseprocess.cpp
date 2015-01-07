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

/************************************************************************
*                      PULSEPROCESS
*
*  File:      pulseprocess.cpp
*  Developer: Beatriz Cobo Martin
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene González Pérez
*             José Ramón Rodón Ortiz
*
*  Revision History:
*
*	version	1.0.0	08/10/07    First version
*	version 1.1.0	19/11/07	Changed the parameters "beforeIndex" and "afterIndex" from double to integer
* 								in the findPulse function
*   version 2.0.0   26/10/09    New findMean and findPulses
*   version 3.0.0	23/07/10	New functions added: "lpf_boxcar", "derMTH", "findTstart" & "findDerPoslength"
*   							Documentation updated
*   version 6.1.0	05/08/10	All initial truncated pulses detected. Modifications in "finDerPosLength" & "findTstart" functions
*								"pi" definition changed in "lpf_boxcar" function
*	version 6.1.1   29/09/10    findMean function documentation
*	version 7.0.0   15/10/10    lpf_boxcar: If status = 1 => tauFALL too small => Cut-off frequency too high => Equivalent to not filter
*	                                        If status = 2 => tauFALL too high => Cut-off frequency too low
*   04/04/11    lpf_box_car: It is going to work with a longer vector to not have strange results
*                            for the last boxLength windows
*   04/04/11    findTstart: tstart modified
*                            - Not the found one by findTstart (from the derivative and spread by the low-pass filter)
*                            - Shifted to the start of the not filtered pulse (+boxLength-10)
*                            - "-10" to be sure any pulse sample being in a future pulse-free interval
*   08/04/11	findTstart: Modifications to the previous "-10" in order to avoid GSL errors
*   04/05/11    findMean: //j_right = j_right-1; and //j_left = j_left-1;
*   24/05/11    findMean: // BE CAREFUL: This function provides different results in a 64-bits or a 32-bits PCs
*   26/05/11    findTstart renamed as findTstart_OLD
*               New functions: findThreshold
*                              kappaClipping
*                              medianKappaClipping
*                              findTstart
*                              findMeanSigma
*  30/05/11    lpf_boxcar: Calculation of 'value' different depending on the sign of 'value'
*              medianKappaClipping: findMeanSigma uses a subvector to not take into account the artifact of the LPF
*  02/06/11    findTstart modified in order to work properly just in case samplesUp=1
*  09/06/11    findTstart modified in order to work if a sample of the event is equal to the threshold
*              derMTH modified to not have 0 values when a pulse is saturated
*  21/07/11    findTstart has new output vectors
*  04/08/11    New functions:
*              derMTHSimple
*              findTstartPrimary
*  24/10/11    Added a new input/output parameter in findTstart
*  14/12/12    New input/output parameter maxDERgsl in findTstartPrimary
*  14/03/13    findPulses renamed as findPulsesOLD
*              New functions:
*              findPulses
*              getB
*              gsl_vector_Sumsubvector
*              findSePrPulses
*              getPulseHeight
*              find_model
*              RS_filter
*              interpolate_model
*  15/02/13   Deleted in/out parameter tend in findPulses
*  18/02/13   Changes in findSePrPulses:
*             	- "if (gsl_vector_get(*newPulses,i) == 1)" instead "if ((gsl_vector_get(*quality,i) != 1.0) && (gsl_vector_get(*quality,i) != 3.0) && (gsl_vector_get(*newPulses,i) == 1))"
*             	- "if ((gsl_vector_get(tstartSecondary,i) == gsl_vector_get(*tstart,j)) || (gsl_vector_get(tstartSecondary,i) == (gsl_vector_get(*tstart,j)+1)))" insted "if (gsl_vector_get(tstartSecondary,i) == (gsl_vector_get(*tstart,j)+1))"
*             Changed getB:
*             	- A pulse-free interval can not be found to apply the running sum filter (B=-999) is not an error
*             	- If any elements in Bgsl is still -999 (after using Bprev) => findSePrPulses is not used (but the task does not stop)
*             Changed findPulses:
*             	- Related to the changes in getB
*             Changed findTstartPrimary:
*               - 'subder' allocation protected
*               - 'ind = ind +1;' in a different line (inside the if)
*  21/02/13   Changes in findPulses to detect saturated pulses
*             Light change in findTstart
*  26/02/13   findTstart: New input/output parameter 'tstartgslOUTNEW'
*			  findPulses: New variable 'tstartNOsmt' (due to the change in findTstart)
*			  findSePrPulses: New input/output parameter 'tstartNOsmt'
*			                  getPulseheight uses 'tstartNOsmt'
*			                  New variable 'tstartPrimaryNOsmt'
*			                  'tstartNOsmt' must also be ordered
*			  getPulseHeight: 'safetymargintstart' not used (deleted as input/output parameter)
*  27/02/13   findSePrPulses: if (gsl_vector_get(tstartSecondary,i) == 0) and all the first elements = 0 => Not a secondary pulse (it belongs to a pulse set to 0 in the adjusted derivative)
*                             'tstartDERPrimary', 'tmaxDERPrimary' and 'tendDERPrimary' modified (added 'tstartDER')
*  28/02/13   findSePrPulses: gsl_vector_set_zero(end0) => In order to detect as saturated the saturated pulses truncated at the end
*                             Not only pulses whose newPulses=1 are taken into account to look for primary pulses (new secondary pulses have newPulses=2)
*  01/03/13   findSePrPulses: 'truePrimary' used
*                             Data in 'tstart', 'SorSeorPr', 'tstartDER'... ordered according 'tstartNOsmt' instead 'tstart'
*                             In order to detect as saturated the saturated pulses truncated at the end:
*                             	if ((i == numSaturated-1) && (gsl_vector_get(end0,i) == 0))	gsl_vector_set(end0,i,vectorin->size);
*  26/04/13   findSePrPulses: Improved the alignment between pulse and model (to build the adjusted derivative)
*  02/05/13   findSePrPulses: Improved the tstart by correcting the calculus with the pulse broadening
*             find_model: Find the pulse broadening by using the pulse models library
*             interpolate_model: Interpolate (in general) between two pulse broadenings
*  08/05/13   findSePrPulses: New input/output parameter 'tstartMOD' (which is the precisely calculated tstart)
*             findPulses: After founding all the pulses, 'tstart' is overwritten with 'tstartPRECISE' (which is the precisely calculated tstart)
*  09/05/13   findPulses: Precise 'tstart' is ordered into ascending order (and consequently, 'energy' and 'quality')
*  13/05/13   findPulses: 'tstart' is only ordered if the number of found pulses is different from 0
*  24/05/13   findPulses: 'tstart' is only ordered if the operation mode is normal (opmode = 1), not calibration (opmode = 0)
*  01/07/13   findSePrPulses: Changes in the alignment between pulse and model
*  02/10/13   findPulses, findSePrPulses, find_model and interpolate_model changed in order to not use PRETRIGS to
*             calculate precisely the tstart
*  08/10/13   Not used functions are deleted
*             Comments of different functions
*  05/08/14   Added some necessary 'gsl_vector_free'
*  10/09/14   findSePrPulses: In 'if (shift==0)', 'for (int k=0;k<modelScaled->size;k++)' instead 'for (int k=0;k<modelShifted->size;k++)'
*             getB: Allocation of some vectors changed in order to not have errors
*  03/10/14   findPulses: The energy of the pulses is also calculated if the operation mode is 0
*  10/10/14   findSePrPulses changed (nNewPulses,lookPrimary...)
*             findPulses changed (if only a primary pulse is found the do while is not run again)
*  16/10/14   findTstartPrimary modified in order to avoid errors if nsamplesUp is 1
*  23/10/14   findTstart: If there is a truncated pulse at the end, put its quality as 1
*             findSePrPulses: If there is a truncated pulse at the end, it is not going to look for Secondary or Primary pulse piled up with it
*    Dec/14   Migrated to CFITSIO (removal of ISDC DAL)
*             Errors processing changed
*             Deleted some non used functions
***********************************************************************/

/******************************************************************************
DESCRIPTION:

The purpose of this package is to support utilities related to find pulses.

MAP OF SECTIONS IN THIS FILE:

 - 1. lpf_boxcar
 - 2. derMTHSimple
 - 3. findMeanSigma
 - 4. medianKappaClipping
 - 5. getB
 - 6. getPulseHeight
 - 7. RS_filter
 - 8. find_model
 - 9. firstSampleModels
 - 10. find_model1stSample
 - 11. interpolate_model
 - 12. findTstart
 - 13. findPulses
 - 14. findSePulses

*******************************************************************************/


/***** SECTION 1 ************************************
*       INCLUDE's
****************************************************/
#include "pulseprocess.h"


/***** SECTION 2 ************************************************************
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
*   sinc(fc)=0, sinc(1)=0 => fc=1
*   fc=kf1 => fc~2f1
*
* - Declare variables
* - Define the LPF (frequency domain) and the box-car function (time domain)
* - It is going to work with a longer vector to not have strange results for the last boxLength windows
* - Apply the box-car window by shifting it along the input vector
* - Free GSL vectors
*
* Parameters:
* - invector: Input/Output vector
* - szVct: Size of the invector
* - tau_fall: Value related to the fall times of the pulses (seconds)
* - sampleRate: Sampling frequency (samples per second)
******************************************************************************/
int lpf_boxcar (gsl_vector **invector, int szVct, double tau_fall, int sampleRate)
{
	int status = EPOK;
	string message = "";

	//Declare variables
	gsl_vector *invectorAux;
	gsl_vector *invectorAux1;
	double cutFreq;	//Frequency domain
	int boxLength;	//Time domain
	double boxSum = 0.0;
	const double pi = 4.0 * atan(1.0);

	//Define the LPF (frequency domain) and the box-car function (time domain)
	cutFreq = 2 * (1/(2*pi*tau_fall));	//According to Jan, sinc(f1)=0.6 where f1=1/(2pi*tau_fall)
						//sinc(0.5)~0.6 => f1~0.5
						//sinc(fc)=0, sinc(1)=0 => fc=1
						//fc=kf1 => fc~2f1
	boxLength =(int) ((1/cutFreq) * sampleRate);

	if (boxLength < 1)boxLength = 1;
	if (boxLength >= szVct)
	{
	    message = "Too  low box car cut-off frequency";
	    EP_PRINT_ERROR(message,4);
	}

	// It is going to work with a longer vector to not have strange results for the last boxLength windows
	invectorAux = gsl_vector_alloc(szVct+boxLength);
	invectorAux1 = gsl_vector_alloc(szVct+boxLength);
	for (int i=0;i<szVct;i++)
	{
		gsl_vector_set(invectorAux,i,gsl_vector_get(*invector,i));
	}
	double value = gsl_vector_get(*invector,szVct-1);
	for (int i=szVct;i<szVct+boxLength;i++)
	{
		if (value > 0) 		value = value-0.01*value;
		else if (value< 0)	value = value+0.01*value;
		gsl_vector_set(invectorAux,i,value);
	}

	//Apply the box-car window by shifting it along the input vector (the longer version of it)
	for (int i=0;i<boxLength;i++)
	{
		boxSum = boxSum + gsl_vector_get(invectorAux,i);
	}
	gsl_vector_set(invectorAux1,0,boxSum/boxLength);

	for (int i=0;i<invectorAux->size-boxLength;i++)
	{
		boxSum = boxSum - gsl_vector_get(invectorAux,i) + gsl_vector_get(invectorAux,i+boxLength);
		gsl_vector_set(invectorAux1,i+1,boxSum/boxLength);
	}

	for (int i=0;i<boxLength-1;i++)
	{
		boxSum = boxSum - gsl_vector_get(invectorAux,invectorAux->size-boxLength+i);
		gsl_vector_set(invectorAux1,invectorAux->size-boxLength+i+1,boxSum/(boxLength-i-1));
	}

	for (int i=0;i<szVct;i++)
	{
		gsl_vector_set(*invector,i,gsl_vector_get(invectorAux1,i));
	}

	// Free GSL vectors
	gsl_vector_free(invectorAux);
	gsl_vector_free(invectorAux1);
	
	if (boxLength == 1)
	{
	    message = "Too high cut-off frequency (no filtering)";
	    EP_PRINT_ERROR(message,3);
	}

	return (EPOK);
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 ************************************************************
* derMTHSimple function: This function applies the derivative method (x_i-x_(i-1)) to the input vector
*
* The derivative method provides more sensitivity to handle with pilot pulses.
* Moreover, little variations of the baseline will not affect.
*
* Arguments:
* 	- invector: Input/Ouput vector (derivative)
* 	- sign: Sign of the derivative
* 	- szVct: Size of invector
******************************************************************************/
int derMTHSimple (gsl_vector **invector,gsl_vector **sign, int szVct)
{
	int status = EPOK;

	for (int i=0; i<szVct-1; i++)
	{
		if ((gsl_vector_get(*invector,i+1)-gsl_vector_get(*invector,i)) > 0)		gsl_vector_set(*sign,i,1);
		else if ((gsl_vector_get(*invector,i+1)-gsl_vector_get(*invector,i)) < 0)	gsl_vector_set(*sign,i,-1);
		else if ((gsl_vector_get(*invector,i+1)-gsl_vector_get(*invector,i)) == 0)	gsl_vector_set(*sign,i,0);
		gsl_vector_set(*invector,i,gsl_vector_get(*invector,i+1)-gsl_vector_get(*invector,i));

	}
	gsl_vector_set(*invector,szVct-1,gsl_vector_get(*invector,szVct-2));
	//gsl_vector_set(*invector,szVct-1,-999);

	return (EPOK);
}
/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
* findMeanSigma function: This function calculates the mean and the standard deviation
******************************************************************************/
int findMeanSigma (gsl_vector *invector, double *mean, double *sigma, FILE *temporalFile)
{
	int status = EPOK;

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
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* medianKappaClipping function: This function calculates a threshold in the first derivative vector by using
*                               a Kappa-clipping method (replacing points beyond mean+-Kappa*sigma with the median).
*
* First, mean and sigma are calculated and invector values out of (mean+Kappa*sigma,mean-Kappa*sigma) are replaced
* with the median (we are trying to look for the baseline). And this process is iteratively repeated until there are
* no points beyond mean+-Kappa *sigma. Finally, the threshold is calculated as mean+nSigmas*sigma ('+' is used because
* if there are pulses in the input invector they are always positive after a previous analysis and the 'necessary' changes).
*
* - Declare variables
* - Calculate the median
* - Iterate until there are no points out of the maximum excursion (kappa*sigma)
* - Establish the threshold as mean+nSigmas*sigma
*
* Parameters:
* - invector: First derivative of the filtered event
* - Kappa: To establish the surroundings of the mean
* - stopCriteria: It is given in %
* - nSigmas: Times sigma to calculate threshold (mean+nSigmas*sigma)
* - boxLPF: Length of the low-pass filtering
* - threshold: Calculated output value
******************************************************************************/
int medianKappaClipping (gsl_vector *invector, double kappa, double stopCriteria, double nSigmas, int boxLPF, double *threshold, FILE * temporalFile)
{
	int status = EPOK;
	string message = "";

	// Declare variables
	int size = invector->size; // Size of the input vector
	double mean1, sg1;
	double mean2, sg2;
	gsl_vector_view temp;
	char val[256];
	// Variables to remove input vector elements higher than the maximum excursion (kappa*sg)
	int i;												// To go through the elements of a vector
	int cnt;											// Number of points inside the excursion (mean+-excursion)
	gsl_vector *invectorNew = gsl_vector_alloc(size);	// Auxiliary vector
	// To calculate the median
	double data[size];									// Intermediate value to use 'gsl_stats_median_from_sorted_data'
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
		temp = gsl_vector_subvector(invectorNew,0,size-boxLPF-1);
		if (findMeanSigma (&temp.vector, &mean1, &sg1, temporalFile))
		{
			message = "Cannot run findMeanSigma routine for kappa-sigma iteration";
			EP_PRINT_ERROR(message,EPFAIL);
		}
		i = 0;
		cnt = 0;

		while (i<invectorNew->size)
		{
			if ((gsl_vector_get(invectorNew,i) >= mean1 + kappa*sg1) || (gsl_vector_get(invectorNew,i) <= mean1 - kappa*sg1))
			// HARDPOINT!!!!!!!!!!!!!!!!!!! (kappa)
			{
				gsl_vector_set(invectorNew,i,median);
				cnt++;
			}
			i++;
		}

		if (cnt != 0)
		// Some points of the invector have been replaced with the median
		{
			temp = gsl_vector_subvector(invectorNew,0,size-boxLPF-1);
			if (findMeanSigma (&temp.vector, &mean2, &sg2, temporalFile))
			{
				message = "Cannot run findMeanSigma routine for kappa-sigma iteration after replacement with the median";
				EP_PRINT_ERROR(message,EPFAIL);
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
	//if (mean2>1e10)	*threshold = 1e10;
	//else	*threshold = mean2+nSigmas*sg2;	// HARDPOINT!!!!!!!!!!!!!!!!!!! (nSigmas)
	*threshold = mean2+nSigmas*sg2;	// HARDPOINT!!!!!!!!!!!!!!!!!!! (nSigmas)

	gsl_vector_free(invectorNew);

	return EPOK;
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
* getB function: This function calculates the sum, B, of the lb digitized data samples of a pulse-free interval immediately
*       before the current pulse.
*       If the pulse-free interval before the current pulse is lower than lb, B is calculated with the available
*       number of samples (although the available number of samples was lower than lb).
*       If there is not a pulse-free interval before the pulse, it is looked for it after the current pulse.
*       The number of samples of the pulse-free interval used to calculate B has to be stored in the lb vector.
*
* First of all, the auxiliary variable 'Baux' is initialized to -999.
* The algorithm is divided into 2 if's:
*  - When the current pulse is the first pulse into the record:
*  		- tstart>=lb => Sum lb samples
*  		- 0<tstart<lb => Sum the available number of samples (although the available number of samples was lower than lb)
*  		- tstart=0 => If there is not a pulse-free interval before the pulse, it is looked for it after the current pulse
*  - When the current pulse is not the first pulse into the record:
*  		- tstart_(i)-tend_(i-1)>=lb => Sum lb samples
*  		- 0<tstart_(i)-tend_(i-1)<lb => Sum the available number of samples (although the available number of samples was lower than lb)
*  		- If there is not a pulse-free interval before the pulse, it is looked for it after the current pulse
*
* If 'Baux' is still -999 => A pulse-free interval can not be found to apply the running sum filter => This has to be taken into account,
* out of the function, to try to get a usable B.
*
* Parameters:
* - vectorin: Input row=record vector
* - tstart: Starting time of the pulses into the row=event
* - nPulses: Number of pulses into the row=event
* - lb: Vector containing the baseline averaging length used for each pulse
*       All its elements are equal to the 'Lb' global variable at the beginning
* - sizePulse:  Size of the pulse in bins, ntaus * tauFALL in bins (equal to 'sizePulse_b' global variable)
* - B: In general, sum of the Lb digitized data samples of a pulse-free interval immediately before
*      the current pulse
****************************************/
int getB(gsl_vector *vectorin, gsl_vector *tstart, int nPulses, gsl_vector **lb, int sizepulse, gsl_vector **B, FILE * temporalFile)
{
	char val[256];

	int status = EPOK;
	string message = "";

	// Declare variables
	*B = gsl_vector_alloc(tstart->size);
	gsl_vector_set_all(*B,-999);
	double Baux = -999;
	double tendprev;

	// Auxiliary variables
	gsl_vector *input;
	gsl_vector_view temp;					// In order to handle with gsl_vector_view (subvectors)

	for (int i=0;i<nPulses;i++)
	{
		if (i == 0)		// First pulse into a row=event
		//  Current pulse
		//      //\\          /\     /\       /\
		// (1) //  \\  (2)   /  \(3)/  \ (4) /  \    (5)
		// ----      --------    ---    -----    ----------
		{
			if (gsl_vector_get(tstart,0)>=gsl_vector_get(*lb,0))
			// tstart>=lb => Sum lb samples
			// length_(1)>=lb
			{
				input = gsl_vector_alloc(gsl_vector_get(*lb,0));
				temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,i)-gsl_vector_get(*lb,0),gsl_vector_get(*lb,0));
				gsl_vector_memcpy(input, &temp.vector);

				// Sum all the elements of 'input'
				if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
				{
					message = "Cannot run gsl_vector_Sumsubvector routine when tstart>=lb";
					EP_PRINT_ERROR(message,EPFAIL);
				}
				gsl_vector_free(input);
			}
			else if ((gsl_vector_get(tstart,0)<gsl_vector_get(*lb,0)) && (gsl_vector_get(tstart,0)>1))
		    // 0<tstart<lb => Sum the available number of samples (although the available number of samples was lower than lb)
			// 0<length_(1)<lb
			{
				input = gsl_vector_alloc(gsl_vector_get(tstart,0));
				temp = gsl_vector_subvector(vectorin,0,gsl_vector_get(tstart,0));
				gsl_vector_memcpy(input, &temp.vector);

				// Sum all the elements of 'input'
				if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
				{
					message = "Cannot run gsl_vector_Sumsubvector routine when tstart<lb";
					EP_PRINT_ERROR(message,EPFAIL);
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
							input = gsl_vector_alloc(gsl_vector_get(*lb,0));
							temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,j)-gsl_vector_get(*lb,0),gsl_vector_get(*lb,0));
							gsl_vector_memcpy(input, &temp.vector);
						}
						else if (gsl_vector_get(tstart,j)-tendprev < gsl_vector_get(*lb,0))
						// 0<length_(j)<lb (j/=nPulses)
						{
							input = gsl_vector_alloc(gsl_vector_get(tstart,j)-tendprev-1);
							temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,j)-tendprev-1);
							gsl_vector_memcpy(input, &temp.vector);
							gsl_vector_set(*lb,0,gsl_vector_get(tstart,j)-tendprev-1);
						}

						// Sum all the elements of input
						if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
						{
							message = "Cannot run gsl_vector_Sumsubvector routine when no pulse free interval before the pulse";
							EP_PRINT_ERROR(message,EPFAIL);
						}
						gsl_vector_free(input);

						break;
					}

					if (j == nPulses-1)	// Last pulse into the row=event
					// (5) is analyzed
					{
						tendprev = gsl_vector_get(tstart,j)+sizepulse-1;
						if (vectorin->size-tendprev >= gsl_vector_get(*lb,0))
						{
							input = gsl_vector_alloc(gsl_vector_get(*lb,0));
							temp = gsl_vector_subvector(vectorin,tendprev,gsl_vector_get(*lb,0));
							gsl_vector_memcpy(input, &temp.vector);

							// Sum all the elements of input
							if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
							{
								message = "Cannot run gsl_vector_Sumsubvector routine when last pulse in row & vectorin>tendprev>lb";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							gsl_vector_free(input);
						}
						else if ((vectorin->size-tendprev < gsl_vector_get(*lb,0)) && (vectorin->size-tendprev > 1))
						{
							input = gsl_vector_alloc(vectorin->size-tendprev-1);
							temp = gsl_vector_subvector(vectorin,tendprev+1,vectorin->size-tendprev-1);
							gsl_vector_memcpy(input, &temp.vector);

							// Sum all the elements of input
							if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
							{
								message = "Cannot run gsl_vector_Sumsubvector routine when last pulse in row";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							gsl_vector_set(*lb,0,vectorin->size-tendprev-1);
							gsl_vector_free(input);
						}
					}
				}
			}
		}
		else	// Not first pulse into a row=event
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
				input = gsl_vector_alloc(gsl_vector_get(*lb,i));
				temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,i)-gsl_vector_get(*lb,i),gsl_vector_get(*lb,i));
				gsl_vector_memcpy(input, &temp.vector);

				// Sum all the elements of input
				if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
				{
					message = "Cannot run gsl_vector_Sumsubvector routine length_(2)>=lb";
					EP_PRINT_ERROR(message,EPFAIL);
				}
				gsl_vector_free(input);
			}
			else if ((gsl_vector_get(tstart,i)-tendprev<gsl_vector_get(*lb,i)) && (gsl_vector_get(tstart,i)-tendprev>1))
			// 0<tstart_(i)-tend_(i-1)<lb => Sum the available number of samples (although the available number of samples was lower than lb)
			// 0<length_(2)<lb
			{
				input = gsl_vector_alloc(gsl_vector_get(tstart,i)-tendprev-1);
				temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,i)-tendprev-1);
				gsl_vector_memcpy(input, &temp.vector);

				// Sum all the elements of input
				if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
				{
					message = "Cannot run gsl_vector_Sumsubvector routine when 0<length_(2)<lb";
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
							input = gsl_vector_alloc(gsl_vector_get(*lb,i));
							temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i),gsl_vector_get(*lb,i));
							gsl_vector_memcpy(input, &temp.vector);

							// Sum all the elements of input
							if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
							{
								message = "Cannot run gsl_vector_Sumsubvector routine when no first pulse in row & tstart-tendprev >= lb";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							gsl_vector_free(input);

							break;
						}
						else if ((gsl_vector_get(tstart,j+1)-tendprev < gsl_vector_get(*lb,i)) && ((gsl_vector_get(tstart,j+1)-tendprev>1)))
						{
							input = gsl_vector_alloc(gsl_vector_get(tstart,j+1)-tendprev-1);
							temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,j+1)-tendprev-1);
							gsl_vector_memcpy(input, &temp.vector);

							// Sum all the elements of input
							if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
							{
								message = "Cannot run gsl_vector_Sumsubvector routine when no first pulse in row & tstart-tendprev < lb";
								EP_PRINT_ERROR(message,EPFAIL);
							}
							gsl_vector_set(*lb,i,gsl_vector_get(tstart,j+1)-tendprev-1);
							gsl_vector_free(input);

							break;
						}
					}

					if (j == nPulses-1)	// Last pulse into the row=event
					// (5) is analyzed
					{
						if (vectorin->size-tendprev >= gsl_vector_get(*lb,i))
						{
							input = gsl_vector_alloc(gsl_vector_get(*lb,i));
							temp = gsl_vector_subvector(vectorin,tendprev,gsl_vector_get(*lb,i));
							gsl_vector_memcpy(input, &temp.vector);
						}
						//else if ((vectorin->size-tendprev < gsl_vector_get(*lb,i)) && (vectorin->size-tendprev != 1))
						else if ((vectorin->size-tendprev < gsl_vector_get(*lb,i)) && (vectorin->size-tendprev > 1))
						{
							input = gsl_vector_alloc(vectorin->size-tendprev-1);
							temp = gsl_vector_subvector(vectorin,tendprev+1,vectorin->size-tendprev-1);
							gsl_vector_memcpy(input, &temp.vector);
							gsl_vector_set(*lb,i,vectorin->size-tendprev-1);
						}
						//else if ((vectorin->size-tendprev < gsl_vector_get(*lb,i)) && (vectorin->size-tendprev == 1))
						else if ((vectorin->size-tendprev < gsl_vector_get(*lb,i)) && (vectorin->size-tendprev <= 1))
						{
							gsl_vector_set(*lb,i,gsl_vector_get(*lb,i-1));

							break;
						}

						// Sum all the elements of input
						if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
						{
							message = "Cannot run gsl_vector_Sumsubvector routine when no pulse-free interval before pulse & last pulse in row";
							EP_PRINT_ERROR(message,EPFAIL);
						}
						gsl_vector_free(input);
					}
				}
			}
		}

		gsl_vector_set(*B,i,Baux);
	}

	return(EPOK);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* getPulseHeight function: This function estimates the pulse height of a pulse by using a running sum filter.
*                 It extracts from the record, vectorin, the pulse whose pulse height is going to be estimated
*                 by using RS_filter.
*
* - Declare variables
* - Extracting from the record the pulse whose pulse height is going to be estimated
* - Apply the running sum filter
*
* Parameters:
* - vectorin: Not filtered row=event
* - tstart: Starting time of the pulse whose pulse height is going to be estimated
* - tstartnext: Starting time of the next pulse whose pulse height is going to be estimated
* - lastPulse: 1 if the pulse is the last one into the row=event or the only one
* - lrs: Running sum length (equal to the 'Lrs' global_variable)
* - lb: Baseline averaging length used for the pulse whose pulse height is going to be estimated
* - B: In general, sum of the Lb digitized data samples of a pulse-free interval immediately before
*      the current pulse
* - sizepulse: Size of the pulse in bins, ntaus * tauFALL in bins (equal to 'sizePulse_b' global variable)
* - pulseheight: Estimated pulse height of the pulse
****************************************/
int getPulseHeight(gsl_vector *vectorin, double tstart, double tstartnext, int lastPulse, double lrs, double lb, double B, int sizepulse, double *pulseheight, FILE * temporalFile)
{
	char val[256];
	char val_aux[256];

	// Declare variables
	int status=EPOK;
	long tend;	// Ending time of the pulse
	double ph;	// Pulse height
	string message = "";
	
	// Auxiliary variables
	gsl_vector *input;		// Segment of 'vectorin' where is the pulse
	gsl_vector_view temp;	// In order to handle with gsl_vector_view (subvectors)

	// Extracting from the record the pulse whose pulse height is going to be estimated
	if (lastPulse == 1)
	{
		tend = tstart+sizepulse-1;
		if (tend >= vectorin->size)	// Truncated pulses at the end of the row
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
		input = gsl_vector_alloc(vectorin->size-tstart-(vectorin->size-tend));	// Only the pulse
		temp = gsl_vector_subvector(vectorin,tstart,input->size);
		gsl_vector_memcpy(input, &temp.vector);

		// Apply the running sum filter
		if (RS_filter (input, lrs, lb, B, &ph, temporalFile))
		{
		    message = "Cannot run RS_filter routine when size-tstart>size-tend";
		    EP_PRINT_ERROR(message,EPFAIL);
		}
		*pulseheight = ph;

		gsl_vector_free(input);
	}

	return(EPOK);
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
* RS_filter function: This function uses the running sum filter to find the pulse height.
*            It always works in time domain.
*
* A running sum filter, RS, is the sum of lrs digitized data samples. It is continuously updated upon the arrival of
* new data point. Simultaneously a baseline filter, B, is the sum of lb digitized data samples without pulses.
*
* The algorithm looks for the time when RS/lrs reaches its maximum. At that time RS is stored, Rs_max, and the baseline
* is scaled with lrs, Bp (Bp=B*lrs/lb). Then, the pulse height related to the pulse pseudoenergy is given by:
*
*                     Rs_max-Bp
*    Pulse height = -------------
*                        lrs
*
* Parameters:
* - vector: Not filtered pulse (extracted from the row=event in 'getPulseHeight')
* - lrs: Running sum length (samples)
* - lb: Baseline averaging length (samples)
* - B: In general, sum of the lb digitized data samples of a pulse-free interval immediately before
*      the current pulse
* - pulseheight: Pulse height of the pulse
*****************************************/
int RS_filter (gsl_vector *vector, double lrs, double lb, double B, double *pulseheight, FILE *temporalFile)
{
	char val[256];

	int status=EPOK;
	string message = "";

	// Declare variables

	double Rs;				// Sum of lrs digitized data samples
	double Rs_max = -1e20;	// Rs is continuously updated upon the arrival of a new data point if the new value is
	                        // higher than the old one
	double Bp;

	if (vector->size<lrs)
	{
		if (gsl_vector_Sumsubvector(vector,0,vector->size,&Rs))
		{
		    message = "Cannot run gsl_vector_Sumsubvector routine when size<lrs";
		    EP_PRINT_ERROR(message,EPFAIL);
		}
		if (Rs > Rs_max)			Rs_max = Rs;
	}
	else
	{
		for (int i=0;i<vector->size;i++)
		{
			if (i+lrs > vector->size)	break;
			if (gsl_vector_Sumsubvector(vector,i,lrs,&Rs))
			{
			    message = "Cannot run gsl_vector_Sumsubvector routine when size>lrs in iteration i=" + i;
			    EP_PRINT_ERROR(message,EPFAIL);
			}
			if (Rs > Rs_max)			{Rs_max = Rs;}
		}
	}

	Bp = B*lrs/lb;
	*pulseheight = (Rs_max-Bp)/lrs;

	return EPOK;
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 9 ************************************************************
* find_model function: This function uses the pulse height (or energy) in order to choose the proper pulse template of
*                     the pulse templates library.
*
* It finds the two energies closer in the pulse models library and interpolates ('interpolate_model')
* 	- If pulse height (energy) is lower than the lower pulse height (energy) in the pulse models library => The model with
*     a lower pulse height (energy) in the pulse models library is chosen
*   - If pulse height (energy) is higher than the higher pulse height (energy) in the pulse models library => The model with
*     a higher pulse height (energy) in the pulse models library is chosen
*
* Parameters:
* - ph: Pulse height (or energy) of the pulse whose pulse template is looking for
* - modelsvalues: Matrix where the values of Energy of each pulse template are stored
* - models: Matrix where all the pulse templates of the pulse templates library are going to be stored
* - modelFound: Found template of the pulse whose pulse height is 'ph'
****************************************/
int find_model(double ph, gsl_matrix *modelsvalues, gsl_matrix *models, gsl_vector **modelFound, FILE * temporalFile)
{
	char val[256];
	char val_aux[256];

	int status=EPOK;
	string message = "";

	long nummodels = modelsvalues->size1;

	if (ph < gsl_matrix_get(modelsvalues,0,0))
	{
		gsl_matrix_get_row(*modelFound,models,0);
	}
	else if (ph > gsl_matrix_get(modelsvalues,nummodels-1,0))
	{
		gsl_matrix_get_row(*modelFound,models,nummodels-1);
	}
	else
	{
		for (int i=0;i<nummodels;i++)
		{
			if (ph == gsl_matrix_get(modelsvalues,i,0))
			{
				gsl_matrix_get_row(*modelFound,models,i);

				break;
			}
			else if ((ph > gsl_matrix_get(modelsvalues,i,0)) && (ph < gsl_matrix_get(modelsvalues,i+1,0)))
			{
				// Interpolate between the two corresponding rows in "models"
				gsl_vector *modelAux = gsl_vector_alloc(models->size2);
				gsl_vector_set_zero(modelAux);
				gsl_vector *model1 = gsl_vector_alloc(models->size2);
				gsl_vector *model2 = gsl_vector_alloc(models->size2);
				gsl_matrix_get_row(model1,models,i);
				gsl_matrix_get_row(model2,models,i+1);
				if (interpolate_model(&modelAux,ph,model1,gsl_matrix_get(modelsvalues,i,0),model2,gsl_matrix_get(modelsvalues,i+1,0), temporalFile))
				{
				    message = "Cannot run interpolate_model with two rows in models";
				    EP_PRINT_ERROR(message,EPFAIL);
				}
				gsl_vector_memcpy(*modelFound,modelAux);
				gsl_vector_free(modelAux);
				gsl_vector_free(model1);
				gsl_vector_free(model2);

				break;
			}
		}
	}

    return(EPOK);
}
/*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 10 ************************************************************
* firstSampleModels function: This function finds the first sample over the threshold of each template of the library
*                             (and the index of that first sample)
*
* Parameters:
* - templates: Matrix which contains all the templates of the library
* - threshold: Threshold
* - firstSamples: Vector which contains the first samples over the threshold of all the templates of the library
* - index_firstSamples: Vector which contains the index of the first samples
****************************************/
int firstSampleModels (gsl_matrix *templates, double threshold, gsl_vector **firstSamples, gsl_vector **index_firstSamples, FILE *temporalFile)
{
  	int status=EPOK;

	int numtemplates = templates->size1;
	gsl_vector *atemplate = gsl_vector_alloc(templates->size2);

	for (int i=0;i<numtemplates;i++)
	{
		gsl_matrix_get_row(atemplate,templates,i);
		for (int j=0;j<atemplate->size;j++)
		{
			if (gsl_vector_get(atemplate,j) > threshold)
			{
				gsl_vector_set(*index_firstSamples,i,j);
				gsl_vector_set(*firstSamples,i,gsl_vector_get(atemplate,j));

				break;
			}
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 11 ************************************************************
* find_model1stSample function: This function uses the value of the first sample of the pulse over the threshold in order
*                               to choose the proper pulse template of the pulse templates library
*
* Parameters:
* - firstSample: Value of the first sample over the threshold of the pulse whose template is looking for
* - firstSamples: Vector with the values of the first sample over the threshold of every template of the library
* - models: Matrix where all the pulse templates of the pulse templates library are stored
* - modelFound: Found template of the pulse
****************************************/
int find_model1stSample(double firstSample, gsl_vector *firstSamples, gsl_matrix *models, gsl_vector **modelFound, FILE * temporalFile)
{
	char val[256];
	char val_aux[256];

  	int status=EPOK;
  	string message = "";

	long nummodels = models->size1;

	if (firstSample < gsl_vector_get(firstSamples,0))
	{
		gsl_matrix_get_row(*modelFound,models,0);
	}
	else if (firstSample > gsl_vector_get(firstSamples,nummodels-1))
	{
		gsl_matrix_get_row(*modelFound,models,nummodels-1);
	}
	else
	{
		for (int i=0;i<nummodels;i++)
		{
			if (firstSample == gsl_vector_get(firstSamples,i))
			{
				gsl_matrix_get_row(*modelFound,models,i);

				break;
			}
			else if ((firstSample > gsl_vector_get(firstSamples,i)) && (firstSample < gsl_vector_get(firstSamples,i+1)))
			{
				// Interpolate between the two corresponding rows in "models"
				gsl_vector *modelAux = gsl_vector_alloc(models->size2);
				gsl_vector_set_zero(modelAux);
				gsl_vector *model1 = gsl_vector_alloc(models->size2);
				gsl_vector *model2 = gsl_vector_alloc(models->size2);
				gsl_matrix_get_row(model1,models,i);
				gsl_matrix_get_row(model2,models,i+1);
				if (interpolate_model(&modelAux,firstSample,model1,gsl_vector_get(firstSamples,i),model2,gsl_vector_get(firstSamples,i+1), temporalFile))
				{
				    message = "Cannot run interpolate_model with two rows in models";
				    EP_PRINT_ERROR(message,EPFAIL);
				}
				gsl_vector_memcpy(*modelFound,modelAux);
				gsl_vector_free(modelAux);
				gsl_vector_free(model1);
				gsl_vector_free(model2);

				break;
			}
		}
	}

    return(EPOK);
}
/*xxxx end of SECTION 11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 12 ************************************************************
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
* Arguments:
*      - modelFound: Found model of the pulse whose pulse height is ph_model
*      - ph_model: Pulse height of the pulse whose model is looking for
*      - modelIn1: Model of the pulse whose pulse height is immediately lower than ph_model in the library FITS file
*      - ph_modelIn1: Pulse height immediately lower than ph_model in the library FITS file
*      - modelIn2: Model of the pulse whose pulse height is immediately greater than ph_model in the library FITS file
*      - ph_modelIn2: Pulse height immediately greater than ph_model in the library FITS file*
****************************************/
int interpolate_model(gsl_vector **modelFound, double ph_model, gsl_vector *modelIn1, double ph_modelIn1, gsl_vector *modelIn2, double ph_modelIn2, FILE * temporalFile)
{
	char val[256];

	int status=EPOK;

	// Declare variables
	double factor1, factor2;
	gsl_vector_set_zero(*modelFound);

	// Method 1: The simplest method
	/*gsl_vector_add(*modelFound,modelIn1);
	gsl_vector_add(*modelFound,modelIn2);
	gsl_vector_scale(*modelFound,0.5);*/

	// Method 2: A bit more intelligent averaging
	factor1 = (ph_modelIn2-ph_model)/(ph_modelIn2-ph_modelIn1);
	factor2 = (ph_model-ph_modelIn1)/(ph_modelIn2-ph_modelIn1);
	gsl_vector_scale(modelIn1,factor1);
	gsl_vector_scale(modelIn2,factor2);
	gsl_vector_add(*modelFound,modelIn1);
	gsl_vector_add(*modelFound,modelIn2);

    return(EPOK);
}
/*xxxx end of SECTION 12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 13 ************************************************************
* findTstart function: This function finds the pulses tstart in the input vector (first derivative of the filtered
* 					   record)
*
* This function scans all values until it finds nSamplesUp consecutive values (due to noise more than 1 value is
* required) in the derivative over the threshold. To look for more pulses, it finds nSamplesUp consecutive values
* (due to noise) in the derivative below the threshold and then, it starts to scan again until it again try to find
* nSamplesUp consecutive values over the threshold.
*
* In general, the tstarts found in the derivative do not coincide with the tstart of the not filtered record
* (by the spreading caused by the low-pass filtered).
*
* - Declare variables
* - Allocate GSL vectors
* - Obtain tstart of each pulse in the derivative:
* 	-If der_i>threshold and prevPulse=false, it looks for nSamplesUp consecutive samples over the threshold
*   	- If not, it looks again for a pulse crossing over the threshold
*       - If yes => Pulse found (truncated if at the beginning)
*   - If der_i>threshold and prevPulse=true, it looks for a sample below the threshold
*     	- If not, it looks again for a sample below the threshold
*       - If yes, it looks for nSamplesUp consecutive samples below the threshold and again it starts to look for a pulse
* - Obtain tstart of each pulse in the not filtered signal (taking into account a safety margin)
*
* Parameters:
* - der: First derivative of the filtered event
* - adaptativeThreshold
* - nSamplesUp: Number of consecutive bins over the threshold to 'find' a pulse
* - safetyMargin: (According to Jan) A safety margin has to be taken into account before the tstart (bins)
* - allPulsesMode: 0-> If it finds a pulse the function returns
*                  1-> It finds all pulses of the event
* - sampling: Samprate (read from the input FITS file)
*             Necessary to recalculate the tstart by establishing a safety margin
* - numberPulses: Number of pulses found in the input (invector)
* - thereIsPulse: 0 -> the algorithm has not found any pulse
*                 1 -> the algorithm has found unless one pulse
* - tstartgslOUT: Pulses tstart (bins)
* - tstartgslOUTNEW: Pulses tstart without taking into account the safetyMargin(bins)
* - flagTruncated: Flag indicating if the pulse is truncated
*                 (inside this function only initial truncated pulses are classified)
* - tstartDERgsl: Point where the first derivative crosses over the threshold
* - tmaxDERgsl: Point where the first derivative is maximum
* - maxDERgsl: Maximum of the first derivative
* - tendDERgsl: Point where the first derivative crosses below the threshold
******************************************************************************/
int findTstart (gsl_vector *der, double adaptativethreshold, int nSamplesUp, double safetyMargin, int allPulsesMode, double sampling, int *numberPulses, int *thereIsPulse, gsl_vector **tstartgslOUT, gsl_vector **tstartgslOUTNEW, gsl_vector **flagTruncated, gsl_vector **tstartDERgsl, gsl_vector **tmaxDERgsl, gsl_vector **maxDERgsl, gsl_vector **tendDERgsl, FILE * temporalFile)
{
	// Provisional => Delete in future
	char val[256];
	char val_aux[256];

	int status=EPOK;
	int verbosity = 3;			// Verbosity level of the output log file

	// Declare variables
	int szRw = der->size;		 // Size of segment of process
	*numberPulses = 0;
	*thereIsPulse = 0;
	int i = 0;					 // To go through the elements of a vector
	bool prevPulse = false;		 // false: It looks for nSamplesUp consecutive samples over the threshold
	                             // true: It looks for nSamplesUp consecutive samples below the threshold
	int cntDown = 0;			 // To taking into account how many consecutive samples are down the threshold
	int cntUp = 0;				 // To taking into account how many consecutive samples are over the threshold
	int possibleTstart;			 // To store the first of the nSamplesUp consecutive samples over the threshold
	int possibleTend;			 // To store the first of the nSamplesUp consecutive samples below the threshold
	int maxDER_index;			 // Maximum of the derivative between tstartDER and tendDER
	gsl_vector_view temp;		 // In order to handle with gsl_vector_view (subvectors)
	bool lastPulse = true;

	// Allocate GSL vectors
	gsl_vector *tstartDER = gsl_vector_alloc(szRw);		// tstarts referred to the derivative
	gsl_vector *tendDER = gsl_vector_alloc(szRw);		// tends referred to the derivative
	gsl_vector_set_zero(tendDER);
	*tstartgslOUT = gsl_vector_alloc(szRw);				// tstarts referred to the not filtered event
	*tstartgslOUTNEW = gsl_vector_alloc(szRw);
	*flagTruncated = gsl_vector_alloc(szRw);
	gsl_vector_set_zero(*flagTruncated);
	*tstartDERgsl = gsl_vector_alloc(szRw);				// Point where the first derivative crosses over the threshold
	*tmaxDERgsl = gsl_vector_alloc(szRw);				// Point where the first derivative is maximum
	*maxDERgsl = gsl_vector_alloc(szRw);				// Maximum of the first derivative
	*tendDERgsl = gsl_vector_alloc(szRw);				// Point where the first derivative crosses below the threshold

	// Obtain tstart of each pulse in the derivative
	while (i < szRw-1)
	{
		if ((gsl_vector_get(der,i) > adaptativethreshold) && (prevPulse == false))
		{
			if (cntUp == 0)
			{
				possibleTstart=i;

				if ((nSamplesUp == 1) || ((nSamplesUp != 1) && (i == szRw-2)))
				{
					if (possibleTstart == 0)
					{
						gsl_vector_set(tstartDER,*numberPulses,possibleTstart);
						gsl_vector_set(*flagTruncated,*numberPulses,1);
						*numberPulses = *numberPulses +1;
						if (allPulsesMode == 0) break;
						prevPulse = true;
					}
					else
					{
						gsl_vector_set(tstartDER,*numberPulses,possibleTstart);
						*numberPulses = *numberPulses +1;
						if (allPulsesMode == 0) break;
						prevPulse = true;
					}
					cntDown = 0;
				}
			}
			else if (cntUp == nSamplesUp-1)
			{
				if (possibleTstart == 0)
				{
					gsl_vector_set(tstartDER,*numberPulses,possibleTstart);
					gsl_vector_set(*flagTruncated,*numberPulses,1);
					*numberPulses = *numberPulses +1;
					if (allPulsesMode == 0) break;
					prevPulse = true;
				}
				else
				{
					gsl_vector_set(tstartDER,*numberPulses,possibleTstart);
					*numberPulses = *numberPulses +1;
					if (allPulsesMode == 0) break;
					prevPulse = true;
				}
				cntDown = 0;
			}

			i++;
			cntUp = cntUp+1;
		}
		else if ((gsl_vector_get(der,i) > adaptativethreshold) && (prevPulse == true))
		{
			i++;
			cntDown = 0;
		}
		else if (gsl_vector_get(der,i) <= adaptativethreshold)
		{
			if (prevPulse == true)
			{
				cntDown = cntDown+1;
				if (cntDown == 1)
				{
					possibleTend = i;
					if (nSamplesUp == 1)
					{
						gsl_vector_set(tendDER,*numberPulses-1,possibleTend);
						prevPulse = false;
					}
				}
				else if (cntDown == nSamplesUp)
				{
					gsl_vector_set(tendDER,*numberPulses-1,possibleTend);
					prevPulse = false;
				}
			}
			cntUp = 0;
			i++;
		}
	}

	if (lastPulse == false) *numberPulses = *numberPulses-1;

	// Just in case there is a truncated pulse at the end of the record whose tend has not been found and it is still 0.0
	// (and protected just in case there are no pulses)
	if ((*numberPulses != 0) && (gsl_vector_get(tendDER,*numberPulses-1) == 0.0))
	{
		gsl_vector_set(tendDER,*numberPulses-1,szRw-1);
		gsl_vector_set(*flagTruncated,*numberPulses-1,1.0);
	}

	// Obtain tstart of each pulse in the not filtered signal
	// Once pulses are located in the first derivative of the filtered event, tstarts in the derivative have to be
	// converted into the tstarts in the not filtered record
	for (int i=0;i<*numberPulses;i++)
	{
		temp = gsl_vector_subvector(der,gsl_vector_get(tstartDER,i),gsl_vector_get(tendDER,i)-gsl_vector_get(tstartDER,i));
		maxDER_index = gsl_vector_max_index(&temp.vector);
		if (maxDER_index+gsl_vector_get(tstartDER,i)-safetyMargin < 0)
		{
			//Truncated at the beginning
			gsl_vector_set(*tstartgslOUT,i,0.0);
			if (i == 0) gsl_vector_set(*flagTruncated,i,1.0);	// In order to not have more than one pulse truncated at the beginning
		}
		else
		{
			gsl_vector_set(*tstartgslOUT,i,gsl_vector_get(tstartDER,i)-safetyMargin);
		}
		gsl_vector_set(*tstartgslOUTNEW,i,gsl_vector_get(tstartDER,i));

		gsl_vector_set(*tmaxDERgsl,i,maxDER_index);
		gsl_vector_set(*maxDERgsl,i,gsl_vector_max(&temp.vector));
	}
	gsl_vector_memcpy(*tstartDERgsl,tstartDER);
	gsl_vector_memcpy(*tendDERgsl,tendDER);

	gsl_vector_free(tstartDER);
	gsl_vector_free(tendDER);

	if (*numberPulses > 0) *thereIsPulse = 1;

	return (EPOK);
}
/*xxxx end of SECTION 13 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 14 ************************************************************
* findPulses: This function is going to find the pulses in a record
*
* - Declare variables
* - First step to look for single pulses
* 	- Call medianKappaClipping
*   - Call findTstart
* - Handling with possible saturated pulses
* - If calibration mode and at least a pulse found
* 	- Get the pulseheight of each found pulse (in order to be used to build the pulse templates library)
* - else if production mode:
* 	- Iterative pulse searching and at least a pulse found
* 		- Call getB (to estimate the pulse height of each pulse in order to find the proper pulse template to build
* 		  the adjusted derivative)
* 		- If any elements in Bgsl is still -999 (after using Bprev) => findSePulses is not used (but the task does not stop)
*   	- Call findSePulses
* - Free allocate of GSL vectors
*
* Parameters:
* - vectorin: Input vector (not filtered row=record)
* - vectorinDER: Derivative of the low-pass filtered 'vectorin' (row=record)
* - tstart: Starting time of the pulses into the row=record
* - quality: Quality of the pulses into the row=record
* - energy: Estimated energy of the pulses in the row=record
*
* - nPulses: Number of found pulses
*
* - opmode: Calibration mode (0) or normal mode (1)
*
* - taufall: Fall time of the pulses (seconds)
* - scalefactor: Scale factor to apply to the fall time of the pulses in order to calculate the LPF box-car length
* - sizepulsebins: Size of the pulse in bins
* - samplingRate: Sampling rate
*
* - samplesup: Number of consecutive samples over the threshold to locate a pulse
* - nsgms: Number of Sigmas to establish the threshold
*
* - lb: Vector containing the baseline averaging length used for each pulse
* - lrs: Running sum length (equal to the 'Lrs' global_variable)
*
* - librarymatrix: Matrix where the values of Energy of each pulse model are stored
* - modelsmatrix: Matrix where all the pulse templates of the pulse templates library are going to be stored
*
* - safetymargintstart: Equal to the 'safetyMarginTstart' global variable
* - stopCriteriamkc: Used in medianKappaClipping (%)
* - kappamkc: Used in medianKappaClipping
* - levelprvpulse: Secondary pulses must be 1/levelPrvPulse times larger than the preceding pulse
****************************************/
int findPulses (
	gsl_vector *vectorin,
	gsl_vector *vectorinDER,
	gsl_vector **tstart,
	gsl_vector **quality,
	gsl_vector **energy,

	int *nPulses,

	int opmode,

	double taufall,
	double scalefactor,
	int sizepulsebins,
	double samplingRate,

	int samplesup,
	double nsgms,

	double lb,
	double lrs,

	gsl_matrix *librarymatrix,
	gsl_matrix *modelsmatrix,

	double safetymargintstart,
	double stopcriteriamkc,
	double kappamkc,
	double levelprvpulse,

	FILE * temporalFile,

	int index)
{
	char val[256];
	char val_aux[256];

 	int status=EPOK;
	string message = "";

	const double pi = 4.0 * atan(1.0);

	// Declare variables
	int pulseFound;
	double thresholdmediankappa;	// Threshold to look for pulses in the first derivative
	gsl_vector *model;				// Pulse which is going to be used as template (or model)
		                        	// (selected row of the PULSE column of the EUR-LIB extension
									// from the pulses templates library file)
									// It will be overwritten with its first derivative
		// To look for single pulses during the first step
	gsl_vector *tstartNOsmt = gsl_vector_alloc(vectorinDER->size);
	gsl_vector *tstartDERgsl = gsl_vector_alloc(vectorinDER->size);
	gsl_vector *tendDERgsl = gsl_vector_alloc(vectorinDER->size);
	gsl_vector *tmaxDERgsl = gsl_vector_alloc(vectorinDER->size);
	gsl_vector *maxDERgsl = gsl_vector_alloc(vectorinDER->size);
		// To look for secondary pulses during the second step
	gsl_vector *vectorDERComposed = gsl_vector_alloc(vectorinDER->size);
	gsl_vector *newPulsesgsl = gsl_vector_alloc(vectorinDER->size); // If a pulse is new => Look again for more pulses
	gsl_vector_set_zero(newPulsesgsl);
	gsl_vector *Lbgsl = gsl_vector_alloc(vectorinDER->size);	    // If there is no free-pulses segments longer than Lb=>
	gsl_vector_set_all(Lbgsl,lb);                                   // segments shorter than Lb will be useed and its length (< Lb)
	                                                                // must be used instead Lb in RS_filter
	gsl_vector *Bgsl;
	double Bprev = -999;
	gsl_vector *Bauxgsl;
	bool flagContinue = true;

	gsl_vector *SorSeorPrgsl = gsl_vector_alloc(vectorinDER->size);	// To know if a pulse is Single(=0), Primary(=1) or Secondary(=2)
	gsl_vector_set_zero(SorSeorPrgsl);								// For example to calculate the threshold in a different way (primary or secondary)

	gsl_vector_set_zero(*quality);
	gsl_vector_set_zero(*energy);									// Estimated energy of the single pulses
																	// In order to choose the proper pulse template to calculate
																	// the adjusted derivative and to fill in the Energy column
		                                                            // in the output FITS file
	// First step to look for single pulses
	if (medianKappaClipping (vectorinDER, kappamkc, stopcriteriamkc, nsgms, (int)(pi*samplingRate*taufall*scalefactor), &thresholdmediankappa, temporalFile))
	{
	    message = "Cannot run medianKappaClipping looking for single pulses";
	    EP_PRINT_ERROR(message,EPFAIL);
	}

	if (findTstart (vectorinDER, thresholdmediankappa, samplesup, safetymargintstart, 1, samplingRate, nPulses, &pulseFound, tstart, &tstartNOsmt, quality, &tstartDERgsl, &tmaxDERgsl, &maxDERgsl,&tendDERgsl, temporalFile))
	{
	  message = "Cannot run findTstart with two rows in models";
	  EP_PRINT_ERROR(message,EPFAIL);
	}

	for (int i=0;i<*nPulses;i++)
	{
		gsl_vector_set(SorSeorPrgsl,i,0);
		gsl_vector_set(newPulsesgsl,i,1);
	}

	// In order to look for saturated pulses
	double maxvectorNOTFIL = gsl_vector_max(vectorin);
	long indexmaxvectorNOTFIL = gsl_vector_max_index(vectorin);
	int cntstart = 0;
	int cntend = 0;
	int possiblestart0 = indexmaxvectorNOTFIL;
	int possibleend0 = indexmaxvectorNOTFIL;
	int prevsaturated = 0;
	gsl_vector *start0 = gsl_vector_alloc(vectorin->size-indexmaxvectorNOTFIL);
	gsl_vector *end0 = gsl_vector_alloc(vectorin->size-indexmaxvectorNOTFIL);
	gsl_vector_set_zero(end0);
	int numSaturated = 0;
	gsl_vector_view temp;					// In order to handle with gsl_vector_view (subvectors)
	gsl_vector *vectorAUX = gsl_vector_alloc(vectorin->size-indexmaxvectorNOTFIL);
	temp = gsl_vector_subvector(vectorin,indexmaxvectorNOTFIL,vectorin->size-indexmaxvectorNOTFIL);
	gsl_vector_memcpy(vectorAUX,&temp.vector);

	for (int i=0;i<vectorAUX->size;i++)
	{
		if ((gsl_vector_get(vectorAUX,i) == maxvectorNOTFIL) && (prevsaturated == 0))
		{
			if (cntstart == 0)
			{
				possiblestart0 = i;
			}

			cntstart = cntstart +1;

			if (cntstart == 2)	// HARDPOINT!!
			{
				gsl_vector_set(start0,numSaturated,possiblestart0+indexmaxvectorNOTFIL);
				prevsaturated = 1;
				numSaturated = numSaturated+1;
			}
			cntend = 0;
		}
		else if ((gsl_vector_get(vectorAUX,i) == maxvectorNOTFIL) && (prevsaturated == 1))
		{
			cntend = 0;
		}
		else
		{
			cntstart = 0;
			if (prevsaturated == 1)
			{
				cntend = cntend +1;
				if (cntend == 1)
				{
					possibleend0 = i;
				//}
				//else if (cntend == 5)	// HARDPOINT!!
				//{
					gsl_vector_set(end0,numSaturated-1,possibleend0+indexmaxvectorNOTFIL);
					prevsaturated = 0;
				}
			}
		}
	}

	for (int i=0;i<numSaturated;i++)
	{
		if ((i == numSaturated-1) && (gsl_vector_get(end0,i) == 0))		gsl_vector_set(end0,i,vectorin->size);

		for (int j=0;j<*nPulses;j++)
		{
			if (j != *nPulses-1)
			{
				if ((gsl_vector_get(start0,i)>gsl_vector_get(*tstart,j)+safetymargintstart) && (gsl_vector_get(end0,i)<gsl_vector_get(*tstart,j+1)+safetymargintstart))
				{
					gsl_vector_set(*quality,j,gsl_vector_get(*quality,j)+2);
				}
			}
			else
			{
				if ((gsl_vector_get(start0,i)>gsl_vector_get(*tstart,j)) && (gsl_vector_get(end0,i)<= vectorin->size))
				{
					gsl_vector_set(*quality,j,gsl_vector_get(*quality,j)+2);
				}
			}

		}
	}

	if ((opmode == 0) && (*nPulses != 0))
	{
		if (getB(vectorin, *tstart, *nPulses, &Lbgsl, sizepulsebins, &Bgsl, temporalFile))
		{
		    message = "Cannot run getB routine with opmode=0 & nPulses != 0";
		    EP_PRINT_ERROR(message,EPFAIL);
		}
		double energyaux = gsl_vector_get(*energy,0);
		for (int i=0;i<*nPulses;i++)
		{
			if (i != *nPulses-1)	// Not last pulse in the record
			{
				if (getPulseHeight(vectorin, gsl_vector_get(tstartNOsmt,i), gsl_vector_get(tstartNOsmt,i+1), 0, lrs, gsl_vector_get(Lbgsl,i), gsl_vector_get(Bgsl,i), sizepulsebins, &energyaux, temporalFile))
				{
				    message = "Cannot run getPulseHeight routine when pulse i=" + boost::lexical_cast<std::string>(i) + " is not the last pulse";
				    EP_PRINT_ERROR(message,EPFAIL);
				}
			}
			else
			{
				if (getPulseHeight(vectorin, gsl_vector_get(tstartNOsmt,i), gsl_vector_get(tstartNOsmt,i+1), 1, lrs, gsl_vector_get(Lbgsl,i), gsl_vector_get(Bgsl,i), sizepulsebins, &energyaux, temporalFile))
				{
				    message = "Cannot run getPulseHeight routine when pulse i=" + boost::lexical_cast<std::string>(i) + " is the last pulse";
				    EP_PRINT_ERROR(message,EPFAIL);
				}
			}
			gsl_vector_set(*energy,i,energyaux);
		}
	}
	else if ((opmode == 1) && (*nPulses != 0))  // Iterative pulse searching
	{
		gsl_vector_memcpy(vectorDERComposed, vectorinDER); 	// Some parts will be overwritten
		model = gsl_vector_alloc(modelsmatrix->size2);

		int nPulsesRowAux;
		int nNewPulses;

		do
		{
			// To estimate the pulse height of each pulse
			// Sum of the Lb digitized data samples of a pulse-free interval immediately before the current pulse, B
			if (getB(vectorin, *tstart, *nPulses, &Lbgsl, sizepulsebins, &Bgsl, temporalFile))
			{
			    message = "Cannot run getB routine with opmode=1 & nPulses != 0";
			    EP_PRINT_ERROR(message,EPFAIL);
			}
			Bauxgsl = gsl_vector_alloc(Bgsl->size);
			gsl_vector_set_all(Bauxgsl,999);
			gsl_vector_add(Bauxgsl,Bgsl);
			if (gsl_vector_isnull(Bauxgsl) == 1) // All the elements of Bgsl are -999
			{
				if (Bprev == -999)	break;	// Out of the do_while
				gsl_vector_set_all(Bgsl,Bprev);
			}
			else
			{
				for (int i=0;i<*nPulses;i++)
				{
				    if (gsl_vector_get(Bgsl,i) == -999)
					{
						gsl_vector_set(Bgsl,i,Bprev);
						if (Bprev == -999)	flagContinue = false;
					}
				}
				if (flagContinue == false) break;	// Out of the do_while

			}
			gsl_vector_free(Bauxgsl);
			Bprev = gsl_vector_get(Bgsl,*nPulses-1);

			nPulsesRowAux = *nPulses;
			if (findSePulses(vectorin, vectorinDER, &vectorDERComposed,
					 thresholdmediankappa,
					 tstart, &tstartNOsmt, quality, energy, &tstartDERgsl, &tmaxDERgsl, &tendDERgsl, &maxDERgsl,
					 &SorSeorPrgsl, &newPulsesgsl,
					 nPulses,
					 start0, end0, numSaturated,
					 taufall, scalefactor, sizepulsebins, samplingRate,
					 samplesup, nsgms,
					 Bgsl, lrs, Lbgsl,
					 librarymatrix, modelsmatrix, model,
					 safetymargintstart, stopcriteriamkc, kappamkc, levelprvpulse,
					 temporalFile,
					 index))
			{
				message = "Cannot run findPulses routine with opmode=1 & nPulses != 0";
				EP_PRINT_ERROR(message,EPFAIL);
			}

			nNewPulses = *nPulses - nPulsesRowAux;
			for (int j=0;j<*nPulses;j++)
			{
				if ((gsl_vector_get(newPulsesgsl,j) == 2.0) && (gsl_vector_get(SorSeorPrgsl,j) == 1.0))
				{
					nNewPulses = nNewPulses-1;
				}
			}

			for (int j=0;j<*nPulses;j++)
			{
				gsl_vector_set(newPulsesgsl,j,gsl_vector_get(newPulsesgsl,j)-1.0);
			}
		} while (nNewPulses > 0);
	}

	// Free allocate of GSL vectors
	gsl_vector_free(tstartDERgsl);
	gsl_vector_free(tendDERgsl);
	gsl_vector_free(tmaxDERgsl);
	gsl_vector_free(maxDERgsl);
	gsl_vector_free(SorSeorPrgsl);
	gsl_vector_free(newPulsesgsl);
	gsl_vector_free(start0);
	gsl_vector_free(end0);
	gsl_vector_free(vectorAUX);
	gsl_vector_free(tstartNOsmt);
	gsl_vector_free(Lbgsl);
	if (opmode == 1)
	{
		gsl_vector_free(model);
		gsl_vector_free(vectorDERComposed);
	}

	return(EPOK);
}
/*xxxx end of SECTION 14 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 15 ************************************************************
* findSePulses: This function looks for the secondary and primary pulses in the input vector.
*
* - Declare variables & allocate GSL vectors
* - The thresholdmediankappaSingle input parameter is applied to the set of templates (from the library) in order to get
	the value of the first sample over the threshold and the index of that sample
* - Build 'vectorDERComposed' to look for secondary pulses into the adjusted first derivative of the whole event.
*   The adjusted first derivative is built by using the differences between each pulse first derivative and
*   the scaled and shifted pulse template first derivative. 'vectorDERComposed' is initialized once at the beginning
*   and iteratively it is updated by overwriting with the difference between a pulse and its template (it is in fact
*   an input/output parameter of the function)
*   - Estimate the amplitude of the first sample over the threshold to choose the proper pulse template ("model")
	  of the pulse templates library ("models")
*   - Find the proper pulse template in the pulse templates library
* 	- Cut the record part corresponding to a pulse
* 	- Difference between each pulse and the corresponding shifted and scale pulse template:
*		- Scale the template according to the maximum of the pulse
*   	- Calculate the shift between the template and the pulse
*       - Allocate the shifted-scaled model and the difference between the pulse and the shifted-scaled model
*       - Shift the scaled template
*   - Overwrite some parts of 'vectorDERComposed' where there is a pulse with 'diff'
*   - Free allocate GSL vectors
* - Look for secondary pulses into 'vectorDERComposed' (call medianKappaClipping and findTstart)
* 	- Secondary pulses must be larger than 1/levelPrvPulse times the preceding pulse to avoid being considered
*     as noise
* 	- Add the single and the secondary pulses found
* 	- Order the found pulses (and all their data) according to 'tstartNOsmt' instead 'tstart'
*
* Parameters:
* - vectorin: Input vector (not filtered row=event)
* - vectorinDER: Derivative of the low-pass filtered 'vectorin' (row=event)
* - vectorinDERComposed: Adjusted first derivative
*
* - thresholdmediankappaSingle: Threshold previously used to find single pulses
*
* - tstart: Starting time of the pulses into the row=record
* - tstartNOsmt: Starting time of the pulses into the row=record without taking into account the safety margin
* - quality: Quality of the pulses into the row=record
* - energy: Estimated energy of the pulses in the row=record
* - tstartDER: Point where the first derivative crosses over the threshold
* - tmaxDER: Point where the first derivative is maximum
* - tendDER: Point where the first derivative crosses below the threshold
* - maxDER: Maximum of the first derivative
*
* - SorSeorPr: Vector whose elements give information about the type of the pulse:
*              Single(0), Secondary(2) or Primary(1)
* - newPulses: Vector whose elements give information about new found pulses(1) or not(0)
*
* - nPulses: Number of single pulses previously found in the input vector ('vectorin')
* - nPulsesNew: Number of new found pulses
*
* - startsaturated: Vector storing the sample where the pulse i starts to be saturated
* - endsaturated: Vector storing the sample where the pulse i ends to be saturated
*
* - taufall: Fall time of the pulses (seconds)
* - scalefactor: Scale factor to apply to the fall time of the pulses in order to calculate the box-car length
* - sizepulse: Size of the pulse in bins
* - samplingRate: Sampling rate
*
* - samplesup: Number of consecutive samples over the threshold to locate a pulse
* - nsgms: Number of Sigmas to establish the threshold

* - B: Vector whose elements are, in general, the sum of the Lb digitized data samples of a pulse-free interval
*      immediately before each pulse
* - lrs: Running sum length (equal to the 'Lrs' global_variable)
* - lb: Vector containing the baseline averaging length used for each pulse
*
* - library: Matrix where the values of Energy of each pulse template are stored
* - models: Matrix where all the pulse templates of the pulse templates library are going to be stored
* - model: Pulse which is going to used as template
*          (selected row of the PULSE column of the EUR-LIB extension from the pulses templates library file)
*
* - safetymargintstart: Equal to the 'safetyMarginTstart' global variable
* - stopCriteriamkc: Used in medianKappaClipping (%)
* - kappamkc: Used in medianKappaClipping
* - levelprvpulse: Secondary pulses must be 1/levelPrvPulse times larger than the preceding pulse
****************************************/
int findSePulses(
	gsl_vector *vectorin,
	gsl_vector *vectorinDER,
	gsl_vector **vectorinDERComposed,

	double thresholdmediankappaSingle,

	gsl_vector **tstart,
	gsl_vector **tstartNOsmt,
	gsl_vector **quality,
	gsl_vector **energy,
	gsl_vector **tstartDER,
	gsl_vector **tmaxDER,
	gsl_vector **tendDER,
	gsl_vector **maxDER,

	gsl_vector **SorSeorPr,
	gsl_vector **newPulses,

	int *nPulses,
	//int *nPulsesNew,

	gsl_vector *startsaturated,
	gsl_vector *endsaturated,
	int nSaturated,

	double taufall,
	double scalefactor,
	int sizepulse,
	double samplingRate,

	int samplesup,
	double nsgms,

	gsl_vector *B,
	double lrs,
	gsl_vector *lb,

	gsl_matrix *library,
	gsl_matrix *models,
	gsl_vector *model,

	double safetymargintstart,
	double stopCriteriamkc,
	double kappamkc,
	double levelprvpulse,

	FILE * temporalFile,

	int indice)
{
	char val[256];
	char val_aux[256];

	int status=EPOK;
	string message = "";

	const double pi = 4.0 * atan(1.0);

	// Declare variables & allocate GSL vectors
	double firstSampleOverThreshold;
	long index_firstSampleOverThreshold;
	double firstSampleOverThreshold_Model;
	long index_firstSampleOverThreshold_Model;
	gsl_vector *firstSamplesgsl = gsl_vector_alloc(models->size1);
	gsl_vector_set_all(firstSamplesgsl,0);
	gsl_vector *index_firstSamplesgsl = gsl_vector_alloc(models->size1);
	gsl_vector_set_all(index_firstSamplesgsl,0);

    int ind;
    int ind1;

    int istherePulse;
    	// To look for secondary pulses
    int lastOne;				// If 1 => The pulse is the last one of the row=record (or the only one)
    gsl_vector *modelScaled;	// Pulse template scaled according to each pulse
    long shift;					// Shift between each pulse and the pulse template
    gsl_vector *modelShifted;	// The scaled pulse template is also shifted to align it with each pulse
    gsl_vector *diff;			// Pulse - Shifted scaled pulse template
    int limInf, limSup;
    gsl_vector *limSup_vector=gsl_vector_alloc(*nPulses);
    gsl_vector_set_zero(limSup_vector);
    double thresholdmediankappaSecondary;
    int nSecondaryPulses = 0;
    gsl_vector *tstartSecondary =  gsl_vector_alloc(vectorinDER->size);
    gsl_vector *tstartSecondaryNOsmt =  gsl_vector_alloc(vectorinDER->size);
    gsl_vector *qualitySecondary =  gsl_vector_alloc(vectorinDER->size);
    gsl_vector *tstartDERSecondary =  gsl_vector_alloc(vectorinDER->size);
    gsl_vector *tmaxDERSecondary =  gsl_vector_alloc(vectorinDER->size);
    gsl_vector *maxDERSecondary =  gsl_vector_alloc(vectorinDER->size);
    gsl_vector *tendDERSecondary =  gsl_vector_alloc(vectorinDER->size);

    int indSaturated = 0;

    double pulseheight;

	// Auxiliary variables
	gsl_vector *pulse;
	gsl_vector_view temp;					// In order to handle with gsl_vector_view (subvectors)

	// The thresholdmediankappaSingle input parameter is applied to the set of templates (from the library) in order to get
	// the value of the first sample over the threshold and the index of that sample
	if (firstSampleModels (models, thresholdmediankappaSingle, &firstSamplesgsl, &index_firstSamplesgsl, temporalFile))
	{
	    message = "Cannot run firstSampleModels routine to get the value of the first sample over the threshold and the index of that sample";	
	    EP_PRINT_ERROR(message,EPFAIL);
	}

	// Build 'vectorinDERComposed'
	modelScaled = gsl_vector_alloc(models->size2);
	for (int i=0;i<*nPulses;i++)
	{
		lastOne = 0;
		if ((*nPulses == 1) || ((*nPulses !=1) && (i == *nPulses-1)))	lastOne = 1;

		if (gsl_vector_get(*newPulses,i) == 1)
		{
			// For the moment, we need to store in the output FITS file the estimated energy
			if (getPulseHeight(vectorin, gsl_vector_get(*tstartNOsmt,i), gsl_vector_get(*tstartNOsmt,i+1), lastOne, lrs, gsl_vector_get(lb,i), gsl_vector_get(B,i), sizepulse, &pulseheight, temporalFile))
			{
			    message = "Cannot run getPulseHeight routine for pulse i=" + boost::lexical_cast<std::string>(i) + " when newPulses = 1";
			    EP_PRINT_ERROR(message,EPFAIL);
			}
			gsl_vector_set(*energy,i,pulseheight);

			// Cut the record part corresponding to a pulse
			if (i != *nPulses-1)		// Neither single pulse nor last pulse
			{
				limInf = gsl_vector_get(*tstartDER,i);
				limSup = gsl_vector_get(*tstartDER,i+1);

				pulse = gsl_vector_alloc(limSup-limInf);
				temp = gsl_vector_subvector(vectorinDER,limInf,limSup-limInf);

				gsl_vector_memcpy(pulse, &temp.vector);
			}
			else						// Either single pulse or last pulse
			{
				limInf = gsl_vector_get(*tstartDER,i);
				limSup = vectorin->size;

				pulse = gsl_vector_alloc(limSup-limInf);
				temp = gsl_vector_subvector(vectorinDER,limInf,limSup-limInf);

				gsl_vector_memcpy(pulse, &temp.vector);
			}

			// Estimate the amplitude of the first sample over the threshold to choose the proper pulse template ("model")
			// of the pulse templates library ("models")
			index_firstSampleOverThreshold = 0; //Because the 'pulse' starts in tstartDER
			firstSampleOverThreshold = gsl_vector_get(pulse,index_firstSampleOverThreshold);

			// Find the proper pulse model in the pulse templates library
			if (find_model1stSample(firstSampleOverThreshold, firstSamplesgsl, models, &model, temporalFile))
			{
			    message = "Cannot run find_model1stSample routine for pulse i=" + boost::lexical_cast<std::string>(i) + " when newPulses = 1";
			    EP_PRINT_ERROR(message,EPFAIL);
			}
			for (int j=0;j<models->size2;j++)
			{
				if (gsl_vector_get(model,j) > thresholdmediankappaSingle)
				{
					index_firstSampleOverThreshold_Model = j;
					firstSampleOverThreshold_Model = gsl_vector_get(model,j);

					break;
				}
			}

			// Scale the template according to the maximum of the pulse
			gsl_vector_memcpy(modelScaled,model);
			gsl_vector_scale(modelScaled,firstSampleOverThreshold/firstSampleOverThreshold_Model);

			// Calculate the shift between the template and the pulse
			shift = index_firstSampleOverThreshold-index_firstSampleOverThreshold_Model;

			// Allocate the shifted-scaled template and the difference between the pulse and the shifted-scaled template
			modelShifted = gsl_vector_alloc(pulse->size);
			diff = gsl_vector_alloc(pulse->size);

			// Shift the scaled template according to 'shift' (load the 'modelShifted' vector)
			if (shift == 0)
			{
				if (modelShifted->size < modelScaled->size)		// 'modelShifted' is a subvector of 'modelScaled'
				{
					temp = gsl_vector_subvector(modelScaled,0,modelShifted->size);
					gsl_vector_memcpy(modelShifted, &temp.vector);
				}
				else											// 'modelShifted' is 'modelScaled' extended
				{
					gsl_vector_set_all(modelShifted,gsl_vector_get(modelScaled,modelScaled->size-1));
					for (int k=0;k<modelScaled->size;k++)
					{
						gsl_vector_set(modelShifted,k,gsl_vector_get(modelScaled,k));
					}
				}
			}
			else if (shift > 0)		// Pulse template must be delayed
			{
				gsl_vector_set_all(modelShifted,gsl_vector_get(modelScaled,0));
				for (int k=0;k<modelShifted->size-shift;k++)
				{
					if (k < modelScaled->size)
					{
						gsl_vector_set(modelShifted,k+shift,gsl_vector_get(modelScaled,k));
					}
				}
			}
			else if (shift < 0)		// Pulse template must be moved forward
			{
				gsl_vector_set_all(modelShifted,gsl_vector_get(modelScaled,modelScaled->size-1));
				for (int k=0;k<modelShifted->size;k++)
				{
					if (k+fabs(shift) < modelScaled->size)
					{
						gsl_vector_set(modelShifted,k,gsl_vector_get(modelScaled,k+fabs(shift)));
					}
				}
			}

			// Difference between each pulse and the corresponding shifted and scaled pulse template
			gsl_vector_memcpy(diff,pulse);
			gsl_vector_sub(diff,modelShifted);

			// Overwrite some parts of 'vectorinDERComposed' where there is a pulse with 'diff'
			// (and with 0 some samples)
			if (i != *nPulses-1)		// Nor single pulse nor last pulse
			{
				limInf = gsl_vector_get(*tstartDER,i);
				limSup = gsl_vector_get(*tstartDER,i+1);

				for (int j=limInf;j<limSup;j++)
				{
					if (gsl_vector_get(diff,j-limInf) < -1e10)
					{
						gsl_vector_set(*vectorinDERComposed,j,-1e10);
					}
					else if (gsl_vector_get(diff,j-limInf) > 1e10)
					{
						gsl_vector_set(*vectorinDERComposed,j,1e10);
					}
					else
					{
						gsl_vector_set(*vectorinDERComposed,j,gsl_vector_get(diff,j-limInf));
					}
				}
			}
			else						// Or single pulse or last pulse
			{
				if (gsl_vector_get(*quality,i) == 1)
				{
					limInf = gsl_vector_get(*tstartDER,i);
					limSup = vectorin->size;	//HARDPOINT!!!

					for (int j=limInf;j<limSup;j++)
					{
						if (gsl_vector_get(diff,j-limInf) < -1e10)
						{
							gsl_vector_set(*vectorinDERComposed,j,-1e10);
						}
						else if (gsl_vector_get(diff,j-limInf) > 1e10)
						{
							gsl_vector_set(*vectorinDERComposed,j,1e10);
						}
						else
						{
							gsl_vector_set(*vectorinDERComposed,j,gsl_vector_get(diff,j-limInf));
						}
					}
				}
				else
				{
					limInf = gsl_vector_get(*tstartDER,i);
					limSup = vectorin->size;

					for (int j=limInf;j<limSup;j++)
					{
						if (gsl_vector_get(diff,j-limInf) < -1e10)
						{
							gsl_vector_set(*vectorinDERComposed,j,-1e10);
						}
						else if (gsl_vector_get(diff,j-limInf) > 1e10)
						{
							gsl_vector_set(*vectorinDERComposed,j,1e10);
						}
						else
						{
							gsl_vector_set(*vectorinDERComposed,j,gsl_vector_get(diff,j-limInf));
						}
					}
				}
			}

			for (int j=0;j<nSaturated;j++)
			{
				for (int k=gsl_vector_get(startsaturated,j);k<gsl_vector_get(endsaturated,j);k++)
				{
					gsl_vector_set(*vectorinDERComposed,k,0.0);
				}
			}

			// Free allocate of GSL vectors
			gsl_vector_free(pulse);
			gsl_vector_free(modelShifted);
			gsl_vector_free(diff);
		}
	}

	gsl_vector_free(modelScaled);

	// Look for secondary pulses into 'vectorinDERComposed'
	if (medianKappaClipping (*vectorinDERComposed, kappamkc, stopCriteriamkc, nsgms, (int)(pi*samplingRate*taufall*scalefactor), &thresholdmediankappaSecondary, temporalFile))
	{
	    message = "Cannot run medianKappaClipping routine to look for secondary pulses into 'vectorinDERComposed'";
	    EP_PRINT_ERROR(message,EPFAIL);
	}

	if (findTstart (*vectorinDERComposed, thresholdmediankappaSecondary, samplesup, safetymargintstart, 1, samplingRate, &nSecondaryPulses, &istherePulse, &tstartSecondary, &tstartSecondaryNOsmt, &qualitySecondary, &tstartDERSecondary, &tmaxDERSecondary, &maxDERSecondary, &tendDERSecondary,temporalFile))
	{
	    message = "Cannot run findTstart routine to look for secondary pulses into 'vectorinDERComposed'";
	    EP_PRINT_ERROR(message,EPFAIL);
	  
	}
	ind = 0;
	ind1 = 0;
	if (nSecondaryPulses != 0)
	{
		gsl_vector *trueSecondary = gsl_vector_alloc(nSecondaryPulses);
		gsl_vector_set_all(trueSecondary,1);

		for (int i=0;i<nSecondaryPulses;i++)
		{
			for (int j=0;j<*nPulses;j++)
			{
				if ((gsl_vector_get(tstartSecondary,i) == gsl_vector_get(*tstart,j)) || (gsl_vector_get(tstartSecondary,i) == (gsl_vector_get(*tstart,j)+1)) || (gsl_vector_get(tstartDERSecondary,i) == gsl_vector_get(limSup_vector,j)))
				{
					gsl_vector_set(trueSecondary,i,0);

					break;
				}
				else
				{
					gsl_vector_set(trueSecondary,i,1);
				}
			}

			if (gsl_vector_get(trueSecondary,i) == 1)
			{
				if (gsl_vector_get(tstartSecondary,i) == 0)
				{
					gsl_vector *pulseaux = gsl_vector_alloc(samplesup);
					temp = gsl_vector_subvector(*vectorinDERComposed,0,pulseaux->size);
					gsl_vector_memcpy(pulseaux, &temp.vector);
					if (gsl_vector_isnull(pulseaux) == 1)
					{
						gsl_vector_set(trueSecondary,i,0);
					}
					gsl_vector_free(pulseaux);
				}
			}

			// Secondary pulses must be larger than 1/levelPrvPulse times the preceding pulse to
			// avoid being considered as noise
			if (gsl_vector_get(trueSecondary,i) == 1)
			{
				if (*nPulses == 1)	// Only a single pulse
				{
					if (gsl_vector_get(maxDERSecondary,i) > gsl_vector_get(*maxDER,0)/levelprvpulse)
					{
						gsl_vector_set(trueSecondary,i,1);
					}
					else
					{
						gsl_vector_set(trueSecondary,i,0);
					}
				}
				else				// More than one a single pulse
				{
					for (int j=0;j<*nPulses;j++)
					{
						if (j != *nPulses-1)
						{
							if ((gsl_vector_get(*tstart,j)<gsl_vector_get(tstartSecondary,i)) && (gsl_vector_get(tstartSecondary,i)<gsl_vector_get(*tstart,j+1)))
							{
								if (gsl_vector_get(maxDERSecondary,i) > gsl_vector_get(*maxDER,j)/levelprvpulse)
								{
									gsl_vector_set(trueSecondary,i,1);
								}
								else
								{
									gsl_vector_set(trueSecondary,i,0);
								}
								break;
							}
						}
						else
						{
							if ((gsl_vector_get(*tstart,j)<gsl_vector_get(tstartSecondary,i)) && (gsl_vector_get(tstartSecondary,i)<vectorin->size))
							{
								if (gsl_vector_get(maxDERSecondary,i) > gsl_vector_get(*maxDER,j)/levelprvpulse)
								{
									gsl_vector_set(trueSecondary,i,1);
								}
								else
								{
									gsl_vector_set(trueSecondary,i,0);
								}
							}
						}
					}
				}
			}
		}

		for (int i=0;i<nSecondaryPulses;i++)
		{
			if (gsl_vector_get(trueSecondary,i) == 1)
			{
				gsl_vector_set(*tstart,ind+*nPulses,gsl_vector_get(tstartSecondary,ind1));
				gsl_vector_set(*tstartNOsmt,ind+*nPulses,gsl_vector_get(tstartSecondaryNOsmt,ind1));
				gsl_vector_set(*tstartDER,ind+*nPulses,gsl_vector_get(tstartDERSecondary,ind1));
				gsl_vector_set(*tmaxDER,ind+*nPulses,gsl_vector_get(tmaxDERSecondary,ind1));
				gsl_vector_set(*maxDER,ind+*nPulses,gsl_vector_get(maxDERSecondary,ind1));
				gsl_vector_set(*tendDER,ind+*nPulses,gsl_vector_get(tendDERSecondary,ind1));
				gsl_vector_set(*SorSeorPr,ind+*nPulses,2);
				gsl_vector_set(*newPulses,ind+*nPulses,2);
				gsl_vector_set(*quality,ind+*nPulses,gsl_vector_get(qualitySecondary,ind1));
				ind = ind+1;
				ind1 = ind1+1;
			}
			else
			{
				ind1 = ind1+1;
			}
		}

		gsl_vector_free(trueSecondary);
	}
	nSecondaryPulses = ind;

	// Add the single and the secondary pulses found
	*nPulses = *nPulses + nSecondaryPulses;

	if (nSecondaryPulses != 0)
	{
		// Order the found pulses (and all their data) according to 'tstartNOsmt'
		gsl_vector *tstartauxNOsmt = gsl_vector_alloc(*nPulses); // 'tstartNOsmt' subvector
		temp = gsl_vector_subvector(*tstartNOsmt,0,*nPulses);
		gsl_vector_memcpy(tstartauxNOsmt, &temp.vector);
		gsl_permutation *perm = gsl_permutation_alloc(tstartauxNOsmt->size);
		// 'gsl_sort_vector_index' indirectly sorts the elements of the vector v into ascending order, storing the resulting
		// permutation in p. The elements of p give the index of the vector element which would have been stored in that position
		// if the vector had been sorted in place. The first element of p gives the index of the least element in v, and the last
		// element of p gives the index of the greatest element in v. The vector v is not changed.
		// Example: tstartaux=(5200 6000 200 3000) tauxsorted=(200 3000 5200 6000) perm=(2 3 0 1)
		gsl_sort_vector_index(perm,tstartauxNOsmt);
		gsl_vector *tstartaux = gsl_vector_alloc(vectorin->size);
		gsl_vector *tstartNOsmtaux = gsl_vector_alloc(vectorin->size);
		gsl_vector *SorSeorPraux = gsl_vector_alloc(vectorin->size);
		gsl_vector *newPulsesaux = gsl_vector_alloc(vectorin->size);
		gsl_vector *energyaux = gsl_vector_alloc(vectorin->size);
		gsl_vector *qualityaux = gsl_vector_alloc(vectorin->size);
		gsl_vector *tstartDERaux = gsl_vector_alloc(vectorin->size);
		gsl_vector *tmaxDERaux = gsl_vector_alloc(vectorin->size);
		gsl_vector *maxDERaux = gsl_vector_alloc(vectorin->size);
		gsl_vector *tendDERaux = gsl_vector_alloc(vectorin->size);
		for (int i=0;i<*nPulses;i++)
		{
			gsl_vector_set(tstartaux,i,gsl_vector_get(*tstart,gsl_permutation_get(perm,i)));
			gsl_vector_set(tstartNOsmtaux,i,gsl_vector_get(*tstartNOsmt,gsl_permutation_get(perm,i)));
			gsl_vector_set(SorSeorPraux,i,gsl_vector_get(*SorSeorPr,gsl_permutation_get(perm,i)));
			gsl_vector_set(newPulsesaux,i,gsl_vector_get(*newPulses,gsl_permutation_get(perm,i)));
			gsl_vector_set(energyaux,i,gsl_vector_get(*energy,gsl_permutation_get(perm,i)));
			gsl_vector_set(qualityaux,i,gsl_vector_get(*quality,gsl_permutation_get(perm,i)));
			gsl_vector_set(tstartDERaux,i,gsl_vector_get(*tstartDER,gsl_permutation_get(perm,i)));
			gsl_vector_set(tmaxDERaux,i,gsl_vector_get(*tmaxDER,gsl_permutation_get(perm,i)));
			gsl_vector_set(maxDERaux,i,gsl_vector_get(*maxDER,gsl_permutation_get(perm,i)));
			gsl_vector_set(tendDERaux,i,gsl_vector_get(*tendDER,gsl_permutation_get(perm,i)));
		}
		gsl_vector_memcpy(*tstart,tstartaux);
		gsl_vector_memcpy(*tstartNOsmt,tstartNOsmtaux);
		gsl_vector_memcpy(*SorSeorPr,SorSeorPraux);
		gsl_vector_memcpy(*newPulses,newPulsesaux);
		gsl_vector_memcpy(*energy,energyaux);
		gsl_vector_memcpy(*quality,qualityaux);
		gsl_vector_memcpy(*tstartDER,tstartDERaux);
		gsl_vector_memcpy(*tmaxDER,tmaxDERaux);
		gsl_vector_memcpy(*maxDER,maxDERaux);
		gsl_vector_memcpy(*tendDER,tendDERaux);
		gsl_vector_free(tstartauxNOsmt);
		gsl_permutation_free(perm);
		gsl_vector_free(tstartaux);
		gsl_vector_free(tstartNOsmtaux);
		gsl_vector_free(SorSeorPraux);
		gsl_vector_free(newPulsesaux);
		gsl_vector_free(energyaux);
		gsl_vector_free(qualityaux);
		gsl_vector_free(tstartDERaux);
		gsl_vector_free(tmaxDERaux);
		gsl_vector_free(maxDERaux);
		gsl_vector_free(tendDERaux);
	}

	gsl_vector_free(limSup_vector);
	gsl_vector_free(firstSamplesgsl);
	gsl_vector_free(index_firstSamplesgsl);

	return(EPOK);
}
/*xxxx end of SECTION 15 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
