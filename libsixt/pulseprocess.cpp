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

   Copyright 2014:  PULSEPROCESS has been developed by the INSTITUTO DE FISICA DE
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
*  26/01/15   In 'getB': tstart -> Convertion to 'int'
*  30/01/15   Changes in 'findSePulses' in order to work properly with close pulses (not piled-up) and truncated pulses at the beginning
*             and at the end of a record
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
 - 15. derivative

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
	//Declare variables
	gsl_vector *invectorAux;
	gsl_vector *invectorAux1;
	double cutFreq;	//Frequency domain
	int boxLength;	//Time domain
	double boxSum = 0.0;

	//Define the LPF (frequency domain) and the box-car function (time domain)
	cutFreq = 2 * (1/(2*pi*tau_fall));	//According to Jan, sinc(f1)=0.6 where f1=1/(2pi*tau_fall)
						//sinc(0.5)~0.6 => f1~0.5
						//sinc(fc)=0, sinc(1)=0 => fc=1
						//fc=kf1 => fc~2f1
	boxLength =(int) ((1/cutFreq) * sampleRate);

	if (boxLength < 1)boxLength = 1;
	if (boxLength >= szVct)return(4);

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
	
	if (boxLength == 1)return(3);

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
int medianKappaClipping (gsl_vector *invector, double kappa, double stopCriteria, double nSigmas, int boxLPF, double *threshold)
{
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
int getB(gsl_vector *vectorin, gsl_vector *tstart, int nPulses, gsl_vector **lb, int sizepulse, gsl_vector **B)
{
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
		gsl_vector_set(tstart,i,(int)(gsl_vector_get(tstart,i)));
	}

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
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
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
							input = gsl_vector_alloc(gsl_vector_get(*lb,0));
							temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,j)-gsl_vector_get(*lb,0),gsl_vector_get(*lb,0));
							gsl_vector_memcpy(input, &temp.vector);

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
							input = gsl_vector_alloc(gsl_vector_get(tstart,j)-tendprev-1);
							temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,j)-tendprev-1);
							gsl_vector_memcpy(input, &temp.vector);
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
										input = gsl_vector_alloc(gsl_vector_get(*lb,i));
										temp = gsl_vector_subvector(vectorin,gsl_vector_get(tstart,j+1)-gsl_vector_get(*lb,i),gsl_vector_get(*lb,i));
										gsl_vector_memcpy(input, &temp.vector);

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
										input = gsl_vector_alloc(gsl_vector_get(tstart,j+1)-tendprev-1);
										temp = gsl_vector_subvector(vectorin,tendprev+1,gsl_vector_get(tstart,j+1)-tendprev-1);
										gsl_vector_memcpy(input, &temp.vector);

										// Sum all the elements of input
										if (gsl_vector_Sumsubvector(input,0,input->size,&Baux))
										{
											message = "Cannot run gsl_vector_Sumsubvector routine when no first pulse in row & tstart-tendprev < lb";
											EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
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
										input = gsl_vector_alloc(vectorin->size-tendprev-1);
										temp = gsl_vector_subvector(vectorin,tendprev+1,vectorin->size-tendprev-1);
										gsl_vector_memcpy(input, &temp.vector);
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
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
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
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
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
								EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
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
								EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
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
							input = gsl_vector_alloc(vectorin->size-tendprev-1);
							temp = gsl_vector_subvector(vectorin,tendprev+1,vectorin->size-tendprev-1);
							gsl_vector_memcpy(input, &temp.vector);
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
							gsl_vector_set(*lb,i,gsl_vector_get(*lb,i-1));

							break;
						}
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
int getPulseHeight(gsl_vector *vectorin, double tstart, double tstartnext, int lastPulse, double lrs, double lb, double B, int sizepulse, double *pulseheight)
{
	// Declare variables
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
int RS_filter (gsl_vector *vector, double lrs, double lb, double B, double *pulseheight)
{
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
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
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
			    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
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
* - modelsvalues: Vector where the values of Energy of each pulse template are stored
* - models: Matrix where all the pulse templates of the pulse templates library are going to be stored
* - modelFound: Found template of the pulse whose pulse height is 'ph'
****************************************/
int find_model_energies(double ph, ReconstructInitSIRENA *reconstruct_init,gsl_vector **modelFound)
{
	string message = "";

	gsl_vector *modelFound_aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);

	long nummodels = reconstruct_init->library_collection->ntemplates;

	if (ph < gsl_vector_get(reconstruct_init->library_collection->energies,0))
	{
		gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_B0[0].ptemplate);
		gsl_vector_scale(modelFound_aux,ph/gsl_vector_get(reconstruct_init->library_collection->energies,0));
	}
	else if (ph > gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1))
	{
		gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_B0[nummodels-1].ptemplate);
		gsl_vector_scale(modelFound_aux,ph/gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1));
	}
	else
	{
		for (int i=0;i<nummodels;i++)
		{
			if (fabs(ph-gsl_vector_get(reconstruct_init->library_collection->energies,i))<1e-6)
			{
				gsl_vector_memcpy(modelFound_aux,reconstruct_init->library_collection->pulse_templates_B0[i].ptemplate);
				gsl_vector_scale(modelFound_aux,ph/gsl_vector_get(reconstruct_init->library_collection->energies,i));

				break;
			}
			else if ((ph > gsl_vector_get(reconstruct_init->library_collection->energies,i)) && (ph < gsl_vector_get(reconstruct_init->library_collection->energies,i+1)))
			{
				// Interpolate between the two corresponding rows in "models"
				gsl_vector *modelAux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);
				gsl_vector_set_zero(modelAux);

				if (interpolate_model(&modelAux,ph,reconstruct_init->library_collection->pulse_templates_B0[i].ptemplate,gsl_vector_get(reconstruct_init->library_collection->energies,i),
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
int find_model_maxDERs(double maxDER, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound)
{
	string message = "";
	gsl_vector *modelFound_aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);

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
		for (int i=0;i<nummodels;i++)
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
				//gsl_vector *modelAux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);
				//gsl_vector_set_zero(modelAux);
				gsl_vector_set_zero(modelFound_aux);

				//if (interpolate_model(&modelAux,maxDER,reconstruct_init->library_collection->pulse_templates_filder[i].ptemplate,gsl_vector_get(reconstruct_init->library_collection->maxDERs,i),
				if (interpolate_model(&modelFound_aux,maxDER,reconstruct_init->library_collection->pulse_templates_filder[i].ptemplate,gsl_vector_get(reconstruct_init->library_collection->maxDERs,i),
						reconstruct_init->library_collection->pulse_templates_filder[i+1].ptemplate,gsl_vector_get(reconstruct_init->library_collection->maxDERs,i+1)))
				{
				    message = "Cannot run interpolate_model with two rows in models";
				    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				//cout<<"Modelo interpolado entre "<<i<<" y "<<i+1<<endl;
				//gsl_vector_memcpy(modelFound_aux,modelAux);
				//gsl_vector_free(modelAux);

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
* firstSampleModels function: This function finds the first sample over the threshold of each template of the library
*                             (and the index of that first sample)
*
* Parameters:
* - templates: Matrix which contains all the templates of the library
* - threshold: Threshold
* - firstSamples: Vector which contains the first samples over the threshold of all the templates of the library
* - index_firstSamples: Vector which contains the index of the first samples
****************************************/
int firstSampleModels (ReconstructInitSIRENA *reconstruct_init, double threshold, gsl_vector **firstSamples, gsl_vector **index_firstSamples)
{
  	int numtemplates = reconstruct_init->library_collection->ntemplates;
  	gsl_vector *atemplate = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);

	for (int i=0;i<numtemplates;i++)
	{
		gsl_vector_memcpy(atemplate,reconstruct_init->library_collection->pulse_templates_filder[i].ptemplate);
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

	gsl_vector_free(atemplate);

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
int find_model1stSample(double firstSample, gsl_vector *firstSamples, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound)
{
  	string message = "";

  	long nummodels = reconstruct_init->library_collection->ntemplates;

	if (firstSample < gsl_vector_get(firstSamples,0))
	{
		gsl_vector_memcpy(*modelFound,reconstruct_init->library_collection->pulse_templates_filder[0].ptemplate);
	}
	else if (firstSample > gsl_vector_get(firstSamples,nummodels-1))
	{
		gsl_vector_memcpy(*modelFound,reconstruct_init->library_collection->pulse_templates_filder[nummodels-1].ptemplate);
	}
	else
	{
		for (int i=0;i<nummodels;i++)
		{
			if (firstSample == gsl_vector_get(firstSamples,i))
			{
				gsl_vector_memcpy(*modelFound,reconstruct_init->library_collection->pulse_templates_filder[i].ptemplate);

				break;
			}
			else if ((firstSample > gsl_vector_get(firstSamples,i)) && (firstSample < gsl_vector_get(firstSamples,i+1)))
			{
				// Interpolate between the two corresponding rows in "models"
				gsl_vector *modelAux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);
				gsl_vector_set_zero(modelAux);

				if (interpolate_model(&modelAux,firstSample,reconstruct_init->library_collection->pulse_templates_filder[i].ptemplate,gsl_vector_get(firstSamples,i),
						reconstruct_init->library_collection->pulse_templates_filder[i+1].ptemplate,gsl_vector_get(firstSamples,i+1)))
				{
				    message = "Cannot run interpolate_model with two rows in models";
				    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}

				//cout<<"Modelo interpolado entre "<<i<<" e "<<i+1<<endl;
				gsl_vector_memcpy(*modelFound,modelAux);
				gsl_vector_free(modelAux);

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
int interpolate_model(gsl_vector **modelFound, double ph_model, gsl_vector *modelIn1, double ph_modelIn1, gsl_vector *modelIn2, double ph_modelIn2)
{
	// Declare variables
	double factor1, factor2;
	gsl_vector_set_zero(*modelFound);
	gsl_vector *modelIn1Aux = gsl_vector_alloc(modelIn1->size);
	gsl_vector *modelIn2Aux = gsl_vector_alloc(modelIn2->size);
	gsl_vector_memcpy(modelIn1Aux,modelIn1);
	gsl_vector_memcpy(modelIn2Aux,modelIn2);

	// Method 1: The simplest method
	/*gsl_vector_add(*modelFound,modelIn1);
	gsl_vector_add(*modelFound,modelIn2);
	gsl_vector_scale(*modelFound,0.5);*/

	// Method 2: A bit more intelligent averaging
	factor1 = (ph_modelIn2-ph_model)/(ph_modelIn2-ph_modelIn1);
	factor2 = (ph_model-ph_modelIn1)/(ph_modelIn2-ph_modelIn1);
	gsl_vector_scale(modelIn1Aux,factor1);
	gsl_vector_scale(modelIn2Aux,factor2);
	gsl_vector_add(*modelFound,modelIn1Aux);
	gsl_vector_add(*modelFound,modelIn2Aux);

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
int findTstartCAL (int maxPulsesPerRecord, gsl_vector *der, double adaptativethreshold, int nSamplesUp, int allPulsesMode, double sampling, int *numberPulses, int *thereIsPulse, gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl)
{
	// Declare variables
	string message="";
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
	bool maxFound = false;

	// Allocate GSL vectors
	gsl_vector *tstartDER = gsl_vector_alloc(maxPulsesPerRecord);		// tstarts referred to the derivative
	gsl_vector *tendDER = gsl_vector_alloc(maxPulsesPerRecord);		// tends referred to the derivative
	gsl_vector_set_zero(tendDER);
	*tstartgsl = gsl_vector_alloc(maxPulsesPerRecord);				// tstarts referred to the not filtered event
	*flagTruncated = gsl_vector_alloc(maxPulsesPerRecord);
	gsl_vector_set_zero(*flagTruncated);
	*maxDERgsl = gsl_vector_alloc(maxPulsesPerRecord);			// Maximum of the first derivative
	gsl_vector_set_all(*maxDERgsl,-1E3);

	// Obtain tstart of each pulse in the derivative
	while (i < szRw-1)
	{
		if ((gsl_vector_get(der,i) > adaptativethreshold) && (prevPulse == false))
		{
			if( i>1 && ((gsl_vector_get(der,i)-gsl_vector_get(der,i-1)) <
			           0.5*(gsl_vector_get(der,i-1)-gsl_vector_get(der,i-2)))  )
			{
				maxFound = true;
			}
		      
			if (cntUp == 0)
			{
				possibleTstart=i;

				if(maxFound==false)
				{
				    gsl_vector_set(*maxDERgsl,*numberPulses,gsl_vector_get(der,i));
				}

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
					if (maxFound==false)
					{
					    gsl_vector_set(*maxDERgsl,*numberPulses,gsl_vector_get(der,i));
					}
					gsl_vector_set(*flagTruncated,*numberPulses,1);
					*numberPulses = *numberPulses +1;
					if (allPulsesMode == 0) break;
					prevPulse = true;
				}
				else
				{
					gsl_vector_set(tstartDER,*numberPulses,possibleTstart);
					if (maxFound==false)
					{
					    gsl_vector_set(*maxDERgsl,*numberPulses,gsl_vector_get(der,i));
					}
					*numberPulses = *numberPulses +1;
					if (allPulsesMode == 0) break;
					prevPulse = true;
				}
				cntDown = 0;
			}
		  
			i++;
			cntUp = cntUp+1;
		}
		else if ((gsl_vector_get(der,i) > adaptativethreshold))
		{
			if (gsl_vector_get(der,i) > gsl_vector_get(*maxDERgsl,*numberPulses-1) &&
					(gsl_vector_get(der,i) > gsl_vector_get(der,i-1)) &&
					(maxFound == false))
			  
			{
				if(((gsl_vector_get(der,i)-gsl_vector_get(der,i-1)) <
				      0.5*(gsl_vector_get(der,i-1)-gsl_vector_get(der,i-2))) && i>1)
				{
					maxFound = true;
				}
				else
				{
					gsl_vector_set(*maxDERgsl,*numberPulses-1,gsl_vector_get(der,i));
				  }
			}
			else if (gsl_vector_get(der,i) <= gsl_vector_get(der,i-1))
			{
				maxFound = true;
			}
				
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
						maxFound = false;
					}
				}
				else if (cntDown == nSamplesUp)
				{
					gsl_vector_set(tendDER,*numberPulses-1,possibleTend);
					prevPulse = false;
					maxFound = false;
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
	/*for (int i=0;i<*numberPulses;i++)
	{
		temp = gsl_vector_subvector(der,gsl_vector_get(tstartDER,i),gsl_vector_get(tendDER,i)-gsl_vector_get(tstartDER,i));
		gsl_vector_set(*maxDERgsl,i,gsl_vector_max(&temp.vector));
	}*/
	gsl_vector_memcpy(*tstartgsl,tstartDER);

	gsl_vector_free(tstartDER);
	gsl_vector_free(tendDER);

	if (*numberPulses > 0) *thereIsPulse = 1;

	/*for (int i=0;i<*numberPulses;i++)
	{
		cout<<"Pulso "<<i<<" en "<<gsl_vector_get(*tstartgsl,i)<<endl;
		cout<<"Pulso "<<i<<" con maxDER de "<<gsl_vector_get(*maxDERgsl,i)<<" en "<<gsl_vector_get(*index_maxDERgsl,i)<<endl;
	}*/

	return (EPOK);
}

int findTstartPROD (int maxPulsesPerRecord, gsl_vector *adjustedDerivative, double adaptativethreshold, int nSamplesUp, ReconstructInitSIRENA *reconstruct_init, int *numberPulses, gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl)
{
	// It only identifies the truncated pulses at the beginning of a record

	string message="";
	bool foundPulse = false;
	int sizeRecord = adjustedDerivative->size;		 // Size of segment of process
	*numberPulses = 0;
	*maxDERgsl = gsl_vector_alloc(maxPulsesPerRecord);			// Maximum of the first derivative
	gsl_vector_set_all(*maxDERgsl,-1E3);
	int cntUp = 0;				 // To taking into account how many consecutive samples are over the threshold
	int possibleTstart;			 // To store the first of the nSamplesUp consecutive samples over the threshold
	bool maxFound = false;
	int i = 0;
	gsl_vector *model = gsl_vector_alloc(reconstruct_init->pulse_length);
	bool validPulse = true;

	do
	{
		cntUp = 0;
		maxFound = false;
		foundPulse = false;

		while (i < sizeRecord-1)
		{
			if (gsl_vector_get(adjustedDerivative,i) > adaptativethreshold)
			{
				if ((i>1) && ((gsl_vector_get(adjustedDerivative,i)-gsl_vector_get(adjustedDerivative,i-1)) <
						0.2*(gsl_vector_get(adjustedDerivative,i-1)-gsl_vector_get(adjustedDerivative,i-2))))
				{
					maxFound = true;
					break;
				}

				if (maxFound == false)	gsl_vector_set(*maxDERgsl,*numberPulses,gsl_vector_get(adjustedDerivative,i));

				if (cntUp == 0)			possibleTstart = i;

				if ((foundPulse == false) && ((nSamplesUp == 1) || ((nSamplesUp != 1) && (i == sizeRecord-2)) || (cntUp == nSamplesUp-1)))
				{
				        // Avoid finding two pulses with the same tstart
					for (int j=0;j<*numberPulses;j++)
					{
					  if (possibleTstart == gsl_vector_get(*tstartgsl,j))
					  {
					    validPulse = false; 
					    break;
					  }
					}
					if (validPulse)
					{
					    gsl_vector_set(*tstartgsl,*numberPulses,possibleTstart);
					    foundPulse = true;
					    //cout<<"Pulso "<<*numberPulses<<" en "<<gsl_vector_get(*tstartgsl,*numberPulses)<<" con maxDER en "<<gsl_vector_get(*maxDERgsl,*numberPulses)<<endl;

					    if (possibleTstart == 0)	gsl_vector_set(*flagTruncated,*numberPulses,1);

					    *numberPulses = *numberPulses +1;
					}
				}

				if (validPulse) cntUp++;
			}

			i++;
		}

		// Subtract the model from the adjusted derivative
		if (foundPulse == true)
		{
			if (find_model_maxDERs(gsl_vector_get(*maxDERgsl,*numberPulses-1), reconstruct_init, &model))
			{
				message = "Cannot run find_model routine for pulse i=" + boost::lexical_cast<std::string>(i) + " when newPulses = 1";
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}

			for (int j=gsl_vector_get(*tstartgsl,*numberPulses-1);j<min(gsl_vector_get(*tstartgsl,*numberPulses-1)+reconstruct_init->pulse_length,(double) sizeRecord);j++)
			{
				gsl_vector_set(adjustedDerivative,j,gsl_vector_get(adjustedDerivative,j)-gsl_vector_get(model,j-gsl_vector_get(*tstartgsl,*numberPulses-1)));
			}
			/*for(int j=0;j<sizeRecord;j++){
			  cout<<j<<" "<<gsl_vector_get(adjustedDerivative,j)<<endl;
			}*/
		}
		i = i - nSamplesUp;

	} while (foundPulse == true);

	return(EPOK);
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
int findPulsesCAL (
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

	int samplesup,
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
	double thresholdmediankappa;	// Threshold to look for pulses in the first derivative

	gsl_vector_set_zero(*quality);
	gsl_vector_set_zero(*pulseheight);								// Pulse height of the single pulses

	// First step to look for single pulses
	if (medianKappaClipping (vectorinDER, kappamkc, stopcriteriamkc, nsgms, (int)(pi*samplingRate*taufall*scalefactor), &thresholdmediankappa))
	{
	    message = "Cannot run medianKappaClipping looking for single pulses";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	*threshold = thresholdmediankappa;
	//cout<<"threshold: "<<*threshold<<endl;

	if (findTstartCAL (reconstruct_init->maxPulsesPerRecord, vectorinDER, thresholdmediankappa, samplesup, 1, samplingRate, nPulses, &pulseFound, tstart, quality, maxDERgsl))
	{
		message = "Cannot run findTstartCAL";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	if (*nPulses != 0)
	{
		gsl_vector *Lbgsl = gsl_vector_alloc(reconstruct_init->maxPulsesPerRecord);		// If there is no free-pulses segments longer than Lb=>
		gsl_vector_set_all(Lbgsl,lb);                                   				// segments shorter than Lb will be useed and its length (< Lb)
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


int findPulsesPROD (
	gsl_vector *vectorinDER,
	gsl_vector **tstart,
	gsl_vector **quality,
	gsl_vector **maxDERgsl,

	int *nPulses,
	double *threshold,

	double taufall,
	double scalefactor,
	double samplingRate,

	int samplesup,
	double nsgms,

	ReconstructInitSIRENA *reconstruct_init,

	double stopcriteriamkc,
	double kappamkc)
{
	string message = "";

	const double pi = 4.0 * atan(1.0);

	// Declare variables
	double thresholdmediankappa;	// Threshold to look for pulses in the first derivative

	gsl_vector_set_zero(*quality);

	// First step to look for single pulses
	if (medianKappaClipping (vectorinDER, kappamkc, stopcriteriamkc, nsgms, (int)(pi*samplingRate*taufall*scalefactor), &thresholdmediankappa))
	{
	    message = "Cannot run medianKappaClipping looking for single pulses";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	*threshold = thresholdmediankappa;
	//cout<<"threshold:"<<*threshold<<endl;
	
	if (findTstartPROD (reconstruct_init->maxPulsesPerRecord, vectorinDER, thresholdmediankappa, samplesup, reconstruct_init, nPulses, tstart, quality, maxDERgsl))
	{
		message = "Cannot run findTstartPROD";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	

	return(EPOK);
}
/*xxxx end of SECTION 14 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/




/***** SECTION X ************************************************************
* derivative function: This function applies the derivative method (x_i-x_(i-1)) to the input vector
*
* The derivative method provides more sensitivity to handle with pilot pulses.
* Moreover, little variations of the baseline will not affect.
*
* Arguments:
* 	- invector: Input/Ouput vector (derivative)
* 	- szVct: Size of invector
******************************************************************************/
int derivative (gsl_vector **invector,int szVct)
{
	for (int i=0; i<szVct-1; i++)
	{
		gsl_vector_set(*invector,i,gsl_vector_get(*invector,i+1)-gsl_vector_get(*invector,i));
	}
	gsl_vector_set(*invector,szVct-1,gsl_vector_get(*invector,szVct-2));
	//gsl_vector_set(*invector,szVct-1,-999);

	return (EPOK);
}
/*xxxx end of SECTION X xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
