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
*                      PULSE PROCESS
*
*  File:      pulseProcess.h
*  Developer: Beatriz Cobo Martin
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene Gonzalez Perez
*             Jose Ramon Rodon Ortiz
*
*  Revision History:
*
*  08/10/08	First version
*  19/01/09	Changed the parameters "beforeIndex" and "afterIndex" from double to integer in the findPulse function
*  26/10/09 Old findMean and findPulses are called findMean_OLD and findPulses_OLD
*           New findMean and findPulses added
*  23/07/10	New functions added: "lpf_boxcar", "derMTH", "findTstart" & "findDerPoslength"
*  29/03/11 Updated .h
*  04/04/11 Input parameters findTstart
*  26/05/11 findTstart renamed as findTstart_OLD
*           New functions: findThreshold
*                          kappaClipping
*                          medianKappaClipping
*                          findTstart
*                          findMeanSigma
*           New GSL libraries used
*           inoutUtils.h used (just in case it was neccessary to write comments in a log file)
*  04/08/11 New functions: derMTHSimple
*                          findTstartPrimary
*  24/10/11 Added a new input/output parameter in findTstart
*  14/12/12 New input/output parameter maxDERgsl in findTstartPrimary
*  14/03/13 findPulses renamed as findPulsesOLD
*           New functions:
*              findPulses
*              getB
*              gsl_vector_Sumsubvector
*              findSePrPulses
*              getPulseHeight
*              find_model
*              RS_filter
*              interpolate_model
*  15/02/13 Deleted in/out parameter tend in findPulses
*  26/02/13 getPulseHeight: 'safetymargintstart' deleted as input/output parameter
*           findTstart: New input/output parameter 'tstartgslOUTNEW'
*           findSePrPulses: New input/output parameter 'tstartNOsmt'
*  02/05/13 findSePrPulses: Improved the tstart by correcting the calculus with the pulse broadening ('ATstart1')
*           find_model: Find the pulse broadening by using the pulse models library ('ATstart2' and 'A')
*           interpolate_model: Interpolate (in general) between two pulse broadenings ('AIn1', 'AIn2' and 'Afound')
*  08/05/13 findSePrPulses: New input/output parameter 'tstartMOD' (which is the precisely calculated tstart)
*  02/10/13 findPulses, findSePrPulses, find_model and interpolate_model changed in order to not use PRETRIGS to
*           calculate precisely the tstart
*  08/10/13 Not used functions are deleted
*  10/10/14 Deleted the input parameter 'nPulsesNew' in findSePrPulses
************************************************************************/

#ifndef PULSEPROCESS_H_
#define PULSEPROCESS_H_

// Utils module

	#include <genutils.h>
	#include <inoututils.h>

	#include <integraSIRENA.h>

	int lpf_boxcar (gsl_vector **invector, int szVct, double tau_fall, int sampleRate);
	int derMTHSimple (gsl_vector **invector,gsl_vector **sign, int szVct);

	int findMeanSigma (gsl_vector *invector, double *mean, double *sigma);
	int medianKappaClipping (gsl_vector *invector, double kappa, double stopCriteria, double nSigmas, int boxLPF,double *threshold);

	int getB(gsl_vector *vectorin, gsl_vector *tstart, int nPulses, gsl_vector **lb, int sizepulse, gsl_vector **B);
	int getPulseHeight(gsl_vector *vectorin, double tstart, double tstartnext, int lastPulse, double lrs, double lb, double B, int sizepulse, double *pulseheight);
	int RS_filter (gsl_vector *vector, double lrs, double lb, double B, double *pulseheight);

	int find_model_energies(double ph, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound);
	int find_model_maxDERs(double ph, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound);
	int find_model1stSample(double firstSample, gsl_vector *firstSamples, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound);
	int firstSampleModels (ReconstructInitSIRENA *reconstruct_init, double threshold, gsl_vector **firstSamples, gsl_vector **index_firstSamples);
	int interpolate_model(gsl_vector **modelFound, double ph_model, gsl_vector *modelIn1, double ph_modelIn1, gsl_vector *modelIn2, double ph_modelIn2);

	int findTstart (gsl_vector *der, double adaptativethreshold, int nSamplesUp,
			int allPulsesMode, double sampling, int *numberPulses, int *thereIsPulse,
			gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl, gsl_vector **index_maxDERgsl);

	int findPulses
	(
		gsl_vector *vectorin,
		gsl_vector *vectorinDER,
		gsl_vector **tstart,
		gsl_vector **quality,
		gsl_vector **energy,
		gsl_vector **maxDERgsl,

		int *nPulses,
		double *threshold,

		int opmode,

		double taufall,
		double scalefactor,
		int sizepulsebins,
		double samplingRate,

		int samplesup,
		double nsgms,

		double lb,
		double lrs,

		ReconstructInitSIRENA *reconstruct_init,

		double stopcriteriamkc,
		double kappamkc,
		double levelprvpulse);

	int findSePulses
	(
		gsl_vector *vectorin,
		gsl_vector *vectorinDER,
		gsl_vector **vectorinDERComposed,

		double thresholdmediankappaSingle,

		gsl_vector **tstart,
		gsl_vector **quality,
		//gsl_vector **energy,
		gsl_vector **maxDER,
		gsl_vector **index_maxDER,

		gsl_vector **newPulses,

		int *nPulses,

		/*gsl_vector *startsaturated,
		gsl_vector *endsaturated,
		int nSaturated,*/

		double taufall,
		double scalefactor,
		int sizepulse,
		double samplingRate,

		int samplesup,
		double nsgms,

		gsl_vector *model,
		ReconstructInitSIRENA *reconstruct_init,

		double stopCriteriamkc,
		double kappamkc,
		double levelprvpulse);

	int derivative (gsl_vector **invector,int szVct);

	using namespace std;

#endif /*PULSEPROCESS_H_*/
