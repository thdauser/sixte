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
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01, 
   ESP2013-48637-C2-1-P and ESP2014-53672-C3-1-P.

/***********************************************************************
*                      PULSEPROCESS
*
*  File:       pulseprocess.h
*  Developers: Beatriz Cobo
* 	       cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef PULSEPROCESS_H_
#define PULSEPROCESS_H_

// Utils module

	#include <genutils.h>
	#include <inoututils.h>

	#include <integraSIRENA.h>

	int lpf_boxcar (gsl_vector **invector, int szVct, double scaleFactor, int sampleRate);
	int differentiate (gsl_vector **invector,int szVct);

	int findMeanSigma (gsl_vector *invector, double *mean, double *sigma);
	int medianKappaClipping (gsl_vector *invector, double kappa, double stopCriteria, double nSigmas, int boxLPF,double *threshold);

	int getB(gsl_vector *vectorin, gsl_vector *tstart, int nPulses, gsl_vector **lb, int sizepulse, gsl_vector **B);
	int getPulseHeight(gsl_vector *vectorin, double tstart, double tstartnext, int lastPulse, double lrs, double lb, double B, int sizepulse, double *pulseheight);
	int RS_filter (gsl_vector *vector, double lrs, double lb, double B, double *pulseheight);

	int find_model_energies(double energy, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound);
	int find_model_maxDERs(double maxDER, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound);
	int find_model_samp1DERs(double samp1DER, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound);
	int interpolate_model(gsl_vector **modelFound, double ph_model, gsl_vector *modelIn1, double ph_modelIn1, gsl_vector *modelIn2, double ph_modelIn2);

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
		double kappamkc);

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
		gsl_vector **maxDERgsl);
	
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

		double *threshold,

		int tstartProvided = 0);
	
	int FindSecondaries
	(
		int maxPulsesPerRecord,

		gsl_vector *adjustedDerivative,
		double threshold,
                double samprate,

		ReconstructInitSIRENA *reconstruct_init,

		int tstartFirstEvent,

		int *numberPulses,
		gsl_vector **tstartgsl,
		gsl_vector **flagTruncated,
		gsl_vector **maxDERgsl,
		gsl_vector **sammp1DERgsl,
                gsl_vector **lagsgsl,
        
                int firstRecord);
        
        int find_model_samp1DERsNoReSCLD(double samp1DER, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound, int *indexMin, int *indexMax);
        int smoothDerivative (gsl_vector **invector, int N);
        
        int FindSecondariesA1
        (       
                int maxPulsesPerRecord,

                gsl_vector *der,
                double adaptativethreshold,
                int nSamplesUp,

                ReconstructInitSIRENA *reconstruct_init,

                int *numberPulses,
                
                gsl_vector **tstartgsl,
                gsl_vector **flagTruncated,
                gsl_vector **maxDERgsl,
                gsl_vector **samp1DERgsl,
                
                int firstRecord);

	using namespace std;

#endif /*PULSEPROCESS_H_*/
