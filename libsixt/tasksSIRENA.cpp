
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

   Copyright 2014:  TASKSSIRENA has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01, 
   ESP2013-48637-C2-1-P and ESP2014-53672-C3-1-P.

***********************************************************************
*                      TASKSSIRENA
*
*  File:       tasksSIRENA.cpp
*  Developers: Beatriz Cobo
* 	       cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************

******************************************************************************
DESCRIPTION:

The purpose of this package is the detection and the reconstruction.

MAP OF SECTIONS IN THIS FILE:

 - A. runDetect 
 - A1. createLibrary
 - A2. createDetectFile
 - A3. filderLibrary
 - A4. loadRecord
 - A5. procRecord
 - A6. writePulses
 - A7. writeTestInfo
 - A8. calculateTemplate
 - A9. createHisto
 - A10. align
 - A11. shiftm
 - A12. shift_m
 - A13. weightMatrix
 - A14. writeLibrary
 - A15. addFirstRow
 - A16. readAddSortParams
 - A17. calculateIntParams
 - A18. matrix2vector
 - A19. vector2matrix
 - A20. convertI2R
 - A21. filterByWavelets
 - B. runEnergy
 - B1. calculus_optimalFilter
 - B2. interpolatePOS
 - B3. find_matchedfilter
 - B4. find_matchedfilterDAB
 - B5. find_optimalfilter
 - B6. find_optimalfilterDAB
 - B7. find_prclwn
 - B8. find_prclofwm
 - B9. find_Esboundary
 - B10. pulseGrading
 - B11. calculateEnergy
 - B12. writeFilterHDU

*******************************************************************************/

#include "tasksSIRENA.h"
#include "log.h"
#include "scheduler.h"

/***** SECTION A ************************************************************
* runDetect: This function is responsible for the detection in SIRENA, record by record.
*            It is used both for library creation ('opmode'=0) and energy reconstruction ('opmode'=1) runnings.
*            In CALIBRATION mode the purpose is the library creation
*
* Conditions:
* 
* - If first record and PRODUCTION ('opmode=1') => 'filderLibrary'
* - If last record and CALIBRATION ('opmode=0') => 'calculateTemplate' and 'writeLibrary'
* - If 'reconstruct_init->intermediate'=1 => 'writeTestInfo' and 'writePulses'
* - If CALIBRATION ('opmode=0') => Find pulses by using 'findPulsesCAL'
* - If PRODUCTION ('opmode=1') => Find pulses by 'InitialTriggering' and 'FindSecondaries' or 'FindSecondariesSTC'
*
* - Declare variables
* - Create library if it is necessary, CALIBRATION ('opmode=0') and 'lastRecord'=1 ('createLibrary')
* - Create intermediate output FITS file if required ('createDetectFile')
* - Filter and differentiate the 'models' of the library (only for the first record in PRODUCTION 'opmode=1') ('filderLibrary')
* - Store the input record in 'invector' ('loadRecord')
* - Convert I into R if 'EnergyMethod' = I2R or I2RALL or I2RNOL or I2RFITTED ('convertI2R')
* - Process each record ('proceRecord')
* 	- Low-pass filter and differentiate
* 	- Find pulses
* 	- Load the found pulses data in the input/output 'foundPulses' structure
* 	- Write test info in intermediate output FITS file if 'reconstruct_init->intermediate'=1 ('writeTestInfo')
* 	- Write pulses info in intermediate output FITS file if 'reconstruct_init->intermediate'=1 ('writePulses')
* - From this point forward, I2R, I2RALL, I2RNOL and I2RFITTED are completely equivalent to OPTFILT
* - If last record in CALIBRATION ('opmode'=0):
* 	- 'calculateTemplate' (and 'weightMatrix')
* 	- 'writeLibrary'
* - If last record and PCA
*       - In order to not have restrictions when providing (*reconstruct_init)->energyPCAx
* 	- Covariance data
* 	- Eigenvalues and eigenvectors
* 	- RSxN (S=2)
* 	- AE straight line: Pto0(x,y) and Pto10(x,y)
*       - Calculus of the rotation angle
*       - Rotation
*       - Histograms of the two clusters (two energies)
*       - Conversion factor from arbitrary unit to eV
*       - Energy calculation
* - Close intermediate output FITS file if it is necessary
* - Free allocated GSL vectors
*
* Parameters:
* - record: Member of TesRecord' structure that contains the input record
* - lastRecord: Integer to verify whether record is the last one (=1) to be read (and thus if library file will be created)
* - pulsesAll: Member of 'PulsesCollection' structure to successively store all the pulses used to create the library. Re-populated after each processed record
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
* - pulsesInRecord: Member of 'PulsesCollection' structure to store all the pulses found in the input record
******************************************************************************/
/*
void runDetect(TesRecord* record, int nRecord, int lastRecord, 
               PulsesCollection *pulsesAll, 
               ReconstructInitSIRENA** reconstruct_init, 
               PulsesCollection** pulsesInRecord)
*/
void runDetect(TesRecord* record, int lastRecord, PulsesCollection *pulsesAll, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord)
{       
	int inputPulseLength = (*reconstruct_init)->pulse_length;
	if ((*reconstruct_init)->opmode == 0)	(*reconstruct_init)->pulse_length = (*reconstruct_init)->largeFilter;
	  
	string message="";
	char valERROR[256];
	int status=EPOK;  
			
	// Declare variables
	fitsfile *inLibObject = NULL;	// Object which contains information of the library FITS file
	bool appendToLibrary = false;	// Calibration library FITS file new (appendToLibrary=false) or not (appendToLibrary=true)

	fitsfile *dtcObject = NULL;	// Object which contains information of the intermediate FITS file ('dtc' comes from 'detectFile')
	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';

	int eventsz = record->trigger_size;
	double tstartRecord;
	// It is not necessary to check the allocation because 'record->trigger_size'='eventsz' has been checked previously
	gsl_vector *invector = gsl_vector_alloc(eventsz);	// Record
        
        gsl_vector *averageRecord = gsl_vector_alloc(eventsz);
	
	if (((*reconstruct_init)->tstartPulse1 != 0) && (((*reconstruct_init)->tstartPulse1 > record->trigger_size) || ((*reconstruct_init)->tstartPulse2 > record->trigger_size) || ((*reconstruct_init)->tstartPulse3 > record->trigger_size)))
	{
		message = "tstartPulsx can not be greater than the record size";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	//assert(((*reconstruct_init)->tstartPulse1 <= record->trigger_size) && ((*reconstruct_init)->tstartPulse2 <= record->trigger_size) && ((*reconstruct_init)->tstartPulse3 <= record->trigger_size));
	
	// Create library if it is necessary
	if (((*reconstruct_init)->opmode == 0) && (lastRecord == 1))
	{
		if (createLibrary(*reconstruct_init, &appendToLibrary, &inLibObject, inputPulseLength))
		{
			message = "Cannot run routine createLibrary to create pulses library";
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	// Create intermediate output FITS file if required (reconstruct_init->intermediate'=1)
	if ((*reconstruct_init)->intermediate == 1)
	{
		if (createDetectFile(*reconstruct_init, 1/record->delta_t, &dtcObject, inputPulseLength))
		{
			message = "Cannot create file " +  string((*reconstruct_init)->detectFile);
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	// (Filter and) differentiate the 'models' of the library (only for the first record)
	if ((*reconstruct_init)->opmode == 1)
	{
		if (filderLibrary(reconstruct_init,1/record->delta_t))
		{
			message = "Cannot run routine filderLibrary to filter & differentiate library if the 1st record";
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	// Store the input record in 'invector'
	if (loadRecord(record, &tstartRecord, &invector))
	{
		message = "Cannot run routine loadRecord";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	eventsz = invector->size;	// Just in case the last record has been filled in with 0's => Re-allocate invector
	
	//gsl_vector_add(averageRecord,invector);

	// Convert I into R if 'EnergyMethod' = I2R or I2RALL or I2RNOL or I2RFITTED
	// It is not necessary to check the allocation because 'invector' size must be > 0
	if ((strcmp((*reconstruct_init)->EnergyMethod,"I2R") == 0) || (strcmp((*reconstruct_init)->EnergyMethod,"I2RALL") == 0) || (strcmp((*reconstruct_init)->EnergyMethod,"I2RNOL") == 0)
		|| (strcmp((*reconstruct_init)->EnergyMethod,"I2RFITTED") == 0))
	{
                //if (convertI2R(reconstruct_init, &record, &invector))
                if (convertI2R((*reconstruct_init)->EnergyMethod,(*reconstruct_init)->i2rdata->R0,(*reconstruct_init)->i2rdata->I0_START,(*reconstruct_init)->i2rdata->IMIN,(*reconstruct_init)->i2rdata->IMAX,(*reconstruct_init)->i2rdata->RPARA,(*reconstruct_init)->i2rdata->TTR,(*reconstruct_init)->i2rdata->LFILTER, 1/record->delta_t, &invector))
		{
			message = "Cannot run routine convertI2R";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		for (int i=0;i<invector->size;i++)		     // Because in 'runEnergy' the record (TesRecord) is used => The I2R, I2RALL, I2RNOL or I2RFITTED transformed record has to be used
                {
                        record->adc_double[i] = gsl_vector_get(invector,i);
                }
	}
	
	// Process each record
	if (procRecord(reconstruct_init, tstartRecord, 1/record->delta_t, dtcObject, invector, *pulsesInRecord))
	{
		message = "Cannot run routine procRecord for record processing";
		EP_EXIT_ERROR(message,EPFAIL);
	}
    
	if ((strcmp((*reconstruct_init)->EnergyMethod,"I2R") == 0) || (strcmp((*reconstruct_init)->EnergyMethod,"I2RALL") == 0) || (strcmp((*reconstruct_init)->EnergyMethod,"I2RNOL") == 0)
		|| (strcmp((*reconstruct_init)->EnergyMethod,"I2RFITTED") == 0))
	{
		strcpy((*reconstruct_init)->EnergyMethod,"OPTFILT"); // From this point forward, I2R, I2RALL, I2RNOL and I2RFITTED are completely equivalent to OPTFILT
	}

	if (((*reconstruct_init)->intermediate == 1) && (lastRecord == 1))
	{		
		// Write output keywords (their values have been previously checked)
		char keyname[10];
		
		char extname[10];
		strncpy(extname,"PULSES",9);
		extname[9]='\0';
		if (fits_movnam_hdu(dtcObject, ANY_HDU,extname, 0, &status))
		{
			message = "Cannot move to HDU " + string(extname) +" in " + string(dtcName);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		long totalpulses;
		if (fits_get_num_rows(dtcObject,&totalpulses, &status))
		{
			message = "Cannot get number of rows in " + string(dtcName);
			EP_EXIT_ERROR(message,EPFAIL);
		}

		int ttpls1 = (int) totalpulses;
		strcpy(keyname,"EVENTCNT");
		if (ttpls1 < 0)
		{
			message = "Legal values for EVENTCNT (PULSES) are integer numbers greater than or equal to 0";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		if(fits_write_key(dtcObject,TINT,keyname,&ttpls1,NULL,&status))
		{
			message = "Cannot write keyword " + string(keyname) +" in " + string(dtcName);
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	if ((lastRecord == 1) && (*reconstruct_init)->opmode == 0 && (pulsesAll->ndetpulses + (*pulsesInRecord)->ndetpulses >0))	// CALIBRATION mode => Calculate the pulse template by averaging some found pulses
	{	
		// Calculate an average record 
                //gsl_vector_scale(averageRecord,1/);

            
                // It is not necessary to check the allocation because 'PulseLength' (input parameter) has been checked previously                    
		gsl_vector *pulsetemplateMaxLengthFixedFilter = gsl_vector_alloc((*reconstruct_init)->largeFilter);
		gsl_vector *pulsetemplate = gsl_vector_alloc((*reconstruct_init)->pulse_length);
		double pulseheighttemplate = 0;
		gsl_matrix *weight = gsl_matrix_alloc(inputPulseLength,inputPulseLength);
		gsl_matrix *covariance = gsl_matrix_alloc(inputPulseLength,inputPulseLength);
		gsl_matrix_set_zero(weight);
		gsl_matrix_set_zero(covariance);
    
                if (calculateTemplate (*reconstruct_init, pulsesAll, *pulsesInRecord, 1/record->delta_t, &pulsetemplate, &pulseheighttemplate, &covariance, &weight, inputPulseLength, &pulsetemplateMaxLengthFixedFilter))
		{
			message = "Cannot run routine calculateTemplate in CALIBRATION mode";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		(*reconstruct_init)->pulse_length = inputPulseLength;
                
		if (writeLibrary(reconstruct_init, 1/record->delta_t, pulseheighttemplate, pulsetemplate, covariance, weight, appendToLibrary, &inLibObject, pulsetemplateMaxLengthFixedFilter))
		{
			message = "Cannot run routine writeLibrary in CALIBRATION mode";
			EP_EXIT_ERROR(message,EPFAIL);
		}

		gsl_vector_free(pulsetemplate);
		gsl_matrix_free(weight);
		gsl_matrix_free(covariance);

		gsl_vector_free(pulsetemplateMaxLengthFixedFilter);
	}
	
	if ((*reconstruct_init)->opmode == 0)	(*reconstruct_init)->pulse_length = inputPulseLength;
	
	// PCA is used  when pulses are farther than 'PulseLength'!!!!!!!!!!
	if ((lastRecord == 1) && (strcmp((*reconstruct_init)->EnergyMethod,"PCA") == 0) && ((pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses))>0)	// PCA
	{
		(*reconstruct_init)->pulse_length = inputPulseLength;
		
		// In order to not have restrictions when providing (*reconstruct_init)->energyPCAx
		double energyPCA1, energyPCA2;
		if ((*reconstruct_init)->energyPCA2 < (*reconstruct_init)->energyPCA1)
		{
			energyPCA1 = (*reconstruct_init)->energyPCA2;
			energyPCA2 = (*reconstruct_init)->energyPCA1;
		}
		else
		{
			energyPCA1 = (*reconstruct_init)->energyPCA1;
			energyPCA2 = (*reconstruct_init)->energyPCA2;
		}
		
		// Covariance data
		// It is not necessary to check the allocation because 'PulseLength' (input parameter) has been checked previously 
		gsl_vector *pulsetemplate = gsl_vector_alloc((*reconstruct_init)->pulse_length);
		gsl_vector_set_zero(pulsetemplate);
		gsl_matrix *covarianceData = gsl_matrix_alloc((*reconstruct_init)->pulse_length,(*reconstruct_init)->pulse_length);
		gsl_matrix *weightData = gsl_matrix_alloc((*reconstruct_init)->pulse_length,(*reconstruct_init)->pulse_length);
		gsl_vector *nonpileup;
		if ((nonpileup = gsl_vector_alloc((pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses))) == 0)
		{  
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);
			message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		gsl_vector_set_all(nonpileup,1);
		if (weightMatrix(*reconstruct_init, false, pulsesAll, *pulsesInRecord, (pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses), nonpileup , pulsetemplate, &covarianceData, &weightData))
		{
			message = "Cannot run weightMatrix routine";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		gsl_vector_free(pulsetemplate);
		gsl_matrix_free(weightData);
		
		// Eigenvalues and eigenvectors
		gsl_vector *eigenvalues;
		gsl_matrix *eigenvectors;	// Eigenvectors in columns
		if (eigenVV(covarianceData,&eigenvectors,&eigenvalues))
		{
			message = "Cannot run eigenVV routine";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		gsl_matrix_free(covarianceData);
		gsl_vector_free(eigenvalues);
		
		// RSxN (S=2)
		gsl_matrix *RowFeatureVectors;		// Eigenvectors in rows
		if ((RowFeatureVectors = gsl_matrix_alloc(eigenvectors->size2,eigenvectors->size1)) == 0)
		{  
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);
			message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		//gsl_matrix *RowFeatureVectors = gsl_matrix_alloc(eigenvectors->size2,eigenvectors->size1);	// Eigenvectors in rows
		gsl_matrix_transpose_memcpy(RowFeatureVectors,eigenvectors);
		gsl_matrix_free(eigenvectors);
		
		// It is not necessary to check the allocation because 'PulseLength'(input parameter) and '(pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulsesc' have been checked previously 
                // S0^1-M0 S0^2-M0...............S0^(nonpileupPulses)-M0
                // S1^1-M1 S1^2-M1...............S1^(nonpileupPulses)-M1
                // ...................................................
                // S1023^1-M1023 S1023^2-M1023...S1023^(nonpileupPulses)-M1023
		gsl_matrix *RowDataAdjust = gsl_matrix_alloc((*reconstruct_init)->pulse_length,(pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses));	 
																				
		for (int m=0;m<(*reconstruct_init)->pulse_length;m++)
		{
			for (int p=0;p<pulsesAll->ndetpulses;p++)
			{
				if (gsl_vector_get(nonpileup,p) == 1)
				{	
					if (m == pulsesAll->pulses_detected[p].pulse_adc->size)
					{
						sprintf(valERROR,"%d",__LINE__+5);
						string str(valERROR);
						message = "PCA should be used with pulses farther than 'PulseLength' => Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
						EP_EXIT_ERROR(message,EPFAIL);
					}
					gsl_matrix_set(RowDataAdjust,m,p,gsl_vector_get(pulsesAll->pulses_detected[p].pulse_adc,m));
				}
			}
			for (int p=0;p<(*pulsesInRecord)->ndetpulses;p++)
			{
				if (gsl_vector_get(nonpileup,pulsesAll->ndetpulses+p) == 1)
				{
					if (m == (*pulsesInRecord)->pulses_detected[p].pulse_adc->size)
					{
						sprintf(valERROR,"%d",__LINE__+5);
						string str(valERROR);
						message = "PCA should be used with pulses farther than 'PulseLength' => Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
						EP_EXIT_ERROR(message,EPFAIL);
					}
					gsl_matrix_set(RowDataAdjust,m,pulsesAll->ndetpulses+p,gsl_vector_get((*pulsesInRecord)->pulses_detected[p].pulse_adc,m));
				}
			}
		}
		gsl_vector_free(nonpileup);
		
		// It is not necessary to check the allocation because the allocation of 'RowFeatureVectors' and '(pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses' have been checked previously
		gsl_matrix *RSrxN = gsl_matrix_alloc(RowFeatureVectors->size1,(pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses));
		if (RowFeatureVectors->size2 != RowDataAdjust->size1)
		{
			sprintf(valERROR,"%d",__LINE__+5);
			string str(valERROR);
			message = "Wrong dimensions to compute matrix-matrix product in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, RowFeatureVectors, RowDataAdjust,0.0, RSrxN);
		gsl_matrix_free(RowFeatureVectors);
		gsl_matrix_free(RowDataAdjust);
		
		// Store RSxN in a FITS file
		FILE * temporalFile;
		char temporalFileName[255];
		char val[256];
		sprintf(temporalFileName,"rsxn.txt");
		temporalFile = fopen (temporalFileName,"w");
		for (int i = 0; i < RSrxN->size2; i++) // All the pulses in a column
		{
			sprintf(val,"%e %e",gsl_matrix_get(RSrxN,0,i),gsl_matrix_get(RSrxN,1,i));
			strcat(val,"\n");
			fputs(val,temporalFile);
		}
		fclose(temporalFile);
		
		// AE straight line: Pto0(x,y) and Pto10(x,y)
		double a,b;	//y=ax+b
		// It is not necessary to check the allocation because '(pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses' has been checked previously
		gsl_vector *xfit = gsl_vector_alloc((pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses));
		gsl_vector *yfit = gsl_vector_alloc((pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses));
		gsl_matrix_get_row(xfit,RSrxN,0);
		gsl_matrix_get_row(yfit,RSrxN,1);
		if (polyFitLinear (xfit, yfit, &a, &b))
		{
			message = "Cannot run routine polyFitLinear";
			EP_EXIT_ERROR(message,EPFAIL);
		}
                gsl_vector *Pto0 = gsl_vector_alloc(2);
		gsl_vector *Pto10 = gsl_vector_alloc(2);
		// It is not necessary to check the allocation because the allocation of 'RSrxN' has been checked previously
		gsl_vector *xRSrxN = gsl_vector_alloc(RSrxN->size2);
		gsl_matrix_get_row(xRSrxN,RSrxN,0);
		gsl_vector_set(Pto0,0,gsl_vector_get(xRSrxN,gsl_vector_min_index(xRSrxN)));
		gsl_vector_set(Pto0,1,a*gsl_vector_get(Pto0,0)+b);
		gsl_vector_set(Pto10,0,gsl_vector_get(xRSrxN,gsl_vector_max_index(xRSrxN))-gsl_vector_get(Pto0,0));
		gsl_vector_set(Pto10,1,a*gsl_vector_get(xRSrxN,gsl_vector_max_index(xRSrxN))+b-gsl_vector_get(Pto0,1));
		
		// Calculus of the rotation angle
		gsl_vector *u = gsl_vector_alloc(2);
		gsl_vector_set(u,0,1);
		gsl_vector_set(u,1,0);
		gsl_vector *v = gsl_vector_alloc(2);
		gsl_vector_memcpy(v,Pto10);
		gsl_vector_free(Pto10);
		double dotResult;
		gsl_blas_ddot (u, v, &dotResult);
		double alfa;
		alfa = acos(dotResult/(gsl_blas_dnrm2(u)*gsl_blas_dnrm2(v)));
		
		// Rotation
		gsl_matrix *RotationMatrix = gsl_matrix_alloc(2,2);
		gsl_matrix_set(RotationMatrix,0,0,cos(alfa));
		gsl_matrix_set(RotationMatrix,0,1,sin(alfa));
		gsl_matrix_set(RotationMatrix,1,0,-sin(alfa));
		gsl_matrix_set(RotationMatrix,1,1,cos(alfa));
		// It is not necessary to check the allocation because the allocation of 'RSrxN' has been checked previously
		gsl_matrix *pointsTranslated = gsl_matrix_alloc(2,RSrxN->size2);
		for (int i = 0; i < RSrxN->size2; i++)
		{
			gsl_matrix_set(pointsTranslated,0,i,gsl_matrix_get(RSrxN,0,i)-gsl_vector_get(Pto0,0));
			gsl_matrix_set(pointsTranslated,1,i,gsl_matrix_get(RSrxN,1,i)-gsl_vector_get(Pto0,1));
		}
		gsl_vector_free(Pto0);
		// It is not necessary to check the allocation because the allocation of 'RSrxN' has been checked previously
		gsl_matrix *pointsTranslatedRotated = gsl_matrix_alloc(2,RSrxN->size2);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, RotationMatrix, pointsTranslated,0.0, pointsTranslatedRotated);
		/*for (int i = 0; i < RSrxN->size2; i++)
		{
			gsl_matrix_set(pointsTranslatedRotated,0,i,gsl_matrix_get(pointsTranslatedRotated,0,i)+gsl_vector_get(Pto0,0));
		}
		sprintf(temporalFileName,"rsxnTransformed.txt");
		temporalFile = fopen (temporalFileName,"w");
		for (int i = 0; i < RSrxN->size2; i++) // All the pulses in a column
		{
			sprintf(val,"%e %e",gsl_matrix_get(pointsTranslatedRotated,0,i),gsl_matrix_get(pointsTranslatedRotated,1,i));
			strcat(val,"\n");
			fputs(val,temporalFile);
		}
		fclose(temporalFile);*/
		gsl_matrix_free(pointsTranslated);
		gsl_matrix_free(RSrxN);
		// It is not necessary to check the allocation because the allocation of 'pointsTranslatedRotated' has been checked previously
		gsl_vector *pointsTranslatedRotated1aux = gsl_vector_alloc(pointsTranslatedRotated->size2);
		gsl_vector *pointsTranslatedRotated2aux = gsl_vector_alloc(pointsTranslatedRotated->size2);
		int num1 = 0;
		int num2 = 0;
		gsl_vector *xpointsTranslatedRotated = gsl_vector_alloc(pointsTranslatedRotated->size2);
		gsl_matrix_get_row(xpointsTranslatedRotated,pointsTranslatedRotated,0);
		double midpoint = gsl_vector_min(xpointsTranslatedRotated)+(gsl_vector_max(xpointsTranslatedRotated)-gsl_vector_min(xpointsTranslatedRotated))/2;
		for (int i = 0; i < xpointsTranslatedRotated->size; i++)
		{
			if (gsl_vector_get(xpointsTranslatedRotated,i) <= midpoint)
			{
				if (num1 >= pointsTranslatedRotated1aux->size)
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
					EP_EXIT_ERROR(message,EPFAIL);
				}
				gsl_vector_set(pointsTranslatedRotated1aux,num1,gsl_vector_get(xpointsTranslatedRotated,i));
				num1++;
			}
			else
			{
				if (num2 >= pointsTranslatedRotated2aux->size)
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
					EP_EXIT_ERROR(message,EPFAIL);
				}
				gsl_vector_set(pointsTranslatedRotated2aux,num2,gsl_vector_get(xpointsTranslatedRotated,i));
				num2++;
			}
		}
		gsl_vector_free(xpointsTranslatedRotated);
		gsl_vector *pointsTranslatedRotated1;
		gsl_vector *pointsTranslatedRotated2;
		if ((pointsTranslatedRotated1 = gsl_vector_alloc(num1)) == 0)
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);
			message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		if ((pointsTranslatedRotated2 = gsl_vector_alloc(num2)) == 0)
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);
			message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		gsl_vector_view temp;
		if ((num1 < 1) || (num1 > pointsTranslatedRotated1aux->size))
		{
			sprintf(valERROR,"%d",__LINE__+5);
			string str(valERROR);
			message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		temp = gsl_vector_subvector(pointsTranslatedRotated1aux,0,num1);
		if (gsl_vector_memcpy(pointsTranslatedRotated1,&temp.vector) != 0)
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);
			message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		if ((num2 < 1) || (num1 > pointsTranslatedRotated2aux->size))
		{
			sprintf(valERROR,"%d",__LINE__+5);
			string str(valERROR);
			message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		temp = gsl_vector_subvector(pointsTranslatedRotated2aux,0,num2);
		if (gsl_vector_memcpy(pointsTranslatedRotated2,&temp.vector) != 0)
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);
			message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		gsl_vector_free(pointsTranslatedRotated1aux);
		gsl_vector_free(pointsTranslatedRotated2aux);
		
		// Histograms of the two clusters (two energies)
		int nbins = 80;
		assert(nbins > 0);
		// It is not necessary to check the allocation because 'nbins' has been checked previously
		gsl_vector *xhisto1 = gsl_vector_alloc(nbins);
		gsl_vector *yhisto1 = gsl_vector_alloc(nbins);
		gsl_vector *xhisto2 = gsl_vector_alloc(nbins);
		gsl_vector *yhisto2 = gsl_vector_alloc(nbins);
		bool histo1neg = false;
		bool histo2neg = false;
		double histo1constant;
		double histo2constant;
		if (gsl_vector_isnonneg(pointsTranslatedRotated1) != 1)
		{
			histo1neg = true;
			histo1constant = fabs(gsl_vector_min(pointsTranslatedRotated1));
			gsl_vector_add_constant(pointsTranslatedRotated1,histo1constant);
		}
		if (createHisto(pointsTranslatedRotated1, nbins, &xhisto1, &yhisto1))
		{
			message = "Cannot run createHisto routine";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		if (histo1neg == true) 
		{
			gsl_vector_add_constant(xhisto1,-histo1constant);
			gsl_vector_add_constant(pointsTranslatedRotated1,-histo1constant);
		}
		sprintf(temporalFileName,"histo1.txt");
		temporalFile = fopen (temporalFileName,"w");
		for (int i = 0; i < xhisto1->size; i++) 
		{
			sprintf(val,"%e %e",gsl_vector_get(xhisto1,i),gsl_vector_get(yhisto1,i));
			strcat(val,"\n");
			fputs(val,temporalFile);
		}
		fclose(temporalFile);
		if (gsl_vector_isnonneg(pointsTranslatedRotated2) != 1)
		{
			histo2neg = true;
			histo2constant = fabs(gsl_vector_min(pointsTranslatedRotated2));
			gsl_vector_add_constant(pointsTranslatedRotated2,histo2constant);
		}
		if (createHisto(pointsTranslatedRotated2, nbins, &xhisto2, &yhisto2))
		{
			message = "Cannot run createHisto routine";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		if (histo2neg == true) 
		{
			gsl_vector_add_constant(xhisto2,-histo2constant);
			gsl_vector_add_constant(pointsTranslatedRotated2,-histo2constant);
		}
		sprintf(temporalFileName,"histo2.txt");
		temporalFile = fopen (temporalFileName,"w");
		for (int i = 0; i < xhisto2->size; i++) 
		{
			sprintf(val,"%e %e",gsl_vector_get(xhisto2,i),gsl_vector_get(yhisto2,i));
			strcat(val,"\n");
			fputs(val,temporalFile);
		}
		fclose(temporalFile);
		double data1[num1];
		double data2[num2];
		if (num1 > pointsTranslatedRotated1->size)
		{
			sprintf(valERROR,"%d",__LINE__+7);
			string str(valERROR);
			message = "Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		for (int i = 0; i < num1; i++)
		{
			data1[i] = gsl_vector_get(pointsTranslatedRotated1,i);
		}
		if (num2 > pointsTranslatedRotated2->size)
		{
			sprintf(valERROR,"%d",__LINE__+7);
			string str(valERROR);
			message = "Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		for (int i = 0; i < num2; i++)
		{
			data2[i] = gsl_vector_get(pointsTranslatedRotated2,i);
		}
		double sigma1 = gsl_stats_sd(data1, 1, num1);
		double sigma2 = gsl_stats_sd(data2, 1, num2);
		double maxcenter1 = gsl_stats_mean(data1, 1, num1);
		double maxcenter2 = gsl_stats_mean(data2, 1, num2);
		
		// Conversion factor from arbitrary unit to eV
		double convertAU2eV = (energyPCA2 - energyPCA1)/(maxcenter2-maxcenter1);	
		
		gsl_vector_free(pointsTranslatedRotated1);
		gsl_vector_free(pointsTranslatedRotated2);
		gsl_vector_free(xhisto1);
		gsl_vector_free(yhisto1);
		gsl_vector_free(xhisto2);
		gsl_vector_free(yhisto2);
		
		// Energy calculation
		for (int i = 0; i < pulsesAll->ndetpulses; i++)
		{
			pulsesAll->pulses_detected[i].energy = gsl_matrix_get(pointsTranslatedRotated,0,i)*convertAU2eV/1e3 + energyPCA1/1e3;		//keV
		}
		for (int i = 0; i < (*pulsesInRecord)->ndetpulses; i++)
		{
			(*pulsesInRecord)->pulses_detected[i].energy = gsl_matrix_get(pointsTranslatedRotated,0,pulsesAll->ndetpulses+i)*convertAU2eV/1e3 + energyPCA1/1e3;	//keV
		}
		  
		gsl_matrix_free(pointsTranslatedRotated);
	}

	// Close intermediate output FITS file if it is necessary
	if ((*reconstruct_init)->intermediate == 1)
	{
		if (fits_close_file(dtcObject,&status))
		{
			message = "Cannot close file " + string(dtcName);
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	// Free allocated GSL vectors
	gsl_vector_free(invector);

	return;
}
/*xxxx end of SECTION A xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/* Run detection for the production mode only in multithread mode*/

std::mutex library_mut;
std::mutex fits_file_mut;

void th_runDetect(TesRecord* record, 
                  //int nRecord, 
                  int lastRecord, 
                  PulsesCollection *pulsesAll, 
                  ReconstructInitSIRENA** reconstruct_init, 
                  PulsesCollection** pulsesInRecord)
{
  scheduler* sc = scheduler::get();
  int inputPulseLength = (*reconstruct_init)->pulse_length;
  
  string message="";
  char valERROR[256];
  int status=EPOK;  
  
  // Declare variables
  fitsfile *inLibObject = NULL;// Object which contains 
                               //information of the library FITS file
  bool appendToLibrary = false;// Calibration library FITS file new 
                               //(appendToLibrary=false) or not 
                               //(appendToLibrary=true)
  
  fitsfile *dtcObject = NULL;  // Object which contains information of 
                               //the intermediate FITS file 
                               //('dtc' comes from 'detectFile')
  
  int eventsz = record->trigger_size;
  double tstartRecord;
  // It is not necessary to check the allocation because 
  //'record->trigger_size'='eventsz' has been checked previously
  gsl_vector *invector = gsl_vector_alloc(eventsz);    // Record
  
  if (((*reconstruct_init)->tstartPulse1 != 0) 
      && (((*reconstruct_init)->tstartPulse1 > record->trigger_size) 
          || ((*reconstruct_init)->tstartPulse2 > record->trigger_size) 
          || ((*reconstruct_init)->tstartPulse3 > record->trigger_size)))
    {
      message = "tstartPulsx can not be greater than the record size";
      EP_EXIT_ERROR(message,EPFAIL);
    }
  
  // Filter and differentiate the 'models' of the library 
  //(only for the first record)
  // thread safe
  // Not treadsafe for the library
  // lock here
  if (filderLibrary(reconstruct_init, 1/record->delta_t))
    {
      message = "Cannot run routine filderLibrary to filter "
        + string("& differentiate library if the 1st record");
      EP_EXIT_ERROR(message,EPFAIL);
    }
  
  // Store the input record in 'invector'
  if (loadRecord(record, &tstartRecord, &invector))
    {
      message = "Cannot run routine loadRecord";
      EP_EXIT_ERROR(message,EPFAIL);
    }
  eventsz = invector->size;// Just in case the last record has been filled 
                           //in with 0's => Re-allocate invector
  
  // Convert I into R if 'EnergyMethod' = I2R or I2RALL or I2RNOL or I2RFITTED
  // It is not necessary to check the allocation because 'invector' 
  // size must be > 0
  if ((strcmp((*reconstruct_init)->EnergyMethod,"I2R") == 0) 
      || (strcmp((*reconstruct_init)->EnergyMethod,"I2RALL") == 0) 
      || (strcmp((*reconstruct_init)->EnergyMethod,"I2RNOL") == 0)
      || (strcmp((*reconstruct_init)->EnergyMethod,"I2RFITTED") == 0))
    {
      // thread safe
      // Not thread safe for the record_file_ptr
      // If the cfitsio is not compiled with -D_REENTRANT, we need to lock
      // here in order to read the file.
      // Also if the reentrant mode is on, the threads should not share the 
      // same pointer.
      if (!sc->is_reentrant()){
        std::unique_lock<std::mutex> lk(fits_file_mut);
        /*if (convertI2R(reconstruct_init, &record, &invector))
          {
            lk.unlock();
            message = "Cannot run routine convertI2R";
            EP_EXIT_ERROR(message,EPFAIL);
          }*/
        if (convertI2R((*reconstruct_init)->EnergyMethod,(*reconstruct_init)->i2rdata->R0,(*reconstruct_init)->i2rdata->I0_START,(*reconstruct_init)->i2rdata->IMIN,(*reconstruct_init)->i2rdata->IMAX,(*reconstruct_init)->i2rdata->RPARA,(*reconstruct_init)->i2rdata->TTR,(*reconstruct_init)->i2rdata->LFILTER, 1/record->delta_t, &invector))
	{
                lk.unlock();
		message = "Cannot run routine convertI2R";
		EP_EXIT_ERROR(message,EPFAIL);
	}
		
	for (int i=0;i<invector->size;i++)		     // Because in 'runEnergy' the record (TesRecord) is used => The I2R, I2RALL, I2RNOL or I2RFITTED transformed record has to be used
        {
                record->adc_double[i] = gsl_vector_get(invector,i);
        }  
        lk.unlock();
      }else{
        if (convertI2R((*reconstruct_init)->EnergyMethod,(*reconstruct_init)->i2rdata->R0,(*reconstruct_init)->i2rdata->I0_START,(*reconstruct_init)->i2rdata->IMIN,(*reconstruct_init)->i2rdata->IMAX,(*reconstruct_init)->i2rdata->RPARA,(*reconstruct_init)->i2rdata->TTR,(*reconstruct_init)->i2rdata->LFILTER, 1/record->delta_t, &invector))
	{
		message = "Cannot run routine convertI2R";
		EP_EXIT_ERROR(message,EPFAIL);
	}
		
	for (int i=0;i<invector->size;i++)		     // Because in 'runEnergy' the record (TesRecord) is used => The I2R, I2RALL, I2RNOL or I2RFITTED transformed record has to be used
        {
                record->adc_double[i] = gsl_vector_get(invector,i);
        }  
        /*if (convertI2R(reconstruct_init, &record, &invector))
        {
          message = "Cannot run routine convertI2R";
          EP_EXIT_ERROR(message,EPFAIL);
        }*/
      }
    }
  
  // Process each record
  // thread safe
  if (procRecord(reconstruct_init, tstartRecord, 1/record->delta_t, dtcObject, 
                 invector, *pulsesInRecord))
    {
      message = "Cannot run routine procRecord for record processing";
      EP_EXIT_ERROR(message,EPFAIL);
    }
  
  if ((strcmp((*reconstruct_init)->EnergyMethod,"I2R") == 0) 
      || (strcmp((*reconstruct_init)->EnergyMethod,"I2RALL") == 0) 
      || (strcmp((*reconstruct_init)->EnergyMethod,"I2RNOL") == 0)
      || (strcmp((*reconstruct_init)->EnergyMethod,"I2RFITTED") == 0))
    {
      strcpy((*reconstruct_init)->EnergyMethod,"OPTFILT"); 
      // From this point forward, I2R, I2RALL, I2RNOL and I2RFITTED 
      //are completely equivalent to OPTFILT
    }

  // PCA is used  when pulses are farther than 'PulseLength'!!!!!!!!!!
  if ((lastRecord == 1) 
      && (strcmp((*reconstruct_init)->EnergyMethod,"PCA") == 0) 
      && ((pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses))>0)// PCA
    {
      (*reconstruct_init)->pulse_length = inputPulseLength;
      
      // In order to not have restrictions when providing 
      //(*reconstruct_init)->energyPCAx
      double energyPCA1, energyPCA2;
      if ((*reconstruct_init)->energyPCA2 < (*reconstruct_init)->energyPCA1)
        {
          energyPCA1 = (*reconstruct_init)->energyPCA2;
          energyPCA2 = (*reconstruct_init)->energyPCA1;
        }
      else
        {
          energyPCA1 = (*reconstruct_init)->energyPCA1;
          energyPCA2 = (*reconstruct_init)->energyPCA2;
        }
      
      // Covariance data
      // It is not necessary to check the allocation because 'PulseLength' 
      //(input parameter) has been checked previously 
      gsl_vector *pulsetemplate = 
        gsl_vector_alloc((*reconstruct_init)->pulse_length);
      
      gsl_vector_set_zero(pulsetemplate);
      
      gsl_matrix *covarianceData = 
        gsl_matrix_alloc((*reconstruct_init)->pulse_length,
                         (*reconstruct_init)->pulse_length);
      
      gsl_matrix *weightData = 
        gsl_matrix_alloc((*reconstruct_init)->pulse_length,
                         (*reconstruct_init)->pulse_length);

      gsl_vector *nonpileup;
      if ((nonpileup = 
           gsl_vector_alloc((pulsesAll->ndetpulses)
                            +((*pulsesInRecord)->ndetpulses))) 
          == 0)
        {  
          sprintf(valERROR,"%d",__LINE__-2);
          string str(valERROR);
          message = "Allocating with <= 0 size in line " + str 
            + " (" + __FILE__ + ")";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      gsl_vector_set_all(nonpileup,1);
      // thread safe
      if (weightMatrix(*reconstruct_init, false, pulsesAll, *pulsesInRecord, 
                       (pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses), 
                       nonpileup , pulsetemplate, &covarianceData, &weightData))
        {
          message = "Cannot run weightMatrix routine";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      gsl_vector_free(pulsetemplate); pulsetemplate = 0;
      gsl_matrix_free(weightData); weightData = 0;
      
      // Eigenvalues and eigenvectors
      gsl_vector *eigenvalues;
      gsl_matrix *eigenvectors;// Eigenvectors in columns
      if (eigenVV(covarianceData,&eigenvectors,&eigenvalues))
        {
          message = "Cannot run eigenVV routine";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      gsl_matrix_free(covarianceData); covarianceData = 0;
      gsl_vector_free(eigenvalues); eigenvalues = 0;
      
      // RSxN (S=2)
      gsl_matrix *RowFeatureVectors;// Eigenvectors in rows
      if ((RowFeatureVectors = gsl_matrix_alloc(eigenvectors->size2,
                                                eigenvectors->size1)) 
          == 0)
        {  
          sprintf(valERROR,"%d",__LINE__-2);
          string str(valERROR);
          message = "Allocating with <= 0 size in line " + str + " (" 
            + __FILE__ + ")";
          EP_EXIT_ERROR(message,EPFAIL);
        }

      gsl_matrix_transpose_memcpy(RowFeatureVectors,eigenvectors);
      gsl_matrix_free(eigenvectors); eigenvectors = 0;
      
      // It is not necessary to check the allocation because 'PulseLength'
      //(input parameter) and 
      //'(pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulsesc' 
      //have been checked previously 
      // S0^1-M0 S0^2-M0...............S0^(nonpileupPulses)-M0
      // S1^1-M1 S1^2-M1...............S1^(nonpileupPulses)-M1
      // ...................................................
      // S1023^1-M1023 S1023^2-M1023...S1023^(nonpileupPulses)-M1023
      gsl_matrix *RowDataAdjust = 
        gsl_matrix_alloc((*reconstruct_init)->pulse_length,
                         (pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses)); 
      
      for (int m=0;m<(*reconstruct_init)->pulse_length;m++)
        {
          for (int p=0;p<pulsesAll->ndetpulses;p++)
            {
              if (gsl_vector_get(nonpileup,p) == 1)
                {
                  if (m == pulsesAll->pulses_detected[p].pulse_adc->size)
                    {
                      sprintf(valERROR,"%d",__LINE__+5);
                      string str(valERROR);
                      message = "PCA should be used with pulses farther than "
                        + string("'PulseLength' => Getting i-th element of vector ")
                        + "out of range in line " + str + " (" + __FILE__ + ")";
                      EP_EXIT_ERROR(message,EPFAIL);
                    }
                  gsl_matrix_set(RowDataAdjust,m,p,
                                 gsl_vector_get(pulsesAll->pulses_detected[p].pulse_adc,m));
                }
            }
          for (int p=0;p<(*pulsesInRecord)->ndetpulses;p++)
            {
              if (gsl_vector_get(nonpileup,pulsesAll->ndetpulses+p) == 1)
                {
                  if (m == (*pulsesInRecord)->pulses_detected[p].pulse_adc->size)
                    {
                      sprintf(valERROR,"%d",__LINE__+5);
                      string str(valERROR);
                      message = "PCA should be used with pulses farther than "
                        + string("'PulseLength' => Getting i-th element of vector ")
                        + "out of range in line " + str + " (" + __FILE__ + ")";
                      EP_EXIT_ERROR(message,EPFAIL);
                    }
                  gsl_matrix_set(RowDataAdjust,m,pulsesAll->ndetpulses+p,
                                 gsl_vector_get((*pulsesInRecord)->pulses_detected[p].pulse_adc,m));
                }
            }
        }
      gsl_vector_free(nonpileup); nonpileup = 0;
      
      // It is not necessary to check the allocation because the allocation 
      //of 'RowFeatureVectors' and 
      //'(pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses' 
      //have been checked previously
      gsl_matrix *RSrxN = 
        gsl_matrix_alloc(RowFeatureVectors->size1,
                         (pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses));
      if (RowFeatureVectors->size2 != RowDataAdjust->size1)
        {
          sprintf(valERROR,"%d",__LINE__+5);
          string str(valERROR);
          message = "Wrong dimensions to compute matrix-matrix product in line "
            + str + " (" + __FILE__ + ")";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, 
                      RowFeatureVectors, RowDataAdjust,0.0, RSrxN);
      gsl_matrix_free(RowFeatureVectors); RowFeatureVectors = 0;
      gsl_matrix_free(RowDataAdjust); RowDataAdjust = 0;
      
      // Store RSxN in a FITS file
      FILE * temporalFile;
      char temporalFileName[255];
      char val[256];
      sprintf(temporalFileName,"%s_rsxn.txt",(*reconstruct_init)->detectFile);
      temporalFile = fopen (temporalFileName,"w");
      for (int i = 0; i < RSrxN->size2; i++) // All the pulses in a column
        {
          sprintf(val,"%e %e",gsl_matrix_get(RSrxN,0,i),gsl_matrix_get(RSrxN,1,i));
          strcat(val,"\n");
          fputs(val,temporalFile);
        }
      fclose(temporalFile);
      
      // AE straight line: Pto0(x,y) and Pto10(x,y)
      double a,b;      //y=ax+b
      // It is not necessary to check the allocation because '(pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses' has been checked previously
      gsl_vector *xfit = gsl_vector_alloc((pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses));
      gsl_vector *yfit = gsl_vector_alloc((pulsesAll->ndetpulses)+((*pulsesInRecord)->ndetpulses));
      gsl_matrix_get_row(xfit,RSrxN,0);
      gsl_matrix_get_row(yfit,RSrxN,1);
      if (polyFitLinear (xfit, yfit, &a, &b))
        {
          message = "Cannot run routine polyFitLinear";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      gsl_vector *Pto0 = gsl_vector_alloc(2);
      gsl_vector *Pto10 = gsl_vector_alloc(2);
      // It is not necessary to check the allocation because the allocation of 'RSrxN' has been checked previously
      gsl_vector *xRSrxN = gsl_vector_alloc(RSrxN->size2);
      gsl_matrix_get_row(xRSrxN,RSrxN,0);
      gsl_vector_set(Pto0,0,gsl_vector_get(xRSrxN,gsl_vector_min_index(xRSrxN)));
      gsl_vector_set(Pto0,1,a*gsl_vector_get(Pto0,0)+b);
      gsl_vector_set(Pto10,0,gsl_vector_get(xRSrxN,gsl_vector_max_index(xRSrxN))-gsl_vector_get(Pto0,0));
      gsl_vector_set(Pto10,1,a*gsl_vector_get(xRSrxN,gsl_vector_max_index(xRSrxN))+b-gsl_vector_get(Pto0,1));
      
      // Calculus of the rotation angle
      gsl_vector *u = gsl_vector_alloc(2);
      gsl_vector_set(u,0,1);
      gsl_vector_set(u,1,0);
      gsl_vector *v = gsl_vector_alloc(2);
      gsl_vector_memcpy(v,Pto10);
      gsl_vector_free(Pto10); Pto10 = 0;
      double dotResult;
      gsl_blas_ddot (u, v, &dotResult);
      double alfa;
      alfa = acos(dotResult/(gsl_blas_dnrm2(u)*gsl_blas_dnrm2(v)));
      
      // Rotation
      gsl_matrix *RotationMatrix = gsl_matrix_alloc(2,2);
      gsl_matrix_set(RotationMatrix,0,0,cos(alfa));
      gsl_matrix_set(RotationMatrix,0,1,sin(alfa));
      gsl_matrix_set(RotationMatrix,1,0,-sin(alfa));
      gsl_matrix_set(RotationMatrix,1,1,cos(alfa));
      // It is not necessary to check the allocation because the allocation of 'RSrxN' has been checked previously
      gsl_matrix *pointsTranslated = gsl_matrix_alloc(2,RSrxN->size2);
      for (int i = 0; i < RSrxN->size2; i++)
        {
          gsl_matrix_set(pointsTranslated,0,i,gsl_matrix_get(RSrxN,0,i)-gsl_vector_get(Pto0,0));
          gsl_matrix_set(pointsTranslated,1,i,gsl_matrix_get(RSrxN,1,i)-gsl_vector_get(Pto0,1));
        }
      gsl_vector_free(Pto0); Pto0 = 0;
      // It is not necessary to check the allocation because the allocation of 'RSrxN' has been checked previously
      gsl_matrix *pointsTranslatedRotated = gsl_matrix_alloc(2,RSrxN->size2);
      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, RotationMatrix, pointsTranslated,0.0, pointsTranslatedRotated);
      /*for (int i = 0; i < RSrxN->size2; i++)
        {
        gsl_matrix_set(pointsTranslatedRotated,0,i,gsl_matrix_get(pointsTranslatedRotated,0,i)+gsl_vector_get(Pto0,0));
        }
        sprintf(temporalFileName,"rsxnTransformed.txt");
        temporalFile = fopen (temporalFileName,"w");
        for (int i = 0; i < RSrxN->size2; i++) // All the pulses in a column
        {
        sprintf(val,"%e %e",gsl_matrix_get(pointsTranslatedRotated,0,i),gsl_matrix_get(pointsTranslatedRotated,1,i));
        strcat(val,"\n");
        fputs(val,temporalFile);
        }
        fclose(temporalFile);*/
      gsl_matrix_free(pointsTranslated); pointsTranslated = 0;
      gsl_matrix_free(RSrxN); RSrxN = 0;
      // It is not necessary to check the allocation because the allocation of 'pointsTranslatedRotated' has been checked previously
      gsl_vector *pointsTranslatedRotated1aux = gsl_vector_alloc(pointsTranslatedRotated->size2);
      gsl_vector *pointsTranslatedRotated2aux = gsl_vector_alloc(pointsTranslatedRotated->size2);
      int num1 = 0;
      int num2 = 0;
      gsl_vector *xpointsTranslatedRotated = gsl_vector_alloc(pointsTranslatedRotated->size2);
      gsl_matrix_get_row(xpointsTranslatedRotated,pointsTranslatedRotated,0);
      double midpoint = gsl_vector_min(xpointsTranslatedRotated)+(gsl_vector_max(xpointsTranslatedRotated)-gsl_vector_min(xpointsTranslatedRotated))/2;
      for (int i = 0; i < xpointsTranslatedRotated->size; i++)
        {
          if (gsl_vector_get(xpointsTranslatedRotated,i) <= midpoint)
            {
              if (num1 >= pointsTranslatedRotated1aux->size)
                {
                  sprintf(valERROR,"%d",__LINE__+5);
                  string str(valERROR);
                  message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                  EP_EXIT_ERROR(message,EPFAIL);
                }
              gsl_vector_set(pointsTranslatedRotated1aux,num1,gsl_vector_get(xpointsTranslatedRotated,i));
              num1++;
            }
          else
            {
              if (num2 >= pointsTranslatedRotated2aux->size)
                {
                  sprintf(valERROR,"%d",__LINE__+5);
                  string str(valERROR);
                  message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                  EP_EXIT_ERROR(message,EPFAIL);
                }
              gsl_vector_set(pointsTranslatedRotated2aux,num2,gsl_vector_get(xpointsTranslatedRotated,i));
              num2++;
            }
        }
      gsl_vector_free(xpointsTranslatedRotated); xpointsTranslatedRotated = 0;
      gsl_vector *pointsTranslatedRotated1;
      gsl_vector *pointsTranslatedRotated2;
      if ((pointsTranslatedRotated1 = gsl_vector_alloc(num1)) == 0)
        {
          sprintf(valERROR,"%d",__LINE__-2);
          string str(valERROR);
          message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      if ((pointsTranslatedRotated2 = gsl_vector_alloc(num2)) == 0)
        {
          sprintf(valERROR,"%d",__LINE__-2);
          string str(valERROR);
          message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      
      gsl_vector_view temp;
      if ((num1 < 1) || (num1 > pointsTranslatedRotated1aux->size))
        {
          sprintf(valERROR,"%d",__LINE__+5);
          string str(valERROR);
          message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      temp = gsl_vector_subvector(pointsTranslatedRotated1aux,0,num1);
      if (gsl_vector_memcpy(pointsTranslatedRotated1,&temp.vector) != 0)
        {
          sprintf(valERROR,"%d",__LINE__-2);
          string str(valERROR);
          message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      if ((num2 < 1) || (num1 > pointsTranslatedRotated2aux->size))
        {
          sprintf(valERROR,"%d",__LINE__+5);
          string str(valERROR);
          message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      temp = gsl_vector_subvector(pointsTranslatedRotated2aux,0,num2);
      if (gsl_vector_memcpy(pointsTranslatedRotated2,&temp.vector) != 0)
        {
          sprintf(valERROR,"%d",__LINE__-2);
          string str(valERROR);
          message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      gsl_vector_free(pointsTranslatedRotated1aux); pointsTranslatedRotated1aux = 0;
      gsl_vector_free(pointsTranslatedRotated2aux); pointsTranslatedRotated2aux = 0;
      
      // Histograms of the two clusters (two energies)
      int nbins = 80;
      assert(nbins > 0);
      // It is not necessary to check the allocation because 'nbins' has been checked previously
      gsl_vector *xhisto1 = gsl_vector_alloc(nbins);
      gsl_vector *yhisto1 = gsl_vector_alloc(nbins);
      gsl_vector *xhisto2 = gsl_vector_alloc(nbins);
      gsl_vector *yhisto2 = gsl_vector_alloc(nbins);
      bool histo1neg = false;
      bool histo2neg = false;
      double histo1constant;
      double histo2constant;
      if (gsl_vector_isnonneg(pointsTranslatedRotated1) != 1)
        {
          histo1neg = true;
          histo1constant = fabs(gsl_vector_min(pointsTranslatedRotated1));
          gsl_vector_add_constant(pointsTranslatedRotated1,histo1constant);
        }
      if (createHisto(pointsTranslatedRotated1, nbins, &xhisto1, &yhisto1))
        {
          message = "Cannot run createHisto routine";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      if (histo1neg == true) 
        {
          gsl_vector_add_constant(xhisto1,-histo1constant);
          gsl_vector_add_constant(pointsTranslatedRotated1,-histo1constant);
        }
      sprintf(temporalFileName,"%s_histo1.txt",(*reconstruct_init)->detectFile);
      temporalFile = fopen (temporalFileName,"w");
      for (int i = 0; i < xhisto1->size; i++) 
        {
          sprintf(val,"%e %e",gsl_vector_get(xhisto1,i),gsl_vector_get(yhisto1,i));
          strcat(val,"\n");
          fputs(val,temporalFile);
        }
      fclose(temporalFile);
      if (gsl_vector_isnonneg(pointsTranslatedRotated2) != 1)
        {
          histo2neg = true;
          histo2constant = fabs(gsl_vector_min(pointsTranslatedRotated2));
          gsl_vector_add_constant(pointsTranslatedRotated2,histo2constant);
        }
      if (createHisto(pointsTranslatedRotated2, nbins, &xhisto2, &yhisto2))
        {
          message = "Cannot run createHisto routine";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      if (histo2neg == true) 
        {
          gsl_vector_add_constant(xhisto2,-histo2constant);
          gsl_vector_add_constant(pointsTranslatedRotated2,-histo2constant);
        }
      sprintf(temporalFileName,"%s_histo2.txt",(*reconstruct_init)->detectFile);
      temporalFile = fopen (temporalFileName,"w");
      for (int i = 0; i < xhisto2->size; i++) 
        {
          sprintf(val,"%e %e",gsl_vector_get(xhisto2,i),gsl_vector_get(yhisto2,i));
          strcat(val,"\n");
          fputs(val,temporalFile);
        }
      fclose(temporalFile);
      double data1[num1];
      double data2[num2];
      if (num1 > pointsTranslatedRotated1->size)
        {
          sprintf(valERROR,"%d",__LINE__+7);
          string str(valERROR);
          message = "Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      for (int i = 0; i < num1; i++)
        {
          data1[i] = gsl_vector_get(pointsTranslatedRotated1,i);
        }
      if (num2 > pointsTranslatedRotated2->size)
        {
          sprintf(valERROR,"%d",__LINE__+7);
          string str(valERROR);
          message = "Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
          EP_EXIT_ERROR(message,EPFAIL);
        }
      for (int i = 0; i < num2; i++)
        {
          data2[i] = gsl_vector_get(pointsTranslatedRotated2,i);
        }
      double sigma1 = gsl_stats_sd(data1, 1, num1);
      double sigma2 = gsl_stats_sd(data2, 1, num2);
      double maxcenter1 = gsl_stats_mean(data1, 1, num1);
      double maxcenter2 = gsl_stats_mean(data2, 1, num2);
      
      // Conversion factor from arbitrary unit to eV
      double convertAU2eV = (energyPCA2 - energyPCA1)/(maxcenter2-maxcenter1);	
      
      gsl_vector_free(pointsTranslatedRotated1); pointsTranslatedRotated1 = 0;
      gsl_vector_free(pointsTranslatedRotated2); pointsTranslatedRotated2 = 0;
      gsl_vector_free(xhisto1); xhisto1 = 0;
      gsl_vector_free(yhisto1); yhisto1 = 0;
      gsl_vector_free(xhisto2); xhisto2 = 0;
      gsl_vector_free(yhisto2); yhisto2 = 0;
      
      // Energy calculation
      for (int i = 0; i < pulsesAll->ndetpulses; i++)
        {
          pulsesAll->pulses_detected[i].energy = 
            gsl_matrix_get(pointsTranslatedRotated,0,i)*convertAU2eV/1e3 
            + energyPCA1/1e3;//keV
        }
      for (int i = 0; i < (*pulsesInRecord)->ndetpulses; i++)
        {
          (*pulsesInRecord)->pulses_detected[i].energy = 
            gsl_matrix_get(pointsTranslatedRotated,0,
                           pulsesAll->ndetpulses+i)*convertAU2eV/1e3 
            + energyPCA1/1e3;//keV
        }
      
      gsl_matrix_free(pointsTranslatedRotated); pointsTranslatedRotated = 0;
    }
  
  // Free allocated GSL vectors
  gsl_vector_free(invector); invector = 0;
  
  return;
}

/***** SECTION A1 ************************************************************
* createLibrary: This function creates the calibration library FITS file, if it does not exist yet. Otherwise, it opens it (to add a new row).
*
* - If it exists => Open it and set 'appendToLibrary' = true
* - If it does not exist => Create it and set 'appendToLibrary' = false
* 	- Write keyword 'EVENTCNT'=1 in the LIBRARY HDU
* 	- Write the whole list of input parameters in "HISTORY" in the Primary HDU
* 
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
*                     In particular, this function uses parameters to call the library ('library_file') and to write input parameters
*                     info in the Primary HDU of the library file
* - appendToLibrary: Used by the function 'writeLibrary'
* - inLibObject: Object which contains information of the library FITS file (used also by 'writeLibrary') 
* - inputPulseLength: PulseLength input parameter
****************************************************************************/
int createLibrary(ReconstructInitSIRENA* reconstruct_init, bool *appendToLibrary, fitsfile **inLibObject, int inputPulseLength)
{
	int status = EPOK;
	string message = "";

	char inLibName[256];
	strncpy(inLibName, reconstruct_init->library_file,255);
	inLibName[255]='\0';
	
	char extname[10];
	strncpy(extname,"LIBRARY",9);
	extname[9]='\0';

	char keyvalstr[1000];
	char *tform[1];
	char *ttype[1];
	char *tunit[1];
	char keyname[10];

	long eventcntLib;

	// If the library exists => Open it
	if (fileExists(string(inLibName)))
	{
		*appendToLibrary = true;

		if (fits_open_file(inLibObject,inLibName,READWRITE,&status))
		{
			message = "Cannot open library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
		{
			message = "Cannot move to HDU " + string(extname);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		
		long numrowsLib;
		if (fits_get_num_rows(*inLibObject,&numrowsLib, &status))
		{
			message = "Cannot get number of rows in " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		for (int i=0;i<reconstruct_init->library_collection->ntemplates;i++)
		{
			if (reconstruct_init->monoenergy == gsl_vector_get(reconstruct_init->library_collection->energies,i))
			{
				message = "Energy which is going to be added to the library is already in the library";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
		}
	}
	else	// If the library does not exist yet => Create it
	{
		*appendToLibrary = false;
		if (fits_create_file(inLibObject, inLibName, &status))
		{
			message = "Cannot create library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		if (fits_open_file(inLibObject,inLibName,READWRITE,&status))
		{
			message = "Cannot open library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Create LIBRARY HDU
		if (fits_create_tbl(*inLibObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
		{
			message = "Cannot create table " + string(extname);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		eventcntLib = 1;
		strncpy(keyname,"EVENTCNT",9);
		keyname[9]='\0';
		assert(eventcntLib > 0);
		if (fits_write_key(*inLibObject,TLONG,keyname,&eventcntLib,NULL,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library";
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		
		strncpy(extname,"FIXFILTT",9);
		extname[9]='\0';
		// Create FIXFILTT HDU
		if (fits_create_tbl(*inLibObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
		{
			message = "Cannot create table " + string(extname);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		
		strncpy(extname,"FIXFILTF",9);
		extname[9]='\0';
		// Create FIXFILTF HDU
		if (fits_create_tbl(*inLibObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
		{
			message = "Cannot create table " + string(extname);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		
		if (reconstruct_init->hduPRECALWN == 1)
		{
			strncpy(extname,"PRECALWN",9);
			extname[9]='\0';
			// Create PRECALWN HDU
			if (fits_create_tbl(*inLibObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
			{
				message = "Cannot create table " + string(extname);
				EP_PRINT_ERROR(message,status); return(EPFAIL);
			}
		}
		
		if (reconstruct_init->hduPRCLOFWM == 1)
		{
			strncpy(extname,"PRCLOFWM",9);
			extname[9]='\0';
			// Create PRCLOFWM HDU
			if (fits_create_tbl(*inLibObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
			{
				message = "Cannot create table " + string(extname);
				EP_PRINT_ERROR(message,status); return(EPFAIL);
			}
		}

		// Primary HDU
		strncpy(extname,"Primary",9);
		extname[9]='\0';
		if (fits_movabs_hdu(*inLibObject, 1, NULL, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strncpy(keyname,"PROC",9);
		keyname[9]='\0';

		string strproc;
		
		char str_largeFilter[125];	snprintf(str_largeFilter,125,"%d",reconstruct_init->largeFilter);
		strproc=string("largeFilter = ") + string(str_largeFilter);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
	               
                HDpar_stamp(*inLibObject, 0, &status);  // Write the whole list of input parameters in HISTORY
		
		if (status != 0)
		{
			message = "Cannot write keyword in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	return (EPOK);
}
/*xxxx end of SECTION A1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A2 ************************************************************
* createDetectFile function: This function creates an intermediate FITS file with some useful info (during the development phase) if the 'intermediate' input parameter is set to 1.
*                            The intermediate FITS file will contain 2 HDUs, PULSES and TESTINFO.
*                            The PULSES HDU will contain some info the found pulses: TSTART, I0 (the pulse itself), TEND, TAURISE, TAUFALL and QUALITY.
*                            The TESTINFO HDU contains FILDER (the low-pass filtered and differentiated records) and THRESHOLD.
*
* If file exists => Check 'clobber' for overwritting it or not. If it does not, then create it.
* 
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
*                     In particular, this function uses some of its values to call the intermediate file ('detectFile'), to write input parameters
*                     info in the Primary HDU of the intermediate file and also 'clobber'
* - samprate: Sampling rate (calculated from input FITS keyword 'DELTAT')
* - dtcObject: Object which contains information of the intermediate FITS file (used also by 'writeTestInfo' and 'writePulses') 
* - inputPulseLength: PulseLength input parameter
***************************************************************************/
int createDetectFile(ReconstructInitSIRENA* reconstruct_init, double samprate, fitsfile **dtcObject, int inputPulseLength)
{
	int status = EPOK;
	string message = "";

	char dtcName[256];
	strncpy(dtcName,reconstruct_init->detectFile,255);
	dtcName[255]='\0'; // enforce zero ending string in case of buffer overflows

	// Create intermediate output FITS file: If it does not exist yet
	// If dtcName does not finish as '.fits' and the file dtcName+'.fits' already exists =>
	// => Data are appended to dtcName file => Must not be allowed
	// '.fits' => 5 characters
	if (strlen(dtcName)<6)
	{
		// dtcName has 5 or less characters => Does not contain '.fits' =>Append '.fits' to dtcName
		strcat(dtcName,".fits");
	}
	else
	{
		// Check if dtcName has '.fits' and if not, append it
		if (strncmp(strndup(dtcName+strlen(dtcName)-5, 5),".fits",5) != 0)
		{
			// dtcName does not finish as '.fits' => Append '.fits' to dtcName
			strncpy(dtcName,strcat(dtcName,".fits"),255);
			dtcName[255]='\0';
		}
	}

	// If firstRecord => Create intermediate file (if file already exists => check clobber)
	if (fileExists(string(dtcName)) && (reconstruct_init->clobber == 1))
	{
		if (remove(dtcName))
		{
			message = "Output intermediate file already exists & cannot be deleted ("+string(strerror(errno))+")";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
	}
	else if (fileExists(string(dtcName)) && (reconstruct_init->clobber == 0))
	{
		message = "Output intermediate file already exists: must not be overwritten";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	
	if (!fileExists(string(dtcName)))
	{
		if (fits_create_file(dtcObject, dtcName, &status))
		{
			message = "Cannot create output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		message = "Create Detect Fits File: " + string(dtcName);
	}

	// Create HDU PULSES
	// To work with tables (HDUs)
	char *tt[1];
	char *tf[1];
	char *tu[1];
	
	// To write keywords
	char keyname[10];
	char keyvalstr[1000];
	char extname[10];

	if (fits_open_file(dtcObject,dtcName,READWRITE,&status))
	{
		message = "Cannot open output intermediate file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (reconstruct_init->clobber == 1)
	{
		strncpy(extname,"PULSES",9);
		extname[9]='\0';
		
		// PULSES HDU
		if (fits_create_tbl(*dtcObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
			message = "Cannot create table " + string(extname) + " in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		if (fits_movnam_hdu(*dtcObject, ANY_HDU,extname, 0, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strncpy(keyname,"OPMODE",9);
		keyname[9]='\0';
		fits_write_key(*dtcObject,TINT,keyname,&(reconstruct_init->opmode),NULL,&status);
		
		strncpy(keyname,"EVENTSZ",9);
		keyname[9]='\0';
		assert(reconstruct_init->pulse_length > 0);
		fits_write_key(*dtcObject,TINT,keyname,&(reconstruct_init->pulse_length),NULL,&status);
		if (reconstruct_init->opmode == 0)
		{
			strcpy(keyname,"ENERGY");
			fits_write_key(*dtcObject,TDOUBLE,keyname,&(reconstruct_init-> monoenergy),NULL,&status);
		}
		strncpy(keyname,"SAMPRATE",9);
		keyname[9]='\0';
		assert(samprate > 0);
		fits_write_key(*dtcObject,TDOUBLE,keyname,&samprate,NULL,&status);
		if (status != 0)
		{
			message = "Cannot write keyword some keyword in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// TESTINFO HDU
		strncpy(extname,"TESTINFO",9);
		extname[9]='\0';
		if (fits_create_tbl(*dtcObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
			message = "Cannot create table " + string(extname) + " in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Primary HDU
		strncpy(extname,"Primary",9);
		extname[9]='\0';
		if (fits_movabs_hdu(*dtcObject, 1, NULL, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		
		strncpy(keyname,"HISTORY",9);
		keyname[9]='\0';
		const char * charhistory= "HISTORY Starting parameter list";
		strncpy(keyvalstr,charhistory,999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		string strhistory;
		char comment[MAXMSG];
		
		sprintf(comment, "RecordFile = %s", reconstruct_init->record_file);
		fits_write_comment(*dtcObject, comment, &status);
		
		sprintf(comment, "TesEventFile = %s", reconstruct_init->event_file);
		fits_write_comment(*dtcObject, comment, &status);
		
		sprintf(comment, "LibraryFile = %s", reconstruct_init->library_file);
		fits_write_comment(*dtcObject, comment, &status);
		
		sprintf(comment, "NoiseFile = %s", reconstruct_init->noise_file);
		fits_write_comment(*dtcObject, comment, &status);

		char str_opmode[125];		snprintf(str_opmode,125,"%d",reconstruct_init->opmode);
		strhistory = string("opmode = ") + string(str_opmode);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strhistory=string("FilterDomain = ") + reconstruct_init->FilterDomain;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strhistory=string("FilterMethod = ") + reconstruct_init->FilterMethod;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strhistory=string("EnergyMethod = ") + reconstruct_init->EnergyMethod;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_LagsOrNot[125];      	snprintf(str_LagsOrNot,125,"%d",reconstruct_init->LagsOrNot);
		strhistory=string("LagsOrNot = ") + string(str_LagsOrNot);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
                
                char str_nLags[125];      	snprintf(str_nLags,125,"%d",reconstruct_init->nLags);
		strhistory=string("nLags = ") + string(str_nLags);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
                
                char str_Parabola3OrFitting5[125];      	snprintf(str_Parabola3OrFitting5,125,"%d",reconstruct_init->Fitting35);
		strhistory=string("Fitting35 = ") + string(str_Parabola3OrFitting5);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_OFIter[125];     	snprintf(str_OFIter,125,"%d",reconstruct_init->OFIter);
		strhistory=string("OFIter = ") + string(str_OFIter);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_OFLib[125];      	snprintf(str_OFLib,125,"%d",reconstruct_init->OFLib);
		strhistory=string("OFLib = ") + string(str_OFLib);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
				
		strhistory=string("OFInterp = ") + reconstruct_init->OFInterp;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
				
		strhistory=string("OFStrategy = ") + reconstruct_init->OFStrategy;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
				
		char str_OFLength[125];		snprintf(str_OFLength,125,"%d",reconstruct_init->OFLength);
		strhistory=string("OFLength = ") + string(str_OFLength);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);

		char str_maxPulsesPerRecord[125];	snprintf(str_maxPulsesPerRecord,125,"%d",reconstruct_init->maxPulsesPerRecord);
		strhistory=string("maxPulsesPerRecord = ") + string(str_maxPulsesPerRecord);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_pulse_length[125];	snprintf(str_pulse_length,125,"%d",inputPulseLength);
		strhistory=string("PulseLength = ") + string(str_pulse_length);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_scaleFactor[125];	snprintf(str_scaleFactor,125,"%f",reconstruct_init->scaleFactor);
		strhistory=string("scaleFactor = ") + string(str_scaleFactor);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_samplesUp[125];	snprintf(str_samplesUp,125,"%f",reconstruct_init->samplesUp);
		strhistory=string("samplesUp = ") + string(str_samplesUp);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_nSgms[125];	    	snprintf(str_nSgms,125,"%f",reconstruct_init->nSgms);
		strhistory=string("nSgms = ") + string(str_nSgms);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_LrsT[125];		snprintf(str_LrsT,125,"%e",reconstruct_init->LrsT);
		strhistory=string("LrsT = ") + string(str_LrsT);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_LbT[125];		snprintf(str_LbT,125,"%e",reconstruct_init->LbT);
		strhistory=string("LbT = ") + string(str_LbT);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_monoenergy[125];	snprintf(str_monoenergy,125,"%f",reconstruct_init->monoenergy);
		strhistory=string("monoenergy = ") + string(str_monoenergy);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_largeFilter[125];	snprintf(str_largeFilter,125,"%d",reconstruct_init->largeFilter);
		strhistory=string("largeFilter = ") + string(str_largeFilter);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_intermediate[125];      snprintf(str_intermediate,125,"%d",reconstruct_init->intermediate);
		strhistory=string("intermediate = ") + string(str_intermediate);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		sprintf(comment, "detectFile = %s", reconstruct_init->detectFile);
		fits_write_comment(*dtcObject, comment, &status);
		
		sprintf(comment, "fiterFile = %s", reconstruct_init->filterFile);
		fits_write_comment(*dtcObject, comment, &status);
		
		char str_clobber[125];      	snprintf(str_clobber,125,"%d",reconstruct_init->clobber);
		strhistory=string("clobber = ") + string(str_clobber);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_SaturationValue[125];	snprintf(str_SaturationValue,125,"%e",reconstruct_init->SaturationValue);
		strhistory=string("SaturationValue = ") + string(str_SaturationValue);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_tstartPulse1[125];	
		if (reconstruct_init->tstartPulse1 == 0)			  	snprintf(str_tstartPulse1,125,"%d",reconstruct_init->tstartPulse1);
		else if (strcmp(reconstruct_init->EnergyMethod,"I2RALL") == 0)		snprintf(str_tstartPulse1,125,"%d",reconstruct_init->tstartPulse1+2);
		else if (strcmp(reconstruct_init->EnergyMethod,"I2RALL") != 0)		snprintf(str_tstartPulse1,125,"%d",reconstruct_init->tstartPulse1+1);
		strhistory=string("tstartPulse1 = ") + string(str_tstartPulse1);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_tstartPulse2[125];	
		if (reconstruct_init->tstartPulse2 == 0)			  	snprintf(str_tstartPulse2,125,"%d",reconstruct_init->tstartPulse2);
		else if (strcmp(reconstruct_init->EnergyMethod,"I2RALL") == 0)		snprintf(str_tstartPulse2,125,"%d",reconstruct_init->tstartPulse2+2);
		else if (strcmp(reconstruct_init->EnergyMethod,"I2RALL") != 0)		snprintf(str_tstartPulse2,125,"%d",reconstruct_init->tstartPulse2+1);
		strhistory=string("tstartPulse2 = ") + string(str_tstartPulse2);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_tstartPulse3[125];	
		if (reconstruct_init->tstartPulse3 == 0)			  	snprintf(str_tstartPulse3,125,"%d",reconstruct_init->tstartPulse3);
		else if (strcmp(reconstruct_init->EnergyMethod,"I2RALL") == 0)		snprintf(str_tstartPulse3,125,"%d",reconstruct_init->tstartPulse3+2);
		else if (strcmp(reconstruct_init->EnergyMethod,"I2RALL") != 0)		snprintf(str_tstartPulse3,125,"%d",reconstruct_init->tstartPulse3+1);
		strhistory=string("tstartPulse3 = ") + string(str_tstartPulse3);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		charhistory= "HISTORY Ending parameter list";
		strncpy(keyvalstr,charhistory,999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		if (status != 0)
		{
			message = "Cannot write some keyword in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	return EPOK;
}
/*xxxx end of SECTION A2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A3 ************************************************************
* filderLibrary: This function calculates the (low-pass filtered and) derivative of the models ('pulse_templates') of the library (only necessary
*                if first record), and it stores the 'pulse_templates_filder' and the 'maxDERs' and 'samp1DERs' in the 'reconstruct_init' structure.
*                The maximum of the (low-pass filtered and) differentiated pulse has to be compared to the 'maxDERs' to select the appropriate model.
                 Or, the 1st sample out of the differentiated pulse has to be compared to the 'samp1DERs' to select the appropriate model.
* 
*
* - Check if it is the first record
* - (Low-pass filtered and) differentiate the models ('pulse_templates') of the library
* - Store the (low-pass filtered) derivatives in 'pulse_templates_filder'
* - Calculate the maximum of the (low-pass filtered and) differentiated models ('maxDERs')
* - Locate the 1st sample of the (low-pass filtered and) differentiated models ('samp1DERs')
* 
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
*                     In particular, this function uses some parameters to filter ('scaleFactor') and others to handle the pulse templates
* 		      and their derivatives ('ntemplates', 'pulse_templates', 'pulse_templates_filder' and 'maxDERs')
* - samprate: Sampling rate
******************************************************************************/
int filderLibrary(ReconstructInitSIRENA** reconstruct_init, double samprate)
{
	int status = EPOK;
	string message = "";
	char valERROR[256];

        scheduler* sc = scheduler::get();

	if ((*reconstruct_init)->library_collection->pulse_templates_filder[0].template_duration == -1)		// First record
	{
		double scaleFactor = (*reconstruct_init)->scaleFactor;

		// Check boxLength
		double cutFreq = 2 * (1/(2*pi*scaleFactor));
		int boxLength = (int) ((1/cutFreq) * samprate);
		if ((boxLength <= 1) && (*reconstruct_init)->tstartPulse1 == 0)
		{
			message = "lpf_boxcar(Model): scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.";
			EP_PRINT_ERROR(message,-999);	// Only a warning
		}

		// 'pulse_templates' are filtered and differentiated
		gsl_vector *model;
		if ((model = gsl_vector_alloc((*reconstruct_init)->library_collection->pulse_templates[0].template_duration)) == 0)
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);
			message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		for (int i=0; i<(*reconstruct_init)->library_collection->ntemplates; i++)
		{
			gsl_vector_memcpy(model, (*reconstruct_init)->library_collection->pulse_templates[i].ptemplate);

			// PULSE TEMPLATE: Low-pass filtering
			status = lpf_boxcar(&model,model->size,scaleFactor,samprate);
			if (status == 1)
			{
				message = "Cannot run routine lpf_boxcar for low-pass filtering";
				EP_PRINT_ERROR(message,status); return(EPFAIL);
			}
			if (status == 3)
			{
				status = EPOK;
			}
			if (status == 4)
			{
				message = "lpf_boxcar: scaleFactor too high => Cut-off frequency too low";
				EP_PRINT_ERROR(message,status); return(EPFAIL);
			}

			// PULSE TEMPLATE: Derivative after filtering
			if (differentiate (&model,model->size))
			{
				message = "Cannot run routine differentiate in filderLibrary";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
                    

                        if(!sc->is_threading()){

                          // Store the low-pass filtered derivatives in 'pulse_templates_filder'
                          gsl_vector_memcpy((*reconstruct_init)->library_collection->pulse_templates_filder[i].ptemplate,model);

                          (*reconstruct_init)->library_collection->pulse_templates_filder[i].template_duration = (*reconstruct_init)->library_collection->pulse_templates[i].template_duration;

                          // Calculate the maximum of the low-pass filtered and differentiated models
                          gsl_vector_set((*reconstruct_init)->library_collection->maxDERs,i,gsl_vector_max(model));

                          // Locate the 1st sample of the (low-pass filtered and) differentiated models ('samp1DERs')
                          gsl_vector_set((*reconstruct_init)->library_collection->samp1DERs,i,gsl_vector_get(model,0));
                        }else{
                          // Since all threads share the same library
                          // we need to lock here to be able to write
                          std::lock_guard<std::mutex> lk(library_mut);
                                                 // Store the low-pass filtered derivatives in 'pulse_templates_filder'
                          gsl_vector_memcpy((*reconstruct_init)->library_collection->pulse_templates_filder[i].ptemplate,model);

                          (*reconstruct_init)->library_collection->pulse_templates_filder[i].template_duration = (*reconstruct_init)->library_collection->pulse_templates[i].template_duration;

                          // Calculate the maximum of the low-pass filtered and differentiated models
                          gsl_vector_set((*reconstruct_init)->library_collection->maxDERs,i,gsl_vector_max(model));

                          // Locate the 1st sample of the (low-pass filtered and) differentiated models ('samp1DERs')
                          gsl_vector_set((*reconstruct_init)->library_collection->samp1DERs,i,gsl_vector_get(model,0));
                        }
			//cout<<i<<" "<<gsl_vector_get((*reconstruct_init)->library_collection->maxDERs,i)<<" "<<gsl_vector_get((*reconstruct_init)->library_collection->samp1DERs,i)<<" "<<gsl_vector_get((*reconstruct_init)->library_collection->energies,i)<<endl;
			//cout<<gsl_vector_get(model,0)<<" "<<gsl_vector_get(model,1)<<" "<<gsl_vector_get(model,2)<<endl;
                        //cout<<gsl_vector_get(model,0)<<" "<<gsl_vector_get(model,1)<<" "<<gsl_vector_get(model,2)<<" "<<gsl_vector_get(model,3)<<" "<<gsl_vector_get(model,4)<<" "<<gsl_vector_get(model,5)<<" "<<gsl_vector_get(model,6)<<" "<<gsl_vector_get(model,7)<<" "<<gsl_vector_get(model,8)<<" "<<gsl_vector_get(model,9)<<" "<<gsl_vector_get(model,10)<<" "<<endl;
		}

		gsl_vector_free(model); model = 0;
	}

	return(EPOK);
}
/*xxxx end of SECTION A3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A4 ************************************************************
* loadRecord: This function loads the struture 'record' in the 'adc_double' GSL vector.
*
* It checks if the record has been filled out with 0's => It only loads the first values (which are different from 0).
* 
* Parameters:
* - record: Member of 'TesRecord' structure that contains the input record
* - time_record: Starting time of the record (output)
* - adc_double: Storage of the record to be processed (input/output) 
******************************************************************************/
int loadRecord(TesRecord* record, double *time_record, gsl_vector **adc_double)
{ 
	string message = "";
	char valERROR[256];
  
	*time_record = record->time;
	for (int i=0;i<record->trigger_size;i++)
	{
		gsl_vector_set(*adc_double,i,record->adc_double[i]);
	}

	// Just in case the last record has been filled out with 0's => Re-allocate 'invector'
	if ((gsl_vector_get(*adc_double,record->trigger_size-1) == 0.0) && (gsl_vector_get(*adc_double,record->trigger_size-2) == 0.0))
	// Baseline could be near to zero, even having negative values of the record
	{
		// Know the new dimension of the last record (elements different from 0)
		long eventszLastRecord = record->trigger_size;
		for (int i=0;i<record->trigger_size;i++)
		{
			if (gsl_vector_get(*adc_double,record->trigger_size-1-i) == 0.0)	eventszLastRecord--;
			else 									break;
		}
		if (eventszLastRecord == 0)
		{
			EP_PRINT_ERROR("Last record size is 0",EPFAIL); return(EPFAIL);
		}

		// Auxiliary vector
		// It is not necessary to check the allocation because 'eventszLastRecord' has been check previously
		gsl_vector *vector_aux = gsl_vector_alloc(eventszLastRecord);
		gsl_vector_view temp;
		if ((eventszLastRecord < 1) || (eventszLastRecord > (*adc_double)->size))
		{
			sprintf(valERROR,"%d",__LINE__+5);
			string str(valERROR);
			message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		temp = gsl_vector_subvector(*adc_double,0,eventszLastRecord);
		gsl_vector_memcpy(vector_aux,&temp.vector);

		// Free invector, allocate with the new dimension and free the auxiliary vector
		gsl_vector_free(*adc_double); 
		// It is not necessary to check the allocation because 'eventszLastRecord' has been check previously
		*adc_double = gsl_vector_alloc(eventszLastRecord);
		gsl_vector_memcpy(*adc_double,vector_aux);
		gsl_vector_free(vector_aux); vector_aux = 0;
	}

	return (EPOK);
}
/*xxxx end of SECTION A4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A5 ************************************************************
* procRecord function:  This function processes the input record.
*
* - Declare and initialize variables
* - Allocate GSL vectors
* - (Low-pass filtering and) differentiation
* - Find the events (pulses) in the record
* - Calculate the tend of the found pulses and check if the pulse is saturated
* - Obtain the approximate rise and fall times of each pulse (to be done)
* - Load the found pulses data in the input/output 'foundPulses' structure
* - Write test info (if 'reconstruct_init->intermediate'=1)
* - Write pulses info in intermediate output FITS file (if 'reconstruct_init->intermediate'=1)
* - Free allocated GSL vectors
*
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
*                     This function uses parameters to filter ('scaleFactor'),
*                     to find pulses ('pulse_length', 'samplesUp', 'nSgms', 'LrsT', 'LbT', 'maxPulsesPerRecord')
*                     and to wrtite info ('detectFile') if 'reconstruct_init->intermediate'=1
* - tstartRecord: Starting time of the record (in order to calculate absolute times)
* - samprate: Sampling rate (in order to low-pass filter)
* - dtcObject: Object which contains information of the intermediate FITS file (to be written if 'intermediate'=1)
* - record: GSL vector with signal values of input record
* - foundPulses: Input/output structure where the found pulses info is stored 
****************************************************************************/
int procRecord(ReconstructInitSIRENA** reconstruct_init, double tstartRecord, double samprate, fitsfile *dtcObject, gsl_vector *record, PulsesCollection *foundPulses)
{
	int status = EPOK;
	string message = "";
	char valERROR[256];
    
	// Declare and initialize variables
	int numPulses = 0;
	double threshold = 0.0;
	double sigmaout = 0.0;

	double stopCriteriaMKC = 1.0;	// Used in medianKappaClipping
	                               	// Given in %
	double kappaMKC = 3.0;		// Used in medianKappaClipping
	double levelPrvPulse = 100.0;  	// Secondary pulses must be 1/levelPrvPulse times larger than the preceding pulse

	gsl_vector_view temp;

	double scaleFactor = (*reconstruct_init)->scaleFactor;
	int sizePulse_b = (*reconstruct_init)->pulse_length;
	double samplesUp = (*reconstruct_init)->samplesUp;
	double nSgms = (*reconstruct_init)->nSgms;
	double Lrs = (int) ((*reconstruct_init)->LrsT*samprate);	// Running sum length (in the RS filter case): 'LrsT' in samples
	double Lb = (int) ((*reconstruct_init)->LbT*samprate); 		// Baseline averaging length (in the RS filter case): 'LbT' in samples

	// Allocate GSL vectors
	// It is not necessary to check the allocation because 'record' size must already be > 0
	gsl_vector *recordNOTFILTERED = gsl_vector_alloc(record->size); // Record without having been filtered
	gsl_vector *recordDERIVATIVE = gsl_vector_alloc(record->size);  // Derivative of 'invectorFILTERED'

	// To look for pulses
	// It is not necessary to check the allocation because '(*reconstruct_init)->maxPulsesPerRecord'='EventListSize'(input parameter) must already be > 0
	gsl_vector *tstartgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *tendgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *qualitygsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *pulseHeightsgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *maxDERgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *samp1DERgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
        gsl_vector *lagsgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector_set_zero(qualitygsl);
	gsl_vector_set_zero(pulseHeightsgsl);		// In order to choose the proper pulse model to calculate
	                                                // the adjusted derivative and to fill in the ESTENRGY column
	                                                // in the output FITS file
	gsl_vector_set_zero(maxDERgsl);
	gsl_vector_set_zero(samp1DERgsl);

	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';
	
	gsl_vector_memcpy(recordNOTFILTERED,record);
	
	// Filtering by using wavelets
	/*int onlyOnce = 0;
	//if (onlyOnce == 0)
	//{
		int length = pow(2,floor(log2(record->size)));
		gsl_vector *recordFilteredByWavelet = gsl_vector_alloc(length);
		temp = gsl_vector_subvector(record,0,length);
		gsl_vector_memcpy(recordFilteredByWavelet,&temp.vector);
		status = filterByWavelets(*reconstruct_init,&recordFilteredByWavelet,length, &onlyOnce);
		if (status == EPFAIL)
		{
			message = "Cannot run routine filterByWavelets";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		cout<<"recordFilteredByWavelet->size: "<<recordFilteredByWavelet->size<<endl;
		for (int i=0;i<recordFilteredByWavelet->size-1;i++)
		{
			//cout<<i<<" "<<gsl_vector_get(record,i)<<" "<<gsl_vector_get(recordFilteredByWavelet,i)<<endl;
			gsl_vector_set(record,i,gsl_vector_get(recordFilteredByWavelet,i));
		}
	
		gsl_vector_free(recordFilteredByWavelet);
	//}*/
	
	// Check boxLength
	double cutFreq = 2 * (1/(2*pi*scaleFactor));
	int boxLength = (int) ((1/cutFreq) * samprate);
	if (boxLength <= 1)
	{
		message = "lpf_boxcar: scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.";
		//?? Too many messages
	}

	// Low-pass filtering
	status = lpf_boxcar(&record,record->size,scaleFactor,samprate);
	if (status == EPFAIL)
	{
		message = "Cannot run routine lpf_boxcar for low pass filtering";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (status == 3)
	{
		status = EPOK;
	}
	if (status == 4)
	{
		message = "lpf_boxcar: scaleFactor too high => Cut-off frequency too low";
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	
	// Differentiate after filtering
	if (differentiate (&record, record->size))
	{
		message = "Cannot run routine differentiate for differentiating after filtering";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	gsl_vector_memcpy(recordDERIVATIVE,record);
        
        // Smooth derivative
        //cout<<"Pulso"<<endl;
        //for (int i=0;i<record->size;i++)  cout<<i<<" "<<gsl_vector_get(record,i)<<endl;
	/*if (smoothDerivative (&record, 4))
	{
		message = "Cannot run routine smoothDerivative";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	for (int i=0;i<4-1;i++)    gsl_vector_set(record,i,0.0);
	gsl_vector_memcpy(recordDERIVATIVE,record);*/
        //cout<<"PulsoDERIVATIVE"<<endl;
        //for (int i=0;i<recordDERIVATIVE->size;i++)  cout<<i<<" "<<gsl_vector_get(recordDERIVATIVE,i)<<endl;
        
	// // It is not necessary to check the allocation because the allocation of 'recordDERIVATIVE' has been checked previously
	gsl_vector *recordDERIVATIVEOriginal = gsl_vector_alloc(recordDERIVATIVE->size);
	gsl_vector_memcpy(recordDERIVATIVEOriginal,recordDERIVATIVE);
        
	// Find the events (pulses) in the record
	if ((*reconstruct_init)->opmode == 1)	// In PRODUCTION mode
	{
		int tstartFirstEvent = 0;
		bool triggerCondition; 
		int flagTruncated;
		int tstartProvided;
		
                if ((*reconstruct_init)->tstartPulse1 != 0)
                {
                        if (noDetect(recordDERIVATIVE, (*reconstruct_init), &numPulses, &tstartgsl, &qualitygsl, &maxDERgsl, &samp1DERgsl))
                        {
                                message = "Cannot run routine noDetect";
                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                        }
                }
                else 
                {
                        if (InitialTriggering (recordDERIVATIVE, nSgms,
                                scaleFactor, samprate, stopCriteriaMKC, kappaMKC,
                                &triggerCondition, &tstartFirstEvent, &flagTruncated,
                                //&threshold,&sigmaout, (*reconstruct_init)->tstartPulse1))
                                &threshold,(*reconstruct_init)->tstartPulse1))
                        {
                                message = "Cannot run routine InitialTriggering";
                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                        }
                        
                        if (strcmp((*reconstruct_init)->detectionMode,"STC") == 0)
                        {
                                if (FindSecondariesSTC ((*reconstruct_init)->maxPulsesPerRecord, recordDERIVATIVE, threshold, (*reconstruct_init), tstartFirstEvent, &numPulses, &tstartgsl, &qualitygsl, &maxDERgsl, &samp1DERgsl))
                                {
                                        message = "Cannot run FindSecondariesSTC";
                                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                        }
                        else if (strcmp((*reconstruct_init)->detectionMode,"AD") == 0)
                        {   
                                if (FindSecondaries ((*reconstruct_init)->maxPulsesPerRecord,
                                        //recordDERIVATIVE, threshold,sigmaout,
                                        recordDERIVATIVE, threshold, samprate,
                                        (*reconstruct_init),
                                        tstartFirstEvent,
                                        &numPulses,&tstartgsl,&qualitygsl, &maxDERgsl,&samp1DERgsl,&lagsgsl))
                                {
                                        message = "Cannot run routine FindSecondaries";
                                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                        }
                }
	}
	else if ((*reconstruct_init)->opmode == 0)	// In CALIBRATION mode
	{
		if (findPulsesCAL (recordNOTFILTERED, recordDERIVATIVE, &tstartgsl, &qualitygsl, &pulseHeightsgsl, &maxDERgsl,
			&numPulses, &threshold, 
		        scaleFactor, samprate,
			samplesUp, nSgms,
			Lb, Lrs,
			(*reconstruct_init),
			stopCriteriaMKC,
			kappaMKC))
		{
			message = "Cannot run routine findPulsesCAL";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
	(*reconstruct_init)->threshold = threshold;
            
	// Write test info
	if ((*reconstruct_init)->intermediate == 1)
	{
		if (writeTestInfo((*reconstruct_init), recordDERIVATIVEOriginal, threshold, dtcObject))
		{
			message = "Cannot run routine writeTestInfo";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
	gsl_vector_free(recordDERIVATIVEOriginal); recordDERIVATIVEOriginal = 0;

	// Calculate the tend of the found pulses and check if the pulse is saturated
	for (int i=0;i<numPulses;i++)
	{
		gsl_vector_set(tendgsl,i,gsl_vector_get(tstartgsl,i)+sizePulse_b);	//tend_i = tstart_i + Pulse_Length

		if (gsl_vector_get(tendgsl,i) >= recordDERIVATIVE->size)		// Truncated pulses at the end of the record
		{
			gsl_vector_set(tendgsl,i,(recordDERIVATIVE->size)-1);
			gsl_vector_set (qualitygsl,i,2);
		}

		if ((numPulses != 1) && (i != numPulses-1)) 				// More than one pulse in the record and not the last one
		{
			if (gsl_vector_get(tendgsl,i) > gsl_vector_get(tstartgsl,i+1))
			{
				gsl_vector_set(tendgsl,i,gsl_vector_get(tstartgsl,i+1));
			}
		}
	}
	
	// Obtain the approximate rise and fall times of each pulse
	// It is not necessary to check the allocation because '(*reconstruct_init)->maxPulsesPerRecord'='EventListSize'(input parameter) must already be > 0
	gsl_vector *tauRisegsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *tauFallgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector_set_zero(tauRisegsl);
	gsl_vector_set_zero(tauFallgsl);
	/*if (obtainTau (recordNOTFILTERED, tstartgsl, tendgsl, *numPulses, &tauRisegsl, &tauFallgsl))
	{
	    message = "Cannot run routine obtainTau to calculate pulses slopes";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}*/

	// Load the found pulses data in the input/output 'foundPulses' structure
	foundPulses->ndetpulses = numPulses;
	foundPulses->pulses_detected = new PulseDetected[numPulses];
	//cout<<"numPulses: "<<numPulses<<endl;
	for (int i=0;i<numPulses;i++)
	{
		//foundPulses->pulses_detected[i].pulse_duration = ceil(gsl_vector_get(tendgsl,i)-gsl_vector_get(tstartgsl,i));
                foundPulses->pulses_detected[i].pulse_duration = floor(gsl_vector_get(tendgsl,i)-gsl_vector_get(tstartgsl,i));
                
                foundPulses->pulses_detected[i].avg_4samplesDerivative = gsl_vector_get(samp1DERgsl,i);
                foundPulses->pulses_detected[i].E_lowres = -999;
                
		// 'grade1' will be known after running 'runEnergy' (but initialize for library creation!)
		foundPulses->pulses_detected[i].grade1 = 0;
                foundPulses->pulses_detected[i].grading = -999;
		if (i!= 0)
		{
                        foundPulses->pulses_detected[i].grade2 = (int)(gsl_vector_get(tstartgsl,i)-gsl_vector_get(tstartgsl,i-1));
		}
		else
		{
			foundPulses->pulses_detected[i].grade2 = (*reconstruct_init)->pulse_length;
		}
		
		// 'phi' and 'lagsShift' will be known after running 'runEnergy' (but initialize for library creation!)
		foundPulses->pulses_detected[i].phi = -999;
                foundPulses->pulses_detected[i].lagsShift = -999;
                
                if (foundPulses->pulses_detected[i].pulse_duration == 0 )
                {
                        foundPulses->pulses_detected[i].pulse_adc = gsl_vector_alloc(1);
                        gsl_vector_set(foundPulses->pulses_detected[i].pulse_adc,0,-999);
                        foundPulses->pulses_detected[i].pulse_duration = 1;
                }
                else        
                {
                        if ((foundPulses->pulses_detected[i].pulse_adc = gsl_vector_alloc(foundPulses->pulses_detected[i].pulse_duration)) == 0)
                        {
                                sprintf(valERROR,"%d",__LINE__-2);
                                string str(valERROR);
                                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                        }
                        temp = gsl_vector_subvector(recordNOTFILTERED,gsl_vector_get(tstartgsl,i),foundPulses->pulses_detected[i].pulse_duration);
                        if (gsl_vector_memcpy(foundPulses->pulses_detected[i].pulse_adc,&temp.vector) != 0)
                        {
                                sprintf(valERROR,"%d",__LINE__-2);
                                string str(valERROR);
                                message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                        }

                }
    
		foundPulses->pulses_detected[i].Tstart = gsl_vector_get(tstartgsl,i)/samprate+tstartRecord;
                foundPulses->pulses_detected[i].TstartSamples = gsl_vector_get(tstartgsl,i);
		foundPulses->pulses_detected[i].Tend = gsl_vector_get(tendgsl,i)/samprate+tstartRecord;
		foundPulses->pulses_detected[i].riseTime = gsl_vector_get(tauRisegsl,i);
		foundPulses->pulses_detected[i].fallTime = gsl_vector_get(tauFallgsl,i);
		foundPulses->pulses_detected[i].pulse_height = gsl_vector_get(pulseHeightsgsl,i);
		foundPulses->pulses_detected[i].maxDER = gsl_vector_get(maxDERgsl,i);
		foundPulses->pulses_detected[i].samp1DER = gsl_vector_get(samp1DERgsl,i);
		// 'energy' will be known after running 'runEnergy'
		foundPulses->pulses_detected[i].quality = gsl_vector_get(qualitygsl,i);
                foundPulses->pulses_detected[i].numLagsUsed = gsl_vector_get(lagsgsl,i);
		//cout<<"Pulse "<<i<<" tstart="<<gsl_vector_get(tstartgsl,i)<<", maxDER= "<<foundPulses->pulses_detected[i].maxDER<<" , samp1DER="<<gsl_vector_get(samp1DERgsl,i)<<", pulse_duration= "<<foundPulses->pulses_detected[i].pulse_duration<<",quality= "<<foundPulses->pulses_detected[i].quality<<" ,lags="<<gsl_vector_get(lagsgsl,i)<<" , Tstart="<<foundPulses->pulses_detected[i].Tstart<<endl;
                //log_debug("Pulse %d", i," tstart=%f",gsl_vector_get(tstartgsl,i), " maxDER=%f",foundPulses->pulses_detected[i].maxDER, " samp1DER=%f",gsl_vector_get(samp1DERgsl,i), " //pulse_duration=%d",foundPulses->pulses_detected[i].pulse_duration," quality=%d",foundPulses->pulses_detected[i].quality," lags=%f",gsl_vector_get(lagsgsl,i));
                //log_debug("Pulse %i tstart=%f maxDER=%f samp1DER=%f pulse_duration=%i quality=%f lags=%f",i,gsl_vector_get(tstartgsl,i),foundPulses->pulses_detected[i].maxDER,gsl_vector_get(samp1DERgsl,i),foundPulses->pulses_detected[i].pulse_duration,foundPulses->pulses_detected[i].quality,gsl_vector_get(lagsgsl,i));
                
	}

	// Write pulses info in intermediate output FITS file
	if ((*reconstruct_init)->intermediate == 1)
	{
		if (writePulses (reconstruct_init, samprate, tstartRecord, recordNOTFILTERED, numPulses, tstartgsl, tendgsl, qualitygsl, tauRisegsl, tauFallgsl, dtcObject))
		{
			message = "Cannot run routine writePulses to write pulses in intermediate output file";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}

	// Free allocated GSL vectors
	gsl_vector_free(recordNOTFILTERED); recordNOTFILTERED = 0;
	gsl_vector_free(recordDERIVATIVE); recordDERIVATIVE = 0;

	gsl_vector_free(tstartgsl); tstartgsl = 0;
	gsl_vector_free(tendgsl); tendgsl = 0;
	gsl_vector_free(qualitygsl); qualitygsl = 0;
	gsl_vector_free(pulseHeightsgsl); pulseHeightsgsl = 0;
	gsl_vector_free(maxDERgsl); maxDERgsl = 0;
	gsl_vector_free(samp1DERgsl); samp1DERgsl = 0;
        gsl_vector_free(lagsgsl); lagsgsl = 0;

	gsl_vector_free(tauRisegsl); tauRisegsl = 0;
	gsl_vector_free(tauFallgsl); tauFallgsl = 0;

	return EPOK;
}
/*xxxx end of SECTION A5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A6 ************************************************************
* writePulses function: This function writes the data of the pulses found in the record in the intermediate FITS file (in the PULSES HDU).
*                       The pulses info given is: TSTART, I0 (the pulse itself), TEND, TAURISE, TAUFALL and QUALITY.
*
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
*                     This function uses 'detectFile', 'pulse_length' and 'clobber'
* - samprate: Sampling rate (to convert samples to seconds)
* - initialtime: Starting time of the record (in order to calculate absolute times)
* - invectorNOTFIL: GSL vector with the original record (neither low-pass filtered nor differentiated)
* - numPulsesRecord: Number of pulses found in the record
* - tstart: GSL vector with the start times of the found pulses
* - tend: GSL vector with the end times of the found pulses
* - quality:  GSL vector with the quality of the found pulses
* - taurise: GSL vector with the rise time constants of the found pulses (to be done)
* - taufall: GSL vector with the fall time constants of the found pulses (to be done)
* - dtcObject: Object which contains information of the intermediate FITS file 
******************************************************************************/
int writePulses(ReconstructInitSIRENA** reconstruct_init, double samprate, double initialtime, gsl_vector *invectorNOTFIL, int numPulsesRecord, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, fitsfile *dtcObject)
{
	int status = EPOK;
	string message = "";
	char valERROR[256];

	// Declare variables
	int t0;		// First value of index of pulse
	gsl_matrix *vgslout2;

	gsl_vector_view temp;

	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';

	// If intermediate=1 => First record, createDetectFile
	//	                Change clobber to 2
	//                   => Not first record, append info to output intermediate file (because clobber is 2)
	long totalpulses;	// It is necessary to know the row where the info is going to be written
	if ((*reconstruct_init)->clobber == 1)
	{
		totalpulses = 0;
		if ((*reconstruct_init)->opmode == 0)	(*reconstruct_init)->clobber = 2;
	}
	else if ((*reconstruct_init)->clobber == 2)
	{
		if (fits_get_num_rows(dtcObject,&totalpulses, &status))
		{
			message = "Cannot get number of rows in " + string(dtcName);
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
	}

	if (numPulsesRecord != 0)
	{
		// It is not necessary to check the allocation because 'numPulsesRecord' and '(*reconstruct_init)->pulse_length' have been checked previously
		vgslout2 = gsl_matrix_alloc(numPulsesRecord,(*reconstruct_init)->pulse_length);

		// Converting samples to seconds
		for (int i=0; i<numPulsesRecord; i++)
		{
			gsl_vector_set(tstart,i,gsl_vector_get(tstart,i)+1);	// To be consistent between the GSL indexes which start from 0 and the real samples starting from 1
			t0 = gsl_vector_get (tstart,i);
			gsl_vector_set(tstart,i,initialtime + (gsl_vector_get (tstart,i) * (1/samprate)));
			gsl_vector_set(tend,i,initialtime + (gsl_vector_get (tend,i) * (1/samprate)));

			if (invectorNOTFIL->size - t0 > (*reconstruct_init)->pulse_length)	//The invectorNOTFIL has more sampless than sizePulse
			{
				if ((t0 < 0) || (t0 > invectorNOTFIL->size-2)
					|| ((*reconstruct_init)->pulse_length < 1) || ((*reconstruct_init)->pulse_length > invectorNOTFIL->size-t0))
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				temp = gsl_vector_subvector(invectorNOTFIL,t0,(*reconstruct_init)->pulse_length);
				gsl_matrix_set_row(vgslout2, i, &temp.vector);
			}
			else 									// The invectorNOTFIL has less bins than sizePulse (truncated)
			{
				for (int j=0; j<(invectorNOTFIL->size) - t0; j++)
				{
					if (t0 == -1) t0 = 0;
					gsl_matrix_set (vgslout2,i,j,gsl_vector_get(invectorNOTFIL,j+t0));
				}

				for (int j=(invectorNOTFIL->size)-t0; j< (*reconstruct_init)->pulse_length; j++) {gsl_matrix_set (vgslout2,i,j,0.0);}
			}
		}

		IOData obj;

		// PULSES HDU
		// Creating TSTART Column
		obj.inObject = dtcObject;
		obj.nameTable = new char [255];
		strcpy(obj.nameTable,"PULSES");
		obj.iniRow = totalpulses+1;
		obj.endRow = totalpulses+numPulsesRecord-1+1;
		obj.iniCol = 0;
		obj.nameCol = new char [255];
		strcpy(obj.nameCol,"TSTART");
		obj.type = TDOUBLE;
		obj.unit = new char [255];
		strcpy(obj.unit,"seconds");

		temp = gsl_vector_subvector(tstart,0,numPulsesRecord);
		if (writeFitsSimple(obj, &temp.vector))
		{
			message = "Cannot run routine writeFitsSimple for column TSTART";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating I0 Column
		strcpy(obj.nameCol,"I0");
		obj.type = TDOUBLE;
		strcpy(obj.unit,"ADC");
		if (writeFitsComplex(obj, vgslout2))
		{
			message = "Cannot run routine writeFitsComplex for column IO";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating TEND Column
		strcpy(obj.nameCol,"TEND");
		obj.type = TDOUBLE;
		strcpy(obj.unit,"seconds");
		temp = gsl_vector_subvector(tend,0,numPulsesRecord);
		if (writeFitsSimple(obj, &temp.vector))
		{
			message = "Cannot run routine writeFitsSimple for column TEND";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating TAURISE Column
		strcpy(obj.nameCol,"TAURISE");
		temp = gsl_vector_subvector(taurise,0,numPulsesRecord);
		if (writeFitsSimple(obj, &temp.vector))
		{
			message = "Cannot run routine writeFitsSimple for column TAURISE";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating TAUFALL Column
		strcpy(obj.nameCol,"TAUFALL");
		temp = gsl_vector_subvector(taufall,0,numPulsesRecord);
		if (writeFitsSimple(obj, &temp.vector))
		{
			message = "Cannot run routine writeFitsSimple for column TAUFALL";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}

		// Creating QUALITY Column
		strcpy(obj.nameCol,"QUALITY");
		obj.type = TSHORT;
		strcpy(obj.unit," ");
		temp = gsl_vector_subvector(quality,0,numPulsesRecord);
		if (writeFitsSimple(obj, &temp.vector))
		{
			message = "Cannot run routine writeFitsSimple for column QUALITY";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		totalpulses = totalpulses + numPulsesRecord;

		// Free allocated GSL vectors
		gsl_matrix_free (vgslout2); vgslout2 = 0;

		delete [] obj.nameTable; obj.nameTable = 0;
		delete [] obj.nameCol; obj.nameCol = 0;
		delete [] obj.unit; obj.unit = 0;
	}

	return (EPOK);
}
/*xxxx end of SECTION A6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A7 ************************************************************
* writeTestInfo function: This function writes the TESTINFO HDU in the intermediate FITS file.
*                         The written columns are FILDER (low-pass filtered and differentiated record) and THRESHOLD.
* 
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
*                     in particular, this function uses 'detectFile'
* - recordDERIVATIVE: GSL vector with input record (low-pass filtered and) differentiated
* - threshold: Threshold used to find pulses
* - dtcObject: Object which contains information of the intermediate FITS file
******************************************************************************/
int writeTestInfo(ReconstructInitSIRENA* reconstruct_init, gsl_vector *recordDERIVATIVE, double threshold, fitsfile *dtcObject)
{
	int status = EPOK;
	string message = "";

	long totalrecords;

	char dtcName[256];
	strncpy(dtcName,reconstruct_init->detectFile,255);
	dtcName[255]='\0'; // enforce zero ending string in case of buffer overflows

	// To work with tables (HDUs)
	char extname[20];

	strcpy(extname,"TESTINFO");
	if (fits_movnam_hdu(dtcObject, ANY_HDU,extname, 0, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in output intermediate file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (fits_get_num_rows(dtcObject,&totalrecords, &status))
	{
		message = "Cannot get number of rows in " + string(dtcName) + " (TESTINFO HDU)";
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	IOData obj;

	// Creating FILDER Column
	obj.inObject = dtcObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"TESTINFO");
	obj.iniRow = totalrecords+1;
	obj.endRow = totalrecords+1;
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit," ");

	// It is not necessary to check the allocation because 'recordDERIVATIVE' size must already be > 0
	gsl_matrix *matrixToWrite = gsl_matrix_alloc(1,recordDERIVATIVE->size);
	strcpy(obj.nameCol,"FILDER");
	gsl_matrix_set_row(matrixToWrite,0,recordDERIVATIVE);
	if (writeFitsComplex(obj, matrixToWrite))
	{
		message = "Cannot run routine writeFitsComplex for recordDERIVATIVE (TESTINFO HDU)";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating THRESHOLD Column
	strcpy(obj.nameCol,"Threshold");
	gsl_vector *vectorToWrite = gsl_vector_alloc(1);
	gsl_vector_set(vectorToWrite,0,threshold);
	if (writeFitsSimple(obj, vectorToWrite))
	{
		message = "Cannot run routine writeFitsSimple for threshold (TESTINFO HDU)";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Free allocated GSL vectors amd matrices
	gsl_matrix_free(matrixToWrite); matrixToWrite = 0;
	gsl_vector_free(vectorToWrite); vectorToWrite = 0;

	delete [] obj.nameTable; obj.nameTable = 0;
	delete [] obj.nameCol; obj.nameCol = 0;
        //cout<<"tend1: "<<gsl_vector_get(tendgsl,i)<<endl;
	delete [] obj.unit; obj.unit = 0;

	return (EPOK);
}
/*xxxx end of SECTION A7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION  A8 ************************************************************
* calculateTemplate function: This function calculates the template (PULSE column in the library) of non piled-up pulses.
*                             Just in case in the detection process some piled-up pulses have not been distinguished as different pulses, a pulseheights histogram is built.
*                             It uses the pulseheights histogram (built by using the PHEIGGHT column of the library), 'Tstart' and 'quality' to select the non piled-up pulses.
*
* - Declare and initialize variables
* - Before building the histogram, select the pulseheihts of the pulseswell separated from other pulses and with 'quality'=0
* - Create the pulseheights histogram
* - Calculate the pulseaverage only taking into account the valid pulses
* 	- Check if the pulse is piled-up or not
* 	- Non piled-up pulses => Average them
* - Calculate covariance and weight matrices
* - Free allocated GSL vectors
*
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
*                     In particular, this function uses 'pulse_length' and 'EnergyMethod'
* - pulsesAll: Collection of pulses found in the previous records
* - pulsesInRecord:  Collection of pulses found in the current record
* - samprate: Sampling rate
* - pulseaverage: GSL vector with the pulseaverage (template) of the non piled-up pulses
* - pulseaverageHeight: Height value of the pulseaverage
* - covariance: GSL matrix with covariance matrix
* - weight: GSL matrix with weight matrix (inverse of covariance matrix)
* - inputPulseLength: Pulse length
* - pulseaverageMaxLengthFixedFilter: GSL vector with the pulseaverage (template) whose length is largeFilter of the non piled-up pulses
******************************************************************************/
int calculateTemplate(ReconstructInitSIRENA *reconstruct_init, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, double samprate, gsl_vector **pulseaverage, double *pulseaverageHeight, gsl_matrix **covariance, gsl_matrix **weight, int inputPulseLength, gsl_vector **pulseaverageMaxLengthFixedFilter)
{
	string message = "";
	char valERROR[256];

	// Declare and initialize variables
	int totalPulses = pulsesAll->ndetpulses + pulsesInRecord->ndetpulses;
	// It is not necessary because 'totalPulses'='pulsesAll->ndetpulses + pulsesInRecord->ndetpulses' has been checked previously
	gsl_vector *tstart = gsl_vector_alloc(totalPulses);
	gsl_vector *pulseheight = gsl_vector_alloc(totalPulses);
	gsl_vector *quality = gsl_vector_alloc(totalPulses);
	for (int i=0;i<pulsesAll->ndetpulses;i++)
	{
		gsl_vector_set(tstart,i,pulsesAll->pulses_detected[i].Tstart);
		gsl_vector_set(pulseheight,i,pulsesAll->pulses_detected[i].pulse_height);
		gsl_vector_set(quality,i,pulsesAll->pulses_detected[i].quality);
	}
	for (int i=0;i<pulsesInRecord->ndetpulses;i++)
	{
		gsl_vector_set(tstart,i+pulsesAll->ndetpulses,pulsesInRecord->pulses_detected[i].Tstart);
		gsl_vector_set(pulseheight,i+pulsesAll->ndetpulses,pulsesInRecord->pulses_detected[i].pulse_height);
		gsl_vector_set(quality,i+pulsesAll->ndetpulses,pulsesInRecord->pulses_detected[i].quality);
	}

	// It is not necessary because 'pulsesAll->ndetpulses + pulsesInRecord->ndetpulses' has been checked previously
	gsl_vector *nonpileup = gsl_vector_alloc(totalPulses);	// Piled-up pulse => Not taken into account to calculate the template
	long nonpileupPulses = totalPulses;			// A priori, all the found pulses are considered as non piled-up
	gsl_vector_set_all(nonpileup,1);

	int nBins;						// Square-root choice (used by Excel and many others)
	gsl_vector *xhisto;					// X-axis of the pulseheights histogram
	gsl_vector *yhisto;					// Y-axis of the pulseheights histogram
	int index_maximumpulseheight;				// Index where the maximum of the pulseheights histogram is
	double maximumpulseheight;				// Maximum of the pulseheights histogram

	bool firstnonpileupPulse = true;
	// It is not necessary because 'reconstruct_init->pulse_length'='PulseLength' (input parameter) has been checked previously
	gsl_vector *pulse = gsl_vector_alloc(reconstruct_init->pulse_length);

	double tstartnext;

	gsl_vector_view temp;

	gsl_vector_scale(tstart,samprate); 			//tstarts not in seconds but in samples

	// Before building the histogram, select the pulseheihts of the pulses which are enough separated from others and whose 'quality' is 0
	// It is not necessary because 'totalPulses'='pulsesAll->ndetpulses + pulsesInRecord->ndetpulses' has been checked previously
	gsl_vector *pulseheightAUX = gsl_vector_alloc(totalPulses);
	int cnt = 0;
	for (int i=0;i<totalPulses;i++)
	{
		if (i == totalPulses-1)		tstartnext = gsl_vector_get(tstart,i)+2*reconstruct_init->pulse_length;
		else				tstartnext = gsl_vector_get(tstart,i+1);

		if ((tstartnext-gsl_vector_get(tstart,i) > reconstruct_init->pulse_length) && ((gsl_vector_get(quality,i) == 0) || (gsl_vector_get(quality,i) == 10)))
		{
			 gsl_vector_set(pulseheightAUX,cnt,gsl_vector_get(pulseheight,i));
			 cnt = cnt +1;
		}
	}
	if (cnt == 0)
	{
		message = "No valid pulses to calculate the template (check as a possibility if PulseLength or largeFilter > Record size)";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (cnt > pulseheightAUX->size)
	{
		sprintf(valERROR,"%d",__LINE__+5);
		string str(valERROR);
		message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	temp = gsl_vector_subvector(pulseheightAUX,0,cnt);
	gsl_vector *pulseheightAUX2;
	if ((pulseheightAUX2 = gsl_vector_alloc(cnt)) == 0)
	{
		sprintf(valERROR,"%d",__LINE__-2);
		string str(valERROR);
	        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	gsl_vector_memcpy(pulseheightAUX2,&temp.vector);
	gsl_vector_free(pulseheightAUX); pulseheightAUX = 0;

	// Create the pulseheights histogram
	nBins = floor(sqrt(cnt)); 
	if ((xhisto = gsl_vector_alloc(nBins)) == 0)	// X-axis of the pulseheights histogram
	{
		sprintf(valERROR,"%d",__LINE__-2);
		string str(valERROR);
	        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	yhisto = gsl_vector_alloc(nBins);	// Y-axis of the pulseheights histogram
	if (createHisto(pulseheightAUX2, nBins, &xhisto, &yhisto))
	{
		message = "Cannot run createHisto routine";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_vector_free(pulseheightAUX2); pulseheightAUX2 = 0;

	index_maximumpulseheight = gsl_vector_max_index(yhisto);
	maximumpulseheight = gsl_vector_get(xhisto,index_maximumpulseheight);

	// Calculate the pulseaverage only taking into account the valid pulses
	gsl_vector_set_all(*pulseaverage,0.0);
	for (int i=0;i<totalPulses;i++)
	{
		if (i == totalPulses-1)		tstartnext = gsl_vector_get(tstart,i)+2*reconstruct_init->pulse_length;
		else 				tstartnext = gsl_vector_get(tstart,i+1);

		// Check if the pulse is piled-up or not
		if ((gsl_vector_get(pulseheight,i) < maximumpulseheight-0.1*maximumpulseheight) || (gsl_vector_get(pulseheight,i) > maximumpulseheight+0.1*maximumpulseheight)
			|| (tstartnext-gsl_vector_get(tstart,i) <= reconstruct_init->pulse_length) || ((gsl_vector_get(quality,i) != 0) && (gsl_vector_get(quality,i) != 10)))
		{
 			gsl_vector_set(nonpileup,i,0);
			nonpileupPulses --;
		}
		else
		{
			// It is note necessary to check 'memcpy' because it is the non pile-up case
			if (i < pulsesAll->ndetpulses)
			{
                          
                          //printf("copy %zu - %zu\n",
                          //       pulse->size, 
                          //       pulsesAll->pulses_detected[i].pulse_adc->size);
                          gsl_vector_memcpy(pulse,pulsesAll->pulses_detected[i].pulse_adc);	
			}
			else
			{
				gsl_vector_memcpy(pulse,pulsesInRecord->pulses_detected[i-pulsesAll->ndetpulses].pulse_adc);
			}

			// Non piled-up pulses => Align and average them
			if (firstnonpileupPulse == true)
			{
				gsl_vector_memcpy(*pulseaverage,pulse);
			}
			else
			{
				/*if (align(samprate, pulseaverage,&pulse))
				{
					message = "Cannot run align for pulse " + boost::lexical_cast<std::string>(i) + " when 1st pulse is piled-up";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}*/
				gsl_vector_add(*pulseaverage,pulse);
			}
			*pulseaverageHeight = *pulseaverageHeight + gsl_vector_get(pulseheight,i);
			if (firstnonpileupPulse == true)	firstnonpileupPulse = false;
		}
	}
	
	gsl_vector_scale(*pulseaverage,1.0/(nonpileupPulses));
 	
	gsl_vector_memcpy(*pulseaverageMaxLengthFixedFilter, *pulseaverage);
	gsl_vector_free(*pulseaverage);
	*pulseaverage = gsl_vector_alloc(inputPulseLength);
	temp = gsl_vector_subvector(*pulseaverageMaxLengthFixedFilter,0,inputPulseLength);
	gsl_vector_memcpy(*pulseaverage,&temp.vector);

        if (reconstruct_init->hduPRECALWN == 1)
        {
                // Calculate covariance and weight matrix
                bool saturatedPulses = false;
                if (pulsesInRecord->pulses_detected[0].quality >= 10)		saturatedPulses = true;
                if (weightMatrix(reconstruct_init, saturatedPulses, pulsesAll, pulsesInRecord, nonpileupPulses, nonpileup, *pulseaverage, covariance, weight))
                {
                        message = "Cannot run weightMatrix routine";
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                }
        }

	*pulseaverageHeight = *pulseaverageHeight/nonpileupPulses;

	// Free allocated GSL vectors
	gsl_vector_free(tstart); tstart = 0;
	gsl_vector_free(pulseheight); pulseheight = 0;
	gsl_vector_free(quality); quality = 0;
	gsl_vector_free(nonpileup); nonpileup = 0;
	gsl_vector_free(xhisto); xhisto = 0;
	gsl_vector_free(yhisto); yhisto = 0;
	gsl_vector_free(pulse); pulse = 0;

	return (EPOK);
}
/*xxxx end of SECTION A8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A9 ************************************************************
* createHisto function: This function builds the histogram of the input vector.
*                       Histogram x-axis values are the different input vector values (pulseheights).
*                       Histogram y-axis values are the the number of cases per unit of the variable on the horizontal axis.
*
* - Declare variables
* - It will work with the positive elements of the input vector -> 'invectoraux2'
* - Check if all the values of 'invector' are the same => Histogram of only one bin
* - Obtain 'invector_max 'and 'invector_min'
* - Obtain 'binSize'
* - Create histogram axis
* - Free allocated GSL vectors
*
* Parameters:
* - invector: GSL nput vector
* - nbins: Number of bins to build the histogram
* - xhistogsl: GSL vector with output histogram x-axis
* - yhistogsl: GSL vector with output histogram y-axis
******************************************************************************/
int createHisto (gsl_vector *invector, int nbins, gsl_vector **xhistogsl, gsl_vector **yhistogsl)
{
	string message = "";
	char valERROR[256];

	// Declare variables
	int size = invector->size;
	double invectormax= 0;      				// Maximum of 'invector'
	double invectormin=1e10;  				// Minimum of 'invector'
	double binSize;						// Size in samples of each bin
	int ind = 0;                				// Index of the bin which contains each 'invector' element
	// It is not necessary to check the allocation because 'invector' size must already be > 0
	gsl_vector *invectoraux = gsl_vector_alloc(size);	// Auxiliary variable
	gsl_vector *invectoraux2;				// Auxiliary variable
	gsl_vector_view temp;					// In order to handle with gsl_vector_view (subvectors)

	// It will work with the positive elements of the input vector -> 'invectoraux2'
	for (int i=0;i<size;i++)
	{
		if (gsl_vector_get(invector,i) > 0)
		{
			gsl_vector_set(invectoraux,ind,gsl_vector_get(invector,i));
			ind = ind+1;
		}
	}
	if ((ind < 1) || (ind > invectoraux->size))
	{
		sprintf(valERROR,"%d",__LINE__+5);
		string str(valERROR);
		message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	temp = gsl_vector_subvector(invectoraux,0,ind);
	if ((invectoraux2 = gsl_vector_alloc(ind)) == 0) 
	{
		sprintf(valERROR,"%d",__LINE__-2);
		string str(valERROR);
	        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	gsl_vector_memcpy(invectoraux2, &temp.vector);
	size = invectoraux2->size;
	gsl_vector_free(invectoraux); invectoraux = 0;

	// To check if all the values of 'invector' are the same
	// For example, high energies and no noise => Saturated pulses
	// It is not necessary to check the allocation because 'invector' size must already be > 0
	invectoraux = gsl_vector_alloc(size);
	gsl_vector *invectoraux1 = gsl_vector_alloc(size);
	gsl_vector_memcpy(invectoraux,invectoraux2);
	gsl_vector_set_all(invectoraux1,gsl_vector_get(invectoraux2,0));
	gsl_vector_scale(invectoraux1,-1.0);
	gsl_vector_add(invectoraux,invectoraux1);
	gsl_vector_free(invectoraux1); invectoraux1 = 0;

	if (gsl_vector_isnull(invectoraux) == 1) // All the values are the same
	{
		gsl_vector_free(*xhistogsl);
		gsl_vector_free(*yhistogsl);
		*xhistogsl = gsl_vector_alloc(1);
		*yhistogsl = gsl_vector_alloc(1);
		gsl_vector_set(*xhistogsl,0,gsl_vector_get(invectoraux2,0));
		gsl_vector_set(*yhistogsl,0,1); //It doesn't matter the value because the maximum is going to look for
	}
	else
	{
		// Obtain 'invector_max'
		for (int i=0; i<size; i++)
		{
			if (invectormax < gsl_vector_get(invectoraux2,i))	invectormax = gsl_vector_get(invectoraux2,i);
		}
		// Obtain 'invector_min'
		for (int i=0; i<size; i++)
		{
			if (invectormin > gsl_vector_get(invectoraux2,i))	invectormin = gsl_vector_get(invectoraux2,i);
		}

		if (invectormax == invectormin)
		{
			message = "Invectormax == invectormin";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		if ((invectormax != 0) && (invectormin != 1e10))
		{
			// Obtain binSize
			binSize = (invectormax - invectormin) / nbins;        // Size of each bin of the histogram

			// Create histogram x-axis
			for (int i=0; i<nbins; i++)	{gsl_vector_set (*xhistogsl,i,binSize*i+invectormin+binSize/2);}

			// Create histogram y-axis
			gsl_vector_set_zero (*yhistogsl);
			for (int i=0; i<size; i++)
			{
				ind = (int) ((gsl_vector_get(invectoraux2,i) - invectormin) / binSize);
				if (ind == nbins) ind--;
				if ((ind < 0) || (ind > (*yhistogsl)->size-1))
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "Setting/Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				gsl_vector_set (*yhistogsl,ind, gsl_vector_get(*yhistogsl,ind) + 1);
			}
		}
	}

	// Free allocated GSL vectors
	gsl_vector_free(invectoraux); invectoraux = 0;
	gsl_vector_free(invectoraux2); invectoraux2 = 0;

	return EPOK;
}
/*xxxx end of SECTION A9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A10 ************************************************************
* align function: This function aligns 'vector1' with 'vector2' (by delaying or moving forward 'vector2') assuming that 'vector1' and 'vector2' are shifted replicas of the same function.
*
* From the discrete function x[n] (n=0,...,N-1 => Length = N) and according to the time shifting property of the Fourier transform:
*
*  x[n]   <------> X[f]
*  x[n-m] <------> X[f]·exp(-j2·pi·m/N)
*
*  Shift = m => Phase due to the shift = 2·pi·m/N => m = Phase due to the shift·N/(2·pi)
*
* - Declare variables
* - FFT of 'vector1'
* - FFT of 'vector2'
* - Phases of the FFT_vector1 and FFT_vector2, *size/(2*pi)
* - Shift between the input vectors
* - 'shiftdouble' into 'shiftint' (because we are working with samples)
* - Move forward or delay 'vector2' depending on positive or negative shift
*
* Parameters:
* - samprate: Sampling rate
* - vector1: GSL vector with input vector
* - vector2: GSL with input vector which is delayed or moved forward to be aligned with 'vector1' 
******************************************************************************/
int align(double samprate, gsl_vector **vector1, gsl_vector **vector2)
{
	const double pi = 4.0 * atan(1.0);
	string message = "";

	// Declare variables
	int size = (*vector1)->size;

	double SelectedTimeDuration = size/samprate;
	// It is not necessary to check the allocation because 'vector1' size must already be > 0
	gsl_vector_complex *vector1fft = gsl_vector_complex_alloc(size);
	gsl_vector_complex *vector2fft = gsl_vector_complex_alloc(size);
	double vector1fft_ph;
	double vector2fft_ph;

	double shiftdouble;
	int shiftint;

	// It is not necessary to check the allocation because 'vector1' size must already be > 0
	gsl_vector *vector2shifted = gsl_vector_alloc(size);

	// FFT of 'vector1'
	if (FFT(*vector1,vector1fft,SelectedTimeDuration))
	{
		message = "Cannot run FFT routine for vector1";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// FFT of 'vector2'
	if (FFT(*vector2,vector2fft,SelectedTimeDuration))
	{
		message = "Cannot run FFT routine for vector2";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Phases of the FFT_vector1 and FFT_vector2, *size/(2*pi)
	vector1fft_ph= gsl_complex_arg(gsl_vector_complex_get(vector1fft,1))*size/(2*pi);
	vector2fft_ph= gsl_complex_arg(gsl_vector_complex_get(vector2fft,1))*size/(2*pi);
	gsl_vector_complex_free(vector1fft); vector1fft = 0;
	gsl_vector_complex_free(vector2fft); vector2fft = 0;

	// Shift between the input vectors
	shiftdouble = vector1fft_ph-vector2fft_ph;

	// 'shiftdouble' into 'shiftint' (because we are working with samples)
	if ((shiftdouble > -1) && (shiftdouble < 1)) 	shiftint = 0;
	else if (shiftdouble > 1)			shiftint = floor(shiftdouble);
	else if (shiftdouble < -1)			shiftint = ceil(shiftdouble);

	// Move forward or delay vector2 depending on positive or negative shift
	if (shiftint > 0)
	{
		if (shift_m(*vector2,vector2shifted,shiftint))
		{
			message = "Cannot run shift_m routine for vector2";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		gsl_vector_memcpy(*vector2,vector2shifted);
	}
	else if (shiftint < 0)
	{
		if (shiftm(*vector2,vector2shifted,fabs(shiftint)))
		{
			message = "Cannot run shiftm routine for vector2";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		gsl_vector_memcpy(*vector2,vector2shifted);
	}

	gsl_vector_free(vector2shifted); vector2shifted = 0;

	return (EPOK);
}
/*xxxx end of SECTION A10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A11 ************************************************************
* shiftm function: This function returns as 'vectorout' the 'vectorin' delayed 'm' samples.
*
* Parameters:
* - vectorin: GSL vector with input vector
* - vectorout: GSL with input vector ('vectorin') is delayed by 'm' samples
* - m: Delay in samples
******************************************************************************/
int shiftm(gsl_vector *vectorin, gsl_vector *vectorout, int m)
{
	string message = "";
	char valERROR[256];
	
	int size = vectorin->size;

	if (size-m > vectorin->size)
	{
		sprintf(valERROR,"%d",__LINE__+14);
		string str(valERROR);
		message = "Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (size-m-1+m > vectorout->size-1)
	{
		sprintf(valERROR,"%d",__LINE__+7);
		string str(valERROR);
		message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	for (int i=0;i<size-m;i++)
	{
		gsl_vector_set(vectorout,i+m,gsl_vector_get(vectorin,i));
	}
	
	if (m-1 > vectorout->size-1)
	{
		sprintf(valERROR,"%d",__LINE__+7);
		string str(valERROR);
		message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	for (int i=0;i<m;i++)
	{
		gsl_vector_set(vectorout,i,gsl_vector_get(vectorin,0));
	}

	return (EPOK);
}
/*xxxx end of SECTION A11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A12 ************************************************************
* shift_m function: This function returns as 'vectorout' the 'vectorin' moved forward 'm' samples.
*
* Parameters:
* - vectorin: GSL vector with input vector
* - vectorout: GSL with input vector (:option:`vectorin`) is moved forward :option:`m` samples
* - m: Advance in samples
******************************************************************************/
int shift_m(gsl_vector *vectorin, gsl_vector *vectorout, int m)
{
	string message = "";
	char valERROR[256];
	
	int size = vectorin->size;

	if (m < 0)
	{
		sprintf(valERROR,"%d",__LINE__+14);
		string str(valERROR);
		message = "Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (size-1-m > vectorout->size-1)
	{
		sprintf(valERROR,"%d",__LINE__+7);
		string str(valERROR);
		message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	for (int i=m;i<size;i++)
	{
		gsl_vector_set(vectorout,i-m,gsl_vector_get(vectorin,i));
	}
	
	if ((size-m < 0) || (m-1 > vectorout->size-1))
	{
		sprintf(valERROR,"%d",__LINE__+7);
		string str(valERROR);
		message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	for (int i=size-m;i<size;i++)
	{
		gsl_vector_set(vectorout,i,gsl_vector_get(vectorin,size-1));
	}

	return (EPOK);
}
/*xxxx end of SECTION A12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A13 ************************************************************
* weightMatrix function: This function calculates the weight matrix by using the non piled-up pulses found in all the records, stored
*                        in 'pulsesAll' (previous records) and 'pulsesInRecord' (current record). The weight matrix of each energy
*                        (and other intermediate values) will be stored in the library by the function 'fillInLibraryData'.
*
* Si^p: Value of the ith-sample of the pulse p
* Mi^p: Value of the ith-sample of the model p (model='pulseaverage')  Mi = <Si> = (1/N)sum(p=1,N){Si^p}
* N: Number of non piled-up pulses
* Di = Si - Mi
* <DiDj> = E[(Si-Mi)(Sj-Mj)] = (1/N)sum(p=1,N){(Si^p-Mi^p)(Sj^p-Mj^p)}
* 		  |<D1D1> <D1D2>...<D1Dn>|
*  Vij = <DiDj> = |<D2D1> <D2D2>...<D2Dn>|	where n is the pulse length     V => (nxn)
*                 |...                   |
*                 |<DnD1> <DnD2>...<DnDn>|
*  W = 1/V
*
* - Calculate the elements of the diagonal of the covariance matrix
* - Calculate the elements out of the diagonal of the covariance matrix
* - If saturated pulses => Covariance matrix is a singular matrix => Non invertible 
*   In order to allow the covariance matrix to be inverted => Replacing 0's (0's are due to the saturated values, equal in the pulse and in the model)
* 	- Elements of the diagonal: Generating a random double f1 between a range (fMin,fMax), (-NoiseStd,NoiseStd), to replace 0's with f1*f1 
*       - Elements out of the diagonal: Generating two random doubles f1 and f2 between a range (fMin,fMax), (-NoiseStd,NoiseStd), to replace 0's with f1*f2 
* - Calculate the weight matrix
*
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
*                     In particular, this function uses 'pulse_length' and 'EnergyMethod'
* - saturatedPulses: If 'true', all the pulses (CALIBRATION => all the pulses have the same energy) are saturated
* - pulsesAll: Collection of pulses found in the previous records
* - pulsesInRecord: Found pulses in the current record
* - nonpileupPulses: Number of non piled-up pulses
* - nonpileup: GSL vector containing info about all the pulses informing if they are piled-up or not
* - pulseaverage: GSL vector with the pulseaverage (= template = model) of the non piled-up pulses
* - covariance: GSL matrix with covariance matrix
* - weight: GSL matrix with weight matrix
******************************************************************************/
int weightMatrix (ReconstructInitSIRENA *reconstruct_init, bool saturatedPulses, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, long nonpileupPulses, gsl_vector *nonpileup, gsl_vector *pulseaverage, gsl_matrix **covariance, gsl_matrix **weight)
{
	string message = "";
	char valERROR[256];
	
	double elementValue = 0.0;
	// It is not necessary to check the allocation because 'reconstruct_init->pulse_length'='PulseLength' (input parameter) has been checked previously
	//gsl_permutation *perm = gsl_permutation_alloc(reconstruct_init->pulse_length);
	gsl_permutation *perm = gsl_permutation_alloc(pulseaverage->size);
	int s=0;

	// Elements of the diagonal of the covariance matrix
	for (int i=0;i<pulseaverage->size;i++)
	{
		for (int p=0;p<pulsesAll->ndetpulses;p++)
		{
			if (gsl_vector_get(nonpileup,p) == 1)
			{
				elementValue = elementValue +
					pow(gsl_vector_get(pulsesAll->pulses_detected[p].pulse_adc,i)-gsl_vector_get(pulseaverage,i),2.0);
			}
		}
		for (int p=0;p<pulsesInRecord->ndetpulses;p++)
		{
			if (gsl_vector_get(nonpileup,pulsesAll->ndetpulses+p) == 1)
			{
				elementValue = elementValue +
					pow(gsl_vector_get(pulsesInRecord->pulses_detected[p].pulse_adc,i)-gsl_vector_get(pulseaverage,i),2.0);
			}
		}
		elementValue = elementValue/nonpileupPulses;
		
		gsl_matrix_set(*covariance,i,i,elementValue);

		elementValue = 0.0;
	}

	// Other elements
	for (int i=0;i<pulseaverage->size;i++)
	{
		for (int j=i+1;j<pulseaverage->size;j++)
		{
			for (int p=0;p<pulsesAll->ndetpulses;p++)
			{
				if (gsl_vector_get(nonpileup,p) == 1)
				{
					elementValue = elementValue +
						(gsl_vector_get(pulsesAll->pulses_detected[p].pulse_adc,i)-gsl_vector_get(pulseaverage,i))*
						(gsl_vector_get(pulsesAll->pulses_detected[p].pulse_adc,j)-gsl_vector_get(pulseaverage,j));	
				}
			}
			for (int p=0;p<pulsesInRecord->ndetpulses;p++)
			{
				if (gsl_vector_get(nonpileup,pulsesAll->ndetpulses+p) == 1)
				{
					elementValue = elementValue +
						(gsl_vector_get(pulsesInRecord->pulses_detected[p].pulse_adc,i)-gsl_vector_get(pulseaverage,i))*
						(gsl_vector_get(pulsesInRecord->pulses_detected[p].pulse_adc,j)-gsl_vector_get(pulseaverage,j));
				}
			}
			elementValue = elementValue/nonpileupPulses;

			gsl_matrix_set(*covariance,i,j,elementValue);
			gsl_matrix_set(*covariance,j,i,elementValue);

			elementValue = 0.0;
		}
	}
	
	if (strcmp(reconstruct_init->EnergyMethod,"PCA") != 0)	// Different from PCA
	{
		// If saturated pulses => Covariance matrix is a singular matrix => Non invertible 
		// In order to allow the covariance matrix to be inverted => Replacing 0's (0's are due to the saturated values, equal in the pulse and in the model)
		//    - Elements of the diagonal: Generating a random double f1 between a range (fMin,fMax), (-NoiseStd,NoiseStd), to replace 0's with f1*f1 
		//    - Elements out of the diagonal: Generating two random doubles f1 and f2 between a range (fMin,fMax), (-NoiseStd,NoiseStd), to replace 0's with f1*f2
		int cnt = 0;
		if (saturatedPulses == true)
		{
			double noiseStd = reconstruct_init->noise_spectrum->noiseStd;
			double f1 = 0;
			double f2 = 0;
			double fMin = (-1.0)*noiseStd;
			double fMax = noiseStd;
				
			for (int i=0;i<pulseaverage->size;i++)
			{
				for (int j=i;j<pulseaverage->size;j++)		// In order to have a symmetric matrix (j=i)
				{
					if (fabs(gsl_matrix_get(*covariance,i,j)) < 1e-24)	// Meaning 'equal to 0'
					{	
						cnt++;
						srand(time(NULL)+i+j+cnt);			// In order to change each run time the seed
						f1 = (double)rand()/RAND_MAX; 
						f1 = fMin + f1 * (fMax - fMin);
						if (i == j)
						{
							srand(time(NULL)+i+j+cnt+1);		// In order to get two different numbers (f1 and f2)
							f2 = (double)rand()/RAND_MAX; 
							f2 = fMin + f2 * (fMax - fMin);
							gsl_matrix_set(*covariance,i,j,f1*f2);
						}
						else	
						{
							gsl_matrix_set(*covariance,i,j,f1*f1);	
							gsl_matrix_set(*covariance,j,i,f1*f1);	// In order to have a symmetric matrix
						}
					}
				}
			}
		}
		
		// Calculate the weight matrix
		// It is not necessary to check the allocation because 'covarianze' size must already be > 0
		gsl_matrix *covarianceaux = gsl_matrix_alloc((*covariance)->size1,(*covariance)->size2);
		gsl_matrix_memcpy(covarianceaux,*covariance);
		gsl_linalg_LU_decomp(covarianceaux, perm, &s);
		if (gsl_linalg_LU_invert(covarianceaux, perm, *weight) != 0) 
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);
			message = "Singular matrix in line " + str + " (" + __FILE__ + ")";
			EP_PRINT_ERROR(message,EPFAIL);	return(EPFAIL);
		}

		gsl_matrix_free(covarianceaux); covarianceaux = 0;
		gsl_permutation_free(perm); perm = 0;
	}
	else	// PCA
	{
		gsl_matrix_memcpy(*weight,*covariance);
	}  

	return (EPOK);
}
/*xxxx end of SECTION A13 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION Axx ************************************************************
* weightMatrixReduced function: This function calculates the weight matrix by using the non piled-up pulses found in all the records which have been transformed
*                               according to a reduced set of eigenvectors.
*
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
*                     In particular, this function uses 'pulse_length'
* - saturatedPulses: If 'true', all the pulses (CALIBRATION => all the pulses have the same energy) are saturated
* - TransformedData: Non piled-up pulses transformed according to a reduced set of eigenvectors
* - TransformedPulseAverage: Average of the non piled-up pulses transformed according to a reduced number of eigenvectors, average of 'TransformedData'
* - covariancer: GSL matrix with reduced covariance matrix
* - weightr: GSL matrix with reduced weight matrix
******************************************************************************/
int weightMatrixReduced (ReconstructInitSIRENA *reconstruct_init, bool saturatedPulses, gsl_matrix *TransformedData, gsl_vector *TransformedPulseAverage , gsl_matrix **covariancer, gsl_matrix **weightr)
{
	string message = "";
	char valERROR[256];
	
	double elementValue = 0.0;
	// It is not necessary to check the allocation because 'covarianzer' size must already be > 0
	gsl_permutation *perm = gsl_permutation_alloc((*covariancer)->size1);
	int s=0;

	// Elements of the diagonal of the covariance matrix
	for (int i=0;i<TransformedData->size1;i++)
	{
		for (int p=0;p<TransformedData->size2;p++)
		{	
			elementValue = elementValue +
					pow(gsl_matrix_get(TransformedData,i,p)-gsl_vector_get(TransformedPulseAverage,i),2.0);	
		}
		
		elementValue = elementValue/TransformedData->size2;

		gsl_matrix_set(*covariancer,i,i,elementValue);

		elementValue = 0.0;
	}

	// Other elements
	for (int i=0;i<TransformedData->size1;i++)
	{
		for (int j=i;j<TransformedData->size1;j++)
		{
			for (int p=0;p<TransformedData->size2;p++)
			{
				elementValue = elementValue +
						(gsl_matrix_get(TransformedData,i,p)-gsl_vector_get(TransformedPulseAverage,i))*
						(gsl_matrix_get(TransformedData,j,p)-gsl_vector_get(TransformedPulseAverage,j));	
			}
			
			elementValue = elementValue/TransformedData->size2;

			gsl_matrix_set(*covariancer,i,j,elementValue);
			gsl_matrix_set(*covariancer,j,i,elementValue);

			elementValue = 0.0;
		}
	}
	
	// If saturated pulses => Covariance matrix is a singular matrix => Non invertible 
	// In order to allow the covariance matrix to be inverted => Replacing 0's (0's are due to the saturated values, equal in the pulse and in the model)
	//    - Elements of the diagonal: Generating a random double f1 between a range (fMin,fMax), (-NoiseStd,NoiseStd), to replace 0's with f1*f1 
	//    - Elements out of the diagonal: Generating two random doubles f1 and f2 between a range (fMin,fMax) ,(-NoiseStd,NoiseStd), to replace 0's with f1*f2
	int cnt = 0;
	if (saturatedPulses == true)
	{
		double noiseStd = reconstruct_init->noise_spectrum->noiseStd;
		double f1 = 0;
		double f2 = 0;
		double fMin = (-1.0)*noiseStd;
		double fMax = noiseStd;
			
		for (int i=0;i<TransformedData->size1;i++)
		{
			for (int j=i;j<TransformedData->size2;j++)			// In order to have a symmetric matrix (j=i)
			{
				if (fabs(gsl_matrix_get(*covariancer,i,j)) < 1e-24)	// Meaning 'equal to 0'
				{	
					cnt++;
					srand(time(NULL)+i+j+cnt);			// In order to change each run time the seed
					f1 = (double)rand()/RAND_MAX; 
					f1 = fMin + f1 * (fMax - fMin);
					if (i == j)
					{
						srand(time(NULL)+i+j+cnt+1);		// In order to get two different numbers (f1 and f2)
						f2 = (double)rand()/RAND_MAX; 
						f2 = fMin + f2 * (fMax - fMin);
						gsl_matrix_set(*covariancer,i,j,f1*f2);
					}
					else	
					{
						gsl_matrix_set(*covariancer,i,j,f1*f1);	
						gsl_matrix_set(*covariancer,j,i,f1*f1);	// In order to have a symmetric matrix
					}
				}
			}
		}
	}
	
	// Calculate the weight matrix
	// It is not necessary to check the allocation because 'covarianzer' size must already be > 0
	gsl_matrix *covarianceaux = gsl_matrix_alloc((*covariancer)->size1,(*covariancer)->size2);
	gsl_matrix_memcpy(covarianceaux,*covariancer);
	gsl_linalg_LU_decomp(covarianceaux, perm, &s);
	if (gsl_linalg_LU_invert(covarianceaux, perm, *weightr) != 0) 
	{
		sprintf(valERROR,"%d",__LINE__-2);
		string str(valERROR);
		message = "Singular matrix in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	gsl_matrix_free(covarianceaux); covarianceaux = 0;

	gsl_permutation_free(perm); perm = 0;

	return (EPOK);
}
/*xxxx end of SECTION Axx xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION Axx ************************************************************
* eigenVV function: This funcion provides the principal eigenvectors and eigenvalues of the input matrix (at the moment, the first two eigenvalues and eigenvectors).
*                   The eigenvalues and eigenvectors are sorted in descending order and only the principal components are provided.
* 
* Parameters:
* - matrixin: Input GSL matrix
* - eigenvectors: Subset of eigenvectors of 'matrixin' chosen by PCA (the first one or the first two ones)
* - eigenvalues: Subset of eigenvalues of 'matrixin' chosen by PCA (the first one or the first two ones)
******************************************************************************/
int eigenVV (gsl_matrix *matrixin, gsl_matrix **eigenvectors, gsl_vector **eigenvalues)
{
	int status = EPOK;
	string message = "";
	char valERROR[256];
	
	// It is not necessary to check the allocation because 'matrixin' size must already be > 0
	gsl_vector *eigenvaluesAll = gsl_vector_alloc(matrixin->size1);
	gsl_matrix *eigenvectorsAll = gsl_matrix_alloc(matrixin->size1,matrixin->size2);
	
	// Calculate the eigenvectors and the eigenvalues
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (matrixin->size1);
	gsl_eigen_symmv (matrixin, eigenvaluesAll, eigenvectorsAll, w);
	gsl_eigen_symmv_free (w); w = 0;

	// Sort the eigenvectors and the eigenvalues in descending order
	gsl_eigen_symmv_sort (eigenvaluesAll, eigenvectorsAll, GSL_EIGEN_SORT_ABS_DESC);
	
	// Choose the main eigenvectors and eigenvalues (the principal components)
	// Choose eigenvectors whose abs(eigenvalues) are grater than 1
	int indexToEndToTakeAccount = eigenvaluesAll->size;
	for (int i = 0; i < eigenvaluesAll->size; i++)
	{
		if (abs(gsl_vector_get(eigenvaluesAll,i)) < 1.0)
		{
			indexToEndToTakeAccount = i;	
			break;
		}
	}
	//Choose the first eigenvector and eigenvalue
	//indexToEndToTakeAccount = 1;
	//Choose the first two eigenvectors and eigenvalues
	indexToEndToTakeAccount = 2;
	gsl_vector_view temp;
	gsl_matrix_view tempm;
	*eigenvalues = gsl_vector_alloc(indexToEndToTakeAccount);
	// It is not necessary to check the allocation because 'eigenvectorsAll' size and 'eigenvalues' size must already be > 0
	*eigenvectors = gsl_matrix_alloc(eigenvectorsAll->size1,(*eigenvalues)->size);	// Eigenvectors in columns
	
	if ((indexToEndToTakeAccount < 1) || (indexToEndToTakeAccount > eigenvaluesAll->size))
	{
		sprintf(valERROR,"%d",__LINE__+5);
		string str(valERROR);
		message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	temp = gsl_vector_subvector(eigenvaluesAll,0,indexToEndToTakeAccount);
	gsl_vector_memcpy(*eigenvalues,&temp.vector);
	if ((indexToEndToTakeAccount < 1) || (indexToEndToTakeAccount > eigenvectorsAll->size2))
	{
		sprintf(valERROR,"%d",__LINE__+5);
		string str(valERROR);
		message = "View goes out of scope the original matrix in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	tempm = gsl_matrix_submatrix(eigenvectorsAll,0,0,eigenvectorsAll->size1,indexToEndToTakeAccount);
	gsl_matrix_memcpy(*eigenvectors,&tempm.matrix);
	
	// Store the eigenvectors and eigenvalues into the corresponding output FITS file
	FILE * temporalFile;
        char temporalFileName[255];
	sprintf(temporalFileName,"eigenvalues.txt");
        temporalFile = fopen (temporalFileName,"w");
	char val[256];
	for (int i = 0; i < eigenvaluesAll->size; i++)	
	{	
		sprintf(val,"%e",gsl_vector_get(eigenvaluesAll,i));
		strcat(val,"\n");
		fputs(val,temporalFile);
	}
        fclose(temporalFile);
	sprintf(temporalFileName,"eigenvectors.txt");
        temporalFile = fopen (temporalFileName,"w");
	for (int i = 0; i < eigenvectorsAll->size1; i++)	
	{	
		sprintf(val,"%e %e %e %e",gsl_matrix_get(eigenvectorsAll,i,0),gsl_matrix_get(eigenvectorsAll,i,1),gsl_matrix_get(eigenvectorsAll,i,2),gsl_matrix_get(eigenvectorsAll,i,3));
		strcat(val,"\n");
		fputs(val,temporalFile);
	}
	fclose(temporalFile);
	
	// Free allocated GSL vectors and matrices
	gsl_vector_free(eigenvaluesAll); eigenvaluesAll = 0;
	gsl_matrix_free(eigenvectorsAll); eigenvectorsAll = 0;
	
	return (EPOK);  
}
/*xxxx end of SECTION Axx xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A14 ************************************************************
* writeLibrary function: This function writes the library (reordering if it is necessary and calculating some intermediate parameters)
*
* - Adding a new row to the library
* - Write the first row of the library
* 
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
* - samprate: Sampling rate
* - estenergy: Pulse height of the template whose energy is going to be added to the library
* - pulsetemplate: GSL vector with the pulse template whose energy is going to be added to the library
* - covariance: GSL matrix with covariance 
* - weight: GSL matrix with weight matrix
* - appendToLibrary: 'true' if adding a new row to the library
*                    'false' if it is the first row to be added
* - inLibObject: FITS object containing information of the library FITS file 
* - pulsetemplateMaxLengthFixedFilter: GSL vector with the largeFilter-length template whose energy is going to be added to the library
******************************************************************************/
int writeLibrary(ReconstructInitSIRENA **reconstruct_init, double samprate, double estenergy, gsl_vector *pulsetemplate, gsl_matrix *covariance, gsl_matrix *weight, bool appendToLibrary, fitsfile **inLibObject, gsl_vector *pulsetemplateMaxLengthFixedFilter)
{
	int status = EPOK;
	string message = "";

	char inLibName[256];
	strncpy(inLibName, (*reconstruct_init)->library_file,255);
	inLibName[255]='\0';
	
	char keyname[10];
	char extname[20];
	
	int runF0orB0val;
	if (strcmp((*reconstruct_init)->FilterMethod,"F0") == 0)		// Deleting the frequency-zero bin
	{
		runF0orB0val = 0;
	}
	else if (strcmp((*reconstruct_init)->FilterMethod,"B0") == 0)	// Working without baseline
	{
		runF0orB0val = 1;
	}

	// Adding a new row to the library
	if (appendToLibrary == true)
	{
		long eventcntLib;
		if (fits_get_num_rows(*inLibObject,&eventcntLib, &status))
		{
			message = "Cannot get number of rows in " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		
		assert(eventcntLib > 0);
		long eventcntLib1 = eventcntLib + 1;
		strcpy(keyname,"EVENTCNT");
		
		if (fits_update_key(*inLibObject,TLONG,keyname, &eventcntLib1,NULL,&status))
		{
			message = "Cannot update keyword " + string(keyname);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		
		if (eventcntLib1-1 > 1)
		{
			strcpy(extname,"PRECALWN");
			if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
			{
				(*reconstruct_init)->hduPRECALWN = 0;
				status = EPOK;
			}
		}
		
		strcpy(extname,"PRCLOFWM");
		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
                {
			(*reconstruct_init)->hduPRCLOFWM = 0;
			status = EPOK;
                }		
		                                
                strcpy(extname,"LIBRARY");
                if (fits_movabs_hdu(*inLibObject, 2, NULL, &status))	// LIBRARY
                {
                        message = "Cannot move to HDU " + string(extname) + " in library file " + string(inLibName);
                        EP_PRINT_ERROR(message,status); return(EPFAIL);
                }

		if (readAddSortParams(*reconstruct_init,inLibObject,samprate,eventcntLib,estenergy,pulsetemplate,covariance,weight,pulsetemplateMaxLengthFixedFilter))
		{
			message = "Cannot run routine readAddSortParams in writeLibrary";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                }
		
		// Primary HDU
                strcpy(extname,"Primary");
                if (fits_movabs_hdu(*inLibObject, 1, NULL, &status))
                {
                        message = "Cannot move to HDU " + string(extname) + " in library file " + string(inLibName);
                        EP_PRINT_ERROR(message,status); return(EPFAIL);
                }

                char keyvalstr[1000];

                char str_procnumber[125];               snprintf(str_procnumber,125,"%ld",eventcntLib);
                string strprocname (string("PROC") + string(str_procnumber));
                strcpy(keyname,strprocname.c_str());
                string strprocval (string("PROC") + string(str_procnumber) + string(" Starting parameter list"));
                strcpy(keyvalstr,strprocval.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		string strproc;
		char comment[MAXMSG];
		
		sprintf(comment, "RecordFile = %s", (*reconstruct_init)->record_file);
		fits_write_comment(*inLibObject, comment, &status);

		sprintf(comment, "TesEventFile = %s", (*reconstruct_init)->event_file);
		fits_write_comment(*inLibObject, comment, &status);
		
		sprintf(comment, "LibraryFile = %s", (*reconstruct_init)->library_file);
		fits_write_comment(*inLibObject, comment, &status);
		
		sprintf(comment, "NoiseFile = %s", (*reconstruct_init)->noise_file);
		fits_write_comment(*inLibObject, comment, &status);
               
                char str_opmode[125];		snprintf(str_opmode,125,"%d",(*reconstruct_init)->opmode);
                strproc = string("opmode = ") + string(str_opmode);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("FilterDomain = ") + (*reconstruct_init)->FilterDomain;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("FilterMethod = ") + (*reconstruct_init)->FilterMethod;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("EnergyMethod = ") + (*reconstruct_init)->EnergyMethod;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_LagsOrNot[125];                snprintf(str_LagsOrNot,125,"%d",(*reconstruct_init)->LagsOrNot);
                strproc=string("LagsOrNot = ") + string(str_LagsOrNot);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
                
                char str_nLags[125];                snprintf(str_nLags,125,"%d",(*reconstruct_init)->nLags);
                strproc=string("nLags = ") + string(str_nLags);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
                
                char str_Parabola3OrFitting5[125];                snprintf(str_Parabola3OrFitting5,125,"%d",(*reconstruct_init)->Fitting35);
                strproc=string("Fitting35 = ") + string(str_Parabola3OrFitting5);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_OFIter[125];                   snprintf(str_OFIter,125,"%d",(*reconstruct_init)->OFIter);
                strproc=string("OFIter = ") + string(str_OFIter);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_OFLib[125];                    snprintf(str_OFLib,125,"%d",(*reconstruct_init)->OFLib);
                strproc=string("OFLib = ") + string(str_OFLib);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
                               
                strproc=string("OFInterp = ") + (*reconstruct_init)->OFInterp;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("OFStrategy = ") + (*reconstruct_init)->OFStrategy;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
                               
                char str_OFLength[125];                 snprintf(str_OFLength,125,"%d",(*reconstruct_init)->OFLength);
                strproc=string("OFLength = ") + string(str_OFLength);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_maxPulsesPerRecord[125];       snprintf(str_maxPulsesPerRecord,125,"%d",(*reconstruct_init)->maxPulsesPerRecord);
                strproc=string("maxPulsesPerRecord = ") + string(str_maxPulsesPerRecord);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_pulse_length[125];             snprintf(str_pulse_length,125,"%d",(*reconstruct_init)->pulse_length);
                strproc=string("PulseLength = ") + string(str_pulse_length);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_scaleFactor[125];              snprintf(str_scaleFactor,125,"%f",(*reconstruct_init)->scaleFactor);
                strproc=string("scaleFactor = ") + string(str_scaleFactor);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_samplesUp[125];                snprintf(str_samplesUp,125,"%f",(*reconstruct_init)->samplesUp);
                strproc=string("samplesUp = ") + string(str_samplesUp);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_nSgms[125];                    snprintf(str_nSgms,125,"%f",(*reconstruct_init)->nSgms);
                strproc=string("nSgms = ") + string(str_nSgms);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_LrsT[125];                     snprintf(str_LrsT,125,"%e",(*reconstruct_init)->LrsT);
                strproc=string("LrsT = ") + string(str_LrsT);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_LbT[125];                      snprintf(str_LbT,125,"%e",(*reconstruct_init)->LbT);
                strproc=string("LbT = ") + string(str_LbT);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_monoenergy[125];               snprintf(str_monoenergy,125,"%f",(*reconstruct_init)->monoenergy);
                strproc=string("monoenergy = ") + string(str_monoenergy);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_hduPRECALWN[125];      	snprintf(str_hduPRECALWN,125,"%d",(*reconstruct_init)->hduPRECALWN);
		strproc=string("hduPRECALWN = ") + string(str_hduPRECALWN);
		strcpy(keyvalstr,strproc.c_str());
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_hduPRCLOFWM[125];      	snprintf(str_hduPRCLOFWM,125,"%d",(*reconstruct_init)->hduPRCLOFWM);
		strproc=string("hduPRCLOFWM = ") + string(str_hduPRCLOFWM);
		strcpy(keyvalstr,strproc.c_str());
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_largeFilter[125];               snprintf(str_largeFilter,125,"%d",(*reconstruct_init)->largeFilter);
                strproc=string("largeFilter = ") + string(str_largeFilter);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_intermediate[125];             snprintf(str_intermediate,125,"%d",(*reconstruct_init)->intermediate);
                strproc=string("intermediate = ") + string(str_intermediate);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               		
		sprintf(comment, "detectFile = %s", (*reconstruct_init)->detectFile);
		fits_write_comment(*inLibObject, comment, &status);
		
		sprintf(comment, "filterFile = %s", (*reconstruct_init)->filterFile);
		fits_write_comment(*inLibObject, comment, &status);
               
                char str_clobber[125];                  snprintf(str_clobber,125,"%d",(*reconstruct_init)->clobber);
                strproc=string("clobber = ") + string(str_clobber);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_SaturationValue[125];          snprintf(str_SaturationValue,125,"%e",(*reconstruct_init)->SaturationValue);
                strproc=string("SaturationValue = ") + string(str_SaturationValue);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_tstartPulse1[125];             
		if ((*reconstruct_init)->tstartPulse1 == 0)			  	snprintf(str_tstartPulse1,125,"%d",(*reconstruct_init)->tstartPulse1);
		else if (strcmp((*reconstruct_init)->EnergyMethod,"I2RALL") == 0)	snprintf(str_tstartPulse1,125,"%d",(*reconstruct_init)->tstartPulse1+2);
		else if (strcmp((*reconstruct_init)->EnergyMethod,"I2RALL") != 0)	snprintf(str_tstartPulse1,125,"%d",(*reconstruct_init)->tstartPulse1+1);
                strproc=string("tstartPulse1 = ") + string(str_tstartPulse1);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_tstartPulse2[125];             
		if ((*reconstruct_init)->tstartPulse2 == 0)			  	snprintf(str_tstartPulse2,125,"%d",(*reconstruct_init)->tstartPulse2);
		else if (strcmp((*reconstruct_init)->EnergyMethod,"I2RALL") == 0)	snprintf(str_tstartPulse2,125,"%d",(*reconstruct_init)->tstartPulse2+2);
		else if (strcmp((*reconstruct_init)->EnergyMethod,"I2RALL") != 0)	snprintf(str_tstartPulse2,125,"%d",(*reconstruct_init)->tstartPulse2+1);
                strproc=string("tstartPulse2 = ") + string(str_tstartPulse2);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_tstartPulse3[125];             
		if ((*reconstruct_init)->tstartPulse3 == 0)			  	snprintf(str_tstartPulse3,125,"%d",(*reconstruct_init)->tstartPulse3);
		else if (strcmp((*reconstruct_init)->EnergyMethod,"I2RALL") == 0)	snprintf(str_tstartPulse3,125,"%d",(*reconstruct_init)->tstartPulse3+2);
		else if (strcmp((*reconstruct_init)->EnergyMethod,"I2RALL") != 0)	snprintf(str_tstartPulse3,125,"%d",(*reconstruct_init)->tstartPulse3+1);
                strproc=string("tstartPulse3 = ") + string(str_tstartPulse3);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strcpy(keyname,strprocname.c_str());
                strprocval = string("PROC") + string(str_procnumber) + string(" Ending parameter list");
                strcpy(keyvalstr,strprocval.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                if (status != 0)
                {
                        message = "Cannot write some keyword in library file " + string(inLibName);
                        EP_PRINT_ERROR(message,status); return(EPFAIL);
                }
	}
	else	// Write the first row of the library
	{
		gsl_vector *energyoutgsl = gsl_vector_alloc(1);
		gsl_vector *estenergyoutgsl = gsl_vector_alloc(1);
		// It is not necessary to check the allocation because 'reconstruct_init->pulse_length'=PulseLength(input parameter) has been checked previously
		gsl_matrix *pulsetemplates_matrix = gsl_matrix_alloc(1,(*reconstruct_init)->pulse_length);
		gsl_matrix *pulsetemplatesMaxLengthFixedFilters_matrix = gsl_matrix_alloc(1,(*reconstruct_init)->largeFilter);
		gsl_matrix *pulsetemplatesb0_matrix = gsl_matrix_alloc(1,(*reconstruct_init)->pulse_length);
		gsl_matrix *matchedfilters_matrix = gsl_matrix_alloc(1,(*reconstruct_init)->pulse_length);
		gsl_matrix *matchedfiltersb0_matrix = gsl_matrix_alloc(1,(*reconstruct_init)->pulse_length);
	
		strcpy(extname,"LIBRARY");
		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
		{
			message = "Cannot move to HDU  " + string(extname) + " in library";
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
		
		if (fits_write_key(*inLibObject,TINT,"EVENTSZ",&(*reconstruct_init)->pulse_length,NULL,&status))
		{
			message = "Cannot write keyword EVENTSZ in library";
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}

		gsl_vector_set (energyoutgsl,0,(*reconstruct_init)->monoenergy);
		gsl_vector_set (estenergyoutgsl,0,estenergy);

		gsl_matrix_set_row(pulsetemplates_matrix,0,pulsetemplate);
		gsl_matrix_set_row(matchedfilters_matrix,0,pulsetemplate);
		gsl_matrix_set_row(pulsetemplatesMaxLengthFixedFilters_matrix,0,pulsetemplateMaxLengthFixedFilter);

		// It is not necessary to check the allocation because 'pulsetemplate' size must already be > 0
		gsl_vector *baselinegsl = gsl_vector_alloc(pulsetemplate->size);
		gsl_vector_set_all(baselinegsl,-1.0*(*reconstruct_init)->noise_spectrum->baseline);
		gsl_vector_add(pulsetemplate,baselinegsl);
		gsl_vector_free(baselinegsl); baselinegsl = 0;
		gsl_matrix_set_row(pulsetemplatesb0_matrix,0,pulsetemplate);
		gsl_matrix_set_row(matchedfiltersb0_matrix,0,pulsetemplate);
		
		gsl_matrix_scale(matchedfilters_matrix,1./(*reconstruct_init)->monoenergy);
		gsl_matrix_scale(matchedfiltersb0_matrix,1./(*reconstruct_init)->monoenergy);
		
		if (addFirstRow(*reconstruct_init, inLibObject, samprate, runF0orB0val, energyoutgsl, estenergyoutgsl, pulsetemplates_matrix, pulsetemplatesb0_matrix, matchedfilters_matrix, matchedfiltersb0_matrix, covariance, weight, pulsetemplatesMaxLengthFixedFilters_matrix))
		{
			message = "Cannot run addFirstRow in writeLibrary";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		
		// Free allocated GSL vectors and matrices
		gsl_vector_free(energyoutgsl); energyoutgsl = 0;
		gsl_vector_free(estenergyoutgsl); estenergyoutgsl = 0;
		gsl_matrix_free(pulsetemplates_matrix); pulsetemplates_matrix = 0;
		gsl_matrix_free(pulsetemplatesMaxLengthFixedFilters_matrix); pulsetemplatesMaxLengthFixedFilters_matrix = 0;
		gsl_matrix_free(pulsetemplatesb0_matrix); pulsetemplatesb0_matrix = 0;
		gsl_matrix_free(matchedfilters_matrix); matchedfilters_matrix = 0;
		gsl_matrix_free(matchedfiltersb0_matrix); matchedfiltersb0_matrix = 0;
	}
	
	if (((*reconstruct_init)->hduPRECALWN == 1) && ((*reconstruct_init)->hduPRCLOFWM == 1))
	{
		message = "Library created with the PRECALWN and PRCLOFWM HDUs (and the COVARM...rE columns in the LIBRARY HDU)";
		EP_PRINT_ERROR(message,-999);	// Only a warning
	}
	else if (((*reconstruct_init)->hduPRECALWN == 0) && ((*reconstruct_init)->hduPRCLOFWM == 0))
	{
		message = "Library created without the PRECALWN and PRCLOFWM HDUs (and without the COVARM...rE columns in the LIBRARY HDU)";
		EP_PRINT_ERROR(message,-999);	// Only a warning
	}
	else if (((*reconstruct_init)->hduPRECALWN == 1) && ((*reconstruct_init)->hduPRCLOFWM == 0))
	{
		message = "Library created with the PRECALWN HDU (and the COVARM...rE columns in the LIBRARY HDU) but no the PRCLOFWM HDU";
		EP_PRINT_ERROR(message,-999);	// Only a warning
	}
	else if (((*reconstruct_init)->hduPRECALWN == 0) && ((*reconstruct_init)->hduPRCLOFWM == 1))
	{
		message = "Library created with the PRCLOFWM HDU but no the PRECALWN HDU (and without the COVARM...rE columns in the LIBRARY HDU)";
		EP_PRINT_ERROR(message,-999);	// Only a warning
	}

	if (fits_close_file(*inLibObject,&status))
	{
		message = "Cannot close file " + string(inLibName);
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
          //*inLibObject = 0;
	
	return (EPOK);
}
/*xxxx end of SECTION A14 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A15 ************************************************************
* addFirstRow function: This function writes the first row of the library (without intermediate AB-related values, because it would be necessary to have at least two rows=energies in the library). It also writes the FIXFILTT and FIXFILTF HDUs with the optimal filters in the time and frequency domain (real and imag parts) with fixed legnth and the PRCLOFWM HDU with the precalculated values for optimal filtering and modeWEIGHTM.
* 
* Parameters: 
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
*                     In particular, this function uses 'opmode' and 'noise_spectrum' in order to run 'calculus_optimalFilter'
* - inLibObject: FITS object containing information of the library FITS file
* - samprate: Sampling rate
* - runF0orB0val: 'FilterMethod' = F0 => 'runF0orB0val' = 1
*                 'FilterMethod' = B0 => 'runF0orB0val' = 0
* - E:	First energy to be included in the library
* - PHEIGHT: Pulse height associated to the first energy to be included in the library
* - PULSE: Pulse template associated to the first energy to be included in the library
* - PULSEB0: Pulse template without baseline associated to the first energy to be included in the library
* - MF: Matched filter associated to the first energy to be included in the library
* - MFB0: Matched filter (baseline subtracted) associated to the first energy to be included in the library
* - COVAR: Covariance matrix associated to the first energy to be included in the library
* - WEIGHT: Weight matrix associated to the first energy to be included in the library
* - PULSEMaxLengthFixedFilter: Pulse template whose length is largeFilter associated to the first energy to be included in the library
******************************************************************************/
int addFirstRow(ReconstructInitSIRENA *reconstruct_init, fitsfile **inLibObject, double samprate, int runF0orB0val, gsl_vector *E, gsl_vector *PHEIGHT, gsl_matrix *PULSE, 
	gsl_matrix *PULSEB0, gsl_matrix *MF, gsl_matrix *MFB0, gsl_matrix *COVAR, gsl_matrix *WEIGHT, gsl_matrix *PULSEMaxLengthFixedFilter)
{ 
	// Declare variables
	int status = EPOK;
	string message = "";
	char valERROR[256];
	
	gsl_vector *optimalfilter = NULL;
	gsl_vector *optimalfilter_f = NULL;
	gsl_vector *optimalfilter_FFT = NULL;
	gsl_vector_complex *optimalfilter_FFT_complex = NULL;
	gsl_vector *optimalfilter_FFT_RI;
	gsl_vector *optimalfilter_x = NULL;
	gsl_vector *optimalfilter_f_x = NULL;
	gsl_vector *optimalfilter_FFT_x = NULL;
	gsl_vector_complex *optimalfilter_FFT_complex_x = NULL;
	
	IOData obj;
	IOData objTIME;
	IOData objFREQ;
	IOData objPRCLOFWM;
	
	char extname[10];
	
	obj.inObject = *inLibObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"LIBRARY");
	obj.iniRow = 1;
	obj.endRow = 1;
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	
	objTIME.inObject = *inLibObject;
	objTIME.nameTable = new char [255];
	strcpy(objTIME.nameTable,"FIXFILTT");
	objTIME.iniRow = 1;
	objTIME.endRow = 1;
	objTIME.iniCol = 0;
	objTIME.nameCol = new char [255];
	
	objFREQ.inObject = *inLibObject;
	objFREQ.nameTable = new char [255];
	strcpy(objFREQ.nameTable,"FIXFILTF");
	objFREQ.iniRow = 1;
	objFREQ.endRow = 1;
	objFREQ.iniCol = 0;
	objFREQ.nameCol = new char [255];
	
	objPRCLOFWM.inObject = *inLibObject;
	objPRCLOFWM.nameTable = new char [255];
	strcpy(objPRCLOFWM.nameTable,"PRCLOFWM");
	objPRCLOFWM.iniRow = 1;
	objPRCLOFWM.endRow = 1;
	objPRCLOFWM.iniCol = 0;
	objPRCLOFWM.nameCol = new char [255];
	
	// Creating ENERGY Column
	strcpy(obj.nameCol,"ENERGY");
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit,"eV");
	if (writeFitsSimple(obj, E))
	{
		message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating PHEIGHT Column
	strcpy(obj.nameCol,"PHEIGHT");
	strcpy(obj.unit,"ADC");
	if (writeFitsSimple(obj, PHEIGHT))
	{
		message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	if (reconstruct_init->largeFilter != reconstruct_init->pulse_length)
	{
		// Creating PLSMXLFF Column
		strcpy(obj.nameCol,"PLSMXLFF");
		strcpy(obj.unit,"ADC");
		if (writeFitsComplex(obj, PULSEMaxLengthFixedFilter))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
	
	// Creating PULSE Column
	strcpy(obj.nameCol,"PULSE");
	strcpy(obj.unit,"ADC");
	if (writeFitsComplex(obj, PULSE))
	{
		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating PULSEB0 Column
	strcpy(obj.nameCol,"PULSEB0");
	strcpy(obj.unit,"ADC");
	if (writeFitsComplex(obj, PULSEB0))
	{
		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating MF Column
	strcpy(obj.nameCol,"MF");
	strcpy(obj.unit," ");
	if (writeFitsComplex(obj, MF))
	{
		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating MFB0 Column
	strcpy(obj.nameCol,"MFB0");
	strcpy(obj.unit," ");
	if (writeFitsComplex(obj, MFB0))
	{
		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	
	// It is not necessary to check the allocation because 'reconstruct_init->pulse_length'=PulseLength(input parameter) has been checked previously
	gsl_vector *matchedfilters_row = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix_get_row(matchedfilters_row,MF,0);	//Matched filter
	// Calculate the optimal filter
	if (calculus_optimalFilter (0, 0, reconstruct_init->opmode, matchedfilters_row, matchedfilters_row->size, samprate, runF0orB0val, reconstruct_init->noise_spectrum->noisefreqs, reconstruct_init->noise_spectrum->noisespec, &optimalfilter, &optimalfilter_f, &optimalfilter_FFT, &optimalfilter_FFT_complex))
	{
		message = "Cannot run routine calculus_optimalFilter in writeLibrary";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	
	if (reconstruct_init->hduPRECALWN == 1)
        {
                // Creating COVARM Column
                strcpy(obj.nameCol,"COVARM");
                strcpy(obj.unit," ");
                if (writeFitsComplex(obj, COVAR))
                {
                        message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                }

                // Creating WEIGHTM Column
                strcpy(obj.nameCol,"WEIGHTM");
                strcpy(obj.unit," ");
                if (writeFitsComplex(obj, WEIGHT))
                {
                        message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                }
        }

	// Writing HDUs with fixed filters in time (FIXFILTT) and frequency (FIXFILTF)
	char str_length[125];
	gsl_vector *fixedlengths;
	if (reconstruct_init->pulse_length != reconstruct_init->largeFilter)
	{
		fixedlengths = gsl_vector_alloc(floor(log2(reconstruct_init->pulse_length))+1); 	//+1 because we are going to add a very long fixed filter (largeFilter)
		
		gsl_vector_set(fixedlengths,0,reconstruct_init->largeFilter);
		for (int i=0;i<floor(log2(reconstruct_init->pulse_length));i++)					
		{
			gsl_vector_set(fixedlengths,i+1,pow(2,floor(log2(reconstruct_init->pulse_length))-i));
		}
	}
	else 	
	{
		fixedlengths = gsl_vector_alloc(floor(log2(reconstruct_init->pulse_length))); 	
		for (int i=0;i<floor(log2(reconstruct_init->pulse_length));i++)					
		{
			gsl_vector_set(fixedlengths,i,pow(2,floor(log2(reconstruct_init->pulse_length))-i)); 
		}
	}
	
	strcpy(extname,"FIXFILTT");
	if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
	{
		message = "Cannot move to HDU  " + string(extname) + " in library";
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	
	// Creating ENERGY Column
	strcpy(objTIME.nameCol,"ENERGY");
	objTIME.type = TDOUBLE;
	objTIME.unit = new char [255];
	strcpy(objTIME.unit,"eV");
	if (writeFitsSimple(objTIME, E))
	{
		message = "Cannot run writeFitsSimple routine for column " + string(objFREQ.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
		
	strcpy(extname,"FIXFILTF");
	if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
	{
		message = "Cannot move to HDU  " + string(extname) + " in library";
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	
	// Creating ENERGY Column
	strcpy(objFREQ.nameTable,"FIXFILTF");
	strcpy(objFREQ.nameCol,"ENERGY");
	objFREQ.type = TDOUBLE;
	objFREQ.unit = new char [255];
	strcpy(objFREQ.unit,"eV");
	if (writeFitsSimple(objFREQ, E))
	{
		message = "Cannot run writeFitsSimple routine for column " + string(objFREQ.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	
	// Creating Ts and Fs Columns (in the FIXFILTT and FIXFILTF HDUs respectively)
	gsl_vector *matchedfiltersSHORT;
	gsl_vector_view(temp);
	gsl_matrix *optimalfiltersT_matrix;
	gsl_matrix *optimalfiltersF_matrix;
	for (int j=0;j<fixedlengths->size;j++)
	{
		if (gsl_vector_get(fixedlengths,j) == optimalfilter_FFT_complex->size)
		{
                        if ((optimalfilter_x = gsl_vector_alloc(optimalfilter_FFT_complex->size)) == 0)
                        {
                                sprintf(valERROR,"%d",__LINE__-2);
                                string str(valERROR);
                                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                EP_PRINT_ERROR(message,EPFAIL);
                        }
			gsl_vector_memcpy(optimalfilter_x,optimalfilter);
		  
			optimalfilter_FFT_complex_x = gsl_vector_complex_alloc(optimalfilter_FFT_complex->size);
			gsl_vector_complex_memcpy(optimalfilter_FFT_complex_x,optimalfilter_FFT_complex);
		}
		else
		{
			// It will enter this 'else' for fixedlengths_i=largeFilter and fixedlengths_i<optimalfilter_FFT_complex->size
                        if ((matchedfiltersSHORT = gsl_vector_alloc(gsl_vector_get(fixedlengths,j))) == 0)
                        {
                                sprintf(valERROR,"%d",__LINE__-2);
                                string str(valERROR);
                                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                EP_PRINT_ERROR(message,EPFAIL);
                        }
			if (gsl_vector_get(fixedlengths,j) == reconstruct_init->largeFilter)
			{
				gsl_vector *matchedfiltersMaxLengthFixedFilter_row = gsl_vector_alloc(reconstruct_init->largeFilter);
				gsl_matrix_get_row(matchedfiltersMaxLengthFixedFilter_row,PULSEMaxLengthFixedFilter,0);	//Matched filter
				gsl_vector_scale(matchedfiltersMaxLengthFixedFilter_row,1.0/reconstruct_init->monoenergy);
                                if ((gsl_vector_get(fixedlengths,j) < 0) || (gsl_vector_get(fixedlengths,j) > matchedfiltersMaxLengthFixedFilter_row->size))
                                {
                                        sprintf(valERROR,"%d",__LINE__+5);
                                        string str(valERROR);
                                        message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                        EP_EXIT_ERROR(message,EPFAIL);
                                }
				temp = gsl_vector_subvector(matchedfiltersMaxLengthFixedFilter_row,0,gsl_vector_get(fixedlengths,j));
                                if (gsl_vector_memcpy(matchedfiltersSHORT,&temp.vector) != 0)
                                {
                                        sprintf(valERROR,"%d",__LINE__-2);
                                        string str(valERROR);
                                        message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                        EP_EXIT_ERROR(message,EPFAIL);
                                }
				gsl_vector_free(matchedfiltersMaxLengthFixedFilter_row); matchedfiltersMaxLengthFixedFilter_row = 0;
			}
			else
			{
                                if ((gsl_vector_get(fixedlengths,j) < 0) || (gsl_vector_get(fixedlengths,j) > matchedfilters_row->size))
                                {
                                        sprintf(valERROR,"%d",__LINE__+5);
                                        string str(valERROR);
                                        message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                        EP_EXIT_ERROR(message,EPFAIL);
                                }
				temp = gsl_vector_subvector(matchedfilters_row,0,gsl_vector_get(fixedlengths,j));
                                if (gsl_vector_memcpy(matchedfiltersSHORT,&temp.vector) != 0)
                                {
                                        sprintf(valERROR,"%d",__LINE__-2);
                                        string str(valERROR);
                                        message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                        EP_EXIT_ERROR(message,EPFAIL);
                                }
				gsl_vector_memcpy(matchedfiltersSHORT,&temp.vector);
			}
			  			
			// Calculate the optimal filter
			if (calculus_optimalFilter (0, 0, reconstruct_init->opmode, matchedfiltersSHORT, matchedfiltersSHORT->size, samprate, runF0orB0val, reconstruct_init->noise_spectrum->noisefreqs, reconstruct_init->noise_spectrum->noisespec, &optimalfilter_x, &optimalfilter_f_x, &optimalfilter_FFT_x, &optimalfilter_FFT_complex_x))
			{
				message = "Cannot run routine calculus_optimalFilter in writeLibrary";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
		
			gsl_vector_free(matchedfiltersSHORT); matchedfiltersSHORT = 0;
			gsl_vector_free(optimalfilter_f_x); optimalfilter_f_x = 0;
			gsl_vector_free(optimalfilter_FFT_x); optimalfilter_FFT_x = 0;
		}
		
		gsl_vector_complex_set(optimalfilter_FFT_complex_x,0,gsl_complex_rect(0.0,0.0));
		for (int i=0;i<optimalfilter_FFT_complex_x->size;i++)
		{
			gsl_vector_complex_set(optimalfilter_FFT_complex_x,i,gsl_complex_conjugate(gsl_vector_complex_get(optimalfilter_FFT_complex_x,i)));
		}
		if ((optimalfilter_FFT_RI = gsl_vector_alloc(optimalfilter_FFT_complex_x->size*2)) == 0)
                {
                        sprintf(valERROR,"%d",__LINE__-2);
                        string str(valERROR);
                        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                        EP_PRINT_ERROR(message,EPFAIL);
                }
		for (int i=0;i<optimalfilter_FFT_complex_x->size;i++)
		{
			gsl_vector_set(optimalfilter_FFT_RI,i,GSL_REAL(gsl_vector_complex_get(optimalfilter_FFT_complex_x,i)));
                        if ((i+optimalfilter_FFT_complex_x->size < 0) || (i+optimalfilter_FFT_complex_x->size > optimalfilter_FFT_RI->size-1))
                        {
                                sprintf(valERROR,"%d",__LINE__+5);
                                string str(valERROR);
                                message = "Setting with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                EP_PRINT_ERROR(message,EPFAIL);
                        }
			gsl_vector_set(optimalfilter_FFT_RI,i+optimalfilter_FFT_complex_x->size,GSL_IMAG(gsl_vector_complex_get(optimalfilter_FFT_complex_x,i)));
		}
		gsl_vector_complex_free(optimalfilter_FFT_complex_x); optimalfilter_FFT_complex_x = 0;
		
		strcpy(extname,"FIXFILTF");
		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
		{
			message = "Cannot move to HDU  " + string(extname) + " in library";
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
		snprintf(str_length,125,"%ld",optimalfilter_FFT_RI->size/2);
		strcpy(objFREQ.nameCol,(string("F")+string(str_length)).c_str());
		strcpy(objFREQ.unit," ");
		optimalfiltersF_matrix = gsl_matrix_alloc(1,optimalfilter_FFT_RI->size);
		gsl_matrix_set_row(optimalfiltersF_matrix,0,optimalfilter_FFT_RI);
		if (writeFitsComplex(objFREQ,optimalfiltersF_matrix))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(objFREQ.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		gsl_vector_free(optimalfilter_FFT_RI); optimalfilter_FFT_RI = 0;
		
		strcpy(extname,"FIXFILTT");
		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
		{
			message = "Cannot move to HDU  " + string(extname) + " in library";
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
		strcpy(objTIME.nameCol,(string("T")+string(str_length)).c_str());
		optimalfiltersT_matrix = gsl_matrix_alloc(1,optimalfilter_x->size);
		gsl_matrix_set_row(optimalfiltersT_matrix,0,optimalfilter_x);
		if (writeFitsComplex(objTIME,optimalfiltersT_matrix))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(objFREQ.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		gsl_vector_free(optimalfilter_x); optimalfilter_x = 0;
				
		gsl_matrix_free(optimalfiltersT_matrix); optimalfiltersT_matrix = 0;
		gsl_matrix_free(optimalfiltersF_matrix); optimalfiltersF_matrix = 0;
	}
	gsl_vector_complex_free(optimalfilter_FFT_complex); optimalfilter_FFT_complex = 0;
		
	gsl_vector_free(optimalfilter); optimalfilter = 0;
	gsl_vector_free(optimalfilter_f); optimalfilter_f = 0;
	gsl_vector_free(optimalfilter_FFT); optimalfilter_FFT = 0;
	gsl_vector_free(matchedfilters_row); matchedfilters_row = 0;
	
	if (reconstruct_init->hduPRCLOFWM == 1)
	{
		strcpy(extname,"PRCLOFWM");
		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
		{
			message = "Cannot move to HDU  " + string(extname) + " in library";
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
		
		// Creating ENERGY Column
		strcpy(objPRCLOFWM.nameTable,"PRCLOFWM");
		strcpy(objPRCLOFWM.nameCol,"ENERGY");
		objPRCLOFWM.type = TDOUBLE;
		objPRCLOFWM.unit = new char [255];
		strcpy(objPRCLOFWM.unit,"eV");
		if (writeFitsSimple(objPRCLOFWM, E))
		{
			message = "Cannot run writeFitsSimple routine for column " + string(objPRCLOFWM.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		
		// Creating ...
		gsl_matrix *R;
		gsl_matrix *R_transW;
		gsl_matrix *R_transWR;
		gsl_vector *PULSENORMshort;	// It is necessary to use the model normalized (divided by the energy)
		gsl_vector *PULSENORM_row;
		gsl_matrix *matrixaux1 = NULL;
		gsl_matrix *aux = gsl_matrix_alloc(2,2);
		gsl_matrix *inv = gsl_matrix_alloc(2,2);
		gsl_vector *vectoraux1_2 = NULL;
		int indexPRCLOFWM;
		gsl_vector *fixedlengthsOFWM = gsl_vector_alloc(reconstruct_init->noise_spectrum->weightMatrixes->size1);
		for (int j=0;j<fixedlengthsOFWM->size;j++)	gsl_vector_set(fixedlengthsOFWM,j,pow(2,reconstruct_init->noise_spectrum->weightMatrixes->size1-j));
		int s=0;
		
		for (int j=0;j<fixedlengths->size;j++)
		{
                        //    | r0 1 | |  .   1|
                        // R =| r1 1 |=|PULSE 1|          R=(PulseLengthx2)
			//    | .    | |  .   1|
			//    | rm 1 | |  .   1|			
			if ((R = gsl_matrix_alloc(gsl_vector_get(fixedlengths,j),2)) == 0)
                        {
                                sprintf(valERROR,"%d",__LINE__-2);
                                string str(valERROR);
                                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                EP_PRINT_ERROR(message,EPFAIL);
                        }
                        //R = gsl_matrix_alloc(gsl_vector_get(fixedlengths,j),2);             	    
			gsl_matrix_set_all(R,1.0);
			PULSENORMshort = gsl_vector_alloc(gsl_vector_get(fixedlengths,j));
			if (gsl_vector_get(fixedlengths,j) == reconstruct_init->largeFilter)
			{
				PULSENORM_row = gsl_vector_alloc(PULSEMaxLengthFixedFilter->size2);
				gsl_matrix_get_row(PULSENORM_row,PULSEMaxLengthFixedFilter,0);
                                if (gsl_vector_memcpy(PULSENORMshort,PULSENORM_row) != 0)
                                {
                                        sprintf(valERROR,"%d",__LINE__-2);
                                        string str(valERROR);
                                        message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                        EP_EXIT_ERROR(message,EPFAIL);
                                }
				//gsl_vector_memcpy(PULSENORMshort,PULSENORM_row);
			}
			else
			{
				PULSENORM_row = gsl_vector_alloc(PULSE->size2);
				gsl_matrix_get_row(PULSENORM_row,PULSE,0);
                                if ((gsl_vector_get(fixedlengths,j) < 0) || (gsl_vector_get(fixedlengths,j) > PULSENORM_row->size))
                                {
                                        sprintf(valERROR,"%d",__LINE__+5);
                                        string str(valERROR);
                                        message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                        EP_EXIT_ERROR(message,EPFAIL);
                                }
				temp = gsl_vector_subvector(PULSENORM_row,0,gsl_vector_get(fixedlengths,j));
				if (gsl_vector_memcpy(PULSENORMshort,&temp.vector) != 0)
                                {
                                        sprintf(valERROR,"%d",__LINE__-2);
                                        string str(valERROR);
                                        message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                        EP_EXIT_ERROR(message,EPFAIL);
                                }   
                                //gsl_vector_memcpy(PULSENORMshort,&temp.vector);
			} 
			gsl_vector *baselinegsl = gsl_vector_alloc(PULSENORMshort->size);
			gsl_vector_set_all(baselinegsl,-1.0*reconstruct_init->noise_spectrum->baseline);
			gsl_vector_add(PULSENORMshort,baselinegsl);
			gsl_vector_free(baselinegsl); baselinegsl = 0;
			gsl_vector_scale(PULSENORMshort,1.0/reconstruct_init->monoenergy);
			gsl_matrix_set_col(R,0,PULSENORMshort);
			gsl_vector_free(PULSENORMshort); PULSENORMshort = 0;
			gsl_vector_free(PULSENORM_row); PULSENORM_row = 0;
			
			gsl_vector *Wi_vector;
			gsl_matrix *Wi_matrix;
			gsl_vector *weightMatrixes_row = gsl_vector_alloc(reconstruct_init->noise_spectrum->weightMatrixes->size2);
			gsl_matrix *PrecalOFWMaux;
			indexPRCLOFWM = 0;
			for (int i=0;i<fixedlengthsOFWM->size;i++)
			{
				if (gsl_vector_get(fixedlengths,j) == gsl_vector_get(fixedlengthsOFWM,i))
				{
					Wi_vector = gsl_vector_alloc(gsl_vector_get(fixedlengthsOFWM,i)*gsl_vector_get(fixedlengthsOFWM,i));
					Wi_matrix = gsl_matrix_alloc(gsl_vector_get(fixedlengthsOFWM,i),gsl_vector_get(fixedlengthsOFWM,i));
					gsl_matrix_get_row(weightMatrixes_row,reconstruct_init->noise_spectrum->weightMatrixes,i);
                                        if ((gsl_vector_get(fixedlengthsOFWM,i) < 0) || (gsl_vector_get(fixedlengthsOFWM,i)*gsl_vector_get(fixedlengthsOFWM,i) > weightMatrixes_row->size))
                                        {
                                                sprintf(valERROR,"%d",__LINE__+5);
                                                string str(valERROR);
                                                message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                                EP_EXIT_ERROR(message,EPFAIL);
                                        }
					temp = gsl_vector_subvector(weightMatrixes_row,0,gsl_vector_get(fixedlengthsOFWM,i)*gsl_vector_get(fixedlengthsOFWM,i));
                                        if (gsl_vector_memcpy(Wi_vector,&temp.vector) != 0)
                                        {
                                                sprintf(valERROR,"%d",__LINE__-2);
                                                string str(valERROR);
                                                message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                                EP_EXIT_ERROR(message,EPFAIL);
                                        }  
					gsl_vector_memcpy(Wi_vector,&temp.vector);
					vector2matrix(Wi_vector,&Wi_matrix);
							
					R_transW = gsl_matrix_alloc(2,gsl_vector_get(fixedlengths,j));      		// R_transW = R'W               R_transW=(2xfixedlengths_j)
					if (R->size1 != Wi_matrix->size1)
                                        {
                                                sprintf(valERROR,"%d",__LINE__+5);
                                                string str(valERROR);
                                                message = "Wrong dimensions to compute matrix-matrix product in line " + str + " (" + __FILE__ + ")";
                                                EP_EXIT_ERROR(message,EPFAIL);
                                        }
                                        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,R,Wi_matrix,0.0,R_transW);
					R_transWR = gsl_matrix_alloc(2,2);                                  		// R_transWR = R'WR             R_transWR=(2x2)
					gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,R_transW,R,0.0,R_transWR);
					matrixaux1 = gsl_matrix_alloc(2,gsl_vector_get(fixedlengths,j));
					gsl_permutation *perm1 = gsl_permutation_alloc(2);
					s=0;                                        
					gsl_matrix_memcpy(aux,R_transWR);
					gsl_linalg_LU_decomp(aux, perm1, &s);
                                        if (gsl_linalg_LU_invert(aux, perm1, inv) != 0) 
                                        {
                                                sprintf(valERROR,"%d",__LINE__-2);
                                                string str(valERROR);
                                                message = "Singular matrix in line " + str + " (" + __FILE__ + ")";
                                                EP_PRINT_ERROR(message,EPFAIL);	return(EPFAIL);
                                        }
					//gsl_linalg_LU_invert(aux, perm1, inv);
					gsl_permutation_free(perm1); perm1 = 0;
					gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,inv,R_transW,0.0,matrixaux1);      // matrixaux1 = [(R'WR)^(-1)]\B7R'W       matrixaux1=(2xfixedlengths_j)
					vectoraux1_2 = gsl_vector_alloc(gsl_vector_get(fixedlengths,j)*2);
					for (int ii=0;ii<2;ii++)
					{
						for (int k=0;k<gsl_vector_get(fixedlengths,j);k++)
						{
                                                        if ((k+ii*gsl_vector_get(fixedlengths,j) < 0) || (k+ii*gsl_vector_get(fixedlengths,j) > vectoraux1_2->size-1))
                                                        {
                                                                sprintf(valERROR,"%d",__LINE__+5);
                                                                string str(valERROR);
                                                                message = "Setting with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                                                EP_PRINT_ERROR(message,EPFAIL);
                                                        }
							gsl_vector_set(vectoraux1_2,k+ii*gsl_vector_get(fixedlengths,j),gsl_matrix_get(matrixaux1,ii,k));
						}
					}
					PrecalOFWMaux = gsl_matrix_alloc(1,gsl_vector_get(fixedlengths,j)*2);
					for (int ii=0;ii<vectoraux1_2->size;ii++)
                                        {		
                                                if ((ii+indexPRCLOFWM < 0) || (ii+indexPRCLOFWM > PrecalOFWMaux->size2-1))
                                                {
                                                        sprintf(valERROR,"%d",__LINE__+5);
                                                        string str(valERROR);
                                                        message = "Setting with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                                        EP_PRINT_ERROR(message,EPFAIL);
                                                }
                                                gsl_matrix_set(PrecalOFWMaux,0,ii+indexPRCLOFWM,gsl_vector_get(vectoraux1_2,ii));
                                            
                                        }
					snprintf(str_length,125,"%d",(int) (gsl_vector_get(fixedlengths,j)));
					strcpy(objPRCLOFWM.nameCol,(string("OFW")+string(str_length)).c_str());
					strcpy(objPRCLOFWM.unit," ");
					if (writeFitsComplex(objPRCLOFWM,PrecalOFWMaux))
					{
						message = "Cannot run writeFitsComplex routine for column " + string(objPRCLOFWM.nameCol);
						EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
					}
					gsl_matrix_free(PrecalOFWMaux); PrecalOFWMaux = 0;
					
					indexPRCLOFWM = indexPRCLOFWM + gsl_vector_get(fixedlengths,j)*2;
					
					gsl_vector_free(Wi_vector); Wi_vector = 0;
					gsl_matrix_free(Wi_matrix); Wi_matrix = 0;
					
					gsl_matrix_free(R_transW); R_transW = 0;
					gsl_matrix_free(R_transWR); R_transWR = 0;
					gsl_matrix_free(matrixaux1); matrixaux1 = 0;
					gsl_vector_free(vectoraux1_2); vectoraux1_2 = 0;
					
					break;
				}
			}
			gsl_vector_free(weightMatrixes_row); weightMatrixes_row = 0;
			gsl_matrix_free(R); R = 0;
		}
		delete [] objPRCLOFWM.unit; objPRCLOFWM.unit = 0;
	}
	
	gsl_vector_free(fixedlengths); fixedlengths = 0;
	
	// Free memory
	delete [] obj.nameTable; obj.nameTable = 0;
	delete [] obj.nameCol; obj.nameCol = 0;
	delete [] obj.unit; obj.unit = 0;
	delete [] objFREQ.nameTable; objFREQ.nameTable = 0;
	delete [] objFREQ.nameCol; objFREQ.nameCol = 0;
	delete [] objFREQ.unit; objFREQ.unit = 0;
	delete [] objTIME.nameTable; objTIME.nameTable = 0;
	delete [] objTIME.nameCol; objTIME.nameCol = 0;
	delete [] objTIME.unit; objTIME.unit = 0;
	delete [] objPRCLOFWM.nameTable; objPRCLOFWM.nameTable = 0;
	delete [] objPRCLOFWM.nameCol; objPRCLOFWM.nameCol = 0;
	
	return (EPOK);
}
/*xxxx end of SECTION A15 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A16 ************************************************************
* readAddSortParams function: This function reads the library data, add new data and sort the data according to an energy-ascending order.
* 
* - Declare variables
* - Load values already in the library
* - Add new values
* - Realign
* - Add intermeadiate values
* - Recalculate intermediate values of some new pairs
* - Write values in the library
* - Free allocated GSL vectors
* 
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
*                     In particular, this function uses 'FilterMethod', 'pulse_length', 'library_collection', 'monoenergy', 'opmode' and 'noise_spectrum'
* - inLibObject: FITS object containing information of the library FITS file
* - samprate: Sampling rate
* - eventcntLib: Number of templates in the library
* - estenergy: Pulse height of the template whose energy is going to be added to the library
* - pulsetemplate: GSL vector with the pulse template whose energy is going to be added to the library
* - covariance: GSL matrix with covariance matrix of the energy which is going to be added to the library
* - weight: GSL matrix with weight matrix of the energy which is going to be added to the library
* - pulsetemplateMaxLengthFixedFilter: GSL vector with the largeFilter-length template whose energy is going to be added to the library
******************************************************************************/
int readAddSortParams(ReconstructInitSIRENA *reconstruct_init,fitsfile **inLibObject,double samprate,int eventcntLib, double estenergy, gsl_vector *pulsetemplate,gsl_matrix *covariance, gsl_matrix *weight, gsl_vector *pulsetemplateMaxLengthFixedFilter)
{	  
	int status = EPOK;
	string message = "";
	char valERROR[256];
	
	IOData obj;
	IOData objFREQ;
	IOData objTIME;
	IOData objWN;
	IOData objOFWM;
	
	int runF0orB0val;
	if (strcmp(reconstruct_init->FilterMethod,"F0") == 0)		// Deleting the frequency-zero bin
	{
		runF0orB0val = 0;
	}
	else if (strcmp(reconstruct_init->FilterMethod,"B0") == 0)	// Working without baseline
	{
		runF0orB0val = 1;
	}
	
	// Declare variables
	// It is not necessary to check the allocation because 'eventcntLib' is already >= 1 and 'reconstruct_init->pulse_length'=PulseLength(input parameter) has been checked previously
	gsl_vector *energycolumn = gsl_vector_alloc(eventcntLib+1);
	gsl_vector *energycolumnORIGINAL = gsl_vector_alloc(eventcntLib);
	gsl_vector *estenergycolumn = gsl_vector_alloc(eventcntLib+1);
	gsl_matrix *modelsMaxLengthFixedFilteraux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->largeFilter);
	gsl_matrix *modelsaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix *modelsb0aux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix *matchedfiltersaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix *matchedfiltersb0aux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix *weightaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix *covarianceaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	
	gsl_matrix *Wabaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix_set_zero(Wabaux);
	gsl_matrix *TVaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(TVaux);
	gsl_vector *tEaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector_set_zero(tEaux);
	gsl_matrix *XMaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix_set_zero(XMaux);
	gsl_matrix *YVaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(YVaux);
	gsl_matrix *ZVaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(ZVaux);
	gsl_vector *rEaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector_set_zero(rEaux);
	gsl_matrix *Pabaux = gsl_matrix_alloc(eventcntLib+1, reconstruct_init->pulse_length);
	gsl_matrix_set_zero(Pabaux);
	gsl_matrix *PabMaxLengthFixedFilteraux = gsl_matrix_alloc(eventcntLib+1, reconstruct_init->largeFilter);
	gsl_matrix_set_zero(PabMaxLengthFixedFilteraux);
	gsl_matrix *Dabaux = gsl_matrix_alloc(eventcntLib+1, reconstruct_init->pulse_length);
	gsl_matrix_set_zero(Dabaux);
	
        gsl_matrix *optimalfiltersFREQaux;
	if ((optimalfiltersFREQaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->library_collection->optimal_filtersFREQ->ofilter_duration)) == 0)
        {
                sprintf(valERROR,"%d",__LINE__-2);
                string str(valERROR);
                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                EP_PRINT_ERROR(message,EPFAIL);
        }
        gsl_matrix *optimalfiltersTIMEaux;
        if ((optimalfiltersTIMEaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->library_collection->optimal_filtersTIME->ofilter_duration)) == 0)
        {
                sprintf(valERROR,"%d",__LINE__-2);
                string str(valERROR);
                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                EP_PRINT_ERROR(message,EPFAIL);
        }
	
	int lengthALL_F = 0;
	int lengthALL_T = 0;
	for (int i=0;i<floor(log2(reconstruct_init->pulse_length));i++)
	{
		lengthALL_F = lengthALL_F + pow(2,floor(log2(reconstruct_init->pulse_length))-i)*2;
		lengthALL_T = lengthALL_T + pow(2,floor(log2(reconstruct_init->pulse_length))-i);
	}
	if (reconstruct_init->largeFilter != reconstruct_init->pulse_length)
	{
		lengthALL_F = lengthALL_F + reconstruct_init->largeFilter*2;
		lengthALL_T = lengthALL_T + reconstruct_init->largeFilter;
	}
	gsl_matrix *optimalfiltersabFREQaux;
	if ((optimalfiltersabFREQaux = gsl_matrix_alloc(eventcntLib+1,lengthALL_F)) == 0)
	{
                sprintf(valERROR,"%d",__LINE__-2);
                string str(valERROR);
                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                EP_PRINT_ERROR(message,EPFAIL);
        }
	gsl_matrix *optimalfiltersabTIMEaux;
	if ((optimalfiltersabTIMEaux = gsl_matrix_alloc(eventcntLib+1,lengthALL_T)) == 0)
	{
                sprintf(valERROR,"%d",__LINE__-2);
                string str(valERROR);
                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                EP_PRINT_ERROR(message,EPFAIL);
        }
		
	int lengthALL_PRCLWN;
	if (reconstruct_init->largeFilter != reconstruct_init->pulse_length)	lengthALL_PRCLWN = lengthALL_F-reconstruct_init->largeFilter*2;
	else									lengthALL_PRCLWN = lengthALL_F;
        gsl_matrix *PRCLWNaux;
        if ((PRCLWNaux = gsl_matrix_alloc(eventcntLib+1,lengthALL_PRCLWN)) == 0)
        {
                sprintf(valERROR,"%d",__LINE__-2);
                string str(valERROR);
                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                EP_PRINT_ERROR(message,EPFAIL);
        }
	gsl_matrix_set_zero(PRCLWNaux);
	
	int lengthALL_PRCLOFWM = lengthALL_F;
	gsl_matrix *PRCLOFWMaux = gsl_matrix_alloc(eventcntLib+1,lengthALL_PRCLOFWM);
	gsl_matrix_set_zero(PRCLOFWMaux);
	
	gsl_vector *vectoraux = gsl_vector_alloc(1);
	gsl_vector *vectoraux1 = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_vector *vectorMaxLengthFixedFilteraux1 = gsl_vector_alloc(reconstruct_init->largeFilter);
	gsl_vector *vectoraux2 = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_vector *vectoraux3 = gsl_vector_alloc(lengthALL_PRCLWN);
	gsl_vector *vectoraux4 = gsl_vector_alloc(lengthALL_PRCLOFWM);
	
	// Load values already in the library
	for (int i=0;i<eventcntLib;i++)
	{
		gsl_vector_set(energycolumn,i,gsl_vector_get(reconstruct_init->library_collection->energies,i));
		gsl_vector_set(energycolumnORIGINAL,i,gsl_vector_get(reconstruct_init->library_collection->energies,i));
		gsl_vector_set(estenergycolumn,i,gsl_vector_get(reconstruct_init->library_collection->pulse_heights,i));
		
		if (reconstruct_init->largeFilter != reconstruct_init->pulse_length)
		{
			gsl_matrix_set_row(modelsMaxLengthFixedFilteraux,i,reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[i].ptemplate);
		}
		else 	gsl_matrix_set_row(modelsMaxLengthFixedFilteraux,i,reconstruct_init->library_collection->pulse_templates[i].ptemplate);
		gsl_matrix_set_row(modelsaux,i,reconstruct_init->library_collection->pulse_templates[i].ptemplate);
		gsl_matrix_set_row(modelsb0aux,i,reconstruct_init->library_collection->pulse_templates_B0[i].ptemplate);
		gsl_matrix_set_row(matchedfiltersaux,i,reconstruct_init->library_collection->matched_filters[i].mfilter);
		gsl_matrix_set_row(matchedfiltersb0aux,i,reconstruct_init->library_collection->matched_filters_B0[i].mfilter);
		 
                if (reconstruct_init->hduPRECALWN == 1)
                {
                        gsl_matrix_get_row(vectoraux2,reconstruct_init->library_collection->W,i);
                        gsl_matrix_set_row(weightaux,i,vectoraux2);

                        gsl_matrix_get_row(vectoraux2,reconstruct_init->library_collection->V,i);
                        gsl_matrix_set_row(covarianceaux,i,vectoraux2);
                }       
		
		if ((eventcntLib > 1) && (i < eventcntLib-1))
		{
			gsl_matrix_get_row(vectoraux1,reconstruct_init->library_collection->PAB,i);
			gsl_matrix_set_row(Pabaux,i,vectoraux1);
                        if (reconstruct_init->largeFilter != reconstruct_init->pulse_length)
                        {
                                gsl_matrix_get_row(vectorMaxLengthFixedFilteraux1,reconstruct_init->library_collection->PABMXLFF,i);
                                gsl_matrix_set_row(PabMaxLengthFixedFilteraux,i,vectorMaxLengthFixedFilteraux1);
                        }
			gsl_matrix_get_row(vectoraux1,reconstruct_init->library_collection->DAB,i);
			gsl_matrix_set_row(Dabaux,i,vectoraux1);
			gsl_matrix_set_row(optimalfiltersabFREQaux,i,reconstruct_init->library_collection->optimal_filtersabFREQ[i].ofilter);
			gsl_matrix_set_row(optimalfiltersabTIMEaux,i,reconstruct_init->library_collection->optimal_filtersabTIME[i].ofilter);
			
			if (reconstruct_init->hduPRECALWN == 1)
			{
                                gsl_matrix_get_row(vectoraux2,reconstruct_init->library_collection->WAB,i);
                                gsl_matrix_set_row(Wabaux,i,vectoraux2);
                                gsl_matrix_get_row(vectoraux1,reconstruct_init->library_collection->T,i);
                                gsl_matrix_set_row(TVaux,i,vectoraux1);
                                gsl_vector_set(tEaux,i,gsl_vector_get(reconstruct_init->library_collection->t,i));  			
                                gsl_matrix_get_row(vectoraux2,reconstruct_init->library_collection->X,i);
                                gsl_matrix_set_row(XMaux,i,vectoraux2);
                                gsl_matrix_get_row(vectoraux1,reconstruct_init->library_collection->Y,i);
                                gsl_matrix_set_row(YVaux,i,vectoraux1);
                                gsl_matrix_get_row(vectoraux1,reconstruct_init->library_collection->Z,i);
                                gsl_matrix_set_row(ZVaux,i,vectoraux1);
                                gsl_vector_set(rEaux,i,gsl_vector_get(reconstruct_init->library_collection->r,i));  
                        
				gsl_matrix_get_row(vectoraux3,reconstruct_init->library_collection->PRECALWN,i);
				gsl_matrix_set_row(PRCLWNaux,i,vectoraux3);
			}
		}
		
		gsl_matrix_set_row(optimalfiltersFREQaux,i,reconstruct_init->library_collection->optimal_filtersFREQ[i].ofilter);
		gsl_matrix_set_row(optimalfiltersTIMEaux,i,reconstruct_init->library_collection->optimal_filtersTIME[i].ofilter);
		
		if (reconstruct_init->hduPRCLOFWM == 1)
		{
			gsl_matrix_get_row(vectoraux4,reconstruct_init->library_collection->PRCLOFWM,i);
			gsl_matrix_set_row(PRCLOFWMaux,i,vectoraux4);
		}
	}
	gsl_vector_free(vectoraux3); vectoraux3 = 0;
	gsl_vector_free(vectoraux4); vectoraux4 = 0;
		
	// Add new values
	gsl_vector_set(energycolumn,eventcntLib,reconstruct_init->monoenergy);	
	gsl_vector_set(estenergycolumn,eventcntLib,estenergy);
	if (reconstruct_init->largeFilter != reconstruct_init->pulse_length)   gsl_matrix_set_row(modelsMaxLengthFixedFilteraux,eventcntLib,pulsetemplateMaxLengthFixedFilter);
	gsl_matrix_set_row(modelsaux,eventcntLib,pulsetemplate);
	
	gsl_vector_memcpy(vectoraux1,pulsetemplate);
	gsl_vector_scale(vectoraux1,1/reconstruct_init->monoenergy);
	gsl_matrix_set_row(matchedfiltersaux,eventcntLib,vectoraux1);
	
	gsl_vector *optimalfilter = NULL;
	gsl_vector *optimalfilter_f = NULL;
	gsl_vector *optimalfilter_FFT = NULL;
	gsl_vector_complex *optimalfilter_FFT_complex = NULL;
	if (calculus_optimalFilter (0, 0, reconstruct_init->opmode, vectoraux1, vectoraux1->size, samprate, runF0orB0val, reconstruct_init->noise_spectrum->noisefreqs, reconstruct_init->noise_spectrum->noisespec, &optimalfilter, &optimalfilter_f, &optimalfilter_FFT, &optimalfilter_FFT_complex))
	{
		message = "Cannot run routine calculus_optimalFilter in writeLibrary";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	gsl_vector *optimalfilter_FFT_RI;
	gsl_vector *optimalfilter_x = NULL;
	gsl_vector *optimalfilter_f_x = NULL;
	gsl_vector *optimalfilter_FFT_x = NULL;
	gsl_vector_complex *optimalfilter_FFT_complex_x = NULL;
	gsl_vector *fixedlengths;
	if (reconstruct_init->largeFilter != reconstruct_init->pulse_length)
	{
		fixedlengths = gsl_vector_alloc(floor(log2(reconstruct_init->pulse_length))+1);
		
		gsl_vector_set(fixedlengths,0,reconstruct_init->largeFilter);
		for (int i=0;i<fixedlengths->size-1;i++)
		{
                        //cout<<i+1<<" "<<pow(2,floor(log2(reconstruct_init->pulse_length))-i)<<endl;
			gsl_vector_set(fixedlengths,i+1,pow(2,floor(log2(reconstruct_init->pulse_length))-i));
		}
	}
	else
	{
		fixedlengths = gsl_vector_alloc(floor(log2(reconstruct_init->pulse_length)));
		
		gsl_vector_set(fixedlengths,0,reconstruct_init->largeFilter);
		for (int i=0;i<fixedlengths->size;i++)
		{
			gsl_vector_set(fixedlengths,i,pow(2,floor(log2(reconstruct_init->pulse_length))-i));
		}
	}
	
	gsl_vector *matchedfiltersSHORT;
	gsl_vector_view(temp);
	int indexT = 0;
	int indexF = 0;
	
	for (int j=0;j<fixedlengths->size;j++)
	{
		if (gsl_vector_get(fixedlengths,j) == optimalfilter_FFT_complex->size)
		{
			if ((optimalfilter_x = gsl_vector_alloc(optimalfilter->size)) == 0) 
                        {
                                sprintf(valERROR,"%d",__LINE__-2);
                                string str(valERROR);
                                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                EP_PRINT_ERROR(message,EPFAIL);
                        }
			gsl_vector_memcpy(optimalfilter_x,optimalfilter);
			
			optimalfilter_FFT_complex_x = gsl_vector_complex_alloc(optimalfilter_FFT_complex->size);
			gsl_vector_complex_memcpy(optimalfilter_FFT_complex_x,optimalfilter_FFT_complex);
			
			gsl_vector_free(optimalfilter_f); optimalfilter_f = 0;
			gsl_vector_free(optimalfilter_FFT); optimalfilter_FFT = 0;
		}
		else
		{
                        if ((matchedfiltersSHORT = gsl_vector_alloc(gsl_vector_get(fixedlengths,j))) == 0)
                        {
                                sprintf(valERROR,"%d",__LINE__-2);
                                string str(valERROR);
                                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                EP_PRINT_ERROR(message,EPFAIL);
                        }
			if (gsl_vector_get(fixedlengths,j) == reconstruct_init->largeFilter)
			{
				gsl_vector_scale(pulsetemplateMaxLengthFixedFilter,1/reconstruct_init->monoenergy);
                                if (gsl_vector_memcpy(matchedfiltersSHORT,pulsetemplateMaxLengthFixedFilter) != 0)
                                {
                                        sprintf(valERROR,"%d",__LINE__-2);
                                        string str(valERROR);
                                        message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                        EP_EXIT_ERROR(message,EPFAIL);
                                }  
			}
			else
			{
				temp = gsl_vector_subvector(vectoraux1,0,gsl_vector_get(fixedlengths,j));
                                if (gsl_vector_memcpy(matchedfiltersSHORT,&temp.vector) != 0)
                                {
                                        sprintf(valERROR,"%d",__LINE__-2);
                                        string str(valERROR);
                                        message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                        EP_EXIT_ERROR(message,EPFAIL);
                                }
			}
			  			
			// Calculate the optimal filter
			if (calculus_optimalFilter (0, 0, reconstruct_init->opmode, matchedfiltersSHORT, matchedfiltersSHORT->size, samprate, runF0orB0val, reconstruct_init->noise_spectrum->noisefreqs, reconstruct_init->noise_spectrum->noisespec, &optimalfilter_x, &optimalfilter_f_x, &optimalfilter_FFT_x, &optimalfilter_FFT_complex_x))
			{
				message = "Cannot run routine calculus_optimalFilter in writeLibrary";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
			gsl_vector_free(matchedfiltersSHORT); matchedfiltersSHORT = 0;
			gsl_vector_free(optimalfilter_f_x); optimalfilter_f_x = 0;
			gsl_vector_free(optimalfilter_FFT_x); optimalfilter_FFT_x = 0;
		}
		
		gsl_vector_complex_set(optimalfilter_FFT_complex_x,0,gsl_complex_rect(0.0,0.0));
		for (int i=0;i<optimalfilter_FFT_complex_x->size;i++)
		{
			gsl_vector_complex_set(optimalfilter_FFT_complex_x,i,gsl_complex_conjugate(gsl_vector_complex_get(optimalfilter_FFT_complex_x,i)));
		}
		optimalfilter_FFT_RI = gsl_vector_alloc(optimalfilter_FFT_complex_x->size*2);
		for (int i=0;i<optimalfilter_FFT_complex_x->size;i++)
		{
			gsl_vector_set(optimalfilter_FFT_RI,i,GSL_REAL(gsl_vector_complex_get(optimalfilter_FFT_complex_x,i)));
                        if ((i+optimalfilter_FFT_complex_x->size < 0) || (i+optimalfilter_FFT_complex_x->size > optimalfilter_FFT_RI->size-1))
                        {
                                sprintf(valERROR,"%d",__LINE__+5);
                                string str(valERROR);
                                message = "Setting with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                EP_PRINT_ERROR(message,EPFAIL);
                        }
			gsl_vector_set(optimalfilter_FFT_RI,i+optimalfilter_FFT_complex_x->size,GSL_IMAG(gsl_vector_complex_get(optimalfilter_FFT_complex_x,i)));
		}
		gsl_vector_complex_free(optimalfilter_FFT_complex_x); optimalfilter_FFT_complex_x = 0;
		
		for (int i=0;i<optimalfilter_FFT_RI->size;i++)
		{
                        if ((i+indexF < 0) || (i+indexF > optimalfiltersFREQaux->size2-1))
                        {
                                sprintf(valERROR,"%d",__LINE__+5);
                                string str(valERROR);
                                message = "Setting with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                EP_PRINT_ERROR(message,EPFAIL);
                        }
			gsl_matrix_set(optimalfiltersFREQaux,eventcntLib,i+indexF,gsl_vector_get(optimalfilter_FFT_RI,i));
		}
		indexF = indexF + optimalfilter_FFT_RI->size;
		
		for (int i=0;i<optimalfilter_x->size;i++)
		{
                        if ((i+indexT < 0) || (i+indexT > optimalfiltersTIMEaux->size2-1))
                        {
                                sprintf(valERROR,"%d",__LINE__+5);
                                string str(valERROR);
                                message = "Setting with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                EP_PRINT_ERROR(message,EPFAIL);
                        }
			gsl_matrix_set(optimalfiltersTIMEaux,eventcntLib,i+indexT,gsl_vector_get(optimalfilter_x,i));
		}
		indexT = indexT + optimalfilter_x->size;
		
		gsl_vector_free(optimalfilter_x); optimalfilter_x = 0;
		gsl_vector_free(optimalfilter_FFT_RI); optimalfilter_FFT_RI = 0;
	}
	gsl_vector_complex_free(optimalfilter_FFT_complex); optimalfilter_FFT_complex = 0;
        
	// It is not necessary to check the allocation because 'pulsetemplate' size must be > 0
	gsl_vector *baselinegsl = gsl_vector_alloc(pulsetemplate->size);
	gsl_vector_set_all(baselinegsl,-1.0*reconstruct_init->noise_spectrum->baseline);
	gsl_vector_add(pulsetemplate,baselinegsl);
	gsl_matrix_set_row(modelsb0aux,eventcntLib,pulsetemplate);
        gsl_vector_free(baselinegsl); baselinegsl = 0;

	gsl_vector_memcpy(vectoraux1,pulsetemplate);
	gsl_vector_scale(vectoraux1,1/reconstruct_init->monoenergy);
	gsl_matrix_set_row(matchedfiltersb0aux,eventcntLib,vectoraux1);

        if (reconstruct_init->hduPRECALWN == 1)
        {
                matrix2vector(weight,&vectoraux2);
                gsl_matrix_set_row(weightaux,eventcntLib,vectoraux2);
                matrix2vector(covariance,&vectoraux2);
                gsl_matrix_set_row(covarianceaux,eventcntLib,vectoraux2);
        }
	
	if (reconstruct_init->hduPRCLOFWM == 1)
	{
		gsl_matrix *R;
		gsl_matrix *R_transW;
		gsl_matrix *R_transWR;
		gsl_matrix *matrixaux1;
		int s = 0;
		gsl_matrix *aux = gsl_matrix_alloc(2,2);
		gsl_matrix *inv = gsl_matrix_alloc(2,2);
		gsl_vector *vectoraux1_2;
		gsl_vector *PULSENORMshort;
		gsl_vector *PULSENORM_row;
		gsl_vector *Wi_vector;
		gsl_matrix *Wi_matrix;
		gsl_vector *weightMatrixes_row = gsl_vector_alloc(reconstruct_init->noise_spectrum->weightMatrixes->size2);
		gsl_vector *fixedlengthsOFWM = gsl_vector_alloc(reconstruct_init->noise_spectrum->weightMatrixes->size1);
		int indexPRCLOFWM = 0;
		for (int j=0;j<fixedlengthsOFWM->size;j++)	gsl_vector_set(fixedlengthsOFWM,j,pow(2,reconstruct_init->noise_spectrum->weightMatrixes->size1-j));
		for (int j=0;j<fixedlengths->size;j++)
		{
                        if ((R = gsl_matrix_alloc(gsl_vector_get(fixedlengths,j),2)) == 0)
                        {
                                sprintf(valERROR,"%d",__LINE__-2);
                                string str(valERROR);
                                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                                EP_PRINT_ERROR(message,EPFAIL);
                        }
                                                                                                    //    | r0 1 | | .    1|
												    // R =| r1 1 |=|PULSE 1|          R=(PulseLengthx2)
												    //    | .    | | .    1|
												    //    | rm 1 | | .    1|			
			gsl_matrix_set_all(R,1.0);
			PULSENORMshort = gsl_vector_alloc(gsl_vector_get(fixedlengths,j));
			if (gsl_vector_get(fixedlengths,j) == reconstruct_init->largeFilter)
			{
				PULSENORM_row = gsl_vector_alloc(modelsMaxLengthFixedFilteraux->size2);
				gsl_matrix_get_row(PULSENORM_row,modelsMaxLengthFixedFilteraux,eventcntLib);
                                if (gsl_vector_memcpy(PULSENORMshort,PULSENORM_row) != 0)
                                {
                                        sprintf(valERROR,"%d",__LINE__-2);
                                        string str(valERROR);
                                        message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                        EP_EXIT_ERROR(message,EPFAIL);
                                } 
			}
			else
			{
				PULSENORM_row = gsl_vector_alloc(modelsaux->size2);
				gsl_matrix_get_row(PULSENORM_row,modelsaux,eventcntLib);
                                if ((gsl_vector_get(fixedlengths,j) < 0) || (gsl_vector_get(fixedlengths,j) > PULSENORM_row->size))
                                        {
                                                sprintf(valERROR,"%d",__LINE__+5);
                                                string str(valERROR);
                                                message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
                                                EP_EXIT_ERROR(message,EPFAIL);
                                        }
				temp = gsl_vector_subvector(PULSENORM_row,0,gsl_vector_get(fixedlengths,j));
                                if (gsl_vector_memcpy(PULSENORMshort,&temp.vector) != 0)
                                {
                                        sprintf(valERROR,"%d",__LINE__-2);
                                        string str(valERROR);
                                        message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
                                        EP_EXIT_ERROR(message,EPFAIL);
                                } 
			}  
                        baselinegsl = gsl_vector_alloc(PULSENORMshort->size);
                        gsl_vector_set_all(baselinegsl,-1.0*reconstruct_init->noise_spectrum->baseline);
			gsl_vector_add(PULSENORMshort,baselinegsl);
                        gsl_vector_free(baselinegsl); baselinegsl = 0;
			gsl_vector_scale(PULSENORMshort,1.0/reconstruct_init->monoenergy);
			gsl_matrix_set_col(R,0,PULSENORMshort);
			gsl_vector_free(PULSENORMshort); PULSENORMshort = 0;
			gsl_vector_free(PULSENORM_row); PULSENORM_row = 0;
				
			for (int i=0;i<fixedlengthsOFWM->size;i++)
			{
				if (gsl_vector_get(fixedlengths,j) == gsl_vector_get(fixedlengthsOFWM,i))
				{
					Wi_vector = gsl_vector_alloc(gsl_vector_get(fixedlengthsOFWM,i)*gsl_vector_get(fixedlengthsOFWM,i));
					Wi_matrix = gsl_matrix_alloc(gsl_vector_get(fixedlengthsOFWM,i),gsl_vector_get(fixedlengthsOFWM,i));
					gsl_matrix_get_row(weightMatrixes_row,reconstruct_init->noise_spectrum->weightMatrixes,i);
					temp = gsl_vector_subvector(weightMatrixes_row,0,gsl_vector_get(fixedlengthsOFWM,i)*gsl_vector_get(fixedlengthsOFWM,i));
					gsl_vector_memcpy(Wi_vector,&temp.vector);
					vector2matrix(Wi_vector,&Wi_matrix);
					
					R_transW = gsl_matrix_alloc(2,gsl_vector_get(fixedlengths,j));      		// R_transW = R'W               R_transW=(2xfixedlengths_j)
					if (R->size1 != Wi_matrix->size1)
                                        {
                                                sprintf(valERROR,"%d",__LINE__+5);
                                                string str(valERROR);
                                                message = "Wrong dimensions to compute matrix-matrix product in line " + str + " (" + __FILE__ + ")";
                                                EP_EXIT_ERROR(message,EPFAIL);
                                        }
					gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,R,Wi_matrix,0.0,R_transW);
					R_transWR = gsl_matrix_alloc(2,2);                                  		// R_transWR = R'WR             R_transWR=(2x2)
					gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,R_transW,R,0.0,R_transWR);
					matrixaux1 = gsl_matrix_alloc(2,gsl_vector_get(fixedlengths,j));
					gsl_permutation *perm1 = gsl_permutation_alloc(2);
					s=0;                                        
					gsl_matrix_memcpy(aux,R_transWR);
					gsl_linalg_LU_decomp(aux, perm1, &s);
					gsl_linalg_LU_invert(aux, perm1, inv);
					gsl_permutation_free(perm1); perm1 = 0;
					gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,inv,R_transW,0.0,matrixaux1);      // matrixaux1 = [(R'WR)^(-1)]\B7R'W       matrixaux1=(2xfixedlengths_j)
					vectoraux1_2 = gsl_vector_alloc(gsl_vector_get(fixedlengths,j)*2);
					for (int i=0;i<2;i++)
					{
						for (int k=0;k<gsl_vector_get(fixedlengths,j);k++)
						{
							gsl_vector_set(vectoraux1_2,k+i*gsl_vector_get(fixedlengths,j),gsl_matrix_get(matrixaux1,i,k));
						}
					}
					for (int i=0;i<vectoraux1_2->size;i++)		gsl_matrix_set(PRCLOFWMaux,eventcntLib,i+indexPRCLOFWM,gsl_vector_get(vectoraux1_2,i));
									
					indexPRCLOFWM = indexPRCLOFWM + gsl_vector_get(fixedlengths,j)*2;
					
					gsl_vector_free(Wi_vector); Wi_vector = 0;
					gsl_matrix_free(Wi_matrix); Wi_matrix = 0;
					
					gsl_matrix_free(R_transW);R_transW = 0;
					gsl_matrix_free(R_transWR); R_transWR = 0;
					gsl_matrix_free(matrixaux1); matrixaux1 = 0;
					gsl_vector_free(vectoraux1_2); vectoraux1_2 = 0;
				}
			}
			gsl_matrix_free(R); R = 0;
		}
		gsl_vector_free(fixedlengthsOFWM); fixedlengthsOFWM  = 0;
		gsl_vector_free(weightMatrixes_row); weightMatrixes_row = 0;
		gsl_matrix_free(aux); aux = 0;
		gsl_matrix_free(inv); inv = 0;
	}
	
	// Realign
	// It is not necessary to check the allocation because 'eventcntLib' is already >= 1 and 'reconstruct_init->pulse_length'=PulseLength(input parameter) 
	// and 'reconstruct_init->library_collection->optimal_filters->ofilter_duration' have been checked previously
	gsl_vector *energycolumnaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector *estenergycolumnaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector *modelsMaxLengthFixedFilterrow = gsl_vector_alloc(reconstruct_init->largeFilter);
	gsl_matrix *modelsMaxLengthFixedFilteraux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->largeFilter);
	gsl_vector *modelsrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *modelsaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_vector *modelsrowb0 = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *modelsb0aux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_vector *matchedfiltersrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *matchedfiltersaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_vector *matchedfiltersrowb0 = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *matchedfiltersb0aux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix *weightaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix *covarianceaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_vector *Wabrow = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix *Wabaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_vector *TVrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *TVaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(TVaux1);
	gsl_vector *tEcolumn = gsl_vector_alloc(eventcntLib+1);
	gsl_vector *tEcolumnaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector_set_zero(tEcolumnaux);
	gsl_vector *XMrow = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix *XMaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_vector *YVrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *YVaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(YVaux1);
	gsl_vector *ZVrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *ZVaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(ZVaux1);
	gsl_vector *rEcolumn = gsl_vector_alloc(eventcntLib+1);
	gsl_vector *rEcolumnaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector_set_zero(rEcolumnaux);
	gsl_vector *Pabrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *Pabaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(Pabaux1);
        gsl_vector *PabMaxLengthFixedFilterrow = gsl_vector_alloc(reconstruct_init->largeFilter);
	gsl_matrix *PabMaxLengthFixedFilteraux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->largeFilter);
	gsl_matrix_set_zero(PabMaxLengthFixedFilteraux1);
	gsl_vector *Dabrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *Dabaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(Dabaux1);
	gsl_vector *optimalfiltersFREQrow = gsl_vector_alloc(reconstruct_init->library_collection->optimal_filtersFREQ->ofilter_duration);
	gsl_matrix *optimalfiltersFREQaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->library_collection->optimal_filtersFREQ->ofilter_duration);
	gsl_vector *optimalfiltersTIMErow = gsl_vector_alloc(reconstruct_init->library_collection->optimal_filtersTIME->ofilter_duration);
	gsl_matrix *optimalfiltersTIMEaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->library_collection->optimal_filtersTIME->ofilter_duration);
	gsl_vector *optimalfiltersabFREQrow = gsl_vector_alloc(lengthALL_F);
	gsl_matrix *optimalfiltersabFREQaux1 = gsl_matrix_alloc(eventcntLib+1,lengthALL_F);
	gsl_vector *optimalfiltersabTIMErow = gsl_vector_alloc(lengthALL_T);
	gsl_matrix *optimalfiltersabTIMEaux1 = gsl_matrix_alloc(eventcntLib+1,lengthALL_T);
	gsl_vector *PRCLWNrow = gsl_vector_alloc(lengthALL_PRCLWN);
	gsl_matrix *PRCLWNaux1 = gsl_matrix_alloc(eventcntLib+1,lengthALL_PRCLWN);
	gsl_vector *PRCLOFWMrow = gsl_vector_alloc(lengthALL_PRCLOFWM);
	gsl_matrix *PRCLOFWMaux1 = gsl_matrix_alloc(eventcntLib+1,lengthALL_PRCLOFWM);
	
	// It is not necessary to check the allocation because 'eventcntLib' is already >= 1
	gsl_permutation *perm = gsl_permutation_alloc(eventcntLib+1);
	// 'gsl_sort_vector_index' indirectly sorts the elements of the vector v into ascending order, storing the resulting
	// permutation in p. The elements of p give the index of the vector element which would have been stored in that position
	// if the vector had been sorted in place. The first element of p gives the index of the least element in v, and the last
	// element of p gives the index of the greatest element in v. The vector v is not changed.
	// Example: tstartaux=(5200 6000 200 3000) tauxsorted=(200 3000 5200 6000) perm=(2 3 0 1)
	gsl_sort_vector_index(perm,energycolumn);
	for (int i=0;i<eventcntLib+1;i++)
	{
		gsl_vector_set(energycolumnaux,i,gsl_vector_get(energycolumn,gsl_permutation_get(perm,i)));
		gsl_vector_set(estenergycolumnaux,i,gsl_vector_get(estenergycolumn,gsl_permutation_get(perm,i)));
		gsl_matrix_get_row(modelsMaxLengthFixedFilterrow,modelsMaxLengthFixedFilteraux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(modelsMaxLengthFixedFilteraux1,i,modelsMaxLengthFixedFilterrow);
		gsl_matrix_get_row(modelsrow,modelsaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(modelsaux1,i,modelsrow);
		gsl_matrix_get_row(modelsrowb0,modelsb0aux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(modelsb0aux1,i,modelsrowb0);
		gsl_matrix_get_row(matchedfiltersrow,matchedfiltersaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(matchedfiltersaux1,i,matchedfiltersrow);
		gsl_matrix_get_row(matchedfiltersrowb0,matchedfiltersb0aux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(matchedfiltersb0aux1,i,matchedfiltersrowb0);
		
		gsl_matrix_get_row(Pabrow,Pabaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(Pabaux1,i,Pabrow);

		gsl_matrix_get_row(PabMaxLengthFixedFilterrow,PabMaxLengthFixedFilteraux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(PabMaxLengthFixedFilteraux1,i,PabMaxLengthFixedFilterrow);
		
		gsl_matrix_get_row(Dabrow,Dabaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(Dabaux1,i,Dabrow);
		
		gsl_matrix_get_row(optimalfiltersFREQrow,optimalfiltersFREQaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(optimalfiltersFREQaux1,i,optimalfiltersFREQrow);
		
		gsl_matrix_get_row(optimalfiltersTIMErow,optimalfiltersTIMEaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(optimalfiltersTIMEaux1,i,optimalfiltersTIMErow);
		
		gsl_matrix_get_row(optimalfiltersabFREQrow,optimalfiltersabFREQaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(optimalfiltersabFREQaux1,i,optimalfiltersabFREQrow);
		
		gsl_matrix_get_row(optimalfiltersabTIMErow,optimalfiltersabTIMEaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(optimalfiltersabTIMEaux1,i,optimalfiltersabTIMErow);
		
                if (reconstruct_init->hduPRECALWN == 1)
		{
                        gsl_matrix_get_row(vectoraux2,weightaux,gsl_permutation_get(perm,i));
                        gsl_matrix_set_row(weightaux1,i,vectoraux2);
                        gsl_matrix_get_row(vectoraux2,covarianceaux,gsl_permutation_get(perm,i));
                        gsl_matrix_set_row(covarianceaux1,i,vectoraux2);
                        
                        gsl_matrix_get_row(Wabrow,Wabaux,gsl_permutation_get(perm,i));
                        gsl_matrix_set_row(Wabaux1,i,Wabrow);
                        
                        gsl_matrix_get_row(TVrow,TVaux,gsl_permutation_get(perm,i));
                        gsl_matrix_set_row(TVaux1,i,TVrow);
                                        
                        gsl_vector_set(tEcolumnaux,i,gsl_vector_get(tEaux,gsl_permutation_get(perm,i)));
                        
                        gsl_matrix_get_row(XMrow,XMaux,gsl_permutation_get(perm,i));
                        gsl_matrix_set_row(XMaux1,i,XMrow);
                        
                        gsl_matrix_get_row(YVrow,YVaux,gsl_permutation_get(perm,i));
                        gsl_matrix_set_row(YVaux1,i,YVrow);
                        
                        gsl_matrix_get_row(ZVrow,ZVaux,gsl_permutation_get(perm,i));
                        gsl_matrix_set_row(ZVaux1,i,ZVrow);
                        
                        gsl_vector_set(rEcolumnaux,i,gsl_vector_get(rEaux,gsl_permutation_get(perm,i)));
                
                        gsl_matrix_get_row(PRCLWNrow,PRCLWNaux,gsl_permutation_get(perm,i));
                        gsl_matrix_set_row(PRCLWNaux1,i,PRCLWNrow);
                }
		
		if (reconstruct_init->hduPRCLOFWM == 1)
		{
			gsl_matrix_get_row(PRCLOFWMrow,PRCLOFWMaux,gsl_permutation_get(perm,i));
			gsl_matrix_set_row(PRCLOFWMaux1,i,PRCLOFWMrow);
		}
	}

	gsl_vector_memcpy(energycolumn,energycolumnaux);
	gsl_vector_memcpy(estenergycolumn,estenergycolumnaux);
	gsl_matrix_memcpy(modelsMaxLengthFixedFilteraux,modelsMaxLengthFixedFilteraux1);
	gsl_matrix_memcpy(modelsaux,modelsaux1);

	gsl_matrix_memcpy(modelsb0aux,modelsb0aux1);
	gsl_matrix_memcpy(matchedfiltersaux,matchedfiltersaux1);
	gsl_matrix_memcpy(matchedfiltersb0aux,matchedfiltersb0aux1);
	
	gsl_matrix_memcpy(weightaux,weightaux1);
	gsl_matrix_memcpy(covarianceaux,covarianceaux1);
	
	gsl_matrix_memcpy(Wabaux,Wabaux1);
	gsl_matrix_memcpy(TVaux,TVaux1);
	gsl_vector_memcpy(tEcolumn,tEcolumnaux);
	gsl_matrix_memcpy(XMaux,XMaux1);
	gsl_matrix_memcpy(YVaux,YVaux1);
	gsl_matrix_memcpy(ZVaux,ZVaux1);
	gsl_vector_memcpy(rEcolumn,rEcolumnaux);
	gsl_matrix_memcpy(Pabaux,Pabaux1);
	gsl_matrix_memcpy(PabMaxLengthFixedFilteraux,PabMaxLengthFixedFilteraux1);
	gsl_matrix_memcpy(Dabaux,Dabaux1);
	
	gsl_matrix_memcpy(optimalfiltersFREQaux,optimalfiltersFREQaux1);
	gsl_matrix_memcpy(optimalfiltersTIMEaux,optimalfiltersTIMEaux1);
	
	gsl_matrix_memcpy(optimalfiltersabFREQaux,optimalfiltersabFREQaux1);
	gsl_matrix_memcpy(optimalfiltersabTIMEaux,optimalfiltersabTIMEaux1);
	
	gsl_matrix_memcpy(PRCLWNaux,PRCLWNaux1);
	
	gsl_matrix_memcpy(PRCLOFWMaux,PRCLOFWMaux1);
			
	gsl_permutation_free(perm); perm = 0;
	gsl_vector_free(energycolumnaux); energycolumnaux = 0;
	gsl_vector_free(estenergycolumnaux); estenergycolumnaux = 0;
	gsl_vector_free(modelsMaxLengthFixedFilterrow); modelsMaxLengthFixedFilterrow = 0;
	gsl_matrix_free(modelsMaxLengthFixedFilteraux1); modelsMaxLengthFixedFilteraux1 = 0;
	gsl_vector_free(modelsrow); modelsrow = 0;
	gsl_matrix_free(modelsaux1); modelsaux1 = 0;
	gsl_vector_free(modelsrowb0); modelsrowb0 = 0;
	gsl_matrix_free(modelsb0aux1); modelsb0aux1 = 0;
	gsl_vector_free(matchedfiltersrow); matchedfiltersrow = 0;
	gsl_matrix_free(matchedfiltersaux1); matchedfiltersaux1 = 0;
	gsl_vector_free(matchedfiltersrowb0); matchedfiltersrowb0 = 0;
	gsl_matrix_free(matchedfiltersb0aux1); matchedfiltersb0aux1 = 0;
	gsl_matrix_free(weightaux1); weightaux1 = 0;
	gsl_matrix_free(covarianceaux1); covarianceaux1 = 0;
	gsl_vector_free(Wabrow); Wabrow = 0;
	gsl_matrix_free(Wabaux1); Wabaux1 = 0;
	gsl_vector_free(TVrow); TVrow = 0;
	gsl_matrix_free(TVaux1); TVaux1 = 0;
	gsl_vector_free(tEaux); tEaux = 0;
	gsl_vector_free(tEcolumnaux); tEcolumnaux = 0;
	gsl_vector_free(XMrow); XMrow = 0;
	gsl_matrix_free(XMaux1); XMaux1 = 0;
	gsl_vector_free(YVrow); YVrow = 0;
	gsl_matrix_free(YVaux1); YVaux1 = 0;
	gsl_vector_free(ZVrow); ZVrow = 0;
	gsl_matrix_free(ZVaux1); ZVaux1 = 0;
	gsl_vector_free(rEaux); rEaux = 0;
	gsl_vector_free(rEcolumnaux); rEcolumnaux = 0;
	gsl_vector_free(Pabrow); Pabrow = 0;
	gsl_matrix_free(Pabaux1); Pabaux1 = 0;
	gsl_vector_free(PabMaxLengthFixedFilterrow); PabMaxLengthFixedFilterrow = 0;
	gsl_matrix_free(PabMaxLengthFixedFilteraux1); PabMaxLengthFixedFilteraux1 = 0;
	gsl_vector_free(Dabrow); Dabrow = 0;
	gsl_matrix_free(Dabaux1); Dabaux1 = 0;
	gsl_matrix_free(optimalfiltersFREQaux1); optimalfiltersFREQaux1 = 0;
	gsl_matrix_free(optimalfiltersTIMEaux1); optimalfiltersTIMEaux1 = 0;
	gsl_matrix_free(optimalfiltersabFREQaux1); optimalfiltersabFREQaux1 = 0;
	gsl_matrix_free(optimalfiltersabTIMEaux1); optimalfiltersabTIMEaux1 = 0;
	gsl_matrix_free(PRCLWNaux1); PRCLWNaux1 = 0;
	gsl_matrix_free(PRCLOFWMaux1); PRCLOFWMaux1 = 0;
			
	// Add intermeadiate values
	int option;
	int indexa=-1, indexb=-1;
	if (reconstruct_init->monoenergy > gsl_vector_get(energycolumnORIGINAL,eventcntLib-1))	
	{	
		// Intermediate values of the last pair of energies are going to be calculated
		option = 2;
		indexa = eventcntLib-1;
		indexb = eventcntLib;
	}
	else if  (reconstruct_init->monoenergy < gsl_vector_get(energycolumnORIGINAL,0))
	{
		// Intermediate values of the first pair of energies are going to be calculated
		option = 0;
		indexa = 0;
		indexb = 1;
	}
	else
	{
		for (int i=0;i<eventcntLib-1;i++)
		{
			if ((reconstruct_init->monoenergy > gsl_vector_get(energycolumnORIGINAL,i)) && (reconstruct_init->monoenergy < gsl_vector_get(energycolumnORIGINAL,i+1)))	
			{
				// Intermediate values of the (i,i+1) pair of energies are going to be calculated
				option = 1;
				indexa = i;
				indexb = i+1;
			      
				break;
			}
		}
	}
	gsl_vector_free(energycolumnORIGINAL); energycolumnORIGINAL = 0;

	// Recalculate intermediate values of some new pairs
	if (calculateIntParams(reconstruct_init,indexa,indexb,samprate,runF0orB0val,modelsaux, covarianceaux,weightaux,energycolumn, 
		&Wabaux,&TVaux,&tEcolumn,&XMaux,&YVaux,&ZVaux,&rEcolumn, &Pabaux, &Dabaux, &PRCLWNaux, &optimalfiltersabFREQaux, &optimalfiltersabTIMEaux, 
		modelsMaxLengthFixedFilteraux, &PabMaxLengthFixedFilteraux))
	{
		message = "Cannot run calculateIntParams routine in readAddSortParams";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	
	if (option == 1)
	{
		if (calculateIntParams(reconstruct_init,indexa+1,indexb+1,samprate,runF0orB0val,modelsaux,covarianceaux,weightaux,energycolumn, 
		  	&Wabaux,&TVaux,&tEcolumn,&XMaux,&YVaux,&ZVaux,&rEcolumn, &Pabaux, &Dabaux, &PRCLWNaux, &optimalfiltersabFREQaux, &optimalfiltersabTIMEaux, 
			 modelsMaxLengthFixedFilteraux, &PabMaxLengthFixedFilteraux))
		{
			message = "Cannot run calculateIntParams routine for option 2 in readAddSortParams";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
	
	// Write values
	// It is not necessary to check the allocation because 'reconstruct_init->pulse_length'=PulseLength(input parameter) has been checked previously
	gsl_matrix *matrixaux = gsl_matrix_alloc(1,reconstruct_init->pulse_length);
	gsl_matrix *matrixaux2 = gsl_matrix_alloc(1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix *matrixaux_2 = gsl_matrix_alloc(1,reconstruct_init->pulse_length*2);
	
	obj.inObject = *inLibObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"LIBRARY");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.unit = new char [255];
	for (int i=0;i<eventcntLib+1;i++)
	{
	  	obj.iniRow = i+1;
		obj.endRow = i+1;
		strcpy(obj.nameCol,"ENERGY");
		obj.type = TDOUBLE;
		strcpy(obj.unit,"eV");
		gsl_vector_set (vectoraux,0,gsl_vector_get(energycolumn,i));
		if (writeFitsSimple(obj, vectoraux))
		{
			message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"PHEIGHT");
		strcpy(obj.unit,"ADC");
		gsl_vector_set (vectoraux,0,gsl_vector_get(estenergycolumn,i));
		if (writeFitsSimple(obj, vectoraux))
		{
			message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		gsl_matrix *matrixauxMaxLengthFixedFilter = gsl_matrix_alloc(1,reconstruct_init->largeFilter);
		if (reconstruct_init->largeFilter != reconstruct_init->pulse_length)
		{
			gsl_vector *vectoraux1MaxLengthFixedFilter = gsl_vector_alloc(reconstruct_init->largeFilter);
			strcpy(obj.nameCol,"PLSMXLFF");
			strcpy(obj.unit,"ADC");
			gsl_matrix_get_row(vectoraux1MaxLengthFixedFilter,modelsMaxLengthFixedFilteraux,i);
			gsl_matrix_set_row(matrixauxMaxLengthFixedFilter,0,vectoraux1MaxLengthFixedFilter);
			if (writeFitsComplex(obj, matrixauxMaxLengthFixedFilter))
			{
				message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			gsl_vector_free(vectoraux1MaxLengthFixedFilter); vectoraux1MaxLengthFixedFilter = 0;
		}
		
		strcpy(obj.nameCol,"PULSE");
		strcpy(obj.unit,"ADC");
		gsl_matrix_get_row(vectoraux1,modelsaux,i);
		gsl_matrix_set_row(matrixaux,0,vectoraux1);
		if (writeFitsComplex(obj, matrixaux))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"PULSEB0");
		strcpy(obj.unit,"ADC");
		gsl_matrix_get_row(vectoraux1,modelsb0aux,i);
		gsl_matrix_set_row(matrixaux,0,vectoraux1);
		if (writeFitsComplex(obj, matrixaux))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"MF");
		strcpy(obj.unit," ");
		gsl_matrix_get_row(vectoraux1,matchedfiltersaux,i);
		gsl_matrix_set_row(matrixaux,0,vectoraux1);
		if (writeFitsComplex(obj, matrixaux))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"MFB0");
		strcpy(obj.unit," ");
		gsl_matrix_get_row(vectoraux1,matchedfiltersb0aux,i);
		gsl_matrix_set_row(matrixaux,0,vectoraux1);
		if (writeFitsComplex(obj, matrixaux))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		if (reconstruct_init->hduPRECALWN == 1)
                {
                        strcpy(obj.nameCol,"COVARM");
                        strcpy(obj.unit," ");
                        gsl_matrix_get_row(vectoraux2,covarianceaux,i);
                        gsl_matrix_set_row(matrixaux2,0,vectoraux2);
                        if (writeFitsComplex(obj, matrixaux2))
                        {
                                message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                        }

                        strcpy(obj.nameCol,"WEIGHTM");
                        strcpy(obj.unit," ");
                        gsl_matrix_get_row(vectoraux2,weightaux,i);
                        gsl_matrix_set_row(matrixaux2,0,vectoraux2);
                        if (writeFitsComplex(obj, matrixaux2))
                        {
                                message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
                                EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                        }
                }
	
		if (i < eventcntLib)
		{
                        if (reconstruct_init->hduPRECALWN == 1)
                        {
                                strcpy(obj.nameCol,"WAB");
                                strcpy(obj.unit," ");
                                gsl_matrix_get_row(vectoraux2,Wabaux,i);
                                gsl_matrix_set_row(matrixaux2,0,vectoraux2);
                                if (writeFitsComplex(obj, matrixaux2))
                                {
                                        message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
                                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                
                                strcpy(obj.nameCol,"TV");
                                strcpy(obj.unit," ");
                                gsl_matrix_get_row(vectoraux1,TVaux,i);
                                gsl_matrix_set_row(matrixaux,0,vectoraux1);
                                if (writeFitsComplex(obj, matrixaux))
                                {
                                        message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
                                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                
                                strcpy(obj.nameCol,"tE");
                                strcpy(obj.unit," ");
                                gsl_vector_set (vectoraux,0,gsl_vector_get(tEcolumn,i));
                                if (writeFitsSimple(obj, vectoraux))
                                {
                                        message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
                                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                
                                strcpy(obj.nameCol,"XM");
                                strcpy(obj.unit," ");
                                gsl_matrix_get_row(vectoraux2,XMaux,i);
                                gsl_matrix_set_row(matrixaux2,0,vectoraux2);
                                if (writeFitsComplex(obj, matrixaux2))
                                {
                                        message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
                                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                
                                strcpy(obj.nameCol,"YV");
                                strcpy(obj.unit," ");
                                gsl_matrix_get_row(vectoraux1,YVaux,i);
                                gsl_matrix_set_row(matrixaux,0,vectoraux1);
                                if (writeFitsComplex(obj, matrixaux))
                                {
                                        message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
                                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                
                                strcpy(obj.nameCol,"ZV");
                                strcpy(obj.unit," ");
                                gsl_matrix_get_row(vectoraux1,ZVaux,i);
                                gsl_matrix_set_row(matrixaux,0,vectoraux1);
                                if (writeFitsComplex(obj, matrixaux))
                                {
                                        message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
                                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                                
                                strcpy(obj.nameCol,"rE");
                                strcpy(obj.unit," ");
                                gsl_vector_set (vectoraux,0,gsl_vector_get(rEcolumn,i));
                                if (writeFitsSimple(obj, vectoraux))
                                {
                                        message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
                                        EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                                }
                        }
			
			strcpy(obj.nameCol,"PAB");
			strcpy(obj.unit," ");
			gsl_matrix_get_row(vectoraux1,Pabaux,i);
			gsl_matrix_set_row(matrixaux,0,vectoraux1);
			if (writeFitsComplex(obj, matrixaux))
			{
				message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			
			if (reconstruct_init->largeFilter != reconstruct_init->pulse_length)
			{
				strcpy(obj.nameCol,"PABMXLFF");
				strcpy(obj.unit," ");
				gsl_matrix_get_row(vectorMaxLengthFixedFilteraux1,PabMaxLengthFixedFilteraux,i);
				gsl_matrix_set_row(matrixauxMaxLengthFixedFilter,0,vectorMaxLengthFixedFilteraux1);
				if (writeFitsComplex(obj, matrixauxMaxLengthFixedFilter))
				{
					message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
			}
			
			strcpy(obj.nameCol,"DAB");
			strcpy(obj.unit," ");
			gsl_matrix_get_row(vectoraux1,Dabaux,i);
			gsl_matrix_set_row(matrixaux,0,vectoraux1);
			if (writeFitsComplex(obj, matrixaux))
			{
				message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
		}
                gsl_matrix_free(matrixauxMaxLengthFixedFilter); matrixauxMaxLengthFixedFilter = 0;
	}
	
	objFREQ.inObject = *inLibObject;
	objFREQ.nameTable = new char [255];
	strcpy(objFREQ.nameTable,"FIXFILTF");
	objFREQ.iniCol = 0;
	objFREQ.nameCol = new char [255];
	objFREQ.unit = new char [255];
	objTIME.inObject = *inLibObject;
	objTIME.nameTable = new char [255];
	strcpy(objTIME.nameTable,"FIXFILTT");
	objTIME.iniCol = 0;
	objTIME.nameCol = new char [255];
	objTIME.unit = new char [255];
	objWN.inObject = *inLibObject;
	objWN.nameTable = new char [255];
	strcpy(objWN.nameTable,"PRECALWN");
	objWN.iniCol = 0;
	objWN.nameCol = new char [255];
	objWN.unit = new char [255];
	objOFWM.inObject = *inLibObject;
	objOFWM.nameTable = new char [255];
	strcpy(objOFWM.nameTable,"PRCLOFWM");
	objOFWM.iniCol = 0;
	objOFWM.nameCol = new char [255];
	objOFWM.unit = new char [255];
	char extname[20];
	char str_length[125];
	for (int i=0;i<eventcntLib+1;i++)
	{
		objFREQ.iniRow = i+1;
		objFREQ.endRow = i+1;
		objTIME.iniRow = i+1;
		objTIME.endRow = i+1;
		objWN.iniRow = i+1;
		objWN.endRow = i+1;
		objOFWM.iniRow = i+1;
		objOFWM.endRow = i+1;
	
		strcpy(extname,"FIXFILTT");
		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
		{
			EP_PRINT_ERROR("Error moving to HDU FIXFILTT in library file",status);
			return(EPFAIL);
		}
	
		strcpy(objTIME.nameCol,"ENERGY");
		objTIME.type = TDOUBLE;
		strcpy(objTIME.unit,"eV");
		gsl_vector_set (vectoraux,0,gsl_vector_get(energycolumn,i));
		if (writeFitsSimple(objTIME, vectoraux))
		{
			message = "Cannot run writeFitsSimple routine for column " + string(objTIME.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		
		strcpy(extname,"FIXFILTF");
		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
		{
			EP_PRINT_ERROR("Error moving to HDU FIXFILTF in library file",status);
			return(EPFAIL);
		}
		
		strcpy(objFREQ.nameCol,"ENERGY");
		objFREQ.type = TDOUBLE;
		strcpy(objFREQ.unit,"eV");
		gsl_vector_set (vectoraux,0,gsl_vector_get(energycolumn,i));
		if (writeFitsSimple(objFREQ, vectoraux))
		{
			message = "Cannot run writeFitsSimple routine for column " + string(objFREQ.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	
		if (reconstruct_init->hduPRECALWN == 1)
		{
			strcpy(extname,"PRECALWN");
			if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
			{
				EP_PRINT_ERROR("Error moving to HDU PRECALWN in library file",status);
				return(EPFAIL);
			}
			
			strcpy(objWN.nameCol,"ENERGY");
			objWN.type = TDOUBLE;
			strcpy(objWN.unit,"eV");
			gsl_vector_set (vectoraux,0,gsl_vector_get(energycolumn,i));
			if (writeFitsSimple(objWN, vectoraux))
			{
				message = "Cannot run writeFitsSimple routine for column " + string(objWN.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
		}


		if (reconstruct_init->hduPRCLOFWM == 1)
		{
			strcpy(extname,"PRCLOFWM");
			if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
			{
				EP_PRINT_ERROR("Error moving to HDU PRCLOFWM in library file",status);
				return(EPFAIL);
			}
			
			strcpy(objOFWM.nameCol,"ENERGY");
			objOFWM.type = TDOUBLE;
			strcpy(objOFWM.unit,"eV");
			gsl_vector_set (vectoraux,0,gsl_vector_get(energycolumn,i));
			if (writeFitsSimple(objOFWM, vectoraux))
			{
				message = "Cannot run writeFitsSimple routine for column " + string(objOFWM.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
		}
		
		indexT = 0;
		indexF = 0;
		int indexPRCLWN = 0;
		int indexPRCLOFWM = 0;
		strcpy(objTIME.unit," ");
		strcpy(objFREQ.unit," ");
		strcpy(objWN.unit," ");
		strcpy(objOFWM.unit," ");
		for (int j=0;j<reconstruct_init->library_collection->nfixedfilters;j++)		  
		{
			strcpy(extname,"FIXFILTT");
			if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
			{
				EP_PRINT_ERROR("Error moving to HDU FIXFILTT in library file",status);
				return(EPFAIL);
			}
		
			snprintf(str_length,125,"%d",(int) gsl_vector_get(fixedlengths,j));
			strcpy(objTIME.nameCol,(string("T")+string(str_length)).c_str());
			gsl_matrix_get_row(optimalfiltersTIMErow,optimalfiltersTIMEaux,i);
			gsl_matrix *optimalfiltersTIME_matrix = gsl_matrix_alloc(1,gsl_vector_get(fixedlengths,j));
			temp = gsl_vector_subvector(optimalfiltersTIMErow,indexT,gsl_vector_get(fixedlengths,j));
			gsl_matrix_set_row(optimalfiltersTIME_matrix,0,&temp.vector);
			if (writeFitsComplex(objTIME,optimalfiltersTIME_matrix))
			{
				message = "Cannot run writeFitsComplex routine for column " + string(objTIME.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			gsl_matrix_free(optimalfiltersTIME_matrix); optimalfiltersTIME_matrix = 0;
			
			if (i < eventcntLib)
			{
				strcpy(objTIME.nameCol,(string("ABT")+string(str_length)).c_str());
				gsl_matrix_get_row(optimalfiltersabTIMErow,optimalfiltersabTIMEaux,i);
				gsl_matrix *optimalfiltersabTIME_matrix = gsl_matrix_alloc(1,gsl_vector_get(fixedlengths,j));
				temp = gsl_vector_subvector(optimalfiltersabTIMErow,indexT,gsl_vector_get(fixedlengths,j));
				gsl_matrix_set_row(optimalfiltersabTIME_matrix,0,&temp.vector);
				if (writeFitsComplex(objTIME, optimalfiltersabTIME_matrix))
				{
					message = "Cannot run writeFitsComplex routine for column " + string(objTIME.nameCol);
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				gsl_matrix_free(optimalfiltersabTIME_matrix); optimalfiltersabTIME_matrix = 0;
			}
			
			indexT = indexT + gsl_vector_get(fixedlengths,j);
		
			strcpy(extname,"FIXFILTF");
			if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
			{
				EP_PRINT_ERROR("Error moving to HDU FIXFILTF in library file",status);
				return(EPFAIL);
			}
		
			strcpy(objFREQ.nameCol,(string("F")+string(str_length)).c_str());
			gsl_matrix_get_row(optimalfiltersFREQrow,optimalfiltersFREQaux,i);
			gsl_matrix *optimalfiltersFREQ_matrix = gsl_matrix_alloc(1,gsl_vector_get(fixedlengths,j)*2);
			temp = gsl_vector_subvector(optimalfiltersFREQrow,indexF,gsl_vector_get(fixedlengths,j)*2);
			gsl_matrix_set_row(optimalfiltersFREQ_matrix,0,&temp.vector);
			if (writeFitsComplex(objFREQ,optimalfiltersFREQ_matrix))
			{
				message = "Cannot run writeFitsComplex routine for column " + string(objFREQ.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			gsl_matrix_free(optimalfiltersFREQ_matrix); optimalfiltersFREQ_matrix = 0;
			
			if (i < eventcntLib)
			{
				strcpy(objFREQ.nameCol,(string("ABF")+string(str_length)).c_str());
				gsl_matrix_get_row(optimalfiltersabFREQrow,optimalfiltersabFREQaux,i);
				gsl_matrix *optimalfiltersabFREQ_matrix = gsl_matrix_alloc(1,gsl_vector_get(fixedlengths,j)*2);
				temp = gsl_vector_subvector(optimalfiltersabFREQrow,indexF,gsl_vector_get(fixedlengths,j)*2);
				gsl_matrix_set_row(optimalfiltersabFREQ_matrix,0,&temp.vector);
				if (writeFitsComplex(objFREQ, optimalfiltersabFREQ_matrix))
				{
					message = "Cannot run writeFitsComplex routine for column " + string(objFREQ.nameCol);
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				gsl_matrix_free(optimalfiltersabFREQ_matrix); optimalfiltersabFREQ_matrix = 0;
			}
			
			indexF = indexF + gsl_vector_get(fixedlengths,j)*2;
			
			if (reconstruct_init->hduPRECALWN == 1)
			{
				if (reconstruct_init->largeFilter != reconstruct_init->pulse_length)
				{
					if (j >= 1)
					{  
						if (i < eventcntLib)
						{
							strcpy(extname,"PRECALWN");
							if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
							{
								EP_PRINT_ERROR("Error moving to HDU PRECALWN in library file",status);
								return(EPFAIL);
							}
						
							strcpy(objWN.nameCol,(string("PCL")+string(str_length)).c_str());
							gsl_matrix_get_row(PRCLWNrow,PRCLWNaux,i);
							gsl_matrix *PRCLWN_matrix;
							PRCLWN_matrix = gsl_matrix_alloc(1,gsl_vector_get(fixedlengths,j)*2);
							temp = gsl_vector_subvector(PRCLWNrow,indexPRCLWN,gsl_vector_get(fixedlengths,j)*2);
							gsl_matrix_set_row(PRCLWN_matrix,0,&temp.vector);
							if (writeFitsComplex(objWN,PRCLWN_matrix))
							{
								message = "Cannot run writeFitsComplex routine for column " + string(objWN.nameCol);
								EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
							}
							gsl_matrix_free(PRCLWN_matrix); PRCLWN_matrix = 0;
						}
						
						indexPRCLWN = indexPRCLWN + gsl_vector_get(fixedlengths,j)*2;
					}
				}
				else
				{
					if (i < eventcntLib)
					{
						strcpy(extname,"PRECALWN");
						if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
						{
							EP_PRINT_ERROR("Error moving to HDU PRECALWN in library file",status);
							return(EPFAIL);
						}
					
						strcpy(objWN.nameCol,(string("PCL")+string(str_length)).c_str());
						gsl_matrix_get_row(PRCLWNrow,PRCLWNaux,i);
						gsl_matrix *PRCLWN_matrix;
						PRCLWN_matrix = gsl_matrix_alloc(1,gsl_vector_get(fixedlengths,j)*2);
						temp = gsl_vector_subvector(PRCLWNrow,indexPRCLWN,gsl_vector_get(fixedlengths,j)*2);
						gsl_matrix_set_row(PRCLWN_matrix,0,&temp.vector);
						if (writeFitsComplex(objWN,PRCLWN_matrix))
						{
							message = "Cannot run writeFitsComplex routine for column " + string(objWN.nameCol);
							EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
						}
						gsl_matrix_free(PRCLWN_matrix); PRCLWN_matrix = 0;
					}
					
					indexPRCLWN = indexPRCLWN + gsl_vector_get(fixedlengths,j)*2;
				}
			}
			
			if (reconstruct_init->hduPRCLOFWM == 1)
			{
				strcpy(extname,"PRCLOFWM");
				if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
				{
					EP_PRINT_ERROR("Error moving to HDU PRCLOFWM in library file",status);
					return(EPFAIL);
				}
			
				strcpy(objOFWM.nameCol,(string("OFW")+string(str_length)).c_str());
				gsl_matrix_get_row(PRCLOFWMrow,PRCLOFWMaux,i);
				gsl_matrix *PRCLOFWM_matrix;
				PRCLOFWM_matrix = gsl_matrix_alloc(1,gsl_vector_get(fixedlengths,j)*2);
				temp = gsl_vector_subvector(PRCLOFWMrow,indexPRCLOFWM,gsl_vector_get(fixedlengths,j)*2);
				gsl_matrix_set_row(PRCLOFWM_matrix,0,&temp.vector);
				if (writeFitsComplex(objOFWM,PRCLOFWM_matrix))
				{
					message = "Cannot run writeFitsComplex routine for column " + string(objOFWM.nameCol);
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				gsl_matrix_free(PRCLOFWM_matrix); PRCLOFWM_matrix = 0;
				
				indexPRCLOFWM = indexPRCLOFWM + gsl_vector_get(fixedlengths,j)*2;
			}
		}
	}
	
	// Free allocated GSL vectors
	gsl_vector_free(optimalfilter); optimalfilter = 0;
		
	gsl_vector_free(vectoraux); vectoraux = 0;
	gsl_vector_free(vectoraux1); vectoraux1 = 0;
	gsl_vector_free(vectorMaxLengthFixedFilteraux1); vectorMaxLengthFixedFilteraux1 = 0;
	gsl_vector_free(vectoraux2); vectoraux2 = 0;
	gsl_matrix_free(matrixaux); matrixaux = 0;
	gsl_matrix_free(matrixaux2); matrixaux2 = 0;
	gsl_matrix_free(matrixaux_2); matrixaux_2 = 0;
		
	gsl_vector_free(energycolumn); energycolumn = 0;
	gsl_vector_free(estenergycolumn); estenergycolumn = 0;
	gsl_matrix_free(modelsMaxLengthFixedFilteraux); modelsMaxLengthFixedFilteraux = 0;
	gsl_matrix_free(modelsaux);modelsaux = 0;
	gsl_matrix_free(modelsb0aux); modelsb0aux = 0;
	gsl_matrix_free(matchedfiltersaux); matchedfiltersaux = 0;
	gsl_matrix_free(matchedfiltersb0aux); matchedfiltersb0aux = 0;
	gsl_matrix_free(weightaux); weightaux = 0;
	gsl_matrix_free(covarianceaux); covarianceaux = 0;
	gsl_matrix_free(Wabaux); Wabaux = 0;
	gsl_matrix_free(TVaux); TVaux = 0;
	gsl_vector_free(tEcolumn); tEcolumn = 0;
	gsl_matrix_free(XMaux); XMaux = 0;
	gsl_matrix_free(YVaux); YVaux = 0;
	gsl_matrix_free(ZVaux); ZVaux = 0;
	gsl_vector_free(rEcolumn); rEcolumn = 0;
	gsl_matrix_free(Pabaux); Pabaux = 0;
        gsl_matrix_free(PabMaxLengthFixedFilteraux); PabMaxLengthFixedFilteraux = 0;
	gsl_matrix_free(Dabaux); Dabaux = 0;
	gsl_matrix_free(optimalfiltersFREQaux);optimalfiltersFREQaux = 0;
	gsl_vector_free(optimalfiltersFREQrow); optimalfiltersFREQrow = 0;
	gsl_matrix_free(optimalfiltersTIMEaux); optimalfiltersTIMEaux = 0;
	gsl_vector_free(optimalfiltersTIMErow); optimalfiltersTIMErow = 0;
	gsl_matrix_free(optimalfiltersabFREQaux); optimalfiltersabFREQaux = 0;
	gsl_vector_free(optimalfiltersabFREQrow); optimalfiltersabFREQrow = 0;
	gsl_matrix_free(optimalfiltersabTIMEaux); optimalfiltersabTIMEaux = 0;
	gsl_vector_free(optimalfiltersabTIMErow); optimalfiltersabTIMErow = 0;
	gsl_matrix_free(PRCLWNaux); PRCLWNaux = 0;
	gsl_vector_free(PRCLWNrow); PRCLWNrow = 0;
	gsl_matrix_free(PRCLOFWMaux); PRCLOFWMaux = 0;
	gsl_vector_free(PRCLOFWMrow); PRCLOFWMrow = 0;
	
	gsl_vector_free(fixedlengths); fixedlengths = 0;
	
	delete [] obj.nameTable; obj.nameTable = 0;
	delete [] obj.nameCol; obj.nameCol = 0;
	delete [] obj.unit; obj.unit = 0;
	
	delete [] objFREQ.nameTable; objFREQ.nameTable = 0;
	delete [] objFREQ.nameCol; objFREQ.nameCol = 0;
	delete [] objFREQ.unit; objFREQ.unit = 0;
	
	delete [] objTIME.nameTable; objTIME.nameTable = 0;
	delete [] objTIME.nameCol; objTIME.nameCol = 0;
	delete [] objTIME.unit; objTIME.unit = 0;
	
	delete [] objWN.nameTable; objWN.nameTable = 0;
	delete [] objWN.nameCol; objWN.nameCol = 0;
	delete [] objWN.unit; objWN.unit = 0;
	
	delete [] objOFWM.nameTable; objOFWM.nameTable = 0;
	delete [] objOFWM.nameCol; objOFWM.nameCol = 0;
	delete [] objOFWM.unit; objOFWM.unit = 0;
		
	return (EPOK);
}
/*xxxx end of SECTION A16 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A17 ************************************************************
* calculateIntParams function: This function calculates some intermediate scalars, vectors and matrices (WAB, TV, tE, XM, YV, ZV, rE, PAB and DAB)
*                              for the interpolation and covariance methods.
*
* - Declare variables and allocate GSL vectors and matrices
* - Calculate intermediate scalars, vectors and matrices 
* - Free allocated GSL vectors and matrices
* 
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
*                     In particular, this function uses 'pulse_length'
* - indexa: Lower index of the library to calculate the intermediate params
* - indexb: Higher index of the library to calculate the intermediate params
* - samprate: Sampling rate
* - runF0orB0val: 'FilterMethod' = F0 => 'runF0orB0val' = 1
*                 'FilterMethod' = B0 => 'runF0orB0val' = 0
* - modelsaux: GSL input matrix with model template
* - covarianceaux: GSL input matrix with covariance matrix  
* - weightaux: GSL input matrix with weight matrix 
* - energycolumn: GSL input vector with list of energies
* - Wabaux: Input/output intermediate parameters
* - TVaux: Input/output intermediate parameters
* - tEcolumn: Input/output intermediate parameters
* - XMaux: Input/output intermediate parameters
* - YVaux: Input/output intermediate parameters
* - ZVaux: Input/output intermediate parameters
* - rEcolumn: Input/output intermediate parameters
* - Pabaux: Input/output intermediate parameters
* - Dabaux: Input/output intermediate parameters
* - precalWMaux: Input/output intermediate parameters 
* - optimalfiltersabFREQaux: Input/output intermediate parameters
* - optimalfiltersabTIMEaux: Input/output intermediate parameters
* - modelsMaxLengthFixedFilteraux: GSL input matrix with model template whose length is largeFilter
* - PabMaxLengthFixedFilteraux: Input/output intermediate parameters
******************************************************************************/
int calculateIntParams(ReconstructInitSIRENA *reconstruct_init, int indexa, int indexb, double samprate, int runF0orB0val, 
	gsl_matrix *modelsaux,gsl_matrix *covarianceaux, gsl_matrix *weightaux,gsl_vector *energycolumn, 
	gsl_matrix **Wabaux,gsl_matrix **TVaux,gsl_vector **tEcolumn,gsl_matrix **XMaux,gsl_matrix **YVaux,gsl_matrix **ZVaux,
	gsl_vector **rEcolumn, gsl_matrix **Pabaux, gsl_matrix **Dabaux, gsl_matrix **PrecalWMaux, gsl_matrix **optimalfiltersabFREQaux, gsl_matrix **optimalfiltersabTIMEaux, gsl_matrix *modelsMaxLengthFixedFilteraux, gsl_matrix **PabMaxLengthFixedFilteraux)
{		
	int status = EPOK;
        string message = "";
	char valERROR[256];
	
	// Declare variables and allocate GSL vectors and matrices
	// It is not necessary to check the allocation because 'reconstruct_init->pulse_length'=PulseLength(input parameter) has been checked previously
	gsl_vector *T;
        gsl_vector *TMaxLengthFixedFilter;
        double t;
        gsl_vector *Y;
        gsl_vector *Z;
        double r;
	gsl_vector *Pab = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_vector *PabMaxLengthFixedFilter = gsl_vector_alloc(reconstruct_init->largeFilter);
	gsl_vector *Dab = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_vector *DabMaxLengthFixedFilter = gsl_vector_alloc(reconstruct_init->largeFilter);
        gsl_vector *Wabvector;
        gsl_matrix *Walpha;
        gsl_vector *Walphavector;
        gsl_vector *Wbetavector;
        gsl_vector *Xvector;
        gsl_matrix *X;
	double Ea,Eb;
	
	gsl_vector *vectoraux = gsl_vector_alloc(1);
	gsl_vector *vectoraux1 = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_vector *vectoraux2 = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix *matrixaux = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);
	
	gsl_permutation *perm;
	int s=0;
	
	gsl_vector *optimalfilter = NULL;
	gsl_vector *optimalfilter_f = NULL;
	gsl_vector *optimalfilter_FFT = NULL;
	gsl_vector_complex *optimalfilter_FFT_complex = NULL;
        
        gsl_vector *vectorMaxLengthFixedFilteraux1 = gsl_vector_alloc(reconstruct_init->largeFilter);
     
	// Calculate
        if (reconstruct_init->hduPRECALWN == 1)
        {
                T = gsl_vector_alloc(reconstruct_init->pulse_length);
                TMaxLengthFixedFilter = gsl_vector_alloc(reconstruct_init->largeFilter);
                Y = gsl_vector_alloc(reconstruct_init->pulse_length);
                Z = gsl_vector_alloc(reconstruct_init->pulse_length);
                Wabvector = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
                Walpha = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);
                Walphavector = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
                Wbetavector = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
                Xvector = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
                X = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);
                
                if ((indexa < 0) || (indexa > weightaux->size1-1))
                {
                        sprintf(valERROR,"%d",__LINE__+5);
                        string str(valERROR);
                        message = "Getting i-th row of matrix out of range in line " + str + " (" + __FILE__ + ")";
                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                }
                gsl_matrix_get_row(Walphavector,weightaux,indexa);
                gsl_vector_memcpy(Wabvector,Walphavector);
                if ((indexb < 0) || (indexb > weightaux->size1-1))
                {
                        sprintf(valERROR,"%d",__LINE__+5);
                        string str(valERROR);
                        message = "Getting i-th row of matrix out of range in line " + str + " (" + __FILE__ + ")";
                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                }
                gsl_matrix_get_row(Wbetavector,weightaux,indexb);
                gsl_vector_add(Wabvector,Wbetavector);
                gsl_vector_scale(Wabvector,1.0/2);
                if ((indexa < 0) || (indexa > (*Wabaux)->size1-1))
                {
                        sprintf(valERROR,"%d",__LINE__+5);
                        string str(valERROR);
                        message = "Setting i-th row of matrix out of range in line " + str + " (" + __FILE__ + ")";
                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                }
                gsl_matrix_set_row(*Wabaux,indexa,Wabvector);
	
                gsl_matrix_get_row(vectoraux1,modelsaux,indexb);
                gsl_vector_memcpy(T,vectoraux1);
                gsl_matrix_get_row(vectoraux1,modelsaux,indexa);
                gsl_vector_sub(T,vectoraux1);
                gsl_matrix_set_row(*TVaux,indexa,T);
        	
                gsl_matrix_get_row(vectorMaxLengthFixedFilteraux1,modelsMaxLengthFixedFilteraux,indexb);
                gsl_vector_memcpy(TMaxLengthFixedFilter,vectorMaxLengthFixedFilteraux1);
                gsl_matrix_get_row(vectorMaxLengthFixedFilteraux1,modelsMaxLengthFixedFilteraux,indexa);
                gsl_vector_sub(TMaxLengthFixedFilter,vectorMaxLengthFixedFilteraux1);
                
                vector2matrix(Walphavector,&Walpha);
                gsl_blas_dgemv(CblasNoTrans, 1.0, Walpha, T, 0.0, vectoraux1);
                gsl_blas_ddot(T,vectoraux1,&t);
                if ((indexa < 0) || (indexa > (*tEcolumn)->size-1))
                {
                        sprintf(valERROR,"%d",__LINE__+5);
                        string str(valERROR);
                        message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                }
                gsl_vector_set(*tEcolumn,indexa,t);
                
                gsl_vector_memcpy(Xvector,Wbetavector);
                gsl_vector_sub(Xvector,Walphavector);
                gsl_vector_scale(Xvector,1/t);
                gsl_matrix_set_row(*XMaux,indexa,Xvector);
                
                gsl_blas_dgemv(CblasNoTrans, 1.0, Walpha, T, 0.0, Y);
                gsl_vector_scale(Y,1/t);
                gsl_matrix_set_row(*YVaux,indexa,Y);
                gsl_matrix_free(Walpha); Walpha = 0;
                
                vector2matrix(Xvector,&X);
                gsl_blas_dgemv(CblasNoTrans, 1.0, X, T, 0.0, Z);
                gsl_matrix_set_row(*ZVaux,indexa,Z);
                gsl_matrix_free(X); X = 0;
                
                gsl_blas_ddot(Z,T,&r);
                r=1/r;
                if ((indexa < 0) || (indexa > (*rEcolumn)->size-1))
                {
                        sprintf(valERROR,"%d",__LINE__+5);
                        string str(valERROR);
                        message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                }
                gsl_vector_set(*rEcolumn,indexa,r);
                
                gsl_vector_free(T); T = 0;
                gsl_vector_free(TMaxLengthFixedFilter); TMaxLengthFixedFilter = 0;
                gsl_vector_free(Y); Y = 0;
                gsl_vector_free(Z); Z = 0;
                //gsl_vector_free(Wabvector);
                gsl_vector_free(Walphavector); Walphavector = 0;
                gsl_vector_free(Wbetavector); Wbetavector = 0;
                gsl_vector_free(Xvector); Xvector = 0;
        }
        	
	if ((indexa < 0) || (indexa > energycolumn->size-1))
	{
		sprintf(valERROR,"%d",__LINE__+5);
		string str(valERROR);
		message = "Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	Ea = gsl_vector_get(energycolumn,indexa);
	if ((indexb < 0) || (indexb > energycolumn->size-1))
	{
		sprintf(valERROR,"%d",__LINE__+5);
		string str(valERROR);
		message = "Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	Eb = gsl_vector_get(energycolumn,indexb);
        gsl_matrix_get_row(vectoraux1,modelsaux,indexb);
        gsl_vector_memcpy(Pab,vectoraux1);
        gsl_matrix_get_row(vectoraux1,modelsaux,indexa);
        gsl_vector_sub(Pab,vectoraux1);
	gsl_vector_scale(Pab,-Ea/(Eb-Ea));
	gsl_matrix_get_row(vectoraux1,modelsaux,indexa);
	gsl_vector_add(Pab,vectoraux1);
	gsl_vector_add_constant(Pab,-1*reconstruct_init->noise_spectrum->baseline);
	gsl_matrix_set_row(*Pabaux,indexa,Pab);
	
        gsl_matrix_get_row(vectorMaxLengthFixedFilteraux1,modelsMaxLengthFixedFilteraux,indexb);
        gsl_vector_memcpy(PabMaxLengthFixedFilter,vectorMaxLengthFixedFilteraux1);
        gsl_matrix_get_row(vectorMaxLengthFixedFilteraux1,modelsMaxLengthFixedFilteraux,indexa);
        gsl_vector_sub(PabMaxLengthFixedFilter,vectorMaxLengthFixedFilteraux1);
	gsl_vector_scale(PabMaxLengthFixedFilter,-Ea/(Eb-Ea));
	gsl_matrix_get_row(vectorMaxLengthFixedFilteraux1,modelsMaxLengthFixedFilteraux,indexa);
	gsl_vector_add(PabMaxLengthFixedFilter,vectorMaxLengthFixedFilteraux1);
	gsl_vector_add_constant(PabMaxLengthFixedFilter,-1*reconstruct_init->noise_spectrum->baseline);
	gsl_matrix_set_row(*PabMaxLengthFixedFilteraux,indexa,PabMaxLengthFixedFilter);
	
        gsl_matrix_get_row(vectorMaxLengthFixedFilteraux1,modelsMaxLengthFixedFilteraux,indexb);
        gsl_vector_memcpy(DabMaxLengthFixedFilter,vectorMaxLengthFixedFilteraux1);
        gsl_matrix_get_row(vectorMaxLengthFixedFilteraux1,modelsMaxLengthFixedFilteraux,indexa);
        gsl_vector_sub(DabMaxLengthFixedFilter,vectorMaxLengthFixedFilteraux1);
        gsl_vector_free(vectorMaxLengthFixedFilteraux1); vectorMaxLengthFixedFilteraux1 = 0;
	//gsl_vector_memcpy(DabMaxLengthFixedFilter,TMaxLengthFixedFilter);
	gsl_vector_scale(DabMaxLengthFixedFilter,1/(Eb-Ea));
        gsl_matrix_get_row(vectoraux1,modelsaux,indexb);
        gsl_vector_memcpy(Dab,vectoraux1);
        gsl_matrix_get_row(vectoraux1,modelsaux,indexa);
        gsl_vector_sub(Dab,vectoraux1);
	gsl_vector_scale(Dab,1/(Eb-Ea));
	gsl_matrix_set_row(*Dabaux,indexa,Dab);
	
        if (reconstruct_init->hduPRECALWN == 1)
        {
                if (calculus_optimalFilter (0, 0, 0, Dab, Dab->size, samprate, runF0orB0val, reconstruct_init->noise_spectrum->noisefreqs, reconstruct_init->noise_spectrum->noisespec, &optimalfilter, &optimalfilter_f, &optimalfilter_FFT, &optimalfilter_FFT_complex))
                {
                        message = "Cannot run routine calculus_optimalFilter in writeLibrary";
                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                }
                
                gsl_vector *optimalfilter_x = NULL;
                gsl_vector *optimalfilter_f_x = NULL;
                gsl_vector *optimalfilter_FFT_x = NULL;
                gsl_vector_complex *optimalfilter_FFT_complex_x = NULL;
                gsl_vector *fixedlengths;
                if (reconstruct_init->largeFilter != reconstruct_init->pulse_length)
                {
                        fixedlengths = gsl_vector_alloc(floor(log2(reconstruct_init->pulse_length))+1);
                        
                        gsl_vector_set(fixedlengths,0,reconstruct_init->largeFilter);
                        for (int i=0;i<fixedlengths->size-1;i++)
                        {
                                gsl_vector_set(fixedlengths,i+1,pow(2,floor(log2(reconstruct_init->pulse_length))-i));
                        }
                }
                else
                {
                        fixedlengths = gsl_vector_alloc(floor(log2(reconstruct_init->pulse_length)));
                        
                        for (int i=0;i<fixedlengths->size;i++)
                        {
                                gsl_vector_set(fixedlengths,i,pow(2,floor(log2(reconstruct_init->pulse_length))-i));
                        }
                }
	
                int startSize = optimalfilter_FFT_complex->size;
                gsl_vector *DabsSHORT;
                gsl_vector_view(temp);
                int indexT = 0;
                int indexF = 0;
                int indexPRCLWN = 0;
                gsl_matrix *R = NULL;
                gsl_matrix *R_transW = NULL;
                gsl_matrix *R_transWR = 0;//gsl_matrix_alloc(2,2);
                gsl_matrix *matrixaux1 = NULL;
                gsl_matrix *aux = gsl_matrix_alloc(2,2);
                gsl_matrix *inv = gsl_matrix_alloc(2,2);
                gsl_vector *vectoraux1_2 = NULL;
                gsl_matrix_view tempm;
                for (int j=0;j<fixedlengths->size;j++)
                {
                        if (gsl_vector_get(fixedlengths,j) == startSize)
                        {
                                optimalfilter_x = gsl_vector_alloc(optimalfilter->size);
                                gsl_vector_memcpy(optimalfilter_x,optimalfilter);
                                
                                optimalfilter_FFT_complex_x = gsl_vector_complex_alloc(optimalfilter_FFT_complex->size);
                                gsl_vector_complex_memcpy(optimalfilter_FFT_complex_x,optimalfilter_FFT_complex);
                                gsl_vector_complex_free(optimalfilter_FFT_complex); optimalfilter_FFT_complex = 0;
                                
                                gsl_vector_free(optimalfilter_f); optimalfilter_f = 0;
                                gsl_vector_free(optimalfilter_FFT); optimalfilter_FFT = 0;
                        }
                        else
                        {
                                DabsSHORT = gsl_vector_alloc(gsl_vector_get(fixedlengths,j));
                                if (gsl_vector_get(fixedlengths,j) == reconstruct_init->largeFilter)
                                {
                                        gsl_vector_memcpy(DabsSHORT,DabMaxLengthFixedFilter);
                                }
                                else
                                {
                                        temp = gsl_vector_subvector(Dab,0,gsl_vector_get(fixedlengths,j));
                                        gsl_vector_memcpy(DabsSHORT,&temp.vector);
                                }

                                // Calculate the optimal filter
                                optimalfilter_x = 0;
                                optimalfilter_f_x = 0;
                                optimalfilter_FFT_x = 0;
                                optimalfilter_FFT_complex_x = 0;
                                if (calculus_optimalFilter (0, 0, reconstruct_init->opmode, DabsSHORT, DabsSHORT->size, samprate, runF0orB0val, reconstruct_init->noise_spectrum->noisefreqs, reconstruct_init->noise_spectrum->noisespec, &optimalfilter_x, &optimalfilter_f_x, &optimalfilter_FFT_x, &optimalfilter_FFT_complex_x))
                                {
                                        message = "Cannot run routine calculus_optimalFilter in writeLibrary";
                                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                }
                                gsl_vector_free(DabsSHORT); DabsSHORT = 0;
                                gsl_vector_free(optimalfilter_f_x); optimalfilter_f_x = 0;
                                gsl_vector_free(optimalfilter_FFT_x); optimalfilter_FFT_x = 0;
                        }
                        
                        gsl_vector_complex_set(optimalfilter_FFT_complex_x,0,gsl_complex_rect(0.0,0.0));
                        for (int i=0;i<optimalfilter_FFT_complex_x->size;i++)
                        {
                                gsl_vector_complex_set(optimalfilter_FFT_complex_x,i,gsl_complex_conjugate(gsl_vector_complex_get(optimalfilter_FFT_complex_x,i)));
                        }
                        
                        for (int i=0;i<optimalfilter_FFT_complex_x->size;i++)
                        {
                                gsl_matrix_set(*optimalfiltersabFREQaux,indexa,i+indexF,GSL_REAL(gsl_vector_complex_get(optimalfilter_FFT_complex_x,i)));
                                
                                gsl_matrix_set(*optimalfiltersabFREQaux,indexa,i+indexF+optimalfilter_FFT_complex_x->size,GSL_IMAG(gsl_vector_complex_get(optimalfilter_FFT_complex_x,i)));
                        }
                
                        gsl_vector_complex_free(optimalfilter_FFT_complex_x); optimalfilter_FFT_complex_x = 0;
                        
                        indexF = indexF + gsl_vector_get(fixedlengths,j)*2;
                        
                        for (int i=0;i<optimalfilter_x->size;i++)
                        {
                                gsl_matrix_set(*optimalfiltersabTIMEaux,indexa,i+indexT,gsl_vector_get(optimalfilter_x,i));
                        }
                
                        indexT = indexT + optimalfilter_x->size;
                        
                        gsl_vector_free(optimalfilter_x); optimalfilter_x = 0;
                        
                        if (((reconstruct_init->largeFilter != reconstruct_init->pulse_length) && (j != 0)) || (reconstruct_init->largeFilter == reconstruct_init->pulse_length))
                        {
                                R = gsl_matrix_alloc(gsl_vector_get(fixedlengths,j),2);             // R=(PulseLengthx2)
                                                                                                    //    | r0 1 | | .  1|
                                                                                                    // R =| r1 1 |=|Dab 1|          
                                                                                                    //    | .    | | .  1|
                                                                                                    //    | rm 1 | | .  1|
            
                                gsl_matrix_set_all(R,1.0);
                                if (matrixaux) gsl_matrix_free(matrixaux); matrixaux = 0;
                                if (gsl_vector_get(fixedlengths,j) == startSize)
                                {
                                        DabsSHORT = gsl_vector_alloc(Dab->size);
                                        gsl_vector_memcpy(DabsSHORT,Dab);
                                        
                                        matrixaux = gsl_matrix_alloc(gsl_vector_get(fixedlengths,j),gsl_vector_get(fixedlengths,j));
                                        vector2matrix(Wabvector,&matrixaux);
                                        gsl_vector_free(Wabvector); Wabvector = 0;
                                }
                                else
                                {
                                        DabsSHORT = gsl_vector_alloc(gsl_vector_get(fixedlengths,j));
                                        temp = gsl_vector_subvector(Dab,0,gsl_vector_get(fixedlengths,j));
                                        gsl_vector_memcpy(DabsSHORT,&temp.vector);
                                        
                                        gsl_vector *covarianceAlphavector = gsl_vector_alloc(covarianceaux->size2);
                                        gsl_vector *covarianceBetavector = gsl_vector_alloc(covarianceaux->size2);
                                        gsl_matrix *covarianceAlphamatrix = gsl_matrix_alloc(sqrt(covarianceaux->size2),sqrt(covarianceaux->size2));
                                        gsl_matrix *covarianceBetamatrix = gsl_matrix_alloc(sqrt(covarianceaux->size2),sqrt(covarianceaux->size2));
                                        gsl_matrix *covarianceAlphamatrixSHORT = gsl_matrix_alloc(gsl_vector_get(fixedlengths,j),gsl_vector_get(fixedlengths,j));
                                        gsl_matrix *covarianceBetamatrixSHORT = gsl_matrix_alloc(gsl_vector_get(fixedlengths,j),gsl_vector_get(fixedlengths,j));
                                        gsl_matrix *weightAlphamatrixSHORT = gsl_matrix_alloc(gsl_vector_get(fixedlengths,j),gsl_vector_get(fixedlengths,j));
                                        gsl_matrix *weightBetamatrixSHORT = gsl_matrix_alloc(gsl_vector_get(fixedlengths,j),gsl_vector_get(fixedlengths,j));
                                        
                                        gsl_matrix_get_row(covarianceAlphavector,covarianceaux,indexa);
                                        vector2matrix(covarianceAlphavector,&covarianceAlphamatrix);
                                        gsl_vector_free(covarianceAlphavector); covarianceAlphavector = 0;
                                        tempm = gsl_matrix_submatrix(covarianceAlphamatrix,0,0,gsl_vector_get(fixedlengths,j),gsl_vector_get(fixedlengths,j));
                                        gsl_matrix_memcpy(covarianceAlphamatrixSHORT,&tempm.matrix);
                                        gsl_matrix_free(covarianceAlphamatrix); covarianceAlphamatrix = 0;
                                        perm = gsl_permutation_alloc(gsl_vector_get(fixedlengths,j));
                                        s = 0;
                                        gsl_linalg_LU_decomp(covarianceAlphamatrixSHORT,perm, &s);
                                        gsl_linalg_LU_invert(covarianceAlphamatrixSHORT,perm, weightAlphamatrixSHORT);
                                        gsl_matrix_free(covarianceAlphamatrixSHORT); covarianceAlphamatrixSHORT = 0;
                                        
                                        gsl_permutation_free(perm); perm = 0;
                                        
                                        gsl_matrix_get_row(covarianceBetavector,covarianceaux,indexb);
                                        vector2matrix(covarianceBetavector,&covarianceBetamatrix);
                                        gsl_vector_free(covarianceBetavector); covarianceBetavector = 0;
                                        tempm = gsl_matrix_submatrix(covarianceBetamatrix,0,0,gsl_vector_get(fixedlengths,j),gsl_vector_get(fixedlengths,j));
                                        gsl_matrix_memcpy(covarianceBetamatrixSHORT,&tempm.matrix);
                                        gsl_matrix_free(covarianceBetamatrix); covarianceBetamatrix = 0;
                                        perm = gsl_permutation_alloc(gsl_vector_get(fixedlengths,j));
                                        s = 0;
                                        gsl_linalg_LU_decomp(covarianceBetamatrixSHORT,perm, &s);
                                        gsl_linalg_LU_invert(covarianceBetamatrixSHORT,perm, weightBetamatrixSHORT);
                                        gsl_matrix_free(covarianceBetamatrixSHORT); covarianceBetamatrixSHORT = 0;
                                        
                                        matrixaux = gsl_matrix_alloc(gsl_vector_get(fixedlengths,j),gsl_vector_get(fixedlengths,j));
                                        gsl_matrix_memcpy(matrixaux,weightAlphamatrixSHORT);
                                        gsl_matrix_add(matrixaux,weightBetamatrixSHORT);
                                        gsl_matrix_scale(matrixaux,1.0/2.0);
                                        gsl_matrix_free(weightAlphamatrixSHORT); weightAlphamatrixSHORT = 0;
                                        gsl_matrix_free(weightBetamatrixSHORT); weightBetamatrixSHORT = 0;
                                        gsl_permutation_free(perm); perm = 0;
                                }
                                gsl_matrix_set_col(R,0,DabsSHORT);
                                
                                R_transW = gsl_matrix_alloc(2,gsl_vector_get(fixedlengths,j));      // R_transW = R'W               R_transW=(2xPL)
                                gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,R,matrixaux,0.0,R_transW);
                                R_transWR = gsl_matrix_alloc(2,2);                                  // R_transWR = R'WR             R_transWR=(2x2)
                                gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,R_transW,R,0.0,R_transWR);
                                matrixaux1 = gsl_matrix_alloc(2,gsl_vector_get(fixedlengths,j));
                                gsl_permutation *perm = gsl_permutation_alloc(2);
                                s=0;                                        
                                gsl_matrix_memcpy(aux,R_transWR);
                                gsl_linalg_LU_decomp(aux, perm, &s);
                                gsl_linalg_LU_invert(aux, perm, inv);
                                gsl_permutation_free(perm); perm = 0;
                                gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,inv,R_transW,0.0,matrixaux1);      // matrixaux1 = [(R'WR)^(-1)]\B7R'W       matrixaux1=(2xPL)
                                vectoraux1_2 = gsl_vector_alloc(gsl_vector_get(fixedlengths,j)*2);
                                for (int i=0;i<2;i++)
                                {
                                        for (int k=0;k<gsl_vector_get(fixedlengths,j);k++)
                                        {
                                                gsl_vector_set(vectoraux1_2,k+i*gsl_vector_get(fixedlengths,j),gsl_matrix_get(matrixaux1,i,k));
                                        }
                                }
                                for (int i=0;i<vectoraux1_2->size;i++)		gsl_matrix_set(*PrecalWMaux,indexa,i+indexPRCLWN,gsl_vector_get(vectoraux1_2,i));
                                
                                indexPRCLWN = indexPRCLWN + gsl_vector_get(fixedlengths,j)*2;
                                
                                gsl_matrix_free(R); R = 0;
                                gsl_matrix_free(R_transW); R_transW = 0;
                                gsl_matrix_free(R_transWR); R_transWR = 0;
                                gsl_matrix_free(matrixaux1); matrixaux1 = 0;
                                gsl_vector_free(vectoraux1_2); vectoraux1_2 = 0;
                                gsl_vector_free(DabsSHORT); DabsSHORT = 0;
                        }
                }

                gsl_matrix_free(aux); aux = 0;
                gsl_matrix_free(inv); inv = 0;
                
                gsl_vector_free(fixedlengths); fixedlengths = 0;
        }
       
    	gsl_vector_free(Pab); Pab = 0;
	gsl_vector_free(Dab); Dab = 0;
	gsl_vector_free(DabMaxLengthFixedFilter);DabMaxLengthFixedFilter = 0;
        
	gsl_vector_free(vectoraux); vectoraux = 0;
	gsl_vector_free(vectoraux1); vectoraux1 = 0;
	gsl_vector_free(vectoraux2); vectoraux2 = 0;
	gsl_matrix_free(matrixaux); matrixaux = 0;
        
	gsl_vector_free(PabMaxLengthFixedFilter); PabMaxLengthFixedFilter = 0;
        
	if (optimalfilter != NULL) {
          gsl_vector_free(optimalfilter); optimalfilter = 0;
        }
		
	return (EPOK);
}
/*xxxx end of SECTION A17 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A18 ************************************************************
* matrix2vector function: This function converts an input [nxn] square matrix into an output n^2 vector.
*                         It puts the first row of the matrix (n elements) in the first n elements of the vector (from 0 to n-1),
*                         the second row of the matrix in the elements from n to 2n-1 of the vector and so on.
*
* Parameters:
* - matrixin: GSL input square matrix whose dimensions are [nxn]
* - vectorout: GSL output vector whose length is n^2
******************************************************************************/
int matrix2vector (gsl_matrix *matrixin, gsl_vector **vectorout)
{
	string message = "";
	char valERROR[256];
  
	// matrixin is a square matrix
	int dim = matrixin->size1;
	// It is not necessary to check the allocation beacuse 'matrixin' size1 must be > 0
	gsl_vector *vectoraux = gsl_vector_alloc(dim);

	for (int i=0;i<dim;i++)
	{
		gsl_matrix_get_row(vectoraux,matrixin,i);
		if (i*dim+dim-1 > (*vectorout)->size-1)
		{
			sprintf(valERROR,"%d",__LINE__+7);
			string str(valERROR);
			message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		for (int j=0;j<dim;j++)
		{
			gsl_vector_set(*vectorout,i*dim+j,gsl_vector_get(vectoraux,j));
		}
	}

	gsl_vector_free(vectoraux); vectoraux = 0;

	return (EPOK);
}
/*xxxx end of SECTION A18 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A19 ************************************************************
* vector2matrix function: This function converts an input n^2 vector into an output [nxn] square matrix.
*                         It puts the first n elements of the vector in the first row of the matrix,
*                         the second group of n elements (from n to 2n-1) of the vector in the second row and so on.
*
* Parameters:
* - vectorin: GSL input vector whose length is n^2
* - matrixout: GSL output matrix whose dimensions are [nxn]
******************************************************************************/
int vector2matrix (gsl_vector *vectorin, gsl_matrix **matrixout)
{
	// matrixin is a square matrix
	double dim = sqrt(vectorin->size);

	for (int i=0;i<dim;i++)
	{
		for (int j=0;j<dim;j++)
		{
			gsl_matrix_set(*matrixout,i,j,gsl_vector_get(vectorin,i*dim+j));
		}
	}

	return (EPOK);
}
/*xxxx end of SECTION A19 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A20 ************************************************************
* convertI2R: This funcion converts the current space into a quasi-resistance space. 
*             The 'invector' filled with currents is filled in with resistances at the output.
* 
* - Read R0, I0_START=Ibias, ADUCNV and IMIN from the 'RecordFile' by using the pointer '(*reconstruct_init)->record_file_fptr'
* - Conversion according to 'EnergyMethod'=I2R
* - Conversion according to 'EnergyMethod'=I2RALL
* 	- Read also RL, LFILTER and TTR from the 'RecordFile'
* - Conversion according to 'EnergyMethod'=I2RNOL
*       - I2RALL neglecting the circuit inductance
* - Conversion according to 'EnergyMethod'=I2RFITTED
* - Because in 'runEnergy' the record (TesRecord) is used => The I2R, I2RALL, I2RNOL or I2RFITTED transformed record has to be used
* 
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
* - record: Structure containing the input record
* - invector: Input current (ADC) vector & output resistance (I2R, I2RALL or I2RNOL) vector
******************************************************************************/
//int convertI2R (ReconstructInitSIRENA** reconstruct_init, TesRecord **record, gsl_vector **invector)
int convertI2R (char* EnergyMethod, double R0, double Ibias, double Imin, double Imax, double TTR, double LFILTER, double RPARA, double samprate, gsl_vector **invector)
{
	int status = EPOK;
	string message="";
  
	//double R0;		// Operating point resistance [Ohm]
	//double Ibias;		// I0_START = Initial bias current [A]
	double aducnv;		// ADU conversion factor [A/ADU]
	//double baselineNew;
	double baseline;
	//double Imin;		// Current corresponding to 0 ADU [A]
	//double Imax;
        //double TTR;                                             // Transformer Turns Ratio
        //double LFILTER;                                         // Circuit filter inductance [H]
        //double RPARA;
        
        /*R0 = (*reconstruct_init)->i2rdata->R0;
        Ibias = (*reconstruct_init)->i2rdata->I0_START;
        Imin = (*reconstruct_init)->i2rdata->IMIN;
        Imax = (*reconstruct_init)->i2rdata->IMAX;
        RPARA = (*reconstruct_init)->i2rdata->RPARA;
        TTR = (*reconstruct_init)->i2rdata->TTR;
        LFILTER = (*reconstruct_init)->i2rdata->LFILTER;*/
        
	aducnv = (Imax-Imin)/65534;    // Quantification levels = 65534

	//if (strcmp((*reconstruct_init)->EnergyMethod,"I2R") == 0)
	if (strcmp(EnergyMethod,"I2R") == 0)
	{	
		// DeltaI <- I-Ibias
		// R <- R0 - R0*(abs(DeltaI)/Ibias)/(1+abs(DeltaI)/Ibias)

		// It is not necessary to check the allocation beacuse 'invector' size must be > 0
		gsl_vector *invector_modified = gsl_vector_alloc((*invector)->size);
		gsl_vector_scale(*invector,aducnv);             	// invector = I(ADC)*ADUCNV
		gsl_vector_add_constant(*invector,Imin);		// invector = I(Amp) = I(ADC)*ADUCNV+Imin
		gsl_vector_add_constant(*invector,-1*Ibias); 		// invector = DeltaI = abs(I(Amp)-Ibias)
		for (int i=0;i<(*invector)->size;i++)
		{
			if (gsl_vector_get(*invector,i)<0) 	gsl_vector_set(*invector,i,(-1.*gsl_vector_get(*invector,i)));
		}
		gsl_vector_scale(*invector,1./Ibias); 			// invector = DeltaI/Ibias = (I(Amp)-Ibias)/Ibias
		gsl_vector_memcpy(invector_modified,*invector);  	// invector_modified = invector = DeltaI/Ibias
		gsl_vector_add_constant(invector_modified,+1.0);	// invector_modified = 1 + DeltaI/Ibias
		gsl_vector_div(*invector,invector_modified);     	// invector = invector/invector_modified = (DeltaI/Ibias)/(1+DeltaI/Ibias)
		gsl_vector_scale(*invector,-1.*R0);			// invector = -R0*(DeltaI/Ibias)/(1+DeltaI/Ibias)
		gsl_vector_add_constant(*invector,R0); 			// invector = R0 - R0*(DeltaI/Ibias)/(1+DeltaI/Ibias)
		gsl_vector_free(invector_modified); invector_modified = 0;
	}
	//else if (strcmp((*reconstruct_init)->EnergyMethod,"I2RALL") == 0)
	else if (strcmp(EnergyMethod,"I2RALL") == 0)
        {
                double RL;                                              // Shunt/load resistor value [oHM]
                double V0;
                double L;
                //double RPARA;
                // It is not necessary to check the allocation beacuse 'invector' size must be > 0
                gsl_vector *I = gsl_vector_alloc((*invector)->size);
                gsl_vector *dI = gsl_vector_alloc((*invector)->size);
                gsl_vector *invector_modified = gsl_vector_alloc((*invector)->size);
               
                RL = RPARA/(pow(TTR,2));                                // RL = RPARA/(TTR)^2
                V0 = Ibias*(R0+RL);                                     // V0 = I0(R0+RL)
                L = LFILTER/(pow(TTR,2));                               // L = LFILTER/(TTR)^2

                // I
                gsl_vector_memcpy(invector_modified,*invector);
                gsl_vector_scale(*invector,aducnv);                     // invector = I(ADC)*ADUCNV
                gsl_vector_add_constant(*invector,Imin);                // invector = I(ADC)*ADUCNV+Imin
                gsl_vector_scale(*invector,-1.0);                       // invector = Ibias-(I(ADC)*ADUCNV+Imin)
                gsl_vector_add_constant(*invector,Ibias);
                gsl_vector_memcpy(I,*invector);

                // dI
                gsl_vector_memcpy(dI,I);
                if (differentiate (&dI, dI->size))
                {
                        message = "Cannot run routine differentiate for convertI2R and I2RALL";
                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                }

                // R = (V0-IRL-LdI/dt)/I
                gsl_vector_memcpy(*invector,dI);
                //gsl_vector_scale(*invector,1/(*record)->delta_t);       // dI/dt
                gsl_vector_scale(*invector,samprate);               // dI/dt
                gsl_vector_scale(*invector,-L);
                gsl_vector_add_constant(*invector,V0);                  // V0-dI/dt
                gsl_vector_memcpy(invector_modified,I);
                gsl_vector_scale(invector_modified,RL);
                gsl_vector_sub(*invector,invector_modified);            // V0-dI/dt-IRL
                gsl_vector_div(*invector,I);                            // (V0-IRL-LdI/dt)/I
               
                gsl_vector_free(I); I = 0;
                gsl_vector_free(dI); dI = 0;
                gsl_vector_free(invector_modified); invector_modified = 0;
        }
        //else if (strcmp((*reconstruct_init)->EnergyMethod,"I2RNOL") == 0)
        else if (strcmp(EnergyMethod,"I2RNOL") == 0)
        {
                double RL;                                      // Shunt/load resistor value [oHM]
                double V0;
                // It is not necessary to check the allocation beacuse 'invector' size must be > 0
                gsl_vector *I = gsl_vector_alloc((*invector)->size);
                gsl_vector *invector_modified = gsl_vector_alloc((*invector)->size);
               
                RL = RPARA/(pow(TTR,2));                        // RL = RPARA/(TTR)^2
                V0 = Ibias*(R0+RL);                             // V0 = I0(R0+RL)

                // I
                gsl_vector_memcpy(invector_modified,*invector);
                gsl_vector_scale(*invector,aducnv);             // invector = I(ADC)*ADUCNV
                gsl_vector_add_constant(*invector,Imin);        // invector = I(ADC)*ADUCNV+Imin
                gsl_vector_scale(*invector,-1.0);               // invector = Ibias-(I(ADC)*ADUCNV+Imin)
                gsl_vector_add_constant(*invector,Ibias);
                gsl_vector_memcpy(I,*invector);

                // R = (V0-IRL)/I
                gsl_vector_memcpy(*invector,I);
                gsl_vector_scale(*invector,-RL);                // IRL
                gsl_vector_add_constant(*invector,V0);          // V0-IRL
                gsl_vector_div(*invector,I);                    // (V0-IRL)/I
               
                gsl_vector_free(I); I = 0;
                gsl_vector_free(invector_modified); invector_modified = 0;
        }
        //else if (strcmp((*reconstruct_init)->EnergyMethod,"I2RFITTED") == 0)
        else if (strcmp(EnergyMethod,"I2RFITTED") == 0)
        {
                double RL;                                      // Shunt/load resistor value [oHM]
                double V0;
                double Ifit = 45.3e-6;                          // Ifit = 45.3 uA
                // It is not necessary to check the allocation beacuse 'invector' size must be > 0
                gsl_vector *I = gsl_vector_alloc((*invector)->size);
                gsl_vector *invector_modified = gsl_vector_alloc((*invector)->size);
                gsl_vector *IfitIgsl = gsl_vector_alloc((*invector)->size);
               
                RL = RPARA/(pow(TTR,2));                        // RL = RPARA/(TTR)^2
                V0 = Ibias*(R0+RL);                             // V0 = I0(R0+RL)

                // I
                gsl_vector_memcpy(invector_modified,*invector);
                gsl_vector_scale(*invector,aducnv);             // invector = I(ADC)*ADUCNV
                gsl_vector_add_constant(*invector,Imin);        // invector = I(ADC)*ADUCNV+Imin
                gsl_vector_scale(*invector,-1.0);               // invector = Ibias-(I(ADC)*ADUCNV+Imin)
                gsl_vector_add_constant(*invector,Ibias);
                gsl_vector_memcpy(I,*invector);

                // R = V0/(Ifit+I)
                gsl_vector_memcpy(IfitIgsl,I);
                gsl_vector_add_constant(IfitIgsl,Ifit);         // Ifit+I
                gsl_vector_set_all(*invector,V0);
                gsl_vector_div(*invector,IfitIgsl);             // V0/(Ifit+I)
               
                gsl_vector_free(I); I = 0;
                gsl_vector_free(invector_modified); invector_modified = 0;
                gsl_vector_free(IfitIgsl); IfitIgsl = 0;
        }

	/*for (int i=0;i<(*invector)->size;i++)		     // Because in 'runEnergy' the record (TesRecord) is used => The I2R, I2RALL, I2RNOL or I2RFITTED transformed record has to be used
	{
		(*record)->adc_double[i] = gsl_vector_get(*invector,i);
	}*/
	  
	return(EPOK);
}
/*xxxx end of SECTION A20 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A21 ************************************************************
* filterByWavelets: This funcion filters the input/output signal 'invector', reducing the noise level.
* 
* - It is only going to work with 'n' elements of 'invector'
* - Discrete Wavelet Treansform 
* - Sorting coefficients
* - Hard thresholding: 'n-nc' coefficeints are deleted (those with low energy)
* - Inverse DWT
* 
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values) 
* - invector: Input/output signal 
* - length: Length of the wavelet transform 
* - onlyOnce: In order to control the times to be executed
******************************************************************************/
int filterByWavelets (ReconstructInitSIRENA* reconstruct_init, gsl_vector **invector, int length, int *onlyOnce)
{
	string message = "";
	char valERROR[256];
	
	*onlyOnce = 1;
	
	// n must be a power of 2
	//int i, n = 256, nc = 20;
	//int i, n = pow(2,floor(log2((*invector)->size))), nc = 2048;
	//int i, n = pow(2,floor(log2((*invector)->size))), nc = n/8;
	//int i, n = length, nc = n/8;
	int i, n = length, nc = n/16;
	//int i, n = pow(2,floor(log2((*invector)->size))), nc = n;
	double data[n];
	double abscoeff[n];
	size_t p[n];
	
	char temporalFileRecordName[255];
	char temporalFiled13Name[255];
	char aux[255];
	char val1[256],val2[256],val3[256],val4[256],val5[256],val6[256],val7[256],val8[256];
	
	sprintf(aux,"_");
	sprintf(val1,"%d",n);
	strcat(aux,val1);
	sprintf(val2,"_");
	strcat(aux,val2);
	sprintf(val3,"%d",nc);
	strcat(aux,val3);
	strcat(aux,val2);
	sprintf(val5,"%d",reconstruct_init->tstartPulse2-reconstruct_init->tstartPulse1);
	strcat(aux,val5);
	strcat(aux,val2);
	sprintf(val6,"%d",(int)(reconstruct_init->monoenergy));
	strcat(aux,val6);
	strcat(aux,".txt");
	sprintf(temporalFileRecordName,"record");
	strcat(temporalFileRecordName,aux);
	sprintf(temporalFiled13Name,"d13");
	strcat(temporalFiled13Name,aux);
	
	FILE * temporalFileRecord;
	temporalFileRecord = fopen (temporalFileRecordName,"w");
	FILE * temporalFiled13;
	temporalFiled13 = fopen (temporalFiled13Name,"w");
		
	gsl_wavelet *w;
	gsl_wavelet_workspace *work;
	
	w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4); 	// Wavelets from the Daubechies family are used as bais fucntions in representing other functions
	// It is not necessary to check the allocation beacuse 'length' must be > 0
	work = gsl_wavelet_workspace_alloc (n);
	
	if (n > (*invector)->size)
	{
		sprintf(valERROR,"%d",__LINE__+7);
		string str(valERROR);
		message = "Getting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	for (i = 0; i < n; i++)
	{
		data[i] = gsl_vector_get(*invector,i); 
	}
	
	// DWT (Discrete Wavelet Transform)
	// An\E1lisis: Genera diferentes subbandas => Distintos niveles de descomposici\F3n se pueden generar de acuerdo a las necesidades de la aplicaci\F3n
	// El an\E1lisis wavelet entrega una serie de coeficientes 'data'. Dado que es un proceso reversible, la se\F1al original puede obtenerse a partir 
	// de los coeficientes obtenidos en el an\E1lisis.
	// Los coeficientes se agrupan en los de baja frecuencia (aproximaciones) y los de alta frecuencia (detalles).
	gsl_wavelet_transform_forward (w, data, 1, n, work);
	
	for (i = 0; i < n; i++)
	{
		abscoeff[i] = fabs (data[i]);
	}
	for (i = 0; i < (*invector)->size/2; i++)
	{
		sprintf(val7,"%d %e",i*2,data[(*invector)->size-(*invector)->size/2+i]);
		strcat(val7,"\n");
		fputs(val7,temporalFiled13);
	}
	
	// Sorting the 'n' elements of 'abscoeff' into ascending order, storing the resulting permutation in 'p'
	gsl_sort_index (p, abscoeff, 1, n);
		
	// 'Wavelet shrinkage' reduce la magnitud de cada coeficiente dependiendo del nivel de ruido que se estima de la se\F1al. Un m\E9todo com\FAn es la designaci\F3n de un umbral. 
	// Al definir un umbral los valores de los coeficientes wavelet que se encuentren por debajo de \E9l son eliminados (hard-threshold) o bien reducidos en magnitud (soft-threshold).
	// Hard thresholding: 'n-nc' coefficeints are deleted (those with low energy)
	for (i = 0; (i + nc) < n; i++)		// If 'n = nc' => No coefficeints are deleted => In = Out
	{
		data[p[i]] = 0;
	}
	
	// IDWT (Inverse Discrete Wavelet Transform) by using the new coefficients
	gsl_wavelet_transform_inverse (w, data, 1, n, work);
	
	for (i = 0; i < n; i++)
	{
		sprintf(val8,"%e %e",gsl_vector_get(*invector,i),data[i]);
		strcat(val8,"\n");
		fputs(val8,temporalFileRecord);
		
		gsl_vector_set(*invector,i,data[i]); 
	}
	
	gsl_wavelet_free (w); w = 0;
	gsl_wavelet_workspace_free (work); work = 0;
	
	fclose(temporalFileRecord);
	fclose(temporalFiled13);
	
	return(EPOK);
}
/*xxxx end of SECTION A21 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B ************************************************************
* runEnergy: This function calculates the pulse energy applying different methods (from 'EnergyMethod' and 'OFNoise').
* 	     It only runs in RECONSTRUCTION mode (not in PCA).
*
* - Declare variables
* - Store the record in 'invector' ('loadRecord')
* - Subtract the baseline if OPTFILT and 'runF0orB0val'= 1 ('FilterMethod'=B0)
* - Subtract the baseline if WEIGHT
* - Check Quality
* - For each pulse:
* 	- Establish the pulse grade (HighRes=1, MidRes=2, LimRes=3, LowRes=4, Rejected=-1, Pileup=-2) and the optimal filter length
* 	- Pulse: Load the proper piece of the record in 'pulse'
* 	- If 'OFIter'=1, in the first iteration ('numiteration'=0) the values of 'maxDER' and 'maxDERs' are used in 'find_matchedfilter', 'find_matchedfilterDAB', 
*         'find_optimalfilter', 'find_optimalfilterDAB' or 'find_Esboundary' getting the values of the 'energies' which straddle the 'maxDER' ('Ealpha' and 'Ebeta'). There will be more
*         iterations if the calculated 'energy' is out of ['Ealpha','Ebeta']. If 'energy' is in ['Ealpha','Ebeta'] the iterative process stops.
* 	  	- If OPTFILT (or I2R, I2RALL, I2RNOL, I2RFITTED) and 'OFLib'=0 and 'OFNOise=NSD':
*	 	    - Filter (find the matched filter and load it in 'filter'):
*		    	- If 'OFInterp'=MF => 'find_matchedfilter'
* 			- If 'OFInterp'=DAB => 'find_matchedfilterDAB'
* 		    - Calculate the optimal filter
* 		- If OPTFILT (or I2R, I2RALL, I2RNOL, I2RFITTED) and 'OFLib'=1 and 'OFNOise=NSD':
*                   - Choose the base-2 system value closest (lower than or equal) to the pulse length
*	 	    - Filter (find the optimal filter and load it in 'optimalfilter'):
*                       - If 'OFInterp'=MF => 'find_optimalfilter'
* 			- If 'OFInterp'=DAB => 'find_optimalfilterDAB'
*	 	    - If 'FilterDomain'=F => Bin0=0 and conjugate
*		- If 'WEIGHT' or 'WEIGHTN':
*		    - Get the indexes of the two energies which straddle the pulse
* 		    - If 'WEIGHTN' and 'OFLib'=1:
*                       - Choose the base-2 system value closest (lower than or equal) to the pulse length
* 		        - 'find_prclwn' to find the appropriate values of the PRECALWN HDU ('PRCLx' columns)
*               - If OPTFILT (or I2R, I2RALL, I2RNOL, I2RFOTTED) and 'OFLib'=1 and 'OFNOise=WEIGHTM':
*                   - Choose the base-2 system value closest (lower than or equal) to the pulse length
* 		    - 'find_prclofwm' to find the appropriate values of the PRCLOFWM HDU ('OFWx' columns)
*		- Calculate the energy of each pulse
*               - If using lags, it is necessary to modify the tstart of the pulse and the length of the filter used
*       - In order to subtract the pulse model, it has to be located in the tstart with jitter and know its values in the digitized samples
*       - Subtract the pulse model from the record
*	- Write info of the pulse in the output intemediate file if 'intermediate'=1
* - Not valid pulse => Its info has to be also stored in the intermediate file (if 'intermediate'=1) and in the structure 'pulsesInRecord'
* - Free allocated GSL vectors
*
* Parameters:
* - record: Structure that contains the input ADC record
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
* - pulsesInRecord: Collection of pulses found in the current record
* - optimalFilter: Optimal filters used in reconstruction
******************************************************************************/
void runEnergy(TesRecord* record,ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord, OptimalFilterSIRENA **optimalFilter)
{
	// Declare variables
	string message="";
	int status = EPOK;
	char valERROR[256];

	fitsfile *dtcObject = NULL;	    // FITS object containing information of the intermediate output FITS file
	if ((*reconstruct_init)->intermediate == 1)
	{
		char dtcName[256];
		strncpy(dtcName,(*reconstruct_init)->detectFile,255);
		dtcName[255]='\0';
	}

	int TorF;
	if (strcmp((*reconstruct_init)->FilterDomain,"T") == 0)		// Time
	{
		TorF=0;
	}
	else if (strcmp((*reconstruct_init)->FilterDomain,"F") == 0)	// Frequency
	{
		TorF=1;
	}

	int runF0orB0val;
	if (strcmp((*reconstruct_init)->FilterMethod,"F0") == 0)	// Deleting the frequency-zero bin
	{
		runF0orB0val = 0;
	}
	else if (strcmp((*reconstruct_init)->FilterMethod,"B0") == 0)	// Working without baseline
	{
		runF0orB0val = 1;
	}

	// I2R, I2RALL, I2RNOL or I2RFITTED methods convert I into R at the beginnig and after that 'I2R', 'I2RALL', 'I2RNOL' or 'I2RFITTED' are equivalent to 'OPTFILT'
	int runEMethod;
	if (strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0)
	{
		runEMethod = 0;
	}
	else if (strcmp((*reconstruct_init)->EnergyMethod,"WEIGHT") == 0)
	{
		runEMethod = 1;
	}
	else if (strcmp((*reconstruct_init)->EnergyMethod,"WEIGHTN") == 0)
	{
		runEMethod = 2;
	}

	int OFlength_strategy;
	if (strcmp((*reconstruct_init)->OFStrategy,"FREE") == 0)
	{
		OFlength_strategy = 0;
	}
	else if (strcmp((*reconstruct_init)->OFStrategy,"BASE2") == 0)
	{
		OFlength_strategy = 1;
	}
	else if (strcmp((*reconstruct_init)->OFStrategy,"BYGRADE") == 0)
	{
		OFlength_strategy = 2;
	}
	else if (strcmp((*reconstruct_init)->OFStrategy,"FIXED") == 0)
	{
		OFlength_strategy = 3;
	}

	long energyInLibrary_row;

	double energy;
	double tstartNewDev = -999.0;    	// Deviation of the starting of the pulses (in samples) respect to the tstart calculated
	//int numlags = 3; 			// Must be odd
	int numlags = (*reconstruct_init)->nLags; 			// Must be odd
	int numlags2 = floor(numlags/2);
        //cout<<"nLags: "<<numlags<<" numlags2: "<<numlags2<<endl;
	int lagsShift = -999;                   // Number of samples shifted to find the maximum of the parabola

	double Ealpha, Ebeta;
	gsl_vector *optimalfilter = NULL;	// Resized optimal filter expressed in the time domain (optimalfilter(t))
	gsl_vector *optimalfilter_f = NULL;	// Resized optimal filter f's when f's are according to [0,...fmax,-fmax,...] (frequency domain)
	gsl_vector *optimalfilter_FFT = NULL;	// Resized optimal filter magnitudes when f's are according to [0,...fmax,-fmax,...] (frequency domain)
	gsl_vector_complex *optimalfilter_FFT_complex = NULL;
	
	gsl_vector *Pab = NULL;
	gsl_matrix *PRCLWN = NULL;
	gsl_matrix *PRCLOFWM = NULL;

	gsl_vector *model;

	gsl_vector *recordAux;

	gsl_vector *pulse = NULL;
	gsl_vector *pulseToCalculateEnergy = NULL;	// Just in case LagsOrNot=1

	gsl_vector *filtergsl = NULL;		// Matched filter values (time domain)

	int iter;
	gsl_vector_view(temp);
        
        double tstartSamplesRecord;		// Tstart of the pulse in samples from the beginning of the record
	double tstartRecord;			// Tstart of the record in seconds
	int tstartRecordSamples = record->time/record->delta_t;	// Tstart of the record in samples (nearest sample)
        //cout<<"tstartRecordSamples: "<<tstartRecordSamples<<endl;
        double tstartJITTER;
        gsl_matrix *optimalfilterMAT;
        gsl_matrix_complex *optimalfilter_FFT_complexMAT;
        int tstartSamplesRecordStartINT;
        double tstartSamplesRecordStartDOUBLE;
        
        double shift;

	int indexEalpha = 0;
	int indexEbeta = 0;

	long resize_mf;
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double energy_lowres;
        long resize_mf_lowres = 4; // In order to get the low resolution energy estimator by filtering with a 4-samples-long filter
        gsl_vector *pulse_lowres = gsl_vector_alloc(4);
        gsl_vector *filtergsl_lowres = NULL;
        if (strcmp((*reconstruct_init)->FilterDomain,"T") == 0)		filtergsl_lowres= gsl_vector_alloc(resize_mf_lowres);
        else if (strcmp((*reconstruct_init)->FilterDomain,"F") == 0)	filtergsl_lowres= gsl_vector_alloc(resize_mf_lowres*2);
        gsl_vector *Pab_lowres = gsl_vector_alloc(resize_mf_lowres);
        gsl_matrix *PRCLWN_lowres = NULL;
	gsl_matrix *PRCLOFWM_lowres = NULL;
        double Ealpha_lowres, Ebeta_lowres;
        gsl_vector *optimalfilter_lowres = gsl_vector_alloc(filtergsl_lowres->size);	// Resized optimal filter expressed in the time domain (optimalfilter(t))
	gsl_vector_complex *optimalfilter_FFT_complex_lowres = NULL;
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
        int pulseGrade; 			// Pileup=-2, Rejected=-1, HighRes=1, MidRes=2, LimRes=3, LowRes=4
	
	// Store the record in 'invector'
	// It is not necessary to check the allocation because 'record->trigger_size' has been checked previously
	recordAux = gsl_vector_alloc(record->trigger_size);
	if (loadRecord(record, &tstartRecord, &recordAux))
	{
		message = "Cannot run routine loadRecord";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	// Subtract the baseline if OPTFILT and 'runF0orB0val'= 1 ('FilterMethod'=B0)
	// Subtract the baseline if WEIGHT
	if (((runF0orB0val == 1) && (runEMethod == 0)) || (runEMethod == 1))
	{
		// It is not necessary to check the allocation because the allocation of 'recordAux' has been checked previously
		gsl_vector *baselinegsl = gsl_vector_alloc(recordAux->size);
		gsl_vector_set_all(baselinegsl,-1.0*(*reconstruct_init)->noise_spectrum->baseline);
		gsl_vector_add(recordAux,baselinegsl);
		gsl_vector_free(baselinegsl); baselinegsl = 0;
	}

	// Check Quality
	iter = 0;
	for (int i=0; i<(*pulsesInRecord)->ndetpulses;i++)
	{
		if ((*pulsesInRecord)->pulses_detected[i].quality >= 10)  // Saturated pulse
		{
			iter++;
		}
	}
	if (iter == (*pulsesInRecord)->ndetpulses)	
	{
		message = "There are no unsaturated pulses in one record";
		EP_PRINT_ERROR(message,-999); // Only a warning
	}

	model =gsl_vector_alloc((*reconstruct_init)->pulse_length);
        
        double valaux;
        int newidx;
        	  
        int extraSizeDueToLags = 0;
	for (int i=0; i<(*pulsesInRecord)->ndetpulses ;i++)
	{
		// Establish the pulse grade (HighRes=1, MidRes=2, LimRes=3, LowRes=4, Rejected=-1, Pileup=-2) and the optimal filter length
		if ((*pulsesInRecord)->pulses_detected[i].quality == 1)		(*pulsesInRecord)->pulses_detected[i].grade1 = -1;
		else								(*pulsesInRecord)->pulses_detected[i].grade1 = (*pulsesInRecord)->pulses_detected[i].pulse_duration;
                
                pulseGrade = 0;
		if (pulseGrading(*reconstruct_init,(*pulsesInRecord)->pulses_detected[i].grade1,(*pulsesInRecord)->pulses_detected[i].grade2,OFlength_strategy,&pulseGrade,&resize_mf))
		{
			message = "Cannot run routine pulseGrading";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		(*pulsesInRecord)->pulses_detected[i].grade1 = resize_mf;
                //cout<<"resize_mf: "<<resize_mf<<endl;

                //cout<<"pulses_detected[i].Tstart0: "<<(*pulsesInRecord)->pulses_detected[i].Tstart<<endl;
		// Pulse: Load the proper piece of the record in 'pulse'
                //cout<<"((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t: "<<((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t<<endl;
                //tstartSamplesRecord = ((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t;
                tstartSamplesRecord = (*pulsesInRecord)->pulses_detected[i].TstartSamples;
                //cout<<"tstartSamplesRecord: "<<tstartSamplesRecord<<endl;
                //tstartSamplesRecordStartINT = ((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t-numlags2;
                //tstartSamplesRecordStartDOUBLE = ((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t-numlags2;
                tstartSamplesRecordStartINT = tstartSamplesRecord-numlags2;
                tstartSamplesRecordStartDOUBLE = tstartSamplesRecord-numlags2;              //He cambiado esto
                //tstartSamplesRecordStartDOUBLE = tstartSamplesRecord-(numlags-1);
                //cout<<"tstartSamplesRecord: "<<tstartSamplesRecord<<endl;
                //cout<<"tstartSamplesRecordStartDOUBLE0: "<<tstartSamplesRecordStartDOUBLE<<endl;
                
                //cout<<"tstart: "<<tstartSamplesRecord<<endl;
                //cout<<"tstartRecord: "<<tstartRecord<<endl;
                //cout<<"tstartRecordSamples: "<<tstartRecordSamples<<endl;
		if ((pulse = gsl_vector_alloc(resize_mf)) == 0)
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);
			message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		if ((tstartSamplesRecord < 0) || (tstartSamplesRecord > recordAux->size-2)
			|| (resize_mf < 1) || (resize_mf > recordAux->size-tstartSamplesRecord))
		{
			sprintf(valERROR,"%d",__LINE__+5);
			string str(valERROR);
			message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL); 
		}
		temp = gsl_vector_subvector(recordAux,tstartSamplesRecord,resize_mf);
		if (gsl_vector_memcpy(pulse, &temp.vector) != 0)
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);	
			message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		//cout<<"tstartRecordSamples: "<<tstartRecordSamples<<endl;
		tstartJITTER = ((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t;
		shift = tstartJITTER - tstartSamplesRecord;
                //cout<<"tstartJITTER: "<<tstartJITTER<<endl;
                //cout<<"shift: "<<shift<<endl;
                
                //cout<<"pulse0: "<<endl;
                //for (int kk=0;kk<10;kk++) cout<<kk<<" "<<gsl_vector_get(pulse,kk)<<endl;
                
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //////////// In order to get the low resolution energy estimator by filtering with a 4-samples-long filter ///////////////////
                // Pulse 
                if (resize_mf_lowres > recordAux->size-tstartSamplesRecord)
		{
			sprintf(valERROR,"%d",__LINE__+5);
			string str(valERROR);
			message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL); 
		}
		temp = gsl_vector_subvector(recordAux,tstartSamplesRecord,resize_mf_lowres);
		if (gsl_vector_memcpy(pulse_lowres, &temp.vector) != 0)
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);	
			message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		// Filter
                if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
                {
                    if (find_optimalfilter((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl_lowres, &Ealpha_lowres, &Ebeta_lowres))
                    {
                        message = "Cannot run routine find_optimalfilter for filter interpolation";
                        EP_EXIT_ERROR(message,EPFAIL);
                    }
                }
                else
                {
                    if (find_optimalfilterDAB((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl_lowres, &Pab_lowres,&Ealpha_lowres, &Ebeta_lowres))
                    {
                        message = "Cannot run routine find_optimalfilterDAB for filter interpolation";
                        EP_EXIT_ERROR(message,EPFAIL);
                    }
                }
                gsl_vector_set_all(optimalfilter_lowres,0);
                if (TorF == 0)     gsl_vector_memcpy(optimalfilter_lowres,filtergsl_lowres);
                else if (TorF == 1)
                {
                    // It is not necessary to check the allocation because 'filtergsl' size has been checked previously
                    optimalfilter_FFT_complex_lowres = gsl_vector_complex_alloc(filtergsl_lowres->size/2);
                    for (int k=0;k<filtergsl_lowres->size/2;k++)
                    {
                        gsl_vector_complex_set(optimalfilter_FFT_complex_lowres,k,gsl_complex_rect(gsl_vector_get(filtergsl_lowres,k),gsl_vector_get(filtergsl_lowres,k+filtergsl_lowres->size/2)));
                    }
                    gsl_vector_complex_set(optimalfilter_FFT_complex_lowres,0,gsl_complex_rect(0.0,0.0));
                    for (int k=0;k<filtergsl_lowres->size/2;k++)
                    {
                        gsl_vector_complex_set(optimalfilter_FFT_complex_lowres,k,gsl_complex_conjugate(gsl_vector_complex_get(optimalfilter_FFT_complex_lowres,k)));
                    }
                }
                // Calculate the low resolution estimator
                if (calculateEnergy(pulse_lowres,1,optimalfilter_lowres,optimalfilter_FFT_complex_lowres,0,0,0,(*reconstruct_init),TorF,1/record->delta_t,Pab_lowres,PRCLWN_lowres,PRCLOFWM_lowres,&energy_lowres,&tstartNewDev,&lagsShift,1))
                {
                    message = "Cannot run calculateEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
                    EP_EXIT_ERROR(message,EPFAIL);
                }
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		if ((*reconstruct_init)->LagsOrNot == 0)	
		{
			pulseToCalculateEnergy = gsl_vector_alloc(pulse->size);
			gsl_vector_memcpy(pulseToCalculateEnergy,pulse);
		}
		else if (((*reconstruct_init)->LagsOrNot == 1) && (strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0))	
		{
                        //cout<<"numlags: "<<numlags<<endl;
			int resize_mfNEW = resize_mf + numlags -1;
                        //cout<<"resize_mfNEW: "<<resize_mfNEW<<endl;
			if ((pulseToCalculateEnergy = gsl_vector_alloc(resize_mfNEW)) == 0)
			{
				sprintf(valERROR,"%d",__LINE__-2);
				string str(valERROR);
				message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
				EP_EXIT_ERROR(message,EPFAIL);
			}
			gsl_vector_set_all(pulseToCalculateEnergy,-999);
			
                        /*cout<<"recordAux->size: "<<recordAux->size<<endl;
                        cout<<"tstartSamplesRecord: "<<tstartSamplesRecord<<endl;*/
                        //cout<<"tstartSamplesRecordStartDOUBLE: "<<tstartSamplesRecordStartDOUBLE<<endl;
			if ((tstartSamplesRecordStartDOUBLE-numlags/2 < 0) || (tstartSamplesRecordStartDOUBLE-numlags/2 > recordAux->size-2)
				|| (resize_mfNEW < 1) || (resize_mfNEW > recordAux->size-tstartSamplesRecordStartDOUBLE+numlags/2))
			{
				sprintf(valERROR,"%d",__LINE__+5);
				string str(valERROR);
				message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
				EP_EXIT_ERROR(message,EPFAIL); 
			}
                        temp = gsl_vector_subvector(recordAux,tstartSamplesRecordStartDOUBLE,resize_mfNEW);
                        //cout<<"recordAux: "<<endl;
                        //for (int kk=0;kk<recordAux->size;kk++) cout<<kk<<" "<<gsl_vector_get(recordAux,kk)<<endl;
                        /*cout<<"995 "<<gsl_vector_get(recordAux,995)<<endl;
                        cout<<"996 "<<gsl_vector_get(recordAux,996)<<endl;
                        cout<<"997 "<<gsl_vector_get(recordAux,997)<<endl;
                        cout<<"998 "<<gsl_vector_get(recordAux,998)<<endl;
                        cout<<"999 "<<gsl_vector_get(recordAux,999)<<endl;
                        cout<<"1000 "<<gsl_vector_get(recordAux,1000)<<endl;
                        cout<<"1001 "<<gsl_vector_get(recordAux,1001)<<endl;
                        cout<<"1002 "<<gsl_vector_get(recordAux,1002)<<endl;
                        cout<<"1003 "<<gsl_vector_get(recordAux,1003)<<endl;
                        cout<<"1004 "<<gsl_vector_get(recordAux,1004)<<endl;
                        cout<<"1005 "<<gsl_vector_get(recordAux,1005)<<endl;*/
			if (gsl_vector_memcpy(pulseToCalculateEnergy, &temp.vector) != 0)
			{
				sprintf(valERROR,"%d",__LINE__-2);
				string str(valERROR);	
				message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
				EP_EXIT_ERROR(message,EPFAIL);
			}
                        
                        //cout<<"pulseToCalculateEnergy->size: "<<pulseToCalculateEnergy->size<<endl;
                        //cout<<"pulseToCalculateEnergy0: "<<endl;
                        //for (int kk=0;kk<pulseToCalculateEnergy->size;kk++) cout<<kk<<" "<<gsl_vector_get(pulseToCalculateEnergy,kk)<<endl;
                        //for (int kk=0;kk<10;kk++) cout<<kk<<" "<<gsl_vector_get(pulseToCalculateEnergy,kk)<<endl;
                        
                        //cout<<"tstartSamplesRecordStartINT: "<<tstartSamplesRecordStartINT<<endl; 
                        //cout<<"tstartSamplesRecordStartDOUBLE: "<<tstartSamplesRecordStartDOUBLE<<endl; 
                        
                        //cout<<"cicloFOR"<<endl;
                        /*for (int kk=0;kk<resize_mfNEW;kk++)
                        {
                            //if ((j+tstartSamplesRecord-numlags/2>=995) && (j+tstartSamplesRecord-numlags/2<1002)) cout<<"recordAux: "<<j+tstartSamplesRecord-numlags/2<<" "<<gsl_vector_get(recordAux,j+tstartSamplesRecord-numlags/2)<<endl;
                            //gsl_vector_set(pulseToCalculateEnergy,j,gsl_vector_get(recordAux,j+tstartSamplesRecord-numlags/2));
                            //if ((j+tstartSamplesRecord-numlags/2>=995) && (j+tstartSamplesRecord-numlags/2<1002)) cout<<"pulseToCalculateEnergy: "<<j<<" "<<gsl_vector_get(pulseToCalculateEnergy,j)<<endl;
                            //if ((kk+tstartSamplesRecordStartDOUBLE>=995) && (kk+tstartSamplesRecordStartDOUBLE<1002)) cout<<kk<<" "<<tstartSamplesRecordStartDOUBLE<<" recordAux: "<<kk+tstartSamplesRecordStartDOUBLE<<" "<<gsl_vector_get(recordAux,kk+tstartSamplesRecordStartDOUBLE)<<endl;
                            //gsl_vector_set(pulseToCalculateEnergy,kk,gsl_vector_get(recordAux,kk+tstartSamplesRecordStartDOUBLE));
                            
                            newidx=kk+tstartSamplesRecordStartDOUBLE;
                            valaux = gsl_vector_get(recordAux,newidx);
                            //cout<<kk<<" "<<tstartSamplesRecordStartDOUBLE<<" "<<newidx<<" valaux: "<<valaux<<" "<<gsl_vector_get(recordAux,1000)<<endl;
                            gsl_vector_set(pulseToCalculateEnergy,kk,valaux);
                            //if (kk == 5)    {cout<<"Hola"<<endl; gsl_vector_set(pulseToCalculateEnergy,5,2305);}
                            //cout<<gsl_vector_get(recordAux,1000)<<" "<<kk<<" "<<tstartSamplesRecordStartDOUBLE<<" "<<" "<<kk+tstartSamplesRecordStartDOUBLE<<" "<<gsl_vector_get(recordAux,kk+tstartSamplesRecordStartDOUBLE)<<" "<<gsl_vector_get(pulseToCalculateEnergy,kk)<<endl;
                            //if ((kk+tstartSamplesRecordStartDOUBLE>=995) && (kk+tstartSamplesRecordStartDOUBLE<1002)) cout<<kk<<" "<<tstartSamplesRecordStartDOUBLE<<" pulseToCalculateEnergy: "<<kk<<" "<<gsl_vector_get(pulseToCalculateEnergy,kk)<<endl;
                        }*/
                        
			extraSizeDueToLags = numlags-1;
                        //cout<<"vector: "<<gsl_vector_get(pulseToCalculateEnergy,995)<<","<<gsl_vector_get(pulseToCalculateEnergy,996)<<","<<gsl_vector_get(pulseToCalculateEnergy,997)<<","<<gsl_vector_get(pulseToCalculateEnergy,998)<<","<<gsl_vector_get(pulseToCalculateEnergy,999)<<","<<gsl_vector_get(pulseToCalculateEnergy,1000)<<endl;
                        //cout<<"pulseToCalculateEnergy: "<<endl;
                        //for (int kk=0;kk<pulseToCalculateEnergy->size;kk++) cout<<kk<<" "<<gsl_vector_get(pulseToCalculateEnergy,kk)<<endl;
                        /*cout<<"0 "<<gsl_vector_get(pulseToCalculateEnergy,0)<<endl;
                        cout<<"1 "<<gsl_vector_get(pulseToCalculateEnergy,1)<<endl;
                        cout<<"2 "<<gsl_vector_get(pulseToCalculateEnergy,2)<<endl;
                        cout<<"3 "<<gsl_vector_get(pulseToCalculateEnergy,3)<<endl;
                        cout<<"4 "<<gsl_vector_get(pulseToCalculateEnergy,4)<<endl;
                        cout<<"5 "<<gsl_vector_get(pulseToCalculateEnergy,5)<<endl;
                        cout<<"6 "<<gsl_vector_get(pulseToCalculateEnergy,6)<<endl;
                        cout<<"7 "<<gsl_vector_get(pulseToCalculateEnergy,7)<<endl;
                        cout<<"8 "<<gsl_vector_get(pulseToCalculateEnergy,8)<<endl;
                        cout<<"9 "<<gsl_vector_get(pulseToCalculateEnergy,9)<<endl;
                        cout<<"10 "<<gsl_vector_get(pulseToCalculateEnergy,10)<<endl;*/
		}
		
		bool iterate = true;
		if ((*reconstruct_init)-> OFIter == 1)	iterate = true;
		else iterate = false;
		// It is not necessary to check the allocation because '(*reconstruct_init)->library_collection->ntemplates' has been check previously
		gsl_matrix *Estraddle = gsl_matrix_alloc(2,(*reconstruct_init)->library_collection->ntemplates);
		gsl_matrix *resultsE = gsl_matrix_alloc(2,(*reconstruct_init)->library_collection->ntemplates);		// Row0 -> calculatedEnergy
															// Row1 -> min[abs(calculatedEnergy-Ealpha),abs(calculatedEnergy-Ebeta)]
		int numiteration = -1;
                
		do
		{
			numiteration++;
			
			if ((strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0) && (strcmp((*reconstruct_init)->OFNoise,"NSD") == 0))
			{
				// Filter (find the matched filter and load it in 'filter')
				if ((*reconstruct_init)->OFLib == 0)
				{	
					// It is not necessary to check the allocation because '(*reconstruct_init)->pulse_length'='PulseLength'(input parameter) has been checked previously
					filtergsl= gsl_vector_alloc(resize_mf);
					Pab = gsl_vector_alloc(resize_mf);
					if (numiteration == 0)
					{
						if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
						{
							if (find_matchedfilter(runF0orB0val, (*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_matchedfilter for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
						else //(*reconstruct_init)->OFInterp = DAB
						{
							if (find_matchedfilterDAB((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl, &Pab, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_matchedfilterDAB for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
					}
					else
					{
						if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
						{
							if (find_matchedfilter(runF0orB0val, energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_matchedfilter for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
						else //(*reconstruct_init)->OFInterp = DAB
						{
							if (find_matchedfilterDAB(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl, &Pab, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_matchedfilterDAB for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
					}
				}
				else if ((*reconstruct_init)->OFLib == 1)
				{
                    
                                        // Choose the base-2 system value closest (lower than or equal) to the pulse length
					if ((strcmp((*reconstruct_init)->OFStrategy,"BASE2") != 0) 
					    && ((((*reconstruct_init)->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != -999) && (resize_mf != (*reconstruct_init)->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration))
					    && (resize_mf != (*reconstruct_init)->library_collection->pulse_templates[0].template_duration)))
					{                                             
                                                resize_mf = pow(2,floor(log2(resize_mf)));
                                                gsl_vector *pulse_aux = gsl_vector_alloc(resize_mf+extraSizeDueToLags);
                                                //gsl_vector_view(temp3);
                                                temp = gsl_vector_subvector(pulseToCalculateEnergy,0,resize_mf+extraSizeDueToLags);
                                                //temp3 = gsl_vector_subvector(pulseToCalculateEnergy,0,resize_mf+extraSizeDueToLags);
                                                gsl_vector_memcpy(pulse_aux,&temp.vector);
                                                //gsl_vector_memcpy(pulse_aux,&temp3.vector);
                                                gsl_vector_free(pulseToCalculateEnergy);
                                                pulseToCalculateEnergy = gsl_vector_alloc(resize_mf+extraSizeDueToLags);
                                                gsl_vector_memcpy(pulseToCalculateEnergy,pulse_aux);
                                                gsl_vector_free(pulse_aux);
                                                //>>>>>>> master
					}
				  
					// It is not necessary to check the allocation because '(*reconstruct_init)->pulse_length'='PulseLength'(input parameter) has been checked previously
					if (strcmp((*reconstruct_init)->FilterDomain,"T") == 0)		filtergsl= gsl_vector_alloc(resize_mf);
					else if (strcmp((*reconstruct_init)->FilterDomain,"F") == 0)	filtergsl= gsl_vector_alloc(resize_mf*2);
					Pab = gsl_vector_alloc(resize_mf);
					if (numiteration == 0)
					{
						if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
						{
							if (find_optimalfilter((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_optimalfilter for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
						else //(*reconstruct_init)->OFInterp = DAB
						{
							if (find_optimalfilterDAB((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl, &Pab,&Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_optimalfilterDAB for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
					}
					else
					{	
						if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
						{
							
							if (find_optimalfilter(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_optimalfilter for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
						else //(*reconstruct_init)->OFInterp = DAB
						{
							if (find_optimalfilterDAB(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl, &Pab, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_optimalfilterDAB for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
					}
				}
			}
			else if ((strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0) && (strcmp((*reconstruct_init)->OFNoise,"WEIGHTM") == 0))
			{
				// Choose the base-2 system value closest (lower than or equal) to the pulse length
				if (strcmp((*reconstruct_init)->OFStrategy,"BASE2") != 0)
				{
					resize_mf = pow(2,floor(log2(resize_mf)));
					gsl_vector *pulse_aux = gsl_vector_alloc(resize_mf);
                                        //gsl_vector_view(temp4);
					temp = gsl_vector_subvector(pulseToCalculateEnergy,0,resize_mf);
                                        //temp4 = gsl_vector_subvector(pulseToCalculateEnergy,0,resize_mf);
					gsl_vector_memcpy(pulse_aux,&temp.vector);
                                        //gsl_vector_memcpy(pulse_aux,&temp4.vector);
					gsl_vector_free(pulseToCalculateEnergy);
					pulseToCalculateEnergy = gsl_vector_alloc(resize_mf);
					gsl_vector_memcpy(pulseToCalculateEnergy,pulse_aux);
					gsl_vector_free(pulse_aux); pulse_aux = 0;
				}
				
				PRCLOFWM = gsl_matrix_alloc(2,resize_mf);
				if (numiteration == 0)
				{
					if (find_prclofwm((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &PRCLOFWM, &Ealpha, &Ebeta))
					{
						message = "Cannot run routine find_prclofwm";
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}
				else
				{
					if (find_prclofwm(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &PRCLOFWM, &Ealpha, &Ebeta))
					{
						message = "Cannot run routine find_prclofwm";
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}
			}
			else if ((strcmp((*reconstruct_init)->EnergyMethod,"WEIGHT") == 0) || (strcmp((*reconstruct_init)->EnergyMethod,"WEIGHTN") == 0))
			{
				if (numiteration == 0)
				{
					// Get the indexes of the two energies which straddle the pulse
					if (find_Esboundary((*pulsesInRecord)->pulses_detected[i].maxDER,(*reconstruct_init)->library_collection->maxDERs,(*reconstruct_init),&indexEalpha,&indexEbeta,&Ealpha, &Ebeta))
					{
						message = "Cannot run routine find_Esboundary for filter interpolation";
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}
				else
				{
					// Get the indexes of the two energies which straddle the pulse
					if (find_Esboundary(energy,(*reconstruct_init)->library_collection->energies,(*reconstruct_init),&indexEalpha,&indexEbeta,&Ealpha, &Ebeta))
					{
						message = "Cannot run routine find_Esboundary for filter interpolation";
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}
				
				if ((indexEalpha == indexEbeta) && ((*reconstruct_init)->library_collection->ntemplates-1 == indexEalpha))
				{
					message = "Not enough info in the library. At least is necessary to add a new last row with energy higher than the energy of the pulses in the input FITS file";
					EP_EXIT_ERROR(message,EPFAIL);
				}
				
				if ((strcmp((*reconstruct_init)->EnergyMethod,"WEIGHTN") == 0) && ((*reconstruct_init)->OFLib == 1))
				{
					// Choose the base-2 system value closest (lower than or equal) to the pulse length
					if (strcmp((*reconstruct_init)->OFStrategy,"BASE2") != 0)
					{
						resize_mf = pow(2,floor(log2(resize_mf)));
						gsl_vector *pulse_aux = gsl_vector_alloc(resize_mf);
                                                //gsl_vector_view(temp5);
						temp = gsl_vector_subvector(pulseToCalculateEnergy,0,resize_mf);
                                                //temp5 = gsl_vector_subvector(pulseToCalculateEnergy,0,resize_mf);
						gsl_vector_memcpy(pulse_aux,&temp.vector);
                                                //gsl_vector_memcpy(pulse_aux,&temp5.vector);
						gsl_vector_free(pulseToCalculateEnergy);
						pulseToCalculateEnergy = gsl_vector_alloc(resize_mf);
						gsl_vector_memcpy(pulseToCalculateEnergy,pulse_aux);
						gsl_vector_free(pulse_aux); pulse_aux = 0;
					}
					
					// Find the appropriate values of the PRECALWN HDU ('PRCLx' columns)
					PRCLWN = gsl_matrix_alloc(2,resize_mf);
					Pab = gsl_vector_alloc(resize_mf);
					if (numiteration == 0)
					{
						if (find_prclwn((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &PRCLWN, &Pab,&Ealpha, &Ebeta))
						{
							message = "Cannot run routine find_prclwn";
							EP_EXIT_ERROR(message,EPFAIL);
						}
					}
					else
					{
						if (find_prclwn(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &PRCLWN, &Pab,&Ealpha, &Ebeta))
						{
							message = "Cannot run routine find_prclwn";
							EP_EXIT_ERROR(message,EPFAIL);
						}
					}
				}
			}
			
			gsl_matrix_set(Estraddle,0,numiteration,Ealpha);
			gsl_matrix_set(Estraddle,1,numiteration,Ebeta);

			if ((strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0) && (strcmp((*reconstruct_init)->OFNoise,"NSD") == 0))
			{
				if ((*reconstruct_init)->OFLib == 0)
				{
					// Calculate the optimal filter
					if (calculus_optimalFilter (TorF, (*reconstruct_init)->intermediate, (*reconstruct_init)->opmode, filtergsl, filtergsl->size, 1.0/record->delta_t, runF0orB0val, (*reconstruct_init)->noise_spectrum->noisefreqs, (*reconstruct_init)->noise_spectrum->noisespec, &optimalfilter, &optimalfilter_f, &optimalfilter_FFT, &optimalfilter_FFT_complex))
					{
						message = "Cannot run routine calculus_optimalFilter";
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}
				else if ((*reconstruct_init)->OFLib == 1)
				{
                         		// It is not necessary to check the allocation because 'filtergsl' size has been checked previously
					optimalfilter = gsl_vector_alloc(filtergsl->size);
					gsl_vector_set_all(optimalfilter,0);
					if (TorF == 0)     gsl_vector_memcpy(optimalfilter,filtergsl);
					else if (TorF == 1)
					{
						// It is not necessary to check the allocation because 'filtergsl' size has been checked previously
						optimalfilter_FFT_complex = gsl_vector_complex_alloc(filtergsl->size/2);
                                         	for (int k=0;k<filtergsl->size/2;k++)
						{
                                                    gsl_vector_complex_set(optimalfilter_FFT_complex,k,gsl_complex_rect(gsl_vector_get(filtergsl,k),gsl_vector_get(filtergsl,k+filtergsl->size/2)));
						}
						gsl_vector_complex_set(optimalfilter_FFT_complex,0,gsl_complex_rect(0.0,0.0));
						for (int k=0;k<filtergsl->size/2;k++)
						{
                                                    gsl_vector_complex_set(optimalfilter_FFT_complex,k,gsl_complex_conjugate(gsl_vector_complex_get(optimalfilter_FFT_complex,k)));
						}
					}
				}
			}
			
			// Calculate the energy of each pulse
			if (calculateEnergy(pulseToCalculateEnergy,pulseGrade,optimalfilter,optimalfilter_FFT_complex,runEMethod,indexEalpha,indexEbeta,(*reconstruct_init),TorF,1/record->delta_t,Pab,PRCLWN,PRCLOFWM,&energy,&tstartNewDev,&lagsShift,0))
			{
				message = "Cannot run calculateEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
				EP_EXIT_ERROR(message,EPFAIL);
			}
			gsl_vector_free(pulseToCalculateEnergy); pulseToCalculateEnergy = 0;
                        
                        //cout<<"Tstart0: "<<endl;
                        //std::cout << std::setprecision(15) << (*pulsesInRecord)->pulses_detected[i].Tstart*156250 << '\n';
                        //cout<<"Tstart0: "<<(*pulsesInRecord)->pulses_detected[i].Tstart<<endl;
			// If using lags, it is necessary to modify the tstart of the pulse and the length of the filter used
			if ((strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0) && tstartNewDev != -999.0)
			{
				(*pulsesInRecord)->pulses_detected[i].Tstart = (*pulsesInRecord)->pulses_detected[i].Tstart + tstartNewDev*record->delta_t; // In seconds
				//if (tstartNewDev != 0) (*pulsesInRecord)->pulses_detected[i].grade1 = (*pulsesInRecord)->pulses_detected[i].grade1 - floor(fabs(tstartNewDev)+0.5);
			}
			//cout<<"Tstart1: "<<(*pulsesInRecord)->pulses_detected[i].Tstart<<endl;
			//cout<<"Tstart1: "<<endl;
                        //std::cout << std::setprecision(15) << (*pulsesInRecord)->pulses_detected[i].Tstart*156250 << '\n';
                        //cout<<"pulses_detected[i].Tstart1: "<<(*pulsesInRecord)->pulses_detected[i].Tstart<<endl;
			
			gsl_matrix_set(resultsE,0,numiteration,energy);
			gsl_matrix_set(resultsE,1,numiteration,min(fabs(energy-Ealpha),fabs(energy-Ebeta)));

			gsl_vector_view(temp);
                        //gsl_vector_view(temp6);
			gsl_vector *subresultsE;
			if ((Ealpha != Ebeta) && ((energy < Ealpha) || (energy > Ebeta)))
			{
				if( numiteration != 0)
				{
					for (int j=0; j<numiteration+1; j++)
					{
						if ((Ealpha == gsl_matrix_get(Estraddle,0,j)) && (Ebeta == gsl_matrix_get(Estraddle,1,j)))
						{
							iterate = false;
							
							subresultsE = gsl_vector_alloc(resultsE->size2);
							gsl_matrix_get_row(subresultsE,resultsE,1);
							temp = gsl_vector_subvector(subresultsE,0,numiteration+1);
                                                        //temp6 = gsl_vector_subvector(subresultsE,0,numiteration+1);
							gsl_vector *SUBsubresultsE = gsl_vector_alloc(numiteration+1);
							gsl_vector_memcpy(SUBsubresultsE,&temp.vector);
                                                        //gsl_vector_memcpy(SUBsubresultsE,&temp6.vector);
							energy = gsl_matrix_get(resultsE,0,gsl_vector_min_index(SUBsubresultsE));
							gsl_vector_free(subresultsE); subresultsE = 0;
							gsl_vector_free(SUBsubresultsE); SUBsubresultsE = 0;
							
							break;
						}
					}
				}
			}
			else 
			{	
				iterate = false;
			
				if (numiteration != 0)
				{
					subresultsE = gsl_vector_alloc(resultsE->size2);
					gsl_matrix_get_row(subresultsE,resultsE,1);
					temp = gsl_vector_subvector(subresultsE,0,numiteration+1);
                                        //temp6 = gsl_vector_subvector(subresultsE,0,numiteration+1);
					gsl_vector *SUBsubresultsE = gsl_vector_alloc(numiteration+1);
					gsl_vector_memcpy(SUBsubresultsE,&temp.vector);
                                        //gsl_vector_memcpy(SUBsubresultsE,&temp6.vector);
					energy = gsl_matrix_get(resultsE,0,gsl_vector_min_index(SUBsubresultsE));
					gsl_vector_free(subresultsE); subresultsE = 0;
					gsl_vector_free(SUBsubresultsE); SUBsubresultsE = 0;
				}
			  
			}
		} while (iterate);
                //gsl_vector_free(pulseToCalculateEnergy);

		gsl_matrix_free(Estraddle); Estraddle = 0;
		gsl_matrix_free(resultsE); resultsE = 0;

		// Subtract the pulse model from the record
		if (find_model_energies(energy, (*reconstruct_init), &model))
		{
			message = "Cannot run find_model_energies routine for pulse i=" + boost::lexical_cast<std::string>(i);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		tstartJITTER = ((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t;
		shift = tstartJITTER - tstartSamplesRecord;
		//cout<<"shiftFIN: "<<shift<<endl;
		// In order to subtract the pulse model, it has to be located in the tstart with jitter and know its values in the digitized samples
		gsl_vector *modelToSubtract = gsl_vector_alloc(model->size);
		for (int j=0;j<model->size;j++)
                {
                    if (shift < 0)
                    {
                        if (j != model->size-1)
                            gsl_vector_set(modelToSubtract,j,(gsl_vector_get(model,j+1)-gsl_vector_get(model,j))*(-shift)+gsl_vector_get(model,j));
                        
                        else 
                            gsl_vector_set(model,j,gsl_vector_get(model,j)); //?????????????????????
                    }
                    else if (shift > 0)
                    {
                        if (j == 0)
                        {
                            gsl_vector_set(modelToSubtract,j,(gsl_vector_get(model,j)-0)*(1-shift)+0);
                        }
                        else 
                        {
                            gsl_vector_set(modelToSubtract,j,(gsl_vector_get(model,j)-gsl_vector_get(model,j-1))*(1-shift)+gsl_vector_get(model,j-1));
                        }
                    }
                    else
                    {
                        gsl_vector_memcpy(modelToSubtract,model);
                    }
                }
                gsl_vector_memcpy(model,modelToSubtract);
                gsl_vector_free(modelToSubtract);

		//for (int j=tstartSamplesRecord;j<min(tstartSamplesRecord+(model->size),(double) record->trigger_size);j++)
                for (int j=tstartSamplesRecord;j<min((double) tstartSamplesRecord+(model->size),(double) record->trigger_size);j++)
		{
			gsl_vector_set(recordAux,j,gsl_vector_get(recordAux,j)-gsl_vector_get(model,j-tstartSamplesRecord));
		}

		// Write info of the pulse in the output intemediate file if 'intermediate'=1
		if ((*reconstruct_init)->intermediate == 1)
		{
			if (writeFilterHDU(reconstruct_init, i,energy, optimalfilter, optimalfilter_f, optimalfilter_FFT, &dtcObject))
			{
				message = "Cannot run writeFilterHDU routine for pulse i=" + boost::lexical_cast<std::string>(i);
				EP_EXIT_ERROR(message,EPFAIL);
			}
		}
		
		(*pulsesInRecord)->pulses_detected[i].energy = energy/1e3;	// In SIXTE, SIGNAL is in keV
		(*pulsesInRecord)->pulses_detected[i].E_lowres = energy_lowres/1e3;
		(*pulsesInRecord)->pulses_detected[i].grading = pulseGrade;	
                ////(*pulsesInRecord)->pulses_detected[i].phi = tstartNewDev - lagsShift;
                //(*pulsesInRecord)->pulses_detected[i].phi = tstartNewDev;
                //(*pulsesInRecord)->pulses_detected[i].lagsShift = lagsShift;
                double intpart;
                (*pulsesInRecord)->pulses_detected[i].phi = modf(tstartNewDev,&intpart);    // fractpart=modf(param,&intpart) Se obtiene la parte entera y                                                      decimal
                (*pulsesInRecord)->pulses_detected[i].lagsShift = lagsShift+intpart;
                //cout<<"phi= "<<modf(tstartNewDev,&intpart)<<endl;
                //cout<<"lagsShift= "<<lagsShift+intpart<<endl;
                //cout<<"phi+n= "<<(*pulsesInRecord)->pulses_detected[i].phi+(*pulsesInRecord)->pulses_detected[i].lagsShift<<endl;
                
		// Free allocated GSL vectors
		gsl_vector_free(optimalfilter); optimalfilter = 0;
		gsl_vector_free(optimalfilter_f); optimalfilter_f = 0;
		gsl_vector_free(optimalfilter_FFT); optimalfilter_FFT = 0;
		if ((*pulsesInRecord)->pulses_detected[i].quality < 10)
		{
                        gsl_vector_free(pulse); pulse = 0;
                        gsl_vector_free(filtergsl); filtergsl = 0;
                        gsl_vector_complex_free(optimalfilter_FFT_complex); optimalfilter_FFT_complex = 0;
                        gsl_vector_free(Pab); Pab = 0;
                        gsl_matrix_free(PRCLWN); PRCLWN = 0;
                        gsl_matrix_free(PRCLOFWM); PRCLOFWM = 0;
		}
	} // End for
	
	gsl_vector_free(recordAux); recordAux = 0;
	gsl_vector_free(model); model = 0;
        
        gsl_vector_free(pulse_lowres); pulse_lowres = 0;
        gsl_vector_free(filtergsl_lowres); filtergsl_lowres = 0;
        if (PRCLWN_lowres != NULL) gsl_matrix_free(PRCLWN_lowres); PRCLWN_lowres = 0;
        if (PRCLOFWM_lowres != NULL) gsl_matrix_free(PRCLOFWM_lowres); PRCLOFWM_lowres = 0;
        gsl_vector_free(optimalfilter_lowres); optimalfilter_lowres = 0;
        if (optimalfilter_FFT_complex_lowres != NULL) gsl_vector_complex_free(optimalfilter_FFT_complex_lowres); optimalfilter_FFT_complex_lowres = 0;

	return;
}
/*xxxx end of SECTION B xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

void th_runEnergy(TesRecord* record, 
                  ReconstructInitSIRENA** reconstruct_init, 
                  PulsesCollection** pulsesInRecord, 
                  OptimalFilterSIRENA **optimalFilter)
{
	// Declare variables
	string message="";
	int status = EPOK;
	char valERROR[256];

	fitsfile *dtcObject = NULL;	    // FITS object containing information of the intermediate output FITS file
	if ((*reconstruct_init)->intermediate == 1)
	{
		char dtcName[256];
		strncpy(dtcName,(*reconstruct_init)->detectFile,255);
		dtcName[255]='\0';
	}

	int TorF;
	if (strcmp((*reconstruct_init)->FilterDomain,"T") == 0)		// Time
	{
		TorF=0;
	}
	else if (strcmp((*reconstruct_init)->FilterDomain,"F") == 0)	// Frequency
	{
		TorF=1;
	}

	int runF0orB0val;
	if (strcmp((*reconstruct_init)->FilterMethod,"F0") == 0)	// Deleting the frequency-zero bin
	{
		runF0orB0val = 0;
	}
	else if (strcmp((*reconstruct_init)->FilterMethod,"B0") == 0)	// Working without baseline
	{
		runF0orB0val = 1;
	}

	// I2R, I2RALL, I2RNOL or I2RFITTED methods convert I into R at the beginnig and after that 'I2R', 'I2RALL', 'I2RNOL' or 'I2RFITTED' are equivalent to 'OPTFILT'
	int runEMethod;
	if (strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0)
	{
		runEMethod = 0;
	}
	else if (strcmp((*reconstruct_init)->EnergyMethod,"WEIGHT") == 0)
	{
		runEMethod = 1;
	}
	else if (strcmp((*reconstruct_init)->EnergyMethod,"WEIGHTN") == 0)
	{
		runEMethod = 2;
	}

	int OFlength_strategy;
	if (strcmp((*reconstruct_init)->OFStrategy,"FREE") == 0)
	{
		OFlength_strategy = 0;
	}
	else if (strcmp((*reconstruct_init)->OFStrategy,"BASE2") == 0)
	{
		OFlength_strategy = 1;
	}
	else if (strcmp((*reconstruct_init)->OFStrategy,"BYGRADE") == 0)
	{
		OFlength_strategy = 2;
	}
	else if (strcmp((*reconstruct_init)->OFStrategy,"FIXED") == 0)
	{
		OFlength_strategy = 3;
	}

	long energyInLibrary_row;

	double energy;
	double tstartNewDev = -999.0;    	// Deviation of the starting of the pulses (in samples) respect to the tstart calculated
	//int numlags = 3; 			// Must be odd
	int numlags = (*reconstruct_init)->nLags; 			// Must be odd
	int numlags2 = floor(numlags/2);
	int lagsShift = -999;                   // Number of samples shifted to find the maximum of the parabola

	double Ealpha, Ebeta;
	gsl_vector *optimalfilter = NULL;	// Resized optimal filter expressed in the time domain (optimalfilter(t))
	gsl_vector *optimalfilter_f = NULL;	// Resized optimal filter f's when f's are according to [0,...fmax,-fmax,...] (frequency domain)
	gsl_vector *optimalfilter_FFT = NULL;	// Resized optimal filter magnitudes when f's are according to [0,...fmax,-fmax,...] (frequency domain)
	gsl_vector_complex *optimalfilter_FFT_complex = NULL;
	
	gsl_vector *Pab = NULL;
	gsl_matrix *PRCLWN = NULL;
	gsl_matrix *PRCLOFWM = NULL;

	gsl_vector *model;

	gsl_vector *recordAux;

	gsl_vector *pulse = NULL;
	gsl_vector *pulseToCalculateEnergy = NULL;	// Just in case LagsOrNot=1

	gsl_vector *filtergsl = NULL;		// Matched filter values (time domain)

	int iter;
	gsl_vector_view(temp);
        
        double tstartSamplesRecord;		// Tstart of the pulse in samples from the beginning of the record
	double tstartRecord;			// Tstart of the record in seconds
	int tstartRecordSamples = record->time/record->delta_t;	// Tstart of the record in samples (nearest sample)
        //cout<<"tstartRecordSamples: "<<tstartRecordSamples<<endl;
        double tstartJITTER;
        gsl_matrix *optimalfilterMAT;
        gsl_matrix_complex *optimalfilter_FFT_complexMAT;
        int tstartSamplesRecordStartINT;
        double tstartSamplesRecordStartDOUBLE;
        
        double shift;

	int indexEalpha = 0;
	int indexEbeta = 0;

	long resize_mf;
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double energy_lowres;
        long resize_mf_lowres = 4; // In order to get the low resolution energy estimator by filtering with a 4-samples-long filter
        gsl_vector *pulse_lowres = gsl_vector_alloc(4);
        gsl_vector *filtergsl_lowres = NULL;
        if (strcmp((*reconstruct_init)->FilterDomain,"T") == 0)		filtergsl_lowres= gsl_vector_alloc(resize_mf_lowres);
        else if (strcmp((*reconstruct_init)->FilterDomain,"F") == 0)	filtergsl_lowres= gsl_vector_alloc(resize_mf_lowres*2);
        gsl_vector *Pab_lowres = gsl_vector_alloc(resize_mf_lowres);
        gsl_matrix *PRCLWN_lowres = NULL;
	gsl_matrix *PRCLOFWM_lowres = NULL;
        double Ealpha_lowres, Ebeta_lowres;
        gsl_vector *optimalfilter_lowres = gsl_vector_alloc(filtergsl_lowres->size);	// Resized optimal filter expressed in the time domain (optimalfilter(t))
	gsl_vector_complex *optimalfilter_FFT_complex_lowres = NULL;
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
	int pulseGrade; 			// Pileup=-2, Rejected=-1, HighRes=1, MidRes=2, LimRes=3, LowRes=4
	
	// Store the record in 'invector'
	// It is not necessary to check the allocation because 'record->trigger_size' has been checked previously
	recordAux = gsl_vector_alloc(record->trigger_size);
	if (loadRecord(record, &tstartRecord, &recordAux))
	{
		message = "Cannot run routine loadRecord";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	// Subtract the baseline if OPTFILT and 'runF0orB0val'= 1 ('FilterMethod'=B0)
	// Subtract the baseline if WEIGHT
	if (((runF0orB0val == 1) && (runEMethod == 0)) || (runEMethod == 1))
	{
		// It is not necessary to check the allocation because the allocation of 'recordAux' has been checked previously
		gsl_vector *baselinegsl = gsl_vector_alloc(recordAux->size);
		gsl_vector_set_all(baselinegsl,-1.0*(*reconstruct_init)->noise_spectrum->baseline);
		gsl_vector_add(recordAux,baselinegsl);
		gsl_vector_free(baselinegsl); baselinegsl = 0;
	}

	// Check Quality
	iter = 0;
	for (int i=0; i<(*pulsesInRecord)->ndetpulses;i++)
	{
		if ((*pulsesInRecord)->pulses_detected[i].quality >= 10)  // Saturated pulse
		{
			iter++;
		}
	}
	if (iter == (*pulsesInRecord)->ndetpulses)	
	{
		message = "There are no unsaturated pulses in one record";
		EP_PRINT_ERROR(message,-999); // Only a warning
	}

	model =gsl_vector_alloc((*reconstruct_init)->pulse_length);
        
        double valaux;
        int newidx;
        	  
        int extraSizeDueToLags = 0;
	for (int i=0; i<(*pulsesInRecord)->ndetpulses ;i++)
	{
		// Establish the pulse grade (HighRes=1, MidRes=2, LimRes=3, LowRes=4, Rejected=-1, Pileup=-2) and the optimal filter length
		if ((*pulsesInRecord)->pulses_detected[i].quality == 1)		(*pulsesInRecord)->pulses_detected[i].grade1 = -1;
		else								(*pulsesInRecord)->pulses_detected[i].grade1 = (*pulsesInRecord)->pulses_detected[i].pulse_duration;      
                pulseGrade = 0;
		if (pulseGrading(*reconstruct_init,(*pulsesInRecord)->pulses_detected[i].grade1,(*pulsesInRecord)->pulses_detected[i].grade2,OFlength_strategy,&pulseGrade,&resize_mf))
		{
			message = "Cannot run routine pulseGrading";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		(*pulsesInRecord)->pulses_detected[i].grade1 = resize_mf;

		// Pulse: Load the proper piece of the record in 'pulse'
                //cout<<"((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t: "<<((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t<<endl;
                //tstartSamplesRecord = ((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t;
                tstartSamplesRecord = (*pulsesInRecord)->pulses_detected[i].TstartSamples;
                //cout<<"tstartSamplesRecord: "<<tstartSamplesRecord<<endl;
                //tstartSamplesRecordStartINT = ((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t-numlags2;
                //tstartSamplesRecordStartDOUBLE = ((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t-numlags2;
                tstartSamplesRecordStartINT = tstartSamplesRecord-numlags2;
                tstartSamplesRecordStartDOUBLE = tstartSamplesRecord-numlags2;
                //tstartSamplesRecordStartDOUBLE = tstartSamplesRecord-(numlags-1);
                
                //cout<<"tstart: "<<tstartSamplesRecord<<endl;
                //cout<<"tstartRecord: "<<tstartRecord<<endl;
                //cout<<"tstartRecordSamples: "<<tstartRecordSamples<<endl;
		if ((pulse = gsl_vector_alloc(resize_mf)) == 0)
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);
			message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		if ((tstartSamplesRecord < 0) || (tstartSamplesRecord > recordAux->size-2)
			|| (resize_mf < 1) || (resize_mf > recordAux->size-tstartSamplesRecord))
		{
			sprintf(valERROR,"%d",__LINE__+5);
			string str(valERROR);
			message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL); 
		}
		//gsl_vector_view(temp1);
		temp = gsl_vector_subvector(recordAux,tstartSamplesRecord,resize_mf);
                //temp1 = gsl_vector_subvector(recordAux,tstartSamplesRecord,resize_mf);
		if (gsl_vector_memcpy(pulse, &temp.vector) != 0)
                //if (gsl_vector_memcpy(pulse, &temp1.vector) != 0)
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);	
			message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		//cout<<"tstartRecordSamples: "<<tstartRecordSamples<<endl;
		//cout<<"pulses_detected[i].Tstart: "<<pulses_detected[i].Tstart<<endl;
		tstartJITTER = ((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t;
		shift = tstartJITTER - tstartSamplesRecord;
                //cout<<"tstartJITTER: "<<tstartJITTER<<endl;
                //cout<<"shift: "<<shift<<endl;
                
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //////////// In order to get the low resolution energy estimator by filtering with a 4-samples-long filter ///////////////////
                // Pulse 
                if (resize_mf_lowres > recordAux->size-tstartSamplesRecord)
		{
			sprintf(valERROR,"%d",__LINE__+5);
			string str(valERROR);
			message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL); 
		}
		temp = gsl_vector_subvector(recordAux,tstartSamplesRecord,resize_mf_lowres);
		if (gsl_vector_memcpy(pulse_lowres, &temp.vector) != 0)
		{
			sprintf(valERROR,"%d",__LINE__-2);
			string str(valERROR);	
			message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		// Filter
		if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
                {
                    if (find_optimalfilter((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl_lowres, &Ealpha_lowres, &Ebeta_lowres))
                    {
                        message = "Cannot run routine find_optimalfilter for filter interpolation";
                        EP_EXIT_ERROR(message,EPFAIL);
                    }
                }
                else
                {
                    if (find_optimalfilterDAB((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl_lowres, &Pab_lowres,&Ealpha_lowres, &Ebeta_lowres))
                    {
                        message = "Cannot run routine find_optimalfilterDAB for filter interpolation";
                        EP_EXIT_ERROR(message,EPFAIL);
                    }
                }
                gsl_vector_set_all(optimalfilter_lowres,0);
                if (TorF == 0)     gsl_vector_memcpy(optimalfilter_lowres,filtergsl_lowres);
                else if (TorF == 1)
                {
                    // It is not necessary to check the allocation because 'filtergsl' size has been checked previously
                    optimalfilter_FFT_complex_lowres = gsl_vector_complex_alloc(filtergsl_lowres->size/2);
                    for (int k=0;k<filtergsl_lowres->size/2;k++)
                    {
                        gsl_vector_complex_set(optimalfilter_FFT_complex_lowres,k,gsl_complex_rect(gsl_vector_get(filtergsl_lowres,k),gsl_vector_get(filtergsl_lowres,k+filtergsl_lowres->size/2)));
                    }
                    gsl_vector_complex_set(optimalfilter_FFT_complex_lowres,0,gsl_complex_rect(0.0,0.0));
                    for (int k=0;k<filtergsl_lowres->size/2;k++)
                    {
                        gsl_vector_complex_set(optimalfilter_FFT_complex_lowres,k,gsl_complex_conjugate(gsl_vector_complex_get(optimalfilter_FFT_complex_lowres,k)));
                    }
                }
                // Calculate the low resolution estimator
                if (calculateEnergy(pulse_lowres,1,optimalfilter_lowres,optimalfilter_FFT_complex_lowres,0,0,0,(*reconstruct_init),TorF,1/record->delta_t,Pab_lowres,PRCLWN_lowres,PRCLOFWM_lowres,&energy_lowres,&tstartNewDev,&lagsShift,1))
                {
                    message = "Cannot run calculateEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
                    EP_EXIT_ERROR(message,EPFAIL);
                }
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		if ((*reconstruct_init)->LagsOrNot == 0)	
		{
			pulseToCalculateEnergy = gsl_vector_alloc(pulse->size);
			gsl_vector_memcpy(pulseToCalculateEnergy,pulse);
		}
		else if (((*reconstruct_init)->LagsOrNot == 1) && (strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0))	
		{
			int resize_mfNEW = resize_mf + numlags -1;
			if ((pulseToCalculateEnergy = gsl_vector_alloc(resize_mfNEW)) == 0)
			{
				sprintf(valERROR,"%d",__LINE__-2);
				string str(valERROR);
				message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
				EP_EXIT_ERROR(message,EPFAIL);
			}
			gsl_vector_set_all(pulseToCalculateEnergy,-999);
			
			//if ((tstartSamplesRecord-numlags/2 < 0) || (tstartSamplesRecord-numlags/2 > recordAux->size-2)
			//	|| (resize_mfNEW < 1) || (resize_mfNEW > recordAux->size-tstartSamplesRecord+numlags/2))
                        if ((tstartSamplesRecordStartDOUBLE < 0) || (tstartSamplesRecordStartDOUBLE > recordAux->size-2)
				|| (resize_mfNEW < 1) || (resize_mfNEW > recordAux->size-tstartSamplesRecordStartDOUBLE))
			{
				sprintf(valERROR,"%d",__LINE__+5);
				string str(valERROR);
				message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";
				EP_EXIT_ERROR(message,EPFAIL); 
			}
			
			//gsl_vector_view(temp1);
			//temp = gsl_vector_subvector(recordAux,tstartSamplesRecord-numlags/2,resize_mfNEW);
                        temp = gsl_vector_subvector(recordAux,tstartSamplesRecordStartDOUBLE,resize_mfNEW);
                        //cout<<"resize_mfNEW: "<<resize_mfNEW<<endl;
                        //cout<<"recordAux: "<<endl;
                        //for (int kk=0;kk<recordAux->size;kk++) cout<<kk<<" "<<gsl_vector_get(recordAux,kk)<<endl;
                        /*cout<<"995 "<<gsl_vector_get(recordAux,995)<<endl;
                        cout<<"996 "<<gsl_vector_get(recordAux,996)<<endl;
                        cout<<"997 "<<gsl_vector_get(recordAux,997)<<endl;
                        cout<<"998 "<<gsl_vector_get(recordAux,998)<<endl;
                        cout<<"999 "<<gsl_vector_get(recordAux,999)<<endl;
                        cout<<"1000 "<<gsl_vector_get(recordAux,1000)<<endl;
                        cout<<"1001 "<<gsl_vector_get(recordAux,1001)<<endl;
                        cout<<"1002 "<<gsl_vector_get(recordAux,1002)<<endl;
                        cout<<"1003 "<<gsl_vector_get(recordAux,1003)<<endl;
                        cout<<"1004 "<<gsl_vector_get(recordAux,1004)<<endl;
                        cout<<"1005 "<<gsl_vector_get(recordAux,1005)<<endl;*/
			if (gsl_vector_memcpy(pulseToCalculateEnergy, &temp.vector) != 0)
			{
				sprintf(valERROR,"%d",__LINE__-2);
				string str(valERROR);	
				message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
				EP_EXIT_ERROR(message,EPFAIL);
			}
                        
                        
                        //cout<<"pulseToCalculateEnergy0: "<<endl;
                        //for (int kk=0;kk<pulseToCalculateEnergy->size;kk++) cout<<kk<<" "<<gsl_vector_get(pulseToCalculateEnergy,kk)<<endl;
                        
                        //cout<<"tstartSamplesRecordStartINT: "<<tstartSamplesRecordStartINT<<endl; 
                        //cout<<"tstartSamplesRecordStartDOUBLE: "<<tstartSamplesRecordStartDOUBLE<<endl; 
                        
                        //cout<<"cicloFOR"<<endl;
                        /*for (int kk=0;kk<resize_mfNEW;kk++)
                        {
                            //if ((j+tstartSamplesRecord-numlags/2>=995) && (j+tstartSamplesRecord-numlags/2<1002)) cout<<"recordAux: "<<j+tstartSamplesRecord-numlags/2<<" "<<gsl_vector_get(recordAux,j+tstartSamplesRecord-numlags/2)<<endl;
                            //gsl_vector_set(pulseToCalculateEnergy,j,gsl_vector_get(recordAux,j+tstartSamplesRecord-numlags/2));
                            //if ((j+tstartSamplesRecord-numlags/2>=995) && (j+tstartSamplesRecord-numlags/2<1002)) cout<<"pulseToCalculateEnergy: "<<j<<" "<<gsl_vector_get(pulseToCalculateEnergy,j)<<endl;
                            //if ((kk+tstartSamplesRecordStartDOUBLE>=995) && (kk+tstartSamplesRecordStartDOUBLE<1002)) cout<<kk<<" "<<tstartSamplesRecordStartDOUBLE<<" recordAux: "<<kk+tstartSamplesRecordStartDOUBLE<<" "<<gsl_vector_get(recordAux,kk+tstartSamplesRecordStartDOUBLE)<<endl;
                            //gsl_vector_set(pulseToCalculateEnergy,kk,gsl_vector_get(recordAux,kk+tstartSamplesRecordStartDOUBLE));
                            
                            newidx=kk+tstartSamplesRecordStartDOUBLE;
                            valaux = gsl_vector_get(recordAux,newidx);
                            //cout<<kk<<" "<<tstartSamplesRecordStartDOUBLE<<" "<<newidx<<" valaux: "<<valaux<<" "<<gsl_vector_get(recordAux,1000)<<endl;
                            gsl_vector_set(pulseToCalculateEnergy,kk,valaux);
                            //if (kk == 5)    {cout<<"Hola"<<endl; gsl_vector_set(pulseToCalculateEnergy,5,2305);}
                            //cout<<gsl_vector_get(recordAux,1000)<<" "<<kk<<" "<<tstartSamplesRecordStartDOUBLE<<" "<<" "<<kk+tstartSamplesRecordStartDOUBLE<<" "<<gsl_vector_get(recordAux,kk+tstartSamplesRecordStartDOUBLE)<<" "<<gsl_vector_get(pulseToCalculateEnergy,kk)<<endl;
                            //if ((kk+tstartSamplesRecordStartDOUBLE>=995) && (kk+tstartSamplesRecordStartDOUBLE<1002)) cout<<kk<<" "<<tstartSamplesRecordStartDOUBLE<<" pulseToCalculateEnergy: "<<kk<<" "<<gsl_vector_get(pulseToCalculateEnergy,kk)<<endl;
                        }*/
                        
			extraSizeDueToLags = numlags-1;
                        //cout<<"vector: "<<gsl_vector_get(pulseToCalculateEnergy,995)<<","<<gsl_vector_get(pulseToCalculateEnergy,996)<<","<<gsl_vector_get(pulseToCalculateEnergy,997)<<","<<gsl_vector_get(pulseToCalculateEnergy,998)<<","<<gsl_vector_get(pulseToCalculateEnergy,999)<<","<<gsl_vector_get(pulseToCalculateEnergy,1000)<<endl;
                        //cout<<"pulseToCalculateEnergy: "<<endl;
                        //for (int kk=0;kk<pulseToCalculateEnergy->size;kk++) cout<<kk<<" "<<gsl_vector_get(pulseToCalculateEnergy,kk)<<endl;
                        /*cout<<"0 "<<gsl_vector_get(pulseToCalculateEnergy,0)<<endl;
                        cout<<"1 "<<gsl_vector_get(pulseToCalculateEnergy,1)<<endl;
                        cout<<"2 "<<gsl_vector_get(pulseToCalculateEnergy,2)<<endl;
                        cout<<"3 "<<gsl_vector_get(pulseToCalculateEnergy,3)<<endl;
                        cout<<"4 "<<gsl_vector_get(pulseToCalculateEnergy,4)<<endl;
                        cout<<"5 "<<gsl_vector_get(pulseToCalculateEnergy,5)<<endl;
                        cout<<"6 "<<gsl_vector_get(pulseToCalculateEnergy,6)<<endl;
                        cout<<"7 "<<gsl_vector_get(pulseToCalculateEnergy,7)<<endl;
                        cout<<"8 "<<gsl_vector_get(pulseToCalculateEnergy,8)<<endl;
                        cout<<"9 "<<gsl_vector_get(pulseToCalculateEnergy,9)<<endl;
                        cout<<"10 "<<gsl_vector_get(pulseToCalculateEnergy,10)<<endl;*/
		}
		
		bool iterate = true;
		if ((*reconstruct_init)-> OFIter == 1)	iterate = true;
		else iterate = false;
		// It is not necessary to check the allocation because '(*reconstruct_init)->library_collection->ntemplates' has been check previously
		gsl_matrix *Estraddle = gsl_matrix_alloc(2,(*reconstruct_init)->library_collection->ntemplates);
		gsl_matrix *resultsE = gsl_matrix_alloc(2,(*reconstruct_init)->library_collection->ntemplates);		// Row0 -> calculatedEnergy
															// Row1 -> min[abs(calculatedEnergy-Ealpha),abs(calculatedEnergy-Ebeta)]
		int numiteration = -1;
                
		do
		{
			numiteration++;
			
			if ((strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0) && (strcmp((*reconstruct_init)->OFNoise,"NSD") == 0))
			{
				// Filter (find the matched filter and load it in 'filter')
				if ((*reconstruct_init)->OFLib == 0)
				{	
					// It is not necessary to check the allocation because '(*reconstruct_init)->pulse_length'='PulseLength'(input parameter) has been checked previously
					filtergsl= gsl_vector_alloc(resize_mf);
					Pab = gsl_vector_alloc(resize_mf);
					if (numiteration == 0)
					{
						if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
						{
							if (find_matchedfilter(runF0orB0val, (*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_matchedfilter for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
						else //(*reconstruct_init)->OFInterp = DAB
						{
							if (find_matchedfilterDAB((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl, &Pab, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_matchedfilterDAB for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
					}
					else
					{
						if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
						{
							if (find_matchedfilter(runF0orB0val, energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_matchedfilter for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
						else //(*reconstruct_init)->OFInterp = DAB
						{
							if (find_matchedfilterDAB(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl, &Pab, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_matchedfilterDAB for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
					}
				}
				else if ((*reconstruct_init)->OFLib == 1)
				{
                    
                                        // Choose the base-2 system value closest (lower than or equal) to the pulse length
					if ((strcmp((*reconstruct_init)->OFStrategy,"BASE2") != 0) 
					    && ((((*reconstruct_init)->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != -999) && (resize_mf != (*reconstruct_init)->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration))
					    && (resize_mf != (*reconstruct_init)->library_collection->pulse_templates[0].template_duration)))
					{                                             
                                                resize_mf = pow(2,floor(log2(resize_mf)));
                                                gsl_vector *pulse_aux = gsl_vector_alloc(resize_mf+extraSizeDueToLags);
                                                //gsl_vector_view(temp3);
                                                temp = gsl_vector_subvector(pulseToCalculateEnergy,0,resize_mf+extraSizeDueToLags);
                                                //temp3 = gsl_vector_subvector(pulseToCalculateEnergy,0,resize_mf+extraSizeDueToLags);
                                                gsl_vector_memcpy(pulse_aux,&temp.vector);
                                                //gsl_vector_memcpy(pulse_aux,&temp3.vector);
                                                gsl_vector_free(pulseToCalculateEnergy);
                                                pulseToCalculateEnergy = gsl_vector_alloc(resize_mf+extraSizeDueToLags);
                                                gsl_vector_memcpy(pulseToCalculateEnergy,pulse_aux);
                                                gsl_vector_free(pulse_aux);
                                                //>>>>>>> master
					}
				  
					// It is not necessary to check the allocation because '(*reconstruct_init)->pulse_length'='PulseLength'(input parameter) has been checked previously
					if (strcmp((*reconstruct_init)->FilterDomain,"T") == 0)		filtergsl= gsl_vector_alloc(resize_mf);
					else if (strcmp((*reconstruct_init)->FilterDomain,"F") == 0)	filtergsl= gsl_vector_alloc(resize_mf*2);
					Pab = gsl_vector_alloc(resize_mf);
					if (numiteration == 0)
					{
						if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
						{
							if (find_optimalfilter((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_optimalfilter for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
						else //(*reconstruct_init)->OFInterp = DAB
						{
							if (find_optimalfilterDAB((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl, &Pab,&Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_optimalfilterDAB for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
					}
					else
					{	
						if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
						{
							
							if (find_optimalfilter(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_optimalfilter for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
						else //(*reconstruct_init)->OFInterp = DAB
						{
							if (find_optimalfilterDAB(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl, &Pab, &Ealpha, &Ebeta))
							{
								message = "Cannot run routine find_optimalfilterDAB for filter interpolation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
						}
					}
				}
			}
			else if ((strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0) && (strcmp((*reconstruct_init)->OFNoise,"WEIGHTM") == 0))
			{
				// Choose the base-2 system value closest (lower than or equal) to the pulse length
				if (strcmp((*reconstruct_init)->OFStrategy,"BASE2") != 0)
				{
					resize_mf = pow(2,floor(log2(resize_mf)));
					gsl_vector *pulse_aux = gsl_vector_alloc(resize_mf);
                                        //gsl_vector_view(temp4);
					temp = gsl_vector_subvector(pulseToCalculateEnergy,0,resize_mf);
                                        //temp4 = gsl_vector_subvector(pulseToCalculateEnergy,0,resize_mf);
					gsl_vector_memcpy(pulse_aux,&temp.vector);
                                        //gsl_vector_memcpy(pulse_aux,&temp4.vector);
					gsl_vector_free(pulseToCalculateEnergy);
					pulseToCalculateEnergy = gsl_vector_alloc(resize_mf);
					gsl_vector_memcpy(pulseToCalculateEnergy,pulse_aux);
					gsl_vector_free(pulse_aux); pulse_aux = 0;
				}
				
				PRCLOFWM = gsl_matrix_alloc(2,resize_mf);
				if (numiteration == 0)
				{
					if (find_prclofwm((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &PRCLOFWM, &Ealpha, &Ebeta))
					{
						message = "Cannot run routine find_prclofwm";
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}
				else
				{
					if (find_prclofwm(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &PRCLOFWM, &Ealpha, &Ebeta))
					{
						message = "Cannot run routine find_prclofwm";
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}
			}
			else if ((strcmp((*reconstruct_init)->EnergyMethod,"WEIGHT") == 0) || (strcmp((*reconstruct_init)->EnergyMethod,"WEIGHTN") == 0))
			{
				if (numiteration == 0)
				{
					// Get the indexes of the two energies which straddle the pulse
					if (find_Esboundary((*pulsesInRecord)->pulses_detected[i].maxDER,(*reconstruct_init)->library_collection->maxDERs,(*reconstruct_init),&indexEalpha,&indexEbeta,&Ealpha, &Ebeta))
					{
						message = "Cannot run routine find_Esboundary for filter interpolation";
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}
				else
				{
					// Get the indexes of the two energies which straddle the pulse
					if (find_Esboundary(energy,(*reconstruct_init)->library_collection->energies,(*reconstruct_init),&indexEalpha,&indexEbeta,&Ealpha, &Ebeta))
					{
						message = "Cannot run routine find_Esboundary for filter interpolation";
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}
				
				if ((indexEalpha == indexEbeta) && ((*reconstruct_init)->library_collection->ntemplates-1 == indexEalpha))
				{
					message = "Not enough info in the library. At least is necessary to add a new last row with energy higher than the energy of the pulses in the input FITS file";
					EP_EXIT_ERROR(message,EPFAIL);
				}
				
				if ((strcmp((*reconstruct_init)->EnergyMethod,"WEIGHTN") == 0) && ((*reconstruct_init)->OFLib == 1))
				{
					// Choose the base-2 system value closest (lower than or equal) to the pulse length
					if (strcmp((*reconstruct_init)->OFStrategy,"BASE2") != 0)
					{
						resize_mf = pow(2,floor(log2(resize_mf)));
						gsl_vector *pulse_aux = gsl_vector_alloc(resize_mf);
                                                //gsl_vector_view(temp5);
						temp = gsl_vector_subvector(pulseToCalculateEnergy,0,resize_mf);
                                                //temp5 = gsl_vector_subvector(pulseToCalculateEnergy,0,resize_mf);
						gsl_vector_memcpy(pulse_aux,&temp.vector);
                                                //gsl_vector_memcpy(pulse_aux,&temp5.vector);
						gsl_vector_free(pulseToCalculateEnergy);
						pulseToCalculateEnergy = gsl_vector_alloc(resize_mf);
						gsl_vector_memcpy(pulseToCalculateEnergy,pulse_aux);
						gsl_vector_free(pulse_aux); pulse_aux = 0;
					}
					
					// Find the appropriate values of the PRECALWN HDU ('PRCLx' columns)
					PRCLWN = gsl_matrix_alloc(2,resize_mf);
					Pab = gsl_vector_alloc(resize_mf);
					if (numiteration == 0)
					{
						if (find_prclwn((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &PRCLWN, &Pab,&Ealpha, &Ebeta))
						{
							message = "Cannot run routine find_prclwn";
							EP_EXIT_ERROR(message,EPFAIL);
						}
					}
					else
					{
						if (find_prclwn(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &PRCLWN, &Pab,&Ealpha, &Ebeta))
						{
							message = "Cannot run routine find_prclwn";
							EP_EXIT_ERROR(message,EPFAIL);
						}
					}
				}
			}
			
			gsl_matrix_set(Estraddle,0,numiteration,Ealpha);
			gsl_matrix_set(Estraddle,1,numiteration,Ebeta);

			if ((strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0) && (strcmp((*reconstruct_init)->OFNoise,"NSD") == 0))
			{
				if ((*reconstruct_init)->OFLib == 0)
				{
					// Calculate the optimal filter
					if (calculus_optimalFilter (TorF, (*reconstruct_init)->intermediate, (*reconstruct_init)->opmode, filtergsl, filtergsl->size, 1.0/record->delta_t, runF0orB0val, (*reconstruct_init)->noise_spectrum->noisefreqs, (*reconstruct_init)->noise_spectrum->noisespec, &optimalfilter, &optimalfilter_f, &optimalfilter_FFT, &optimalfilter_FFT_complex))
					{
						message = "Cannot run routine calculus_optimalFilter";
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}
				else if ((*reconstruct_init)->OFLib == 1)
				{
                         		// It is not necessary to check the allocation because 'filtergsl' size has been checked previously
					optimalfilter = gsl_vector_alloc(filtergsl->size);
					gsl_vector_set_all(optimalfilter,0);
					if (TorF == 0)     gsl_vector_memcpy(optimalfilter,filtergsl);
					else if (TorF == 1)
					{
						// It is not necessary to check the allocation because 'filtergsl' size has been checked previously
						optimalfilter_FFT_complex = gsl_vector_complex_alloc(filtergsl->size/2);
                                         	for (int k=0;k<filtergsl->size/2;k++)
						{
                                                    gsl_vector_complex_set(optimalfilter_FFT_complex,k,gsl_complex_rect(gsl_vector_get(filtergsl,k),gsl_vector_get(filtergsl,k+filtergsl->size/2)));
						}
						gsl_vector_complex_set(optimalfilter_FFT_complex,0,gsl_complex_rect(0.0,0.0));
						for (int k=0;k<filtergsl->size/2;k++)
						{
                                                    gsl_vector_complex_set(optimalfilter_FFT_complex,k,gsl_complex_conjugate(gsl_vector_complex_get(optimalfilter_FFT_complex,k)));
						}
					}
				}
			}
			
			// Calculate the energy of each pulse
			if (calculateEnergy(pulseToCalculateEnergy,pulseGrade,optimalfilter,optimalfilter_FFT_complex,runEMethod,indexEalpha,indexEbeta,(*reconstruct_init),TorF,1/record->delta_t,Pab,PRCLWN,PRCLOFWM,&energy,&tstartNewDev,&lagsShift,0))
			{
				message = "Cannot run calculateEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
				EP_EXIT_ERROR(message,EPFAIL);
			}
			gsl_vector_free(pulseToCalculateEnergy); pulseToCalculateEnergy = 0;
                        
                        //cout<<"Tstart0: "<<endl;
                        //std::cout << std::setprecision(15) << (*pulsesInRecord)->pulses_detected[i].Tstart*156250 << '\n';
                        //cout<<"Tstart0: "<<(*pulsesInRecord)->pulses_detected[i].Tstart<<endl;
			// If using lags, it is necessary to modify the tstart of the pulse and the length of the filter used
			if ((strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0) && tstartNewDev != -999.0)
			{
				(*pulsesInRecord)->pulses_detected[i].Tstart = (*pulsesInRecord)->pulses_detected[i].Tstart + tstartNewDev*record->delta_t; // In seconds
				//if (tstartNewDev != 0) (*pulsesInRecord)->pulses_detected[i].grade1 = (*pulsesInRecord)->pulses_detected[i].grade1 - floor(fabs(tstartNewDev)+0.5);
			}
			//cout<<"pulses_detected[i].Tstart1: "<<pulses_detected[i].Tstart<<endl;
			//cout<<"Tstart1: "<<(*pulsesInRecord)->pulses_detected[i].Tstart<<endl;
			//cout<<"Tstart1: "<<endl;
                        //std::cout << std::setprecision(15) << (*pulsesInRecord)->pulses_detected[i].Tstart*156250 << '\n';
			
			gsl_matrix_set(resultsE,0,numiteration,energy);
			gsl_matrix_set(resultsE,1,numiteration,min(fabs(energy-Ealpha),fabs(energy-Ebeta)));

			gsl_vector_view(temp);
                        //gsl_vector_view(temp6);
			gsl_vector *subresultsE;
			if ((Ealpha != Ebeta) && ((energy < Ealpha) || (energy > Ebeta)))
			{
				if (numiteration != 0)
				{
					for (int j=0; j<numiteration+1; j++)
					{
						if ((Ealpha == gsl_matrix_get(Estraddle,0,j)) && (Ebeta == gsl_matrix_get(Estraddle,1,j)))
						{
							iterate = false;
							
							subresultsE = gsl_vector_alloc(resultsE->size2);
							gsl_matrix_get_row(subresultsE,resultsE,1);
							temp = gsl_vector_subvector(subresultsE,0,numiteration+1);
                                                        //temp6 = gsl_vector_subvector(subresultsE,0,numiteration+1);
							gsl_vector *SUBsubresultsE = gsl_vector_alloc(numiteration+1);
							gsl_vector_memcpy(SUBsubresultsE,&temp.vector);
                                                        //gsl_vector_memcpy(SUBsubresultsE,&temp6.vector);
							energy = gsl_matrix_get(resultsE,0,gsl_vector_min_index(SUBsubresultsE));
							gsl_vector_free(subresultsE); subresultsE = 0;
							gsl_vector_free(SUBsubresultsE); SUBsubresultsE = 0;
							
							break;
						}
					}
				}
			}
			else 
			{	
				iterate = false;
			
				if (numiteration != 0)
				{
					subresultsE = gsl_vector_alloc(resultsE->size2);
					gsl_matrix_get_row(subresultsE,resultsE,1);
					temp = gsl_vector_subvector(subresultsE,0,numiteration+1);
                                        //temp6 = gsl_vector_subvector(subresultsE,0,numiteration+1);
					gsl_vector *SUBsubresultsE = gsl_vector_alloc(numiteration+1);
					gsl_vector_memcpy(SUBsubresultsE,&temp.vector);
                                        //gsl_vector_memcpy(SUBsubresultsE,&temp6.vector);
					energy = gsl_matrix_get(resultsE,0,gsl_vector_min_index(SUBsubresultsE));
					gsl_vector_free(subresultsE); subresultsE = 0;
					gsl_vector_free(SUBsubresultsE); SUBsubresultsE = 0;
				}
			  
			}
		} while (iterate);
                //gsl_vector_free(pulseToCalculateEnergy);

		gsl_matrix_free(Estraddle); Estraddle = 0;
		gsl_matrix_free(resultsE); resultsE = 0;
		

		// Subtract the pulse model from the record
		if (find_model_energies(energy, (*reconstruct_init), &model))
		{
			message = "Cannot run find_model_energies routine for pulse i=" + boost::lexical_cast<std::string>(i);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		tstartJITTER = ((*pulsesInRecord)->pulses_detected[i].Tstart-record->time)/record->delta_t;
		shift = tstartJITTER - tstartSamplesRecord;
		//cout<<"shiftFIN: "<<shift<<endl;
		// In order to subtract the pulse model, it has to be located in the tstart with jitter and know its values in the digitized samples
		gsl_vector *modelToSubtract = gsl_vector_alloc(model->size);
		for (int j=0;j<model->size;j++)
                {
                    if (shift < 0)
                    {
                        if (j != model->size-1)
                            gsl_vector_set(modelToSubtract,j,(gsl_vector_get(model,j+1)-gsl_vector_get(model,j))*(-shift)+gsl_vector_get(model,j));
                        
                        else 
                            gsl_vector_set(model,j,gsl_vector_get(model,j)); //?????????????????????
                    }
                    else if (shift > 0)
                    {
                        if (j == 0)
                        {
                            gsl_vector_set(modelToSubtract,j,(gsl_vector_get(model,j)-0)*(1-shift)+0);
                        }
                        else 
                        {
                            gsl_vector_set(modelToSubtract,j,(gsl_vector_get(model,j)-gsl_vector_get(model,j-1))*(1-shift)+gsl_vector_get(model,j-1));
                        }
                    }
                    else
                    {
                        gsl_vector_memcpy(modelToSubtract,model);
                    }
                }
                gsl_vector_memcpy(model,modelToSubtract);
                gsl_vector_free(modelToSubtract);

		//for (int j=tstartSamplesRecord;j<min(tstartSamplesRecord+(model->size),(double) record->trigger_size);j++)
                for (int j=tstartSamplesRecord;j<min((double) tstartSamplesRecord+(model->size),(double) record->trigger_size);j++)
		{
			gsl_vector_set(recordAux,j,gsl_vector_get(recordAux,j)-gsl_vector_get(model,j-tstartSamplesRecord));
		}

		// Write info of the pulse in the output intemediate file if 'intermediate'=1
		if ((*reconstruct_init)->intermediate == 1)
		{
			if (writeFilterHDU(reconstruct_init, i,energy, optimalfilter, optimalfilter_f, optimalfilter_FFT, &dtcObject))
			{
				message = "Cannot run writeFilterHDU routine for pulse i=" + boost::lexical_cast<std::string>(i);
				EP_EXIT_ERROR(message,EPFAIL);
			}
		}
		
		(*pulsesInRecord)->pulses_detected[i].energy = energy/1e3;	// In SIXTE, SIGNAL is in keV
		(*pulsesInRecord)->pulses_detected[i].grading = pulseGrade;	
                ////(*pulsesInRecord)->pulses_detected[i].phi = tstartNewDev - lagsShift;
                //(*pulsesInRecord)->pulses_detected[i].phi = tstartNewDev;
                //(*pulsesInRecord)->pulses_detected[i].lagsShift = lagsShift;
                double intpart;
                (*pulsesInRecord)->pulses_detected[i].phi = modf(tstartNewDev,&intpart);
                (*pulsesInRecord)->pulses_detected[i].lagsShift = lagsShift+intpart;
                
		// Free allocated GSL vectors
		gsl_vector_free(optimalfilter); optimalfilter = 0;
		gsl_vector_free(optimalfilter_f); optimalfilter_f = 0;
		gsl_vector_free(optimalfilter_FFT); optimalfilter_FFT = 0;
		if ((*pulsesInRecord)->pulses_detected[i].quality < 10)
		{
                  gsl_vector_free(pulse); pulse = 0;
                  gsl_vector_free(filtergsl); filtergsl = 0;
                  gsl_vector_complex_free(optimalfilter_FFT_complex); optimalfilter_FFT_complex = 0;
                  gsl_vector_free(Pab); Pab = 0;
                  gsl_matrix_free(PRCLWN); PRCLWN = 0;
                  gsl_matrix_free(PRCLOFWM); PRCLOFWM = 0;
		}
	} // End for
	
	gsl_vector_free(recordAux); recordAux = 0;
	gsl_vector_free(model); model = 0;
        
        gsl_vector_free(pulse_lowres); pulse_lowres = 0;
        gsl_vector_free(filtergsl_lowres); filtergsl_lowres = 0;
        if (PRCLWN_lowres != NULL) gsl_matrix_free(PRCLWN_lowres); PRCLWN_lowres = 0;
        if (PRCLOFWM_lowres != NULL) gsl_matrix_free(PRCLOFWM_lowres); PRCLOFWM_lowres = 0;
        gsl_vector_free(optimalfilter_lowres); optimalfilter_lowres = 0;
        if (optimalfilter_FFT_complex_lowres != NULL) gsl_vector_complex_free(optimalfilter_FFT_complex_lowres); optimalfilter_FFT_complex_lowres = 0;

	return;
}

/***** SECTION B1 ************************************************************
* calculus_optimalFilter: This function calculates the optimal filter for a pulse whose matched filter is provided as input parameter.
*                         An optimal filter is just a matched filter that has been adjusted based on the noise spectrum of the system.
*
* It is assumed that all pulses are scaled versions of a template. In the frequency domain (as noise can be frequency dependent), the raw data
* can be expressed as P(f)=E·S(f)+N(f), where S(f) is the normalized model pulse shape in the frequency domain, N(f) is the power spectrum of the noise and
* E is the scalar amplitude for the photon energy.
*
* The second assumption is that the noise is stationary, i.e. it does not vary with time. The amplitude of each pulse can then be estimated by 
* minimizing (weighted least-squares sense) the difference between the noisy data and the model pulse shape, being the X^2 condition 
* to be minimized:
* 
*            (P(f)-E·S(f))^2
*  X^2 = SUM ---------------- 
*                N(f)^2
*
* In the time domain, the amplitude is the best weighted (optimally filtered) sum of the values in the pulse
* 
* E = k·SUM p(t)*op(t)
* 
* where of(t) is the time domain expression of optimal filter in frequency domain
*                          
*                    S*(f)               1                               |S(f)|^2
* OptimalFilter(f)= -------  		--- = NormalizationFactor = SUM ---------
*                    N(f)^2              k                               |N(f)|^2
*
* 
* - FFT calculus of the matched filter (filter template)
* 	- Declare variables
* 	- Complex FFT values for positive and negative frequencies
* 	- FFT calculus
* 	- Generation of the frequencies (positive and negative)
* 	- Magnitude and argument for positive and negative frequencies
*       - Free allocated GSL vectors
* - N(f)
* - To divide MatchedFilter(f)/N^2(f) => MatchedFilter(f) and N(f) must have the same number of points
* 	- 'if (mf_size' < freqgsl->size)' 
* 		- 'if ((freqgsl->size)%mf_size == 0)' => Decimate noise samples
* 		- 'else' => It is necessary to work only with the positive frequencies in order to not handle the f=0
* 		            N(f) interpolation ('interpolatePOS')
* 	- 'else if (mf_size > freqgsl->size)' => Error: Noise spectrum must have more samples than pulse spectrum
*	- 'else if (mf_size == freqgsl->size)' => It is not necessary to do anything
* - OptimalFilter = MatchedFilter'(f)/N^2(f)
* - Calculus of the normalization factor
* - Inverse FFT (to get the expression of the optimal filter in time domain)
*	- Complex OptimalFilter(f) => Taking into account magnitude (MatchedFilter(f)/N^2(f)) and phase (given by MatchedFilter(f))
* - Free allocated GSL vectors
*
* Parameters:
* - TorF: 'FilterDomain' = T => 'TorF' = 0
*         'FilterDomain' = F => 'TorF' = 1
* - intermediate: 'intermediate' = 0 => Not write an intermediate file
*                 'intermediate' = 1 => Write an intermediate file
* - opmode: 'opmode' = 0 => CALIBRATION
*           'opmode' = 1 => PRODUCTION
* - matchedfiltergsl: Matched filter associated to the pulse (in general, from the interpolation between two matched filters of the library)
* - mf_size: Matched filter size
* - samprate: Sampling rate
* - runF0orB0val: 'FilterMethod' = F0 => 'runF0orB0val' = 1
*                 'FilterMethod' = B0 => 'runF0orB0val' = 0
* - freqgsl: Frequency axis of the current noise spectral density (input)
* - csdgsl: Current noise spectral density (input)
* - optimal_filtergsl: Optimal filter in time domain (output)
* - of_f: Frequency axis of the optimal filter spectrum (output)
* - of_FFT: Optimal filter spectrum (absolute values) (output)
* - of_FFT_complex: Optimal filter spectrum (complex values) (output)
****************************************/
int calculus_optimalFilter(int TorF, int intermediate, int opmode, gsl_vector *matchedfiltergsl, long mf_size, double samprate, int runF0orB0val, gsl_vector *freqgsl, 
	gsl_vector *csdgsl, gsl_vector **optimal_filtergsl,gsl_vector **of_f, gsl_vector **of_FFT, gsl_vector_complex **of_FFT_complex)
{  
	// FFT calculus of the matched filter (filter template)
	// Declare variables
	double SelectedTimeDuration = mf_size/samprate;
	string message = "";
	char valERROR[256];
	
	// Complex FFT values for positive and negative frequencies
	// It is not necessary to check the allocation because 'mf_size' must already be > 0
	gsl_vector_complex *mfFFTcomp = gsl_vector_complex_alloc(mf_size);
	gsl_vector_complex *mfFFTcomp_conj = gsl_vector_complex_alloc(mf_size);
	gsl_vector *mf_arg = gsl_vector_alloc(mf_size);				// Argument for positive and negative frequencies
	gsl_vector *mf_f = gsl_vector_alloc(mf_size);				// Sorted frequencies according [-fmax,...,0,...,fmax]
	gsl_vector *mf_FFT = gsl_vector_alloc(mf_size);				// Sorted magnitude according [-fmax,...,0,...,fmax]

	// FFT calculus
	if (FFT(matchedfiltergsl,mfFFTcomp,SelectedTimeDuration))
	{
		message = "Cannot run FFT routine to calculate filter template";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (runF0orB0val == 0)		gsl_vector_complex_set(mfFFTcomp,0,gsl_complex_rect(0.0,0.0));

	// Generation of the frequencies (positive and negative)
	// Here is a table which shows the layout of the array data, and the correspondence between the time-domain data z,
	// and the frequency-domain data x

	// index    z               x = FFT(z)

	// 0        z(t = 0)        x(f = 0)
	// 1        z(t = 1)        x(f = 1/(n Delta))
	// 2        z(t = 2)        x(f = 2/(n Delta))
	// .        ........        ..................
	// n/2      z(t = n/2)      x(f = +1/(2 Delta),
	//                                -1/(2 Delta))
	// .        ........        ..................
	// n-3      z(t = n-3)      x(f = -3/(n Delta))
	// n-2      z(t = n-2)      x(f = -2/(n Delta))
	// n-1      z(t = n-1)      x(f = -1/(n Delta))

	if (mf_size == 1)
	{
		gsl_vector_set(mf_f,0,0);
	}
	else if (mf_size%2 == 0)	//Even
	{
		for (int i=0; i<=(mf_size)/2; i++)
		{
			gsl_vector_set(mf_f,i,i/SelectedTimeDuration);
		}
		for (int i=1; i<(mf_size/2); i++)
		{
			gsl_vector_set(mf_f,i+mf_size/2,(-1.0)*(i+mf_size/2-i*2)/SelectedTimeDuration);
		}
	}
	else	//Odd
	{
		for (int i=0; i<=(mf_size)/2; i++)
		{
			gsl_vector_set(mf_f,i,i/SelectedTimeDuration);
		}
		gsl_vector_set(mf_f,mf_size/2+1,(-1.0)*(mf_size/2)/SelectedTimeDuration);
		for (int i=2; i<=(mf_size/2); i++)
		{
			gsl_vector_set(mf_f,i+mf_size/2,(-1.0)*(1+mf_size/2-i)/SelectedTimeDuration);
		}
	}

	// Magnitude and argument for positive and negative frequencies
	gsl_vector_complex_absIFCA(mf_FFT,mfFFTcomp);	// Magnitude
	for (int i=0;i<mf_size;i++)
	{
		gsl_vector_complex_set(mfFFTcomp_conj,i,gsl_complex_conjugate(gsl_vector_complex_get(mfFFTcomp,i)));
	}
	gsl_vector_complex_argIFCA(mf_arg,mfFFTcomp);	// Argument

	// Free allocated GSL vectors
	gsl_vector_complex_free(mfFFTcomp); mfFFTcomp = 0;

	// N(f)
	gsl_vector *n_f;
	gsl_vector *n_FFT;

	// To divide MatchedFilter(f)/N^2(f) => MatchedFilter(f) and N(f) must have the same number of points
	gsl_vector_view temp;
	if (mf_size < freqgsl->size)	
	{
		if ((freqgsl->size)%mf_size == 0)	// Decimate noise samples
		{
			int timesNoverMF = freqgsl->size/mf_size;
			// It is not necessary to check the allocation because 'mf_size' must already be > 0
			n_f = gsl_vector_alloc(mf_size);
			n_FFT = gsl_vector_alloc(mf_size);
			for (int i=0;i<n_f->size;i++)
			{
				gsl_vector_set(n_f,i,gsl_vector_get(freqgsl,i*timesNoverMF));
				gsl_vector_set(n_FFT,i,gsl_vector_get(csdgsl,i*timesNoverMF));
			}
		}
		else
		{
			// It is necessary to work only with the positive frequencies in order to not handle the f=0
			int noisePOS_size = floor(freqgsl->size/2);
			// It is not necessary to check the allocation because 'freqgsl' size must already be > 0
			gsl_vector *freqgsl_POS = gsl_vector_alloc(noisePOS_size);
			temp = gsl_vector_subvector(freqgsl,1,noisePOS_size);
			gsl_vector_memcpy(freqgsl_POS,&temp.vector);
			gsl_vector *csdgsl_POS = gsl_vector_alloc(noisePOS_size);
			temp = gsl_vector_subvector(csdgsl,1,noisePOS_size);
			gsl_vector_memcpy(csdgsl_POS,&temp.vector);

			// N(f) interpolation
			int mfPOS_size = floor(mf_size/2);
			// It is not necessary to check the allocation because 'mf_size' must already be > 0
			gsl_vector *n_f_interp = gsl_vector_alloc(mfPOS_size);
			gsl_vector *n_FFT_interp = gsl_vector_alloc(mfPOS_size);
		   	if (interpolatePOS (freqgsl_POS, csdgsl_POS, mfPOS_size, gsl_vector_get(mf_f,1)-gsl_vector_get(mf_f,0), &n_f_interp, &n_FFT_interp))
		   	{
				message = "Cannot run routine interpolate for matched filter interpolation";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}

			// It is not necessary to check the allocation because 'mf_size' must already be > 0
		   	n_f = gsl_vector_alloc(mf_size);
		   	n_FFT = gsl_vector_alloc(mf_size);
		   	gsl_vector_set(n_f,0,0);
		   	gsl_vector_set(n_FFT,0,gsl_vector_get(csdgsl,0));
			for (int i=0;i<n_f_interp->size;i++)
		   	{
		   		gsl_vector_set(n_f,i+1,gsl_vector_get(n_f_interp,i));
		   		gsl_vector_set(n_FFT,i+1,gsl_vector_get(n_FFT_interp,i));
		   	}
			if (mf_size%2 == 0)
			{
				if (n_f->size-n_f_interp->size+1 > n_f->size-1)
				{
					sprintf(valERROR,"%d",__LINE__+7);
					string str(valERROR);
					message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				for (int i=0;i<n_f_interp->size-1;i++)
				{
					gsl_vector_set(n_f,n_f->size-i-1,-1.0*gsl_vector_get(n_f_interp,i));
					gsl_vector_set(n_FFT,n_f->size-i-1,gsl_vector_get(n_FFT_interp,i));
				}
			}
			else
			{
				if (n_f->size-n_f_interp->size-1-1 > n_f->size-1)
				{
					sprintf(valERROR,"%d",__LINE__+7);
					string str(valERROR);
					message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				for (int i=0;i<n_f_interp->size;i++)
				{
					gsl_vector_set(n_f,n_f->size-i-1,-1.0*gsl_vector_get(n_f_interp,i));
					gsl_vector_set(n_FFT,n_f->size-i-1,gsl_vector_get(n_FFT_interp,i));
				}
			}

			gsl_vector_free(n_f_interp); n_f_interp = 0;
		   	gsl_vector_free(n_FFT_interp); n_FFT_interp = 0;
		   	gsl_vector_free(freqgsl_POS); freqgsl_POS = 0;
		   	gsl_vector_free(csdgsl_POS); csdgsl_POS = 0;
		}
	}
	else if (mf_size > freqgsl->size)	// Error
	{
		message = "Noise must have more samples than pulse";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	else if (mf_size == freqgsl->size)	// It is not necessary to do anything
	{
		// It is not necessary to check the allocation because 'freqgsl' size must already be > 0
		n_f = gsl_vector_alloc(freqgsl->size);
		n_FFT = gsl_vector_alloc(freqgsl->size);
		gsl_vector_memcpy(n_f,freqgsl);
		gsl_vector_memcpy(n_FFT,csdgsl);
	}

	// It is not necessary to check the allocation because 'freqgsl' size must already be > 0
	// OptimalFilter = MatchedFilter'(f)/N^2(f)
	// It is not necessary to check the allocation because 'mf_f' size ('mf_size') must already be > 0
	*of_f = gsl_vector_alloc(mf_f->size);
	*of_FFT = gsl_vector_alloc(mf_f->size);
	gsl_vector_memcpy(*of_f,mf_f);

	// It is not necessary to check the allocation because 'mf_f' size ('mf_size') must already be > 0
	gsl_vector *mf_FFT_2 = gsl_vector_alloc(mf_f->size);
	gsl_vector *n_FFT_2 = gsl_vector_alloc(mf_f->size);
	gsl_vector_memcpy(mf_FFT_2,mf_FFT);
	gsl_vector_mul(mf_FFT_2,mf_FFT_2);
	if (gsl_vector_memcpy(n_FFT_2,n_FFT) != 0)
	{
		sprintf(valERROR,"%d",__LINE__-2);
		string str(valERROR);
		message = "Copying vectors of different length in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	gsl_vector_mul(n_FFT_2,n_FFT_2);

	// It is not necessary to check the allocation because 'mf_size' must already be > 0
	*of_FFT_complex = gsl_vector_complex_alloc(mf_size);
	gsl_vector_complex_set(*of_FFT_complex,0,gsl_complex_rect(0.0,0.0));
	for (int i=1;i<mf_size;i++)
	{
		gsl_vector_complex_set(*of_FFT_complex,i,gsl_complex_div_real(gsl_vector_complex_get(mfFFTcomp_conj,i),gsl_vector_get(n_FFT_2,i)));
	}

	// Calculus of the normalization factor
	double normalizationFactor = 0;
	for (int i=1; i<mf_f->size; i++)
	{
		normalizationFactor = normalizationFactor + gsl_vector_get(mf_FFT_2,i)/gsl_vector_get(n_FFT_2,i);
	}
	
	// Apply the normalization factor
	gsl_vector_complex_scale(*of_FFT_complex,gsl_complex_rect(1.0/normalizationFactor,0.0));
	
	gsl_vector_complex_absIFCA(*of_FFT,*of_FFT_complex);

	gsl_vector_free(mf_FFT_2); mf_FFT_2 = 0;
	gsl_vector_free(n_FFT_2); n_FFT_2 = 0;
	
	if ((TorF == 0) || (intermediate == 1) || (opmode == 0))
	{
		// Inverse FFT (to get the expression of the optimal filter in time domain)
		// Complex OptimalFilter(f) => Taking into account magnitude (MatchedFilter(f)/N^2(f)) and phase (given by MatchedFilter(f))
		// It is not necessary to check the allocation because 'mf_size' must already be > 0
		gsl_vector_complex *of_FFTcomp = gsl_vector_complex_alloc(mf_size);
		*optimal_filtergsl = gsl_vector_alloc(mf_size);
		for (int i=0;i<mf_size;i++)
		{
			gsl_vector_complex_set(of_FFTcomp,i,gsl_complex_polar(gsl_complex_abs(gsl_vector_complex_get(*of_FFT_complex,i)),gsl_vector_get(mf_arg,i)));
		}
		if (FFTinverse(of_FFTcomp,*optimal_filtergsl,SelectedTimeDuration))
		{
			message = "Cannot run routine FFTinverse to get optimal filter in time domain";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}

		gsl_vector_complex_free(of_FFTcomp); of_FFTcomp = 0;
	}
	else
	{
		// It is not necessary to check the allocation because 'mf_size' must already be > 0
		*optimal_filtergsl = gsl_vector_alloc(mf_size);
		gsl_vector_set_zero(*optimal_filtergsl);
	}
	
	// Free allocated GSL vectors
	gsl_vector_complex_free(mfFFTcomp_conj); mfFFTcomp_conj = 0;
	gsl_vector_free(mf_f); mf_f = 0;
	gsl_vector_free(mf_FFT); mf_FFT = 0;
	gsl_vector_free(mf_arg); mf_arg = 0;
	gsl_vector_free(n_f); n_f = 0;
	gsl_vector_free(n_FFT); n_FFT = 0;

	return(EPOK);
}
/*xxxx end of SECTION B1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B2 ************************************************************
* interpolatePOS: This function interpolates an input vector ('x_in', 'y_in'), creating an output vector ('x_out', 'y_out') with the size and
*                 frequency step given.
*                 POS is due to the fact that the input spectrum only has positive frequencies (in order to not handle the f=0).
*
* - Declare and initialize variables
* - Method applied to interpolate
* - Generate the interpolated output vector
* - Free memory
*
* Parameters:
* - x_in: GSL input vector with the abscissas of the vector which is going to be interpolated 
* - y_in: GSL input vector with the ordinates of the vector which is going to be interpolated
* - size: Size of the interpolated output vector
* - step: Frequency step of the interpolated output vector
* - x_out: GSL output vector with the abscissas of the interpolated vector
* - y_out: GSL output vector with the ordinates of the interpolated vector
****************************************/
int interpolatePOS (gsl_vector *x_in, gsl_vector *y_in, long size, double step, gsl_vector **x_out, gsl_vector **y_out)
{
	string message = "";
	char valERROR[256];
	
	// Declare variables
	gsl_vector_set_zero(*x_out);
	gsl_vector_set_zero(*y_out);
	int N = x_in->size;
	double x[N], y[N];
	for (int i=0; i<N; i++)
	{
		x[i] = gsl_vector_get(x_in,i);
		y[i] = gsl_vector_get(y_in,i);
	}
	double xi, yi;

	// Method applied to interpolate
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	const gsl_interp_type *t = gsl_interp_cspline;

	// It is not necessary to check the allocation because 'x_in' size ('N') must already be > 0
	gsl_spline *spline = gsl_spline_alloc (t, N);

	gsl_spline_init (spline, x, y, N);

	// Generate the interpolated output vector
	xi = step;
	if (size-1 > (*x_out)->size-1)
	{
		sprintf(valERROR,"%d",__LINE__+12);
		string str(valERROR);
		message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	for (int i=0; i<size; i++)
	{
		if ((i == size-1) && (xi <= x[N-1]+x[N-1]/1000)) xi = x[N-1];
		if ((xi >= x[0]) && (xi <= x[N-1]))
		{
			yi = gsl_spline_eval (spline, xi, acc);

			gsl_vector_set(*x_out,i,xi);
			gsl_vector_set(*y_out,i,yi);
		}

	    xi = xi+step;
	}

	// Free memory
	gsl_spline_free (spline); spline = 0;
	gsl_interp_accel_free (acc); acc = 0;

	return EPOK;
}
/*xxxx end of SECTION B2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B3 ************************************************************
* find_matchedfilter: This function chooses the proper matched filter of the calibration library ('MF' or 'MFB0')
*                     by comparing the maximum of their derivatives ('maxDER' versus 'maxDERs').
*
* It finds the two embracing 'maxDERs' in the calibration library and interpolates ('interpolate_model')
*   - If 'maxDER' is lower than the lowest 'maxDERs' in the library => The matched filter with the lowest 'maxDERs' (first row) in the library is chosen
*   - If 'maxDER' is higher than the highest 'maxDERs' in the library => The matched filter with the highest 'maxDERs' (last row) in the library is chosen
* 
* Parameters:
* - runF0orB0val: 'FilterMethod' = F0 => 'runF0orB0val' = 1
*                 'FilterMethod' = B0 => 'runF0orB0val' = 0
* - maxDER: Max value of the derivative of the (filtered) pulse whose matched filter is being sought
* - maxDERs: GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values).
*                     In particular, this function uses the info in the library ('matched_filters')
* - matchedfilterFound: GSL vector with the matched filter selected
* - Ealpha: Energy (in eV) which straddle the 'maxDER' in the lower limit
* - Ebeta: Energy (in eV) which straddle the 'maxDER' in the higher limit
****************************************/
int find_matchedfilter(int runF0orB0val, double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **matchedfilterFound, double *Ealpha, double *Ebeta)
{
	string message = "";
	char valERROR[256];
        
        gsl_vector *maxDERs_LIB1row;
        long nummodels;
        int index_E_LIB1row;
        if (reconstruct_init->filtEev == 0)
        {
                nummodels = maxDERs->size;
                maxDERs_LIB1row = gsl_vector_alloc(nummodels);
                gsl_vector_memcpy(maxDERs_LIB1row,maxDERs);
        }
        else
        {
                nummodels = 1;
                maxDERs_LIB1row = gsl_vector_alloc(1);
                for (int i=0;i<maxDERs->size;i++)
                {
                        if (gsl_vector_get(reconstruct_init->library_collection->energies,i) == reconstruct_init->filtEev)
                        {
                                gsl_vector_set(maxDERs_LIB1row,0,gsl_vector_get(maxDERs,i));
                                index_E_LIB1row = i;
                                break;
                        }
                }
        }
       		
	gsl_vector *matchedfilterFound_aux;
	if ((matchedfilterFound_aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration)) == 0)
	{
		sprintf(valERROR,"%d",__LINE__-2);
		string str(valERROR);
	        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	if (maxDER < gsl_vector_get(maxDERs_LIB1row,0))
	{
                if (reconstruct_init->filtEev == 0)
                {
                        if (runF0orB0val == 0)	gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection->matched_filters[0].mfilter);
                        else    		gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection->matched_filters_B0[0].mfilter);
                        *Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
                }
                else
                {
                        if (runF0orB0val == 0)	gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection->matched_filters[index_E_LIB1row].mfilter);
                        else    		gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection->matched_filters_B0[index_E_LIB1row].mfilter);
                        *Ebeta = reconstruct_init->filtEev;
                }
		
		*Ealpha = 0.0;
	}
	else if (maxDER > gsl_vector_get(maxDERs_LIB1row,nummodels-1))
	{
                if (reconstruct_init->filtEev == 0)
                {
                        if (runF0orB0val == 0)		gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection
->matched_filters[nummodels-1].mfilter);
                        *Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1);
                }
                else
                {
                        if (runF0orB0val == 0)		gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection
->matched_filters[index_E_LIB1row].mfilter);
                        *Ealpha = reconstruct_init->filtEev;
                }
		
		*Ebeta = 0.0;
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
                        if (maxDER == gsl_vector_get(maxDERs_LIB1row,i))
			{
				if (runF0orB0val == 0)	gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection->matched_filters[i].mfilter);
				else			gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection->matched_filters_B0[i].mfilter);

				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i);

				break;
			}
			else if ((maxDER > gsl_vector_get(maxDERs_LIB1row,i)) && (maxDER < gsl_vector_get(maxDERs_LIB1row,i+1)))
			{
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i+1);

				// Interpolate between the two corresponding rows in "models"
				gsl_vector *matchedfilterAux;
				if ((matchedfilterAux = gsl_vector_alloc(reconstruct_init->library_collection->matched_filters[0].mfilter_duration)) == 0)
				{
					sprintf(valERROR,"%d",__LINE__-2);
					string str(valERROR);
					message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				gsl_vector_set_zero(matchedfilterAux);
				if (runF0orB0val == 0)
				{
					if (interpolate_model(&matchedfilterAux,maxDER,reconstruct_init->library_collection->matched_filters[i].mfilter,gsl_vector_get(maxDERs_LIB1row,i),
						reconstruct_init->library_collection->matched_filters[i+1].mfilter,gsl_vector_get(maxDERs_LIB1row,i+1)))
					{
						message = "Cannot run interpolate_model routine for model interpolation";
						EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
					}
				}
				else
				{
                                        if (interpolate_model(&matchedfilterAux,maxDER,reconstruct_init->library_collection->matched_filters_B0[i].mfilter,gsl_vector_get(maxDERs_LIB1row,i),
						reconstruct_init->library_collection->matched_filters_B0[i+1].mfilter,gsl_vector_get(maxDERs_LIB1row,i+1)))
					{
						message = "Cannot run interpolate_model routine for model interpolation";
						EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
					}
				}
				gsl_vector_memcpy(matchedfilterFound_aux,matchedfilterAux);

				gsl_vector_free(matchedfilterAux); matchedfilterAux = 0;
				
				break;
			}
		}
	}
	
	gsl_vector_view temp;
	temp = gsl_vector_subvector(matchedfilterFound_aux,0,(*matchedfilterFound)->size);
	gsl_vector_memcpy(*matchedfilterFound,&temp.vector);

	gsl_vector_free(matchedfilterFound_aux); matchedfilterFound_aux = 0;
        
        gsl_vector_free(maxDERs_LIB1row); maxDERs_LIB1row = 0;

	return(EPOK);
}
/*xxxx end of SECTION B3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B4 ************************************************************
* find_matchedfilterDAB: This function selects the proper matched filter (normalized template) from the calibration library (from column 'DAB')
*                        by comparing the maximum value of the pulse derivative ('maxDER') to the list of maximums in the library ('maxDERs') for the DAB interpolation method.
*                        It also selects the proper row from the column 'PAB'.
*
* It finds the two embracing 'maxDERs' in the calibration library
*   - If 'maxDER' is lower than the lowest 'maxDERs' in the library => The data with the lowest 'maxDERs' (first row) in the library are chosen
*   - If 'maxDER' is higher than the highest 'maxDERs' in the library => The data of the penultimate row in the library are chosen
*
* Parameters:
* - runF0orB0val: 'FilterMethod' = F0 => 'runF0orB0val' = 1
*                 'FilterMethod' = B0 => 'runF0orB0val' = 0
* - maxDER: Max value of the derivative of the (filtered) pulse whose matched filter is being sought
* - maxDERs: GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values).
*                     In particular, this function uses the info in the library ('matched_filters')
* - matchedfilterFound: GSL vector with the matched filter selected
* - PabFound: PAB column from the library
* - Ealpha: Energy (in eV) which straddle the 'maxDER' in the lower limit
* - Ebeta: Energy (in eV) which straddle the 'maxDER' in the higher limit
****************************************/
int find_matchedfilterDAB(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **matchedfilterFound, gsl_vector **PabFound, double *Ealpha, double *Ebeta)
{
	string message = "";
	char valERROR[256];
	
	long nummodels = maxDERs->size;
	
	gsl_vector *matchedfilterFound_aux;
	if ((matchedfilterFound_aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration)) == 0)
	{
		sprintf(valERROR,"%d",__LINE__-2);
		string str(valERROR);
	        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	// It is not necessary to check the allocation because the allocation of 'matchedfilterFound_aux' has been checked previously 
	gsl_vector *PabFound_aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);

	if (maxDER < gsl_vector_get(maxDERs,0))
	{
		gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection->matched_filters[0].mfilter);
	
		gsl_matrix_get_row(PabFound_aux,reconstruct_init->library_collection->PAB,0);
				
		*Ealpha = 0.0;
		*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
	}
	else if (maxDER > gsl_vector_get(maxDERs,nummodels-1))
	{
                if (nummodels == 1)
                {
                        gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection->matched_filters[0].mfilter);
                                        
                        gsl_matrix_get_row(PabFound_aux,reconstruct_init->library_collection->PAB,0);
                        
                        *Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,0);
                        *Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
                }
                else
                {
                        gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection->matched_filters[nummodels-2].mfilter);
                                        
                        gsl_matrix_get_row(PabFound_aux,reconstruct_init->library_collection->PAB,nummodels-2);	//!!!!!! -2 !!!!!!!
                                                                                                                // -1 because of the GSL indexes start in 0
                                                                                                                // -1 because the info AB between two energies A<B is stored in the row of the energy A 
                        
                        *Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-2);
                        *Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1);
                }
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
			if (maxDER == gsl_vector_get(maxDERs,i))
			{
				gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection->matched_filters[i].mfilter);
				
				gsl_matrix_get_row(PabFound_aux,reconstruct_init->library_collection->PAB,i);

				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i);

				break;
			}
			else if ((maxDER > gsl_vector_get(maxDERs,i)) && (maxDER < gsl_vector_get(maxDERs,i+1)))
			{
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i+1);

				gsl_vector_memcpy(matchedfilterFound_aux,reconstruct_init->library_collection->matched_filters[i].mfilter);
				
				gsl_matrix_get_row(PabFound_aux,reconstruct_init->library_collection->PAB,i);
	
				break;
			}
		}
	}
	
	gsl_vector_view temp;
	temp = gsl_vector_subvector(matchedfilterFound_aux,0,(*matchedfilterFound)->size);
	gsl_vector_memcpy(*matchedfilterFound,&temp.vector);
	temp = gsl_vector_subvector(PabFound_aux,0,(*matchedfilterFound)->size);
	gsl_vector_memcpy(*PabFound,&temp.vector);
	
	gsl_vector_free(matchedfilterFound_aux); matchedfilterFound_aux = 0;
	gsl_vector_free(PabFound_aux); PabFound_aux = 0;

	return(EPOK);
}
/*xxxx end of SECTION B4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B5 ************************************************************
* find_optimalfilter: This function selects the proper optimal filter from the calibration library ('Tx', 'ABTx', 'Fx' or 'ABFx')
*                     by comparing the maximum value of the pulse derivative ('maxDER') to the list of maximums in the library ('maxDERs').
*
* It finds the two embracing 'maxDERs' in the calibration library and interpolates ('interpolate_model')
*   - If 'maxDER' is lower than the lowest 'maxDERs' in the  library => The optimal filter with the lowest 'maxDERs' (first row) in the library is chosen
*   - If 'maxDER' is higher than the highest 'maxDERs' in the library => The optimal filter with the highest 'maxDERs' (last row) in the library is chosen
*
* Parameters:
* - maxDER: Max value of the derivative of the filtered pulse whose matched filter is being sought.
* - maxDERs: GSL vector with the maximum values of the derivatives of the templates in the library to be compare with the pulse being analysed
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values).
*                     In particular, this function uses the info in the library ('optimal_filters')
* - optimalfilterFound: GSL vector with the optimal filter selected
* - Ealpha: Energy (in eV) which straddle the 'maxDER' in the lower limit
* - Ebeta: Energy (in eV) which straddle the 'maxDER' in the higher limit
****************************************/
int find_optimalfilter(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **optimalfilterFound, double *Ealpha, double *Ebeta)
{
	string message = "";
	char valERROR[256];
	
        gsl_vector *maxDERs_LIB1row;
        long nummodels;
        int index_E_LIB1row;
        if (reconstruct_init->filtEev == 0)
        {
                nummodels = maxDERs->size;
                maxDERs_LIB1row = gsl_vector_alloc(nummodels);
                gsl_vector_memcpy(maxDERs_LIB1row,maxDERs);
        }
        else
        {
                nummodels = 1;
                maxDERs_LIB1row = gsl_vector_alloc(1);
                for (int i=0;i<maxDERs->size;i++)
                {
                        if (gsl_vector_get(reconstruct_init->library_collection->energies,i) == reconstruct_init->filtEev)
                        {
                                gsl_vector_set(maxDERs_LIB1row,0,gsl_vector_get(maxDERs,i));
                                index_E_LIB1row = i;
                                break;
                        }
                }
        }
        
	int sizeFilter;
	if (strcmp(reconstruct_init->FilterDomain,"F") == 0)		sizeFilter = (*optimalfilterFound)->size/2.0;
	else if (strcmp(reconstruct_init->FilterDomain,"T") == 0)	sizeFilter = (*optimalfilterFound)->size;
	
	gsl_vector *optimalfilterFound_Aux;
        if ((optimalfilterFound_Aux = gsl_vector_alloc(reconstruct_init->library_collection->optimal_filters[0].ofilter_duration)) == 0)
        {
                sprintf(valERROR,"%d",__LINE__-2);
                string str(valERROR);
                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
        }

	if (maxDER < gsl_vector_get(maxDERs_LIB1row,0))
	{  
                if (reconstruct_init->filtEev == 0)
                {
                        gsl_vector_memcpy(optimalfilterFound_Aux,reconstruct_init->library_collection->optimal_filters[0].ofilter);
                        *Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
                }
                else
                {
                        gsl_vector_memcpy(optimalfilterFound_Aux,reconstruct_init->library_collection->optimal_filters[index_E_LIB1row].ofilter);
                        *Ebeta = reconstruct_init->filtEev;
                }

		*Ealpha = 0.0;
	}
	else if (maxDER > gsl_vector_get(maxDERs_LIB1row,nummodels-1))
	{
                if (reconstruct_init->filtEev == 0)
                {
                        gsl_vector_memcpy(optimalfilterFound_Aux,reconstruct_init->library_collection->optimal_filters[nummodels-1].ofilter);
                        *Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1);
                }
                else
                {
                        gsl_vector_memcpy(optimalfilterFound_Aux,reconstruct_init->library_collection->optimal_filters[index_E_LIB1row].ofilter);
                        *Ealpha = reconstruct_init->filtEev;
                }
		
		*Ebeta = 0.0;
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
			if (maxDER == gsl_vector_get(maxDERs,i))
			{
				gsl_vector_memcpy(optimalfilterFound_Aux,reconstruct_init->library_collection->optimal_filters[i].ofilter);

				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i);

				break;
			}
			else if ((maxDER > gsl_vector_get(maxDERs,i)) && (maxDER < gsl_vector_get(maxDERs,i+1)))
			{
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i+1);

				// Interpolate between the two corresponding rows in "models"
				gsl_vector *optimalfilter_Aux2;
				if ((optimalfilter_Aux2 = gsl_vector_alloc(reconstruct_init->library_collection->optimal_filters[0].ofilter_duration)) == 0)
				{
					sprintf(valERROR,"%d",__LINE__-2);
					string str(valERROR);
					message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				gsl_vector_set_zero(optimalfilter_Aux2);
				if (interpolate_model(&optimalfilter_Aux2,maxDER,reconstruct_init->library_collection->optimal_filters[i].ofilter,gsl_vector_get(maxDERs,i),
					reconstruct_init->library_collection->optimal_filters[i+1].ofilter,gsl_vector_get(maxDERs,i+1)))
				{
					message = "Cannot run interpolate_model routine for model interpolation";
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}

				gsl_vector_memcpy(optimalfilterFound_Aux,optimalfilter_Aux2);

				gsl_vector_free(optimalfilter_Aux2); optimalfilter_Aux2 = 0;
							
				break;
			}
		}
	}
	
	gsl_vector *fixedlengths = gsl_vector_alloc(reconstruct_init->library_collection->nfixedfilters);
	if ((reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != reconstruct_init->library_collection->pulse_templates[0].template_duration) &&
	   (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != -999)
	   && (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != 999))
	   
	{
		for (int i=0;i<reconstruct_init->library_collection->nfixedfilters-1;i++)
		{
			gsl_vector_set(fixedlengths,reconstruct_init->library_collection->nfixedfilters-1-i,pow(2,1+i));
		}
		gsl_vector_set(fixedlengths,0,reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration);
	}
	else
	{
		for (int i=0;i<reconstruct_init->library_collection->nfixedfilters;i++)
		{
			gsl_vector_set(fixedlengths,reconstruct_init->library_collection->nfixedfilters-1-i,pow(2,1+i));
		}
	}
		
	int index = 0;
	gsl_vector_view temp;
	for (int i=0;i<reconstruct_init->library_collection->nfixedfilters;i++)
	{
		if (gsl_vector_get(fixedlengths,i) == sizeFilter)
		{
			if (strcmp(reconstruct_init->FilterDomain,"F") == 0)		temp = gsl_vector_subvector(optimalfilterFound_Aux,index,sizeFilter*2);
			else if (strcmp(reconstruct_init->FilterDomain,"T") == 0)	temp = gsl_vector_subvector(optimalfilterFound_Aux,index,sizeFilter);
		    
			gsl_vector_memcpy(*optimalfilterFound,&temp.vector);
			
			break;
		}
		
		if (strcmp(reconstruct_init->FilterDomain,"F") == 0) 	        index = index + gsl_vector_get(fixedlengths,i)*2;
		else if (strcmp(reconstruct_init->FilterDomain,"T") == 0) 	index = index + gsl_vector_get(fixedlengths,i);
	}
	
	gsl_vector_free(optimalfilterFound_Aux); optimalfilterFound_Aux = 0;
        
        gsl_vector_free(maxDERs_LIB1row); maxDERs_LIB1row = 0;

	return(EPOK);
}
/*xxxx end of SECTION B5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B6 ************************************************************
* find_optimalfilterDAB: This function selects the proper optimal filter from the calibration library ('Tx', 'ABTx', 'Fx' or 'ABFx')
*                        by comparing the maximum value of the pulse derivative ('maxDER') to the list of maximums in the library ('maxDERs') for the DAB interpolation method.
* 			 It also selects the proper row from the column 'PAB'.
* 
* It finds the embracing  'maxDERs' in the calibration library
*   - If 'maxDER' is lower than the lowest 'maxDERs' in the ibrary => The data with the lowest 'maxDERs' (first row) in the library are chosen
*   - If 'maxDER' is higher than the highest 'maxDERs' in the library => The data of the penultimate row in the library are chosen
*
* Parameters:
* - maxDER: Max value of the derivative of the filtered pulse whose matched filter is being sought
* - maxDERs: GSL vector with the maximum values of the derivatives of the templates in the library to be compare with the pulse being analysed
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values).
*                     In particular, this function uses the info in the library ('optimal_filters')
* - optimalfilterFound: GSL vector with the optimal filter selected
* - PabFound: PAB column from the library
* - Ealpha: Energy (in eV) which straddle the 'maxDER' in the lower limit
* - Ebeta: Energy (in eV) which straddle the 'maxDER' in the higher limit
****************************************/
int find_optimalfilterDAB(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **optimalfilterFound, gsl_vector **PabFound, double *Ealpha, double *Ebeta)
{
	string message = "";
	char valERROR[256];
        
        long nummodels = maxDERs->size;
	int sizeFilter = (*optimalfilterFound)->size;
	
	gsl_vector *optimalfilterFound_Aux;
        if ((optimalfilterFound_Aux = gsl_vector_alloc(reconstruct_init->library_collection->optimal_filters[0].ofilter_duration)) == 0)
        {
                sprintf(valERROR,"%d",__LINE__-2);
                string str(valERROR);
                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
        }
        gsl_vector *PabFound_Aux;
        // It is not necessary to check the allocation because the allocation of 'optimalfilterFound_aux' has been checked previously
        if (((*PabFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration) 
	   && (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != -999))
		PabFound_Aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration);
	else 
		PabFound_Aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);

	if (maxDER < gsl_vector_get(maxDERs,0))
	{
		gsl_vector_memcpy(optimalfilterFound_Aux,reconstruct_init->library_collection->optimal_filters[0].ofilter);

		if (((*PabFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration) 
		   && (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != -999))
                        gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PABMXLFF,0);
		else
                        gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PAB,0);
		
		*Ealpha = 0.0;
		*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
	}
	else if (maxDER > gsl_vector_get(maxDERs,nummodels-1))
	{
                if (nummodels == 1)
                {
                        gsl_vector_memcpy(optimalfilterFound_Aux,reconstruct_init->library_collection->optimal_filters[0].ofilter);
                        if (((*PabFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration) 
                            && (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != -999))
                                gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PABMXLFF,0);
                        else
                                gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PAB,0);
		
                        *Ealpha = 0.0;
                        *Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
                }
                else
                {
                        gsl_vector_memcpy(optimalfilterFound_Aux,reconstruct_init->library_collection->optimal_filters[nummodels-2].ofilter);
                        
                        if (((*PabFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration) 
                        && (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != -999))
                                gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PABMXLFF,nummodels-2);	//!!!!!! -2 !!!!!!!
                                                                                                                                // -1 because of the GSL indexes start in 0
                                                                                                                                // -1 because the info AB between two energies A<B is stored in the row of the energy A
                        else
                                gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PAB,nummodels-2);	//!!!!!! -2 !!!!!!!
                                                                                                                        // -1 because of the GSL indexes start in 0
                                                                                                                        // -1 because the info AB between two energies A<B is stored in the row of the energy A
                        *Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-2);
                        *Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1);
                }
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
			if (maxDER == gsl_vector_get(maxDERs,i))
			{
				gsl_vector_memcpy(optimalfilterFound_Aux,reconstruct_init->library_collection->optimal_filters[i].ofilter);
				
				if (((*PabFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration)
				   && (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != -999))
					gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PABMXLFF,i);
				else
					gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PAB,i);
				
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i);

				break;
			}
			else if ((maxDER > gsl_vector_get(maxDERs,i)) && (maxDER < gsl_vector_get(maxDERs,i+1)))
			{
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i+1);
		
				gsl_vector_memcpy(optimalfilterFound_Aux,reconstruct_init->library_collection->optimal_filters[i].ofilter);

				if (((*PabFound)->size == reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration)
				   && (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != -999))
				{
					gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PABMXLFF,i);
				}
				else
					gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PAB,i);
								
				break;
			}
		}
	}
	
	gsl_vector *fixedlengths = gsl_vector_alloc(reconstruct_init->library_collection->nfixedfilters);
	if ((reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != reconstruct_init->library_collection->pulse_templates[0].template_duration)
	   && (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != -999)
	   && (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != 999))
	{
		for (int i=0;i<reconstruct_init->library_collection->nfixedfilters-1;i++)
		{
			gsl_vector_set(fixedlengths,reconstruct_init->library_collection->nfixedfilters-1-i,pow(2,1+i));
		}
		gsl_vector_set(fixedlengths,0,reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration);
	}
	else
	{
		for (int i=0;i<reconstruct_init->library_collection->nfixedfilters;i++)
		{
			gsl_vector_set(fixedlengths,reconstruct_init->library_collection->nfixedfilters-1-i,pow(2,1+i));
		}
	}
	
	int index = 0;
	gsl_vector_view temp;
	for (int i=0;i<reconstruct_init->library_collection->nfixedfilters;i++)
	{
		if (gsl_vector_get(fixedlengths,i) == (*PabFound)->size)
                {
			temp = gsl_vector_subvector(optimalfilterFound_Aux,index,sizeFilter);
			gsl_vector_memcpy(*optimalfilterFound,&temp.vector);
		    
			if (((*PabFound)->size != reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration)
			   && (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != -999))
			{
				if (strcmp(reconstruct_init->FilterDomain,"F") == 0)	temp = gsl_vector_subvector(PabFound_Aux,0,sizeFilter/2);
				else 							temp = gsl_vector_subvector(PabFound_Aux,0,sizeFilter);
                                
                                gsl_vector_memcpy(*PabFound,&temp.vector);
			}
			else	gsl_vector_memcpy(*PabFound,PabFound_Aux);
		    
			break;
		}
		
		if (strcmp(reconstruct_init->FilterDomain,"F") == 0) 	        index = index + gsl_vector_get(fixedlengths,i)*2;
		else if (strcmp(reconstruct_init->FilterDomain,"T") == 0) 	index = index + gsl_vector_get(fixedlengths,i);
	}
		
	gsl_vector_free(optimalfilterFound_Aux); optimalfilterFound_Aux = 0;
	gsl_vector_free(PabFound_Aux); PabFound_Aux = 0;
	gsl_vector_free(fixedlengths); fixedlengths = 0;
        
	return(EPOK);
}
/*xxxx end of SECTION B6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B7 ************************************************************
* find_prclwn: When EnergyMethod=WEIGHTN this function selects the proper precalculated values ('PRCLx') from the calibration library 
*              by comparing the maximum value of the pulse derivative ('maxDER') to the list of maximums in the library ('maxDERs') for the OFLib=yes.
* 	       It also selects the proper row from the column 'PAB'.
* 
* It finds the embracing  'maxDERs' in the calibration library
*   - If 'maxDER' is lower than the lowest 'maxDERs' in the ibrary => The data with the lowest 'maxDERs' (first row) in the library are chosen
*   - If 'maxDER' is higher than the highest 'maxDERs' in the library => The data of the penultimate row in the library are chosen
*
* Parameters:
* - maxDER: Max value of the derivative of the filtered pulse whose matched filter is being sought
* - maxDERs: GSL vector with the maximum values of the derivatives of the templates in the library to be compare with the pulse being analysed
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values).
*                     In particular, this function uses the info in the library ('optimal_filters')
* - PRCLWNFound: GSL vector with the precalculated matrix 'PRCLWNx' selected
* - PabFound: PAB column from the library
* - Ealpha: Energy (in eV) which straddle the 'maxDER' in the lower limit
* - Ebeta: Energy (in eV) which straddle the 'maxDER' in the higher limit
****************************************/
int find_prclwn(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_matrix **PRCLWNFound, gsl_vector **PabFound, double *Ealpha, double *Ebeta)
{	  
	string message = "";
	char valERROR[256];
	
	long nummodels = maxDERs->size;
	int sizePRCLWN = (*PRCLWNFound)->size2;
	
	gsl_vector *PRCLWNFound_Aux;
        if ((PRCLWNFound_Aux = gsl_vector_alloc(reconstruct_init->library_collection->PRECALWN->size2)) == 0)
        {
                sprintf(valERROR,"%d",__LINE__-2);
                string str(valERROR);
                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
        }
        gsl_vector *PabFound_Aux = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);

	if (maxDER < gsl_vector_get(maxDERs,0))
	{
		gsl_matrix_get_row(PRCLWNFound_Aux,reconstruct_init->library_collection->PRECALWN,0);
		
		gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PAB,0);
		
		*Ealpha = 0.0;
		*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
	}
	else if (maxDER > gsl_vector_get(maxDERs,nummodels-1))
	{
		gsl_matrix_get_row(PRCLWNFound_Aux,reconstruct_init->library_collection->PRECALWN,nummodels-2);
		
		gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PAB,nummodels-2);
		
		*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-2);
		*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1);
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
			if (maxDER == gsl_vector_get(maxDERs,i))
			{
				gsl_matrix_get_row(PRCLWNFound_Aux,reconstruct_init->library_collection->PRECALWN,i);
				
				gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PAB,i);
				
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				
				break;
			}
			else if ((maxDER > gsl_vector_get(maxDERs,i)) && (maxDER < gsl_vector_get(maxDERs,i+1)))
			{
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i+1);
				
				gsl_matrix_get_row(PRCLWNFound_Aux,reconstruct_init->library_collection->PRECALWN,i);
				
				gsl_matrix_get_row(PabFound_Aux,reconstruct_init->library_collection->PAB,i);
								
				break;
			}
		}
	}
	
	gsl_vector *fixedlengths = gsl_vector_alloc(reconstruct_init->library_collection->nfixedfilters);
	for (int i=0;i<reconstruct_init->library_collection->nfixedfilters;i++)
	{
		gsl_vector_set(fixedlengths,reconstruct_init->library_collection->nfixedfilters-1-i,pow(2,1+i));
	}
	
	int index = 0;
	gsl_vector_view temp;
	for (int i=0;i<reconstruct_init->library_collection->nfixedfilters;i++)
	{	
		if (gsl_vector_get(fixedlengths,i) == sizePRCLWN)
		{
			temp = gsl_vector_subvector(PRCLWNFound_Aux,index,sizePRCLWN);
			gsl_matrix_set_row(*PRCLWNFound,0,&temp.vector);
			temp = gsl_vector_subvector(PRCLWNFound_Aux,index+sizePRCLWN,sizePRCLWN);
			gsl_matrix_set_row(*PRCLWNFound,1,&temp.vector);
			
			temp = gsl_vector_subvector(PabFound_Aux,0,sizePRCLWN);
			gsl_vector_memcpy(*PabFound,&temp.vector);
			
			break;
		}
		index = index + gsl_vector_get(fixedlengths,i)*2; 
	}
	
	gsl_vector_free(PRCLWNFound_Aux); PRCLWNFound_Aux = 0;
	gsl_vector_free(PabFound_Aux); PabFound_Aux = 0;
	gsl_vector_free(fixedlengths); fixedlengths = 0;

	return(EPOK);
}
/*xxxx end of SECTION B7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B8 ************************************************************
* find_prclofwm: When EnergyMethod=OPTFILT and OFNoise=WEIGHTM this function selects the proper precalculated values ('OFWx') from the calibration library 
*              by comparing the maximum value of the pulse derivative ('maxDER') to the list of maximums in the library ('maxDERs') for the OFLib=yes.
* 	      
* It finds the embracing  'maxDERs' in the calibration library
*   - If 'maxDER' is lower than the lowest 'maxDERs' in the ibrary => The data with the lowest 'maxDERs' (first row) in the library are chosen
*   - If 'maxDER' is higher than the highest 'maxDERs' in the library => The data of the penultimate row in the library are chosen
*
* Parameters:
* - maxDER: Max value of the derivative of the filtered pulse whose matched filter is being sought
* - maxDERs: GSL vector with the maximum values of the derivatives of the templates in the library to be compare with the pulse being analysed
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values).
*                     In particular, this function uses the info in the library ('optimal_filters')
* - PRCLOFWMFound: GSL vector with the precalculated matrix 'PRCLOFWMx' selected
* - Ealpha: Energy (in eV) which straddle the 'maxDER' in the lower limit
* - Ebeta: Energy (in eV) which straddle the 'maxDER' in the higher limit
****************************************/
int find_prclofwm(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_matrix **PRCLOFWMFound, double *Ealpha, double *Ebeta)
{
	string message = "";
	char valERROR[256];
	
        gsl_vector *maxDERs_LIB1row;
        long nummodels;
        int index_E_LIB1row;
        if (reconstruct_init->filtEev == 0)
        {
                nummodels = maxDERs->size;
                maxDERs_LIB1row = gsl_vector_alloc(nummodels);
                gsl_vector_memcpy(maxDERs_LIB1row,maxDERs);
        }
        else
        {
                nummodels = 1;
                maxDERs_LIB1row = gsl_vector_alloc(1);
                for (int i=0;i<maxDERs->size;i++)
                {
                        if (gsl_vector_get(reconstruct_init->library_collection->energies,i) == reconstruct_init->filtEev)
                        {
                                gsl_vector_set(maxDERs_LIB1row,0,gsl_vector_get(maxDERs,i));
                                index_E_LIB1row = i;
                                break;
                        }
                }
        }
        
	int sizePRCLOFWM = (*PRCLOFWMFound)->size2;
	
	gsl_vector *PRCLOFWMFound_Aux;
        if ((PRCLOFWMFound_Aux = gsl_vector_alloc(reconstruct_init->library_collection->PRCLOFWM->size2)) == 0)
        {
                sprintf(valERROR,"%d",__LINE__-2);
                string str(valERROR);
                message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
        }

	if (maxDER < gsl_vector_get(maxDERs_LIB1row,0))
	{
                if (reconstruct_init->filtEev == 0)
                {
                        gsl_matrix_get_row(PRCLOFWMFound_Aux,reconstruct_init->library_collection->PRCLOFWM,0);
                        *Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
                }
                else
                {
                        gsl_matrix_get_row(PRCLOFWMFound_Aux,reconstruct_init->library_collection->PRCLOFWM,index_E_LIB1row);
                        *Ebeta = reconstruct_init->filtEev;
                }
		
		*Ealpha = 0.0;
	}
	else if (maxDER > gsl_vector_get(maxDERs_LIB1row,nummodels-1))
	{
                if (reconstruct_init->filtEev == 0)
                {
                        gsl_matrix_get_row(PRCLOFWMFound_Aux,reconstruct_init->library_collection->PRCLOFWM,nummodels-1);
                        *Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,0);
                }
                else
                {
                        gsl_matrix_get_row(PRCLOFWMFound_Aux,reconstruct_init->library_collection->PRCLOFWM,index_E_LIB1row);
                        *Ealpha = reconstruct_init->filtEev;
                }
                
                *Ebeta = 0.0;
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
			if (maxDER == gsl_vector_get(maxDERs_LIB1row,i))
			{
				gsl_matrix_get_row(PRCLOFWMFound_Aux,reconstruct_init->library_collection->PRCLOFWM,i);
				
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i);
			}
			else if ((maxDER > gsl_vector_get(maxDERs_LIB1row,i)) && (maxDER < gsl_vector_get(maxDERs_LIB1row,i+1)))
			{		
				gsl_matrix_get_row(PRCLOFWMFound_Aux,reconstruct_init->library_collection->PRCLOFWM,i);
				
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i+1);
								
				break;
			}
		}
	}
	
	gsl_vector *fixedlengths = gsl_vector_alloc(reconstruct_init->library_collection->nfixedfilters);
	if ((reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != reconstruct_init->library_collection->pulse_templates[0].template_duration)
	   && (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != -999)
	   && (reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration != 999))
	{
		for (int i=0;i<reconstruct_init->library_collection->nfixedfilters-1;i++)
		{
			gsl_vector_set(fixedlengths,reconstruct_init->library_collection->nfixedfilters-1-i,pow(2,1+i));
		}
		gsl_vector_set(fixedlengths,0,reconstruct_init->library_collection->pulse_templatesMaxLengthFixedFilter[0].template_duration);
	}
	else
	{
		for (int i=0;i<reconstruct_init->library_collection->nfixedfilters;i++)
		{
			gsl_vector_set(fixedlengths,reconstruct_init->library_collection->nfixedfilters-1-i,pow(2,1+i));
		}
	}
		
	int index = 0;
	gsl_vector_view temp;
	for (int i=0;i<reconstruct_init->library_collection->nfixedfilters;i++)
	{
		if (gsl_vector_get(fixedlengths,i) == sizePRCLOFWM)
		{
			temp = gsl_vector_subvector(PRCLOFWMFound_Aux,index,sizePRCLOFWM);
			gsl_matrix_set_row(*PRCLOFWMFound,0,&temp.vector);
			temp = gsl_vector_subvector(PRCLOFWMFound_Aux,index+sizePRCLOFWM,sizePRCLOFWM);
			gsl_matrix_set_row(*PRCLOFWMFound,1,&temp.vector);
			
			break;
		}
		
		index = index + gsl_vector_get(fixedlengths,i)*2;
	}
	
	gsl_vector_free(PRCLOFWMFound_Aux); PRCLOFWMFound_Aux = 0;
	gsl_vector_free(fixedlengths); fixedlengths = 0;
        
        gsl_vector_free(maxDERs_LIB1row); maxDERs_LIB1row = 0;

	return(EPOK);
}
/*xxxx end of SECTION B8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B9 ************************************************************
* find_Esboundary: This function provides the indexes of the two energies which straddle the pulse energy, by
*                  comparing the maximum value of the pulse derivative ('maxDER') to the list of maximums in the library ('maxDERs').
*
* It finds the two embracing 'maxDERs' in the calibration library
*   - If 'maxDER' is lower than the lowest 'maxDERs' in the library => `indexEalpha`=`indexEbeta`= 0
*   - If 'maxDER' is higher than the higher 'maxDERs' in the pulse models library => `indexEalpha`=`indexEbeta`= Number of templates -1
*
* Parameters:
* - maxDER: Max value of the derivative of the filtered pulse whose matched filter is being sought
* - maxDERs: GSL vector with the maximum values of the derivatives of the templates in the library to be compare with the pulse being analysed
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values).
*                     In particular, this function uses the info in the library ('energies')
* - indexEalpha: Index of the energy lower than the energy of the pulse which is being analyzed
* - indexEbeta: Index of the energy higher than the energy of the pulse which is being analyzed
* - Ealpha: Energy (in eV) which straddle the 'maxDER' in the lower limit
* - Ebeta: Energy (in eV) which straddle the 'maxDER' in the higher limit
****************************************/
int find_Esboundary(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, int *indexEalpha, int *indexEbeta, double *Ealpha, double *Ebeta)
{
	string message = "";
	int nummodels = maxDERs->size;

	if (maxDER < gsl_vector_get(maxDERs,0))
	{
		*indexEalpha = 0;
		*indexEbeta = 0;
		
		*Ealpha = 0.0;
		*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
	}
	else if (maxDER > gsl_vector_get(maxDERs,nummodels-1))
	{
		*indexEalpha = nummodels - 1;
		*indexEbeta = nummodels - 1;
		
		*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-2);
		*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1);
	}
	else
	{
		for (int i=0;i<nummodels-1;i++)
		{
			if (maxDER == gsl_vector_get(maxDERs,i))
			{
				*indexEalpha = i;
				//*indexEbeta = i;
				*indexEbeta = i+1;
				
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				//*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i+1);

				break;
			}
			else if ((maxDER > gsl_vector_get(maxDERs,i)) && (maxDER < gsl_vector_get(maxDERs,i+1)))
			{
				*indexEalpha = i;
				*indexEbeta = i+1;
				
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i+1);
				
				break;
			}
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION B9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B10 ************************************************************
* pulseGrading: This function provides the pulse grade (Pileup=-2, Rejected=-1, HighRes=1, MidRes=2, Limres= 3, LowRes=4) and the optimal filter length
*               by taking into account the info read from the XML file and the 'OFStrategy' (FREE, BASE2, BYGRADE or FIXED).
*
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values).
*                     In particular, this function uses the 'OFLength' and 'grading'.
* - grade1: Pulse duration (length of the optimal filter applied)
* - grade2: Difference between the start time of the pulse and the start time of the previous pulse
* - OFlength_strategy: 'OFStrategy' (input)
* - pulseGrade: Pulse grade (output)
* - OFlength: Optimal filter length (='OFLength' only if 'OFStrategy'=FIXED and 'OFLength' <= grade1) (output)
****************************************/
int pulseGrading (ReconstructInitSIRENA *reconstruct_init, int grade1, int grade2, int OFlength_strategy, int *pulseGrade, long *OFlength)
{
	string message = "";
	char valERROR[256];
        	
	// For the moment, in the 'xifu_detector_hex_baseline.xml' there is not info about all the pixels (SPA, LPA1, LPA2 or LPA3)
        // It works with the current structure in the .xml file (maybe NOT with a different structure)
        int L2 = gsl_matrix_get(reconstruct_init->grading->gradeData,3,1);	// 'gradelim_post' if 'value'(grading num) = 4
        int LIMITED = gsl_matrix_get(reconstruct_init->grading->gradeData,2,1);	// 'gradelim_post' if 'value'(grading num) = 3
	int M1 = gsl_matrix_get(reconstruct_init->grading->gradeData,1,1);	// 'gradelim_post' if 'value'(grading num) = 2
	int H1 = gsl_matrix_get(reconstruct_init->grading->gradeData,0,1);	// 'gradelim_post' if 'value'(grading num) = 1
        /*cout<<"H1: "<<H1<<endl;
        cout<<"M1: "<<M1<<endl;
        cout<<"LIMITED: "<<LIMITED<<endl;
        cout<<"L2: "<<L2<<endl;*/
	gsl_vector *gradelim;
	if ((gradelim = gsl_vector_alloc(reconstruct_init->grading->ngrades)) == 0)
	{
		sprintf(valERROR,"%d",__LINE__-2);
		string str(valERROR);
	        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
        if (gradelim->size < 4)
	{
		sprintf(valERROR,"%d",__LINE__+7);
		string str(valERROR);
		message = "Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
        for (int i=0;i<reconstruct_init->grading->ngrades;i++)
                gsl_vector_set(gradelim,i,gsl_matrix_get(reconstruct_init->grading->gradeData,i,0));	
	int gradelim_pre = gsl_vector_max(gradelim);
	gsl_vector_free(gradelim);
	
	// pulseGrade and OF length
        if (OFlength_strategy == 0)		*OFlength = grade1;
        else if (OFlength_strategy == 1) 	*OFlength = pow(2,floor(log2(grade1)));
        else if (OFlength_strategy == 3)
        {
                *OFlength = min(reconstruct_init->OFLength,grade1);
                if (reconstruct_init->OFLength > grade1)
                {
                        message = "OFLength provided as input parameter > Pulse duration (there can be a pulse in its tail) => OFLength=Pulse duration";
                        EP_PRINT_ERROR(message,-999);	// Only a warning
                }
        }
        
	if (grade1 >= H1)	// High res
        {
                *pulseGrade = 1;	
                if (OFlength_strategy == 2) 	*OFlength = H1;
        }
        else if (grade1 >= M1)	// Medium res
        {
                *pulseGrade = 2;
                if (OFlength_strategy == 2)	*OFlength = M1;	
        }
        else if (grade1 >= LIMITED)	// Limited res
        {
                *pulseGrade = 3;
                if (OFlength_strategy == 2)	*OFlength = LIMITED;
        }
        else if (grade1 > L2)	// Low res
        {
                *pulseGrade = 4;
                if (OFlength_strategy == 2) 	*OFlength = L2;
        }
        else if (grade1 <= L2)	//Pileup
        {
                *pulseGrade = -2;
                if (OFlength_strategy == 2) 	*OFlength = grade1; 
        }  
        
        if ((grade2 < gradelim_pre) || (grade1 == -1))	// Rejected: It is distinguished by the grade but currently its energy is also calculated (by using the info of the next pulse to establish the OFLength)
	{
		*pulseGrade = -1;	
		////// It has no sense to provide the '*OFlength' because its energy will not be calculated
		//// It does not matter '*OFlength' because its energy will not be calculated
		////*OFlength = 1;
		//*OFlength = grade1;
	}

	return EPOK;
}
/*xxxx end of SECTION B10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B11 ************************************************************
* calculateEnergy function: This function calculates the energy of a pulse ('vector') depending on
*                           the 'EnergyMethod', 'OFNoise' and the 'FilterDomain' basically.
*
* OPTFILT (= I2R =I2RALL =I2RNOL =I2RFITTED) and NSD: Optimal filter = Wiener filter
*
*   Once the filter template has been created ('filter' or 'filterFFT'), pulse height analysis is performed by aligning the template
*   with a pulse and multiplying each point in the template by the corresponding point in the pulse. The sum of these products is the energy.
*
* 	Time domain: E = SUM p(t)·of(t)
*
* 	Frequency domain: E = SUM P(f)·OF(f)
*
* 	In practice, the alignment of the pulse relative to the trigger is not completely accurate, so a number of n lags could be used in
* 	order to find the peak value of the energy. The n peak values are fitted to a parabola to find the most accurate energy and a corrected starting time.
*
* 	(*) IXO Onboard processing Trade-Off
* 
* 	If 'OFInterp'=DAB, E = SUM {(p(t)-PAB(t))·of(t)} or E = SUM {(P(f)-PAB(f))·OF(f)}
* 
* OPTFILT and WEIGHTM: 
* 
* 	PRCLOFWM is =(X'.W.X)^(-1) .X'.W) being W the noise weight matrix calculated from noise intervals.
*
* WEIGHT:
*
* 	It is an algorithm for the least-squares measure of pulse size in the presence of non-stationary system noise through the signal event, which
* 	is optimal in both the linear and nonlinear regime. In the case of a linear detector with stationary noise, the algorithm reduces to the Wiener filter.
*
* 	On the transition the TES response is approximately linear, but large changes in resistance can result in a significant change in noise level.
* 	So, modelling the response is not enough. To get the optimum performace, we must also model the system noise. We describe the noise by its covariance
* 	matrix <didj>, which for an ensemble of time series of length n, is a nxn matrix.
*
* 	Consider a set of measurements of the current in the dectector element, Si. An average of these is used as a model template, Mi=<Si>=(1/N)SUM(p=1,N){Sip}
* 	(Mi is the i-sample of the pulseaverage). And the deviationsfrom this mean, Di=Si-Mi, are used to construct a covariance matrix, Vij=<DiDj> (weight
* 	matrix W=1/V). (W is the noise weight matrix calculated from the subtraction of the pulse model from the pulses)
*
*	Given a number of models, M, with their associated weight matrices W, only a crude estimation of signal size is sufficient to determine which two
*	calibration points alpha and beta straddle the unkonwn signal U (the pulse whose energy we want to calculate).
*
*	With a linear interpolation of the signal and weight matrix the best energy estimate is:
*
*   E = Ealpha + (Ebeta-Ealpha)(r/3){(2DZ-1)+sqrt[(2DZ-1)^2+3(2DY-DXD)/r]}
*
*   where D=U-Salpha. U and Salpha are signals without baseline, i.e., we are assuming that the baseline is known or that the baseline is constant (it is the same in 
*   calibration and during the measurement).
* 
*   The terms T=Sbeta-Salpha, t=TWalphaT, X=(Wbeta-Walpha)/t, Y=WalphaT/t, Z=XT and r=1/(ZT) can be precalculated with the calibration
*   data alone.
*
* 	(*) Fixsen et all, "Pulse estimation in nonllinear detectors with nonstationary noise"
*
* WEIGHTN:
*
*	The starting idea is the same as in WEIGHT, i.e. minimizing (S-M)W(S-M) but without interpolating W (W is the noise weight matrix calculated from the subtraction of the pulse model from the pulses). To do that, a first order development of the pulse
*	shape between 2 calibration points:
*
*	p(E) = p(Ea) + (E-Ea)/(Eb-Ea) . (p(Eb)-p(Ea)) (p without the baseline) => P(E)- bm - P(Ea) - b0 + Ea/(Eb-Ea) . (P(Eb)-P(Ea)) = E . (P(Eb) - P(Ea))/(Eb - Ea)
* 
*       bm -> Baseline of the measured data
*       b0 -> Baseline in the calibration
*
*   => P(E) - Pab = E . Dab + (bm-b0) => Datos = E · Modelo + Baseline => y = E·x + B (condition equation)
*
*   where Pab = P(Ea) - Ea/(Eb-Ea) * (P(Eb)-P(Ea))  and Dab = (P(Eb) - P(Ea))/(Eb - Ea) can be pre-calculated. That way, you see that if you subtract Pab
*   to your data, you end up again in an optimal filter like situation with your data modeled by something that is proportional to a template
*   (that here would be more accurately called a "differential" template).
* 
*   Data (P(E)) and the model P(Ea) are on the top of a baseline.
* 
*   Anyway, this is now a linear chi^2 problem which has an exact solution:
*
*   y = E·X + B => y0 = E·x0 + B => SUM{xiyi} = E·SUM{xi^2} + mB => Matricial expression => chi^2 taking into account the errors
*                  y1 = E·x1 + B                                   X'Y = |E|·X'X           chi^2 = SUM{(data-model)^2}/sigma^2
*                  ...                                                   |B|
*                  yN = E·xN + B                                   |E| = (1/(X'X))X'Y      y/sigma = E·x/sigma
*                                                                  |B|
*                                                            
*                                                                  |E| = (1/(X'WX))X'WY
*                                                                  |B|
*
*   |E| = (X'.W.X)^(-1) .X'.W.Y with ' being the transposition operation and Y=P-Pab=|y0| and X=Dab=|x0 1|
*   |B|                                                                              |. |           |.  1|
*                                                                                    |ym|           |xm 1|
*   If OFLib=no, PAB, DAB and WAB columns from the library are used to calculate (X'.W.X)^(-1) .X'.W.Y. In particular, if pulse size is different from the 
*   template length in the library, instead of wrongly cutting WAB (because of the inverse to convert V into W) it is necessary to work with COVARM column and recalculate *   WAB for the appropiate length.
*   If OFLib=yes, information in the PRECALWN HDU (PRCLx) of the library is used (=(X'.W.X)^(-1) .X'.W).
* 
* Parameters:
* - vector: Pulse whose energy has to be determined (if LagsOrNot=1 => Pulse is numlags-1 samples longer but only filterFFT->size samples will be used each lag)
* - pulseGrade: Grade of the input pulse (to decide whether a full or only a rough estimation of energy is required)
* - filter: Optimal filter in time domain
* - filterFFT: Optimal filter in frequency domain
* - runEMethod: 'EnergyMethod' = OPTFILT => 'runEMethod' = 0
* 		'EnergyMethod' = I2R => 'runEMethod' = 0
* 		'EnergyMethod' = I2RALL => 'runEMethod' = 0
* 		'EnergyMethod' = I2RNOL => 'runEMethod' = 0
* 		'EnergyMethod' = I2RFITTED => 'runEMethod' = 0
* 		'EnergyMethod' = WEIGHT => 'runEMethod' = 1
* 		'EnergyMethod' = WEIGHTN => 'runEMethod' = 2
* - indexEalpha: Index of the energy lower than the energy of the pulse which is being analyzed
* - indexEbeta: Index of the energy higher than the energy of the pulse which is being analyzed
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values)
* - domain: 'FilterDomain' = T => 'domain' = 0
*           'FilterDomain' = F => 'domain' = 1
* - samprate: Sampling rate
* - Pab: PAB column in the library
* - PRCLWN: Appropriate PRCLx column in the library
* - PRCLOFWM: Appropriate OFWx column in the library
* - calculatedEnergy: Calculated energy in eV
*                     If 'pulseGrade'=-1 (rejected) the provided calculated energy is -1
* - tstartNewDev: Additional deviation of the tstart (if lags)
* - lagsShift: Number of samples shifted to find the maximum of the parabola
* - LowRes: 1 if the low resolution energy estimator (without lags) is going to be calculated
****************************************************************************/
int calculateEnergy (gsl_vector *vector, int pulseGrade, gsl_vector *filter, gsl_vector_complex *filterFFT,int runEMethod, int indexEalpha, int indexEbeta, ReconstructInitSIRENA *reconstruct_init, int domain, double samprate, gsl_vector *Pab, gsl_matrix *PRCLWN, gsl_matrix *PRCLOFWM, double *calculatedEnergy, double *tstartNewDev, int *lagsShift, int LowRes)
{
        //cout<<"filter->size: "<<filter->size<<endl;
        //cout<<"filterFFT->size: "<<filterFFT->size<<endl;
        //cout<<"vector->size: "<<vector->size<<endl;
        
        string message = "";
	char valERROR[256];

        // To calculate Elowres (optimal filter length = 4 samples)
	double LagsOrNot = reconstruct_init->LagsOrNot;
        char *OFNoise = reconstruct_init->OFNoise;
        if (LowRes == 1)  
        {
            LagsOrNot = 0;
            runEMethod = 0;  
            strncpy(OFNoise,"NSD",4);
        }
        
        *tstartNewDev = 0;
    
        int numlags; // Different from numlags = nLags/2
        if (reconstruct_init->Fitting35 == 3)         numlags = 3;
        else if (reconstruct_init->Fitting35 == 5)    numlags = 5;
        *lagsShift = 0;
        
        double calculatedEnergy2 = 0.0;
	
        //if ((vector->size <= numlags) && (runEMethod == 0) && (strcmp(reconstruct_init->OFNoise,"NSD") == 0))
        if ((vector->size <= numlags) && (runEMethod == 0) && (strcmp(OFNoise,"NSD") == 0))
        {
                *calculatedEnergy = -1.0;
        }
	/*if (pulseGrade == 4)
	{
	    *calculatedEnergy = -1;
	}
	else
	{*/
        else
        {
		//if ((runEMethod == 0) && (strcmp(reconstruct_init->OFNoise,"NSD") == 0))
                if ((runEMethod == 0) && (strcmp(OFNoise,"NSD") == 0))
		// OPTFILT	I2R or I2RALL or I2RNOL or I2RFITTED => OPTFILT
		{
                        //cout<<"NSD"<<endl;
                        gsl_vector_view temp;
                        
			if (strcmp(reconstruct_init->OFInterp,"DAB") == 0)	// DAB
			{
				gsl_vector_scale(Pab,-1.0);
				gsl_vector_add(vector,Pab);
			}
			
			gsl_vector *lags_vector;
			//if ((LagsOrNot == 0) || (pulseGrade == 4))	//NOLAGS or LowRes
                        if (LagsOrNot == 0)	//NOLAGS
			{
				numlags = 1;
				lags_vector = gsl_vector_alloc(numlags);
				gsl_vector_set(lags_vector,0,0);
			}
			else
			{
				lags_vector = gsl_vector_alloc(numlags);
				for (int i=0;i<numlags;i++)
				{
					gsl_vector_set(lags_vector,i,-numlags/2+i);
				}
			}
			
			int index_vector;
			// It is not necessary to check the allocation because 'numlags' has been fixed > 0 (1 or 5)
			gsl_vector *calculatedEnergy_vector = gsl_vector_alloc(numlags);
			gsl_vector_set_zero(calculatedEnergy_vector);
			double a,b,c;
			double xmax;

			double SelectedTimeDuration;
			if (domain == 0)	SelectedTimeDuration = filter->size/samprate;
			else 			SelectedTimeDuration = filterFFT->size/samprate;
			*calculatedEnergy = 0.0;
                    
			if (domain == 0)	// Time domain filtering
			{
                                //cout<<"Time Domain"<<endl;
				if ((numlags == 0) && (vector->size != filter->size)) *calculatedEnergy = 0.0;
				else
				{
					//if ((pulseGrade == 4) || (LagsOrNot == 0))	*calculatedEnergy = gsl_vector_get(calculatedEnergy_vector,0);
					if (LagsOrNot == 0)   
                                        {
                                                
                                                for (int i=0;i<filter->size;i++)
                                                {
                                                    gsl_vector_set(calculatedEnergy_vector,0,gsl_vector_get(calculatedEnergy_vector,0)+gsl_vector_get(vector,i+0)*gsl_vector_get(filter,i));
                                                    //if (j==1) cout<<i<<" "<<gsl_vector_get(vector,i+numlagsOUT/2+j-1)<<" "<<gsl_vector_get(filter,i)<<" "<<gsl_vector_get(calculatedEnergy_vector,j)<<endl;
                                                }
                                                gsl_vector_set(calculatedEnergy_vector,0,fabs(gsl_vector_get(calculatedEnergy_vector,0)*2*SelectedTimeDuration));
                                                    
                                                *calculatedEnergy = gsl_vector_get(calculatedEnergy_vector,0);
                                                //cout<<"energy: "<<*calculatedEnergy<<endl;
                                                
                                        }
					else
					{
                                                //cout<<"LagsOrNot=1"<<endl;
                                                    
                                                int indexmax; 
                                                int indexLags = 0;
                                                int newLag = 0;
                                                bool exitLags = false;
                                                double newEnergy;
                                                                                        
                                                if ((reconstruct_init->Fitting35) == 3)   // Parabola by using 3 points
                                                {   
                                                    //cout<<"Parabola"<<endl;
                                                    for (int j=0;j<numlags;j++)
                                                    {
                                                            for (int i=0;i<filter->size;i++)
                                                            {
                                                                    gsl_vector_set(calculatedEnergy_vector,j,gsl_vector_get(calculatedEnergy_vector,j)+gsl_vector_get(vector,i+(reconstruct_init->nLags)/2+j-1)*gsl_vector_get(filter,i));
                                                                    //if (i<10 || i>8185)   cout<<"j="<<j<<" i="<<i<<" "<<gsl_vector_get(vector,i+(reconstruct_init->nLags)/2+j-1)<<" "<<gsl_vector_get(filter,i)<<endl;
                                                            }
                                                            gsl_vector_set(calculatedEnergy_vector,j,fabs(gsl_vector_get(calculatedEnergy_vector,j)*2*SelectedTimeDuration));
                                                            //cout<<"calculatedEnergy_vector("<<j<<"): "<<gsl_vector_get(calculatedEnergy_vector,j)<<endl;
                                                    }
                                                    indexmax = gsl_vector_max_index(calculatedEnergy_vector);
                                                    //cout<<"indexmax= "<<indexmax<<endl;
                                                    
                                                    if (parabola3Pts (lags_vector, calculatedEnergy_vector, &a, &b, &c))
                                                    {
                                                        message = "Cannot run routine parabola3Pts";
                                                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                    }
                                                    xmax = -b/(2*a);
                                                    //if ((xmax < -1) || (xmax > 1))
                                                    if (((xmax < -1) || (xmax > 1)) && (reconstruct_init->nLags > 3))
                                                    {
                                                        do
                                                        {  
                                                            indexLags = indexLags + 1;
                                                            //cout<<"indexLags= "<<indexLags<<endl;
                                                            if (indexmax == 0)  
                                                            {      
                                                                //newLag = gsl_vector_get(lags_vector,0)-indexLags;
                                                                gsl_vector_set(calculatedEnergy_vector,2,gsl_vector_get(calculatedEnergy_vector,1));
                                                                gsl_vector_set(calculatedEnergy_vector,1,gsl_vector_get(calculatedEnergy_vector,0));
                                                                
                                                                *lagsShift = *lagsShift - 1;
                                                            }
                                                            else    
                                                            {
                                                                //newLag = gsl_vector_get(lags_vector,2)+indexLags;
                                                                gsl_vector_set(calculatedEnergy_vector,0,gsl_vector_get(calculatedEnergy_vector,1));
                                                                gsl_vector_set(calculatedEnergy_vector,1,gsl_vector_get(calculatedEnergy_vector,2));
                                                                
                                                                *lagsShift = *lagsShift + 1;
                                                            }
                                                            //cout<<"indexmax= "<<indexmax<<endl;
                                                            //cout<<"newLag= "<<newLag<<endl;
                                                            
                                                            newEnergy = 0.0;
                                                            for (int k=0;k<filter->size;k++)
                                                            {
                                                                newEnergy = newEnergy + gsl_vector_get(vector,(reconstruct_init->nLags)/2+newLag+k)*gsl_vector_get(filter,k);
                                                                //cout<<k<<" "<<gsl_vector_get(vector,k+5)<<" "<<gsl_vector_get(filter,k)<<" "<<newEnergy<<endl;
                                                            }
                                                            newEnergy = fabs(newEnergy*2*SelectedTimeDuration);
                                                            //cout<<"newEnergy: "<<newEnergy<<endl;
                                                            
                                                            if (indexmax == 0)
                                                            {
                                                                /*gsl_vector_set(lags_vector,0,newLag);
                                                                gsl_vector_set(lags_vector,1,newLag+1);
                                                                gsl_vector_set(lags_vector,2,newLag+2);*/
                                                                gsl_vector_set(calculatedEnergy_vector,0,newEnergy);
                                                            }
                                                            else if (indexmax == 2)
                                                            {
                                                                /*gsl_vector_set(lags_vector,2,newLag);
                                                                gsl_vector_set(lags_vector,1,newLag-1);
                                                                gsl_vector_set(lags_vector,0,newLag-2);*/
                                                                gsl_vector_set(calculatedEnergy_vector,2,newEnergy);
                                                            }
                                                            
                                                            /*cout<<"lAGS_vector: "<<gsl_vector_get(lags_vector,0)<<" "<<gsl_vector_get(calculatedEnergy_vector,0)<<endl;
                                                            cout<<"lAGS_vector: "<<gsl_vector_get(lags_vector,1)<<" "<<gsl_vector_get(calculatedEnergy_vector,1)<<endl;
                                                            cout<<"lAGS_vector: "<<gsl_vector_get(lags_vector,2)<<" "<<gsl_vector_get(calculatedEnergy_vector,2)<<endl;*/
                                                            
                                                            if (parabola3Pts (lags_vector, calculatedEnergy_vector, &a, &b, &c))
                                                            {
                                                                message = "Cannot run routine parabola3Pts";
                                                                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                            }
                                                            xmax = -b/(2*a);
                                                            if ((xmax >= -1) && (xmax <= 1))      exitLags = true;
                                                            else                                indexmax = gsl_vector_max_index(calculatedEnergy_vector); 
                                                            /*cout<<"xmax: "<<xmax<<endl;
                                                            cout<<"exitLags: "<<exitLags<<endl;
                                                            cout<<"indexLags: "<<indexLags<<" limite="<<(reconstruct_init->nLags)/2-1<<endl;*/
                                                        } while ((exitLags == false) && (indexLags < (reconstruct_init->nLags)/2-1));
                                                    }
                                                }
                                                else if ((reconstruct_init->Fitting35) == 5) //Fitting by using 5 points
                                                {
                                                    //cout<<"Fitting"<<endl;
                                                    for (int j=0;j<numlags;j++)
                                                    {
                                                            for (int i=0;i<filter->size;i++)
                                                            {
                                                                    gsl_vector_set(calculatedEnergy_vector,j,gsl_vector_get(calculatedEnergy_vector,j)+gsl_vector_get(vector,i+(reconstruct_init->nLags)/2+j-2)*gsl_vector_get(filter,i));
                                                            }
                                                            gsl_vector_set(calculatedEnergy_vector,j,fabs(gsl_vector_get(calculatedEnergy_vector,j)*2*SelectedTimeDuration));
                                                            //cout<<"calculatedEnergy_vector("<<j<<"): "<<gsl_vector_get(calculatedEnergy_vector,j)<<endl;
                                                    }
                                                    indexmax = gsl_vector_max_index(calculatedEnergy_vector);
                                                    //cout<<"indexmax= "<<indexmax<<endl;
                                                    
                                                    if (polyFit(lags_vector, calculatedEnergy_vector, &a, &b, &c))
                                                    {
                                                        message = "Cannot run routine polyFit";
                                                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                    }
                                                    xmax = -b/(2*a);
                                                    //if ((xmax < -2) || (xmax > 2))
                                                    if (((xmax < -2) || (xmax > 2)) && (reconstruct_init->nLags > 5))
                                                    {
                                                        do
                                                        {  
                                                            indexLags = indexLags + 1;
                                                            if (indexmax == 0)  
                                                            {      
                                                                newLag = gsl_vector_get(lags_vector,0)-indexLags;
                                                                gsl_vector_set(calculatedEnergy_vector,4,gsl_vector_get(calculatedEnergy_vector,3));
                                                                gsl_vector_set(calculatedEnergy_vector,3,gsl_vector_get(calculatedEnergy_vector,2));
                                                                gsl_vector_set(calculatedEnergy_vector,2,gsl_vector_get(calculatedEnergy_vector,1));
                                                                gsl_vector_set(calculatedEnergy_vector,1,gsl_vector_get(calculatedEnergy_vector,0));
                                                                
                                                                *lagsShift = *lagsShift - 1;
                                                            }
                                                            else if (indexmax == 4)    
                                                            {
                                                                newLag = gsl_vector_get(lags_vector,2)+indexLags;
                                                                gsl_vector_set(calculatedEnergy_vector,0,gsl_vector_get(calculatedEnergy_vector,1));
                                                                gsl_vector_set(calculatedEnergy_vector,1,gsl_vector_get(calculatedEnergy_vector,2));
                                                                gsl_vector_set(calculatedEnergy_vector,2,gsl_vector_get(calculatedEnergy_vector,3));
                                                                gsl_vector_set(calculatedEnergy_vector,3,gsl_vector_get(calculatedEnergy_vector,4));
                                                                
                                                                *lagsShift = *lagsShift + 1;
                                                            }
                                                            //cout<<"indexmax= "<<indexmax<<endl;
                                                            //cout<<"newLag= "<<newLag<<endl;
                                                            
                                                            newEnergy = 0.0;
                                                            for (int k=0;k<filter->size;k++)
                                                            {
                                                                newEnergy = newEnergy + gsl_vector_get(vector,(reconstruct_init->nLags)/2+newLag+k)*gsl_vector_get(filter,k);
                                                                //cout<<k<<" "<<gsl_vector_get(vector,k+5)<<" "<<gsl_vector_get(filter,k)<<" "<<newEnergy<<endl;
                                                            }
                                                            newEnergy = fabs(newEnergy*2*SelectedTimeDuration);
                                                            
                                                            if (indexmax == 0)
                                                            {
                                                                /*gsl_vector_set(lags_vector,0,newLag);
                                                                gsl_vector_set(lags_vector,1,newLag+1);
                                                                gsl_vector_set(lags_vector,2,newLag+2);
                                                                gsl_vector_set(lags_vector,3,newLag+3);
                                                                gsl_vector_set(lags_vector,4,newLag+4);*/
                                                                gsl_vector_set(calculatedEnergy_vector,0,newEnergy);
                                                            }
                                                            else if (indexmax == 4)
                                                            {
                                                                /*gsl_vector_set(lags_vector,4,newLag);
                                                                gsl_vector_set(lags_vector,3,newLag-1);
                                                                gsl_vector_set(lags_vector,2,newLag-2);
                                                                gsl_vector_set(lags_vector,1,newLag-3);
                                                                gsl_vector_set(lags_vector,0,newLag-4);*/
                                                                gsl_vector_set(calculatedEnergy_vector,4,newEnergy);
                                                            }
                                                            
                                                            //cout<<"lAGS_vector: "<<gsl_vector_get(lags_vector,0)<<" "<<gsl_vector_get(calculatedEnergy_vector,0)<<endl;
                                                            //           c out<<"lAGS_vector: "<<gsl_vector_get(lags_vector,1)<<*" "<<gsl_vector_ge*t(calculatedEnergy_vector,1)<<endl;
                                                            //           cout<<"lAGS_vector: "<<gsl_vector_get(lags_vector,2)<<" "<<gsl_vector_get(calculatedEnergy_vector,2)<<endl;
                                                            
                                                            if (polyFit(lags_vector, calculatedEnergy_vector, &a, &b, &c))
                                                            {
                                                                message = "Cannot run routine polyFit";
                                                                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                            }
                                                            xmax = -b/(2*a);
                                                            if ((xmax >= -2) && (xmax <= 2))    exitLags = true;
                                                            else                                indexmax = gsl_vector_max_index(calculatedEnergy_vector); 
                                                            
                                                        } while ((exitLags == false) && (indexLags < (reconstruct_init->nLags)/2-2));
                                                    }
                                                }
                                                    
                                                //xmax = -b/(2*a);
                                                //cout<<"xmax= "<<xmax<<endl;
                                                *tstartNewDev = xmax;
                                                //cout<<"*tstartNewDev= "<<*tstartNewDev<<endl;
                                                *calculatedEnergy = a*pow(xmax,2.0) + b*xmax +c;
                                                //cout<<"*calculatedEnergy= "<<*calculatedEnergy<<endl;
                                                //cout<<"phi= "<<xmax-(*lagsShift)<<endl;
                                                //cout<<"lagsShift= "<<*lagsShift<<endl;                                   
                                                //gsl_vector_free(calculatedEnergy_vectorABS);
                                               
					}
				}
			}
			else if (domain == 1)	// Frequency domain filtering (multiply vectorFFT and filterFFT)
			{
				if ((numlags == 0) && (vector->size != filterFFT->size)) *calculatedEnergy = 0.0;
				else
				{
					// It is not necessary to check the allocation because 'numlags' has been fixed > 0
					gsl_vector_complex *calculatedEnergy_vectorcomplex = gsl_vector_complex_alloc(numlags);
					gsl_vector_complex_set_zero(calculatedEnergy_vectorcomplex);

					// It is not necessary to check the allocation because 'vector' size must already be > 0
					gsl_vector_complex *vectorFFT = gsl_vector_complex_alloc(filterFFT->size);

					*calculatedEnergy = 0.0;
										
					//if ((pulseGrade == 4) || (LagsOrNot == 0))
					if (LagsOrNot == 0)
                                        {
                                                gsl_vector  *vectorSHORT = gsl_vector_alloc(filterFFT->size);
                                                temp = gsl_vector_subvector(vector,0,filterFFT->size);
                                                gsl_vector_memcpy(vectorSHORT,&temp.vector);
                                                if (FFT(vectorSHORT,vectorFFT,SelectedTimeDuration))
                                                {
                                                    message = "Cannot run routine FFT";
                                                    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                }
                                                gsl_vector_free(vectorSHORT);
                                                
                                                for (int i=0; i<filterFFT->size; i++)
                                                {
                                                    gsl_vector_complex_set(calculatedEnergy_vectorcomplex,0,gsl_complex_add(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,0),gsl_complex_mul(gsl_vector_complex_get(vectorFFT,i),gsl_vector_complex_get(filterFFT,i))));
                                                }
                                                gsl_vector_set(calculatedEnergy_vector,0,fabs(GSL_REAL(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,0))));
                                                //cout<<gsl_vector_get(lags_vector,0)<<" "<<gsl_vector_get(calculatedEnergy_vector,0)<<endl;
                                                					
                                                *calculatedEnergy = gsl_vector_get(calculatedEnergy_vector,0);
                                        }
					else
					{  
                                                //cout<<"numlags: "<<numlags<<endl;
                                                //cout<<"vector: "<<gsl_vector_get(vector,0)<<","<<gsl_vector_get(vector,1)<<","<<gsl_vector_get(vector,2)<<","<<gsl_vector_get(vector,3)<<","<<gsl_vector_get(vector,4)<<","<<gsl_vector_get(vector,5)<<endl;
                                                //cout<<"numlagsOUT: "<<numlagsOUT<<endl;
                                                int indexmax;
                                                
                                                bool exitLags = false;
                                                int indexLags = 0;
                                                int newLag = 0;
                                                double newEnergy;
                                                
                                                if ((reconstruct_init->Fitting35) == 3)   // Parabola by using 3 points
                                                {
                                                    for (int j=0;j<numlags;j++)
                                                    {
                                                            gsl_vector  *vectorSHORT = gsl_vector_alloc(filterFFT->size);
                                                            temp = gsl_vector_subvector(vector,(reconstruct_init->nLags)/2+j-1,filterFFT->size);
                                                            gsl_vector_memcpy(vectorSHORT,&temp.vector);
                                                            if (FFT(vectorSHORT,vectorFFT,SelectedTimeDuration))
                                                            {
                                                                    message = "Cannot run routine FFT";
                                                                    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                            }
                                                            gsl_vector_free(vectorSHORT);
                                                    
                                                            for (int i=0; i<filterFFT->size; i++)
                                                            {
                                                                    gsl_vector_complex_set(calculatedEnergy_vectorcomplex,j,gsl_complex_add(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j),gsl_complex_mul(gsl_vector_complex_get(vectorFFT,i),gsl_vector_complex_get(filterFFT,i))));
                                                            }
                                                            gsl_vector_set(calculatedEnergy_vector,j,fabs(GSL_REAL(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j))));
                                                            //cout<<gsl_vector_get(lags_vector,j)<<" "<<gsl_vector_get(calculatedEnergy_vector,j)<<endl;
                                                    }
                                                    indexmax = gsl_vector_max_index(calculatedEnergy_vector);
                                                    //cout<<"indexmax: "<<indexmax<<endl;
                                                
                                                    if (parabola3Pts (lags_vector, calculatedEnergy_vector, &a, &b, &c))
                                                    {
                                                        message = "Cannot run routine parabola3Pts";
                                                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                    }
                                                    xmax = -b/(2*a);
                                                    //cout<<"xmax0: "<<xmax<<endl;
                                                    //if ((xmax < -1) || (xmax > 1))
                                                    if (((xmax < -1) || (xmax > 1)) && (reconstruct_init->nLags > 3))
                                                    {
                                                        do
                                                        {  
                                                            indexLags = indexLags + 1;
                                                            if (indexmax == 0)  
                                                            {      
                                                                newLag = gsl_vector_get(lags_vector,0)-indexLags;
                                                                gsl_vector_set(calculatedEnergy_vector,2,gsl_vector_get(calculatedEnergy_vector,1));
                                                                gsl_vector_set(calculatedEnergy_vector,1,gsl_vector_get(calculatedEnergy_vector,0));
                                                                
                                                                *lagsShift = *lagsShift - 1;
                                                            }
                                                            else    
                                                            {
                                                                newLag = gsl_vector_get(lags_vector,2)+indexLags;
                                                                gsl_vector_set(calculatedEnergy_vector,0,gsl_vector_get(calculatedEnergy_vector,1));
                                                                gsl_vector_set(calculatedEnergy_vector,1,gsl_vector_get(calculatedEnergy_vector,2));
                                                                
                                                                *lagsShift = *lagsShift + 1;
                                                            }
                                                            //cout<<"newLag= "<<newLag<<endl;
                                                            //cout<<"(reconstruct_init->nLags)/2+newLag= "<<(reconstruct_init->nLags)/2+newLag<<endl;
                                                            
                                                            gsl_vector  *vectorSHORT = gsl_vector_alloc(filterFFT->size);
                                                            temp = gsl_vector_subvector(vector,(reconstruct_init->nLags)/2+newLag,filterFFT->size);
                                                            gsl_vector_memcpy(vectorSHORT,&temp.vector);
                                                            if (FFT(vectorSHORT,vectorFFT,SelectedTimeDuration))
                                                            {
                                                                message = "Cannot run routine FFT";
                                                                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                            }
                                                            gsl_vector_free(vectorSHORT);
                                                            
                                                            gsl_complex newEnergyComplex = gsl_complex_rect(0.0,0.0);
                                                            for (int i=0; i<filterFFT->size; i++)
                                                            {
                                                                newEnergyComplex = gsl_complex_add(newEnergyComplex,gsl_complex_mul(gsl_vector_complex_get(vectorFFT,i),gsl_vector_complex_get(filterFFT,i))); 
                                                            }
                                                            newEnergy = fabs(GSL_REAL(newEnergyComplex));
                                                            //cout<<"newEnergy: "<<newEnergy<<endl;
                                                            //gsl_vector_set(calculatedEnergy_vector,j,fabs(GSL_REAL(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j))));
                                                            
                                                            if (indexmax == 0)
                                                            {
                                                                gsl_vector_set(calculatedEnergy_vector,0,newEnergy);
                                                            }
                                                            else if (indexmax == 2)
                                                            {
                                                                gsl_vector_set(calculatedEnergy_vector,2,newEnergy);
                                                            }
                                                            
                                                            /*cout<<"lAGS_vector: "<<gsl_vector_get(lags_vector,0)<<" "<<gsl_vector_get(calculatedEnergy_vector,0)<<endl;
                                                            cout<<"lAGS_vector: "<<gsl_vector_get(lags_vector,1)<<*" "<<gsl_vector_get(calculatedEnergy_vector,1)<<endl;
                                                            cout<<"lAGS_vector: "<<gsl_vector_get(lags_vector,2)<<" "<<gsl_vector_get(calculatedEnergy_vector,2)<<endl;*/
                                                            
                                                            if (parabola3Pts (lags_vector, calculatedEnergy_vector, &a, &b, &c))
                                                            {
                                                                    message = "Cannot run routine parabola3Pts";
                                                                    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                            }
                                                            xmax = -b/(2*a);
                                                            //cout<<"xmaxi: "<<xmax<<endl;
                                                            if ((xmax >= -1) && (xmax <= 1))      exitLags = true;
                                                            else                                indexmax = gsl_vector_max_index(calculatedEnergy_vector); 
                                                                
                                                        } while ((exitLags == false) && (indexLags < (reconstruct_init->nLags)/2-1));
                                                    }
                                                }
                                                else if ((reconstruct_init->Fitting35) == 5) //Fitting by using 5 points
                                                {
                                                    for (int j=0;j<numlags;j++)
                                                    {
                                                            gsl_vector  *vectorSHORT = gsl_vector_alloc(filterFFT->size);
                                                            temp = gsl_vector_subvector(vector,(reconstruct_init->nLags)/2+j-2,filterFFT->size);
                                                            gsl_vector_memcpy(vectorSHORT,&temp.vector);
                                                            if (FFT(vectorSHORT,vectorFFT,SelectedTimeDuration))
                                                            {
                                                                    message = "Cannot run routine FFT";
                                                                    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                            }
                                                            gsl_vector_free(vectorSHORT);
                                                    
                                                            for (int i=0; i<filterFFT->size; i++)
                                                            {
                                                                    gsl_vector_complex_set(calculatedEnergy_vectorcomplex,j,gsl_complex_add(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j),gsl_complex_mul(gsl_vector_complex_get(vectorFFT,i),gsl_vector_complex_get(filterFFT,i))));
                                                            }
                                                            gsl_vector_set(calculatedEnergy_vector,j,fabs(GSL_REAL(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j))));
                                                            //cout<<gsl_vector_get(lags_vector,j)<<" "<<gsl_vector_get(calculatedEnergy_vector,j)<<endl;
                                                    }
                                                    indexmax = gsl_vector_max_index(calculatedEnergy_vector);
                                                    //cout<<"indexmax: "<<indexmax<<endl;
                                                    
                                                    if (polyFit(lags_vector, calculatedEnergy_vector, &a, &b, &c))
                                                    {
                                                        message = "Cannot run routine polyFit";
                                                        EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                    }
                                                    xmax = -b/(2*a);
                                                    //cout<<"xmax0: "<<xmax<<endl;
                                                    //if ((xmax < -2) || (xmax > 2))
                                                    if (((xmax < -1) || (xmax > 1)) && (reconstruct_init->nLags > 5))
                                                    {
                                                        do
                                                        {  
                                                            indexLags = indexLags + 1;
                                                            if (indexmax == 0)  
                                                            {      
                                                                newLag = gsl_vector_get(lags_vector,0)-indexLags;
                                                                gsl_vector_set(calculatedEnergy_vector,4,gsl_vector_get(calculatedEnergy_vector,3));
                                                                gsl_vector_set(calculatedEnergy_vector,3,gsl_vector_get(calculatedEnergy_vector,2));
                                                                gsl_vector_set(calculatedEnergy_vector,2,gsl_vector_get(calculatedEnergy_vector,1));
                                                                gsl_vector_set(calculatedEnergy_vector,1,gsl_vector_get(calculatedEnergy_vector,0));
                                                                
                                                                *lagsShift = *lagsShift - 1;
                                                            }
                                                            else if (indexmax == 4)    
                                                            {
                                                                newLag = gsl_vector_get(lags_vector,2)+indexLags;
                                                                gsl_vector_set(calculatedEnergy_vector,0,gsl_vector_get(calculatedEnergy_vector,1));
                                                                gsl_vector_set(calculatedEnergy_vector,1,gsl_vector_get(calculatedEnergy_vector,2));
                                                                gsl_vector_set(calculatedEnergy_vector,2,gsl_vector_get(calculatedEnergy_vector,3));
                                                                gsl_vector_set(calculatedEnergy_vector,3,gsl_vector_get(calculatedEnergy_vector,4));
                                                                
                                                                *lagsShift = *lagsShift + 1;
                                                            }
                                                            //cout<<"indexmax= "<<indexmax<<endl;
                                                            //cout<<"newLag= "<<newLag<<endl;
                                                            
                                                            gsl_vector  *vectorSHORT = gsl_vector_alloc(filterFFT->size);
                                                            temp = gsl_vector_subvector(vector,(reconstruct_init->nLags)/2+newLag,filterFFT->size);
                                                            gsl_vector_memcpy(vectorSHORT,&temp.vector);
                                                            if (FFT(vectorSHORT,vectorFFT,SelectedTimeDuration))
                                                            {
                                                                message = "Cannot run routine FFT";
                                                                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                            }
                                                            gsl_vector_free(vectorSHORT);
                                                            
                                                            gsl_complex newEnergyComplex = gsl_complex_rect(0.0,0.0);
                                                            for (int i=0; i<filterFFT->size; i++)
                                                            {
                                                                newEnergyComplex = gsl_complex_add(newEnergyComplex,gsl_complex_mul(gsl_vector_complex_get(vectorFFT,i),gsl_vector_complex_get(filterFFT,i))); 
                                                                //gsl_vector_complex_set(calculatedEnergy_vectorcomplex,j,gsl_complex_add(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j),gsl_complex_mul(gsl_vector_complex_get(vectorFFT,i),gsl_vector_complex_get(filterFFT,i))));
                                                            }
                                                            newEnergy = fabs(GSL_REAL(newEnergyComplex));
                                                            //cout<<"newEnergy: "<<newEnergy<<endl;
                                                            
                                                            if (indexmax == 0)
                                                            {
                                                                /*gsl_vector_set(lags_vector,0,newLag);
                                                                gsl_vector_set(lags_vector,1,newLag+1);
                                                                gsl_vector_set(lags_vector,2,newLag+2);
                                                                gsl_vector_set(lags_vector,3,newLag+3);
                                                                gsl_vector_set(lags_vector,4,newLag+4);*/
                                                                gsl_vector_set(calculatedEnergy_vector,0,newEnergy);
                                                            }
                                                            else if (indexmax == 4)
                                                            {
                                                                /*gsl_vector_set(lags_vector,4,newLag);
                                                                gsl_vector_set(lags_vector,3,newLag-1);
                                                                gsl_vector_set(lags_vector,2,newLag-2);
                                                                gsl_vector_set(lags_vector,1,newLag-3);
                                                                gsl_vector_set(lags_vector,0,newLag-4);*/
                                                                gsl_vector_set(calculatedEnergy_vector,4,newEnergy);
                                                            }
                                                            
                                                            /*cout<<"lAGS_vector: "<<gsl_vector_get(lags_vector,0)<<" "<<gsl_vector_get(calculatedEnergy_vector,0)<<endl;
                                                            cout<<"lAGS_vector: "<<gsl_vector_get(lags_vector,1)<<" "<<gsl_vector_get(calculatedEnergy_vector,1)<<endl;
                                                            cout<<"lAGS_vector: "<<gsl_vector_get(lags_vector,2)<<" "<<gsl_vector_get(calculatedEnergy_vector,2)<<endl;*/
                                                            
                                                            if (polyFit(lags_vector, calculatedEnergy_vector, &a, &b, &c))
                                                            {
                                                                message = "Cannot run routine polyFit";
                                                                EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                                                            }
                                                            xmax = -b/(2*a);
                                                            //cout<<"xmaxi: "<<xmax<<endl;
                                                            if ((xmax >= -2) && (xmax <= 2))    exitLags = true;
                                                            else                                indexmax = gsl_vector_max_index(calculatedEnergy_vector); 
                                                            
                                                        } while ((exitLags == false) && (indexLags < (reconstruct_init->nLags)/2-2));
                                                    }
                                                }
                                                
                                                //xmax = -b/(2*a);
                                                //cout<<"xmax= "<<xmax<<endl;
                                                //cout<<xmax<<endl;
                                                //if (*lagsShift==0) cout<<xmax<<endl;
                                                *tstartNewDev = xmax;
                                                //cout<<"*tstartNewDev= "<<*tstartNewDev<<endl;
						*calculatedEnergy = a*pow(xmax,2.0) + b*xmax +c;
                                                //cout<<"*calculatedEnergy= "<<*calculatedEnergy<<endl;
                                                //cout<<"phi= "<<xmax-(*lagsShift)<<endl;
                                                //cout<<"phi= "<<xmax<<endl;
                                                //cout<<"lagsShift= "<<*lagsShift<<endl;
                                                //gsl_vector_free(calculatedEnergy_vectorABS);
					}
					
					gsl_vector_complex_free(calculatedEnergy_vectorcomplex); calculatedEnergy_vectorcomplex = 0;
					gsl_vector_complex_free(vectorFFT); vectorFFT = 0;
				}
			}

			gsl_vector_free(lags_vector); lags_vector = 0;
			gsl_vector_free(calculatedEnergy_vector); calculatedEnergy_vector = 0;
		}
		//else if ((runEMethod == 0) && (strcmp(reconstruct_init->OFNoise,"WEIGHTM") == 0))
		else if ((runEMethod == 0) && (strcmp(OFNoise,"WEIGHTM") == 0))
		{
			gsl_vector *EB = gsl_vector_alloc(2);
		
			gsl_blas_dgemv(CblasNoTrans,1.0,PRCLOFWM,vector,0.0,EB);                           // [(R'WR)^(-1)]R'W\B7pulse
			
			*calculatedEnergy = gsl_vector_get(EB,0);
		
			//gsl_vector_free(P_Pab);
			gsl_vector_free(EB); EB = 0;
		}
		else if (runEMethod == 1) //WEIGHT
		{
			int pulselength = reconstruct_init->pulse_length;
			// It is not necessary to check the allocation because 'vector' size must already be > 0 and 'reconstruct_init->pulse_length'=PulseLength(input parameter) has been checked previously
			gsl_vector *D = gsl_vector_alloc(vector->size);
			gsl_vector *Salpha = gsl_vector_alloc(vector->size);
			gsl_matrix *X = gsl_matrix_alloc(vector->size,vector->size);
			gsl_vector *Y = gsl_vector_alloc(vector->size);
			gsl_vector *Z = gsl_vector_alloc(vector->size);
			double r;
			gsl_vector_view tempv;
			double scalar_aux1;
			double scalar_aux2;
			double scalar_aux3;
			gsl_vector *vector_aux = gsl_vector_alloc(vector->size);
			
			if ((indexEalpha == indexEbeta) && (indexEalpha == 0))		indexEbeta = 1;
			else if ((indexEalpha == indexEbeta) && (indexEalpha != 0))	indexEbeta = reconstruct_init->library_collection->ntemplates-1;
			
			tempv = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_B0[indexEalpha].ptemplate,0,vector->size);
			gsl_vector_memcpy(Salpha,&tempv.vector);	// Salpha
			gsl_vector_memcpy(D,vector);
			gsl_vector_sub(D,Salpha);			// D = U - Salpha
			
			if (vector->size == reconstruct_init->library_collection->pulse_templates[0].template_duration)
			{
				gsl_matrix_get_row(Y,reconstruct_init->library_collection->Y,indexEalpha);	// Y
				
				gsl_vector *vector_auxlong = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
				gsl_matrix_get_row(vector_auxlong,reconstruct_init->library_collection->X,indexEalpha);
				vector2matrix(vector_auxlong,&X);						// X
				gsl_vector_free(vector_auxlong); vector_auxlong = 0;
				
				gsl_matrix_get_row(Z,reconstruct_init->library_collection->Z,indexEalpha);	// Z
				
				r = gsl_vector_get(reconstruct_init->library_collection->r,indexEalpha);	// r
			}
			else
			{
				gsl_permutation *perm = gsl_permutation_alloc(vector->size);
				int s = 0;
				gsl_matrix_view tempm;
				double t;
				gsl_vector *Sbeta = gsl_vector_alloc(vector->size);
				gsl_vector *WalphaT = gsl_vector_alloc(vector->size);
				
				tempv = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_B0[indexEalpha].ptemplate,0,vector->size);
				gsl_vector_memcpy(Salpha,&tempv.vector);					// Salpha
				tempv = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates_B0[indexEbeta].ptemplate,0,vector->size);
				gsl_vector_memcpy(Sbeta,&tempv.vector);						// Sbeta
				gsl_vector *T_short = gsl_vector_alloc(vector->size);
				gsl_vector_memcpy(T_short,Sbeta);
				gsl_vector_sub(T_short,Salpha);							// T_short
				gsl_vector_free(Sbeta); Sbeta = 0;
				
				gsl_vector *covarianceAlphavector = gsl_vector_alloc(reconstruct_init->library_collection->V->size2);
				gsl_vector *covarianceBetavector = gsl_vector_alloc(reconstruct_init->library_collection->V->size2);
				gsl_matrix *covarianceAlphamatrix = gsl_matrix_alloc(sqrt(reconstruct_init->library_collection->V->size2),sqrt(reconstruct_init->library_collection->V->size2));
				gsl_matrix *covarianceBetamatrix = gsl_matrix_alloc(sqrt(reconstruct_init->library_collection->V->size2),sqrt(reconstruct_init->library_collection->V->size2));
				gsl_matrix *covarianceAlphamatrixSHORT = gsl_matrix_alloc(vector->size,vector->size);
				gsl_matrix *covarianceBetamatrixSHORT = gsl_matrix_alloc(vector->size,vector->size);
				gsl_matrix *weightAlphamatrixSHORT = gsl_matrix_alloc(vector->size,vector->size);
				gsl_matrix *weightBetamatrixSHORT = gsl_matrix_alloc(vector->size,vector->size);
				
				gsl_matrix_get_row(covarianceAlphavector,reconstruct_init->library_collection->V,indexEalpha);
				vector2matrix(covarianceAlphavector,&covarianceAlphamatrix);
				gsl_vector_free(covarianceAlphavector); covarianceAlphavector = 0;
				tempm = gsl_matrix_submatrix(covarianceAlphamatrix,0,0,vector->size,vector->size);
				gsl_matrix_memcpy(covarianceAlphamatrixSHORT,&tempm.matrix);
				gsl_matrix_free(covarianceAlphamatrix); covarianceAlphamatrix = 0;
				gsl_linalg_LU_decomp(covarianceAlphamatrixSHORT,perm, &s);
				gsl_linalg_LU_invert(covarianceAlphamatrixSHORT,perm, weightAlphamatrixSHORT);	// Walpha_short
				gsl_matrix_free(covarianceAlphamatrixSHORT); covarianceAlphamatrixSHORT = 0;
				gsl_permutation_free(perm); perm = 0;
				
				gsl_matrix_get_row(covarianceBetavector,reconstruct_init->library_collection->V,indexEbeta);
				vector2matrix(covarianceBetavector,&covarianceBetamatrix);
				gsl_vector_free(covarianceBetavector); covarianceBetavector = 0;
				tempm = gsl_matrix_submatrix(covarianceBetamatrix,0,0,vector->size,vector->size);
				gsl_matrix_memcpy(covarianceBetamatrixSHORT,&tempm.matrix);
				gsl_matrix_free(covarianceBetamatrix); covarianceBetamatrix = 0;
				perm = gsl_permutation_alloc(vector->size);
				s = 0;
				gsl_linalg_LU_decomp(covarianceBetamatrixSHORT,perm, &s);
				gsl_linalg_LU_invert(covarianceBetamatrixSHORT,perm, weightBetamatrixSHORT);	// Wbeta_short
				gsl_matrix_free(covarianceBetamatrixSHORT); covarianceBetamatrixSHORT = 0;
				gsl_permutation_free(perm); perm = 0;
				
				gsl_blas_dgemv(CblasNoTrans, 1.0, weightAlphamatrixSHORT, T_short, 0.0, WalphaT);
				gsl_blas_ddot(T_short,WalphaT,&t);						// t
				
				gsl_vector_memcpy(Y,WalphaT);
				gsl_vector_scale(Y,1/t);							// Y
				gsl_vector_free(WalphaT); WalphaT = 0;
				
				gsl_matrix_memcpy(X,weightBetamatrixSHORT);
				gsl_matrix_sub(X,weightAlphamatrixSHORT);
				gsl_matrix_scale(X,1/t);							// X
	
				gsl_blas_dgemv(CblasNoTrans, 1.0, X, T_short, 0.0, Z);				// Z
	
				gsl_blas_ddot(Z,T_short,&r);
				r=1/r;										// r
								
				gsl_vector_free(T_short); T_short = 0;
				gsl_matrix_free(weightAlphamatrixSHORT); weightAlphamatrixSHORT = 0;
				gsl_matrix_free(weightBetamatrixSHORT); weightBetamatrixSHORT = 0;
			}
			
			gsl_blas_ddot(D,Y,&scalar_aux1);		// DY
			scalar_aux1 = 2*scalar_aux1;			// 2DY

			gsl_blas_dgemv(CblasNoTrans, 1.0, X, D, 0.0, vector_aux);	// XD
			gsl_blas_ddot(D,vector_aux,&scalar_aux2);	// DXD

			scalar_aux1 = scalar_aux1-scalar_aux2;		// 2DY-DXD
			scalar_aux1 = 3*scalar_aux1/r;			// 3(2DY-DXD)/r

			gsl_blas_ddot(D,Z,&scalar_aux3);		// DZ
			scalar_aux2 = pow(2*scalar_aux3-1,2.0);		// (2DZ-1)\B2

			scalar_aux1 = sqrt(scalar_aux1+scalar_aux2);	// sqrt[(2DZ-1)\B2 + 3(2DY-DXD)/r]

			scalar_aux3 = 2*scalar_aux3-1;			// (2DZ-1)

			// (Ebeta-Ealpha)*(r/3)*{(2DZ-1) + sqrt[(2DZ-1)\B2 + 3(2DY-DXD)/r]}
			scalar_aux1 = (scalar_aux3 + scalar_aux1)*(r/3)*(gsl_vector_get(reconstruct_init->library_collection->energies,indexEbeta)-gsl_vector_get(reconstruct_init->library_collection->energies,indexEalpha));

			// Ealpha + (Ebeta-Ealpha)*(r/3)*{(2DZ-1) + sqrt[(2DZ-1)\B2 + 3(2DY-DXD)/r]}
			*calculatedEnergy = gsl_vector_get(reconstruct_init->library_collection->energies,indexEalpha) + scalar_aux1;

			gsl_vector_free(D); D = 0;
			gsl_vector_free(Salpha); Salpha = 0;
			gsl_matrix_free(X); X = 0;
			gsl_vector_free(Y); Y = 0;
			gsl_vector_free(Z); Z = 0;
			gsl_vector_free(vector_aux); vector_aux = 0;
		}
		else if (runEMethod == 2) //WEIGHTN
		{
			if (reconstruct_init->OFLib == 0)
			{
				// It is not necessary to check the allocation because 'vector' size must already be > 0 and 'reconstruct_init->pulse_length'=PulseLength(input parameter) has been checked previously
				gsl_vector *Pab = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);
				gsl_vector *P_Pab = gsl_vector_alloc(vector->size);
				gsl_vector *Dab = gsl_vector_alloc(reconstruct_init->library_collection->pulse_templates[0].template_duration);
				gsl_matrix *Wm_short = gsl_matrix_alloc(vector->size,vector->size);
				//gsl_permutation *perm = 0;
                                gsl_permutation *perm = gsl_permutation_alloc(vector->size);
				int s = 0;
			      
				gsl_matrix_get_row(Dab,reconstruct_init->library_collection->DAB,indexEalpha);  // Dab
				// It is not necessary to check the allocation because 'vector' size must already be > 0
				gsl_vector *Dab_short = gsl_vector_alloc(vector->size);
				gsl_vector_view temp;
				temp = gsl_vector_subvector(Dab,0,vector->size);
				gsl_vector_memcpy(Dab_short,&temp.vector);

				gsl_vector_memcpy(P_Pab,vector);
				//if (vector->size == 1) cout<<gsl_vector_get(vector,0)<<endl;
				gsl_matrix_get_row(Pab,reconstruct_init->library_collection->PAB,indexEalpha);  // Pab
				// It is not necessary to check the allocation because 'vector' size must already be > 0
				gsl_vector *Pab_short = gsl_vector_alloc(vector->size);
				temp = gsl_vector_subvector(Pab,0,vector->size);
				gsl_vector_memcpy(Pab_short,&temp.vector);
				gsl_vector_sub(P_Pab,Pab_short);                                                // P-Pab
				gsl_vector_free(Pab_short); Pab_short = 0;
			      
				if (vector->size == reconstruct_init->library_collection->pulse_templates[0].template_duration)
				{
					gsl_vector *Wabv = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
					gsl_matrix_get_row(Wabv,reconstruct_init->library_collection->WAB,indexEalpha); // WAB vector => (Walpha + Wbeta)/2
					vector2matrix(Wabv,&Wm_short);    
					gsl_vector_free(Wabv); Wabv = 0;
				}
				else if (vector->size != reconstruct_init->library_collection->pulse_templates[0].template_duration)
				{
					if ((indexEalpha == indexEbeta) && (indexEalpha == 0))		indexEbeta = 1;
					else if ((indexEalpha == indexEbeta) && (indexEalpha != 0))	indexEbeta = reconstruct_init->library_collection->ntemplates-1;
			
					gsl_vector *covarianceAlphavector = gsl_vector_alloc(reconstruct_init->library_collection->V->size2);
					gsl_vector *covarianceBetavector = gsl_vector_alloc(reconstruct_init->library_collection->V->size2);
					gsl_matrix *covarianceAlphamatrix = gsl_matrix_alloc(sqrt(reconstruct_init->library_collection->V->size2),sqrt(reconstruct_init->library_collection->V->size2));
					gsl_matrix *covarianceBetamatrix = gsl_matrix_alloc(sqrt(reconstruct_init->library_collection->V->size2),sqrt(reconstruct_init->library_collection->V->size2));
					gsl_matrix *covarianceAlphamatrixSHORT = gsl_matrix_alloc(vector->size,vector->size);
					gsl_matrix *covarianceBetamatrixSHORT = gsl_matrix_alloc(vector->size,vector->size);
					gsl_matrix *weightAlphamatrixSHORT = gsl_matrix_alloc(vector->size,vector->size);
					gsl_matrix *weightBetamatrixSHORT = gsl_matrix_alloc(vector->size,vector->size);
					gsl_matrix_view tempm;
                                        
                                        perm  = gsl_permutation_alloc(vector->size);
					
					gsl_matrix_get_row(covarianceAlphavector,reconstruct_init->library_collection->V,indexEalpha);
					vector2matrix(covarianceAlphavector,&covarianceAlphamatrix);
					gsl_vector_free(covarianceAlphavector); covarianceAlphavector = 0;
					tempm = gsl_matrix_submatrix(covarianceAlphamatrix,0,0,vector->size,vector->size);
					gsl_matrix_memcpy(covarianceAlphamatrixSHORT,&tempm.matrix);
					gsl_matrix_free(covarianceAlphamatrix); covarianceAlphamatrix = 0;
					gsl_linalg_LU_decomp(covarianceAlphamatrixSHORT,perm, &s);
					gsl_linalg_LU_invert(covarianceAlphamatrixSHORT,perm, weightAlphamatrixSHORT);
					gsl_matrix_free(covarianceAlphamatrixSHORT); covarianceAlphamatrixSHORT = 0;
					gsl_permutation_free(perm); perm = 0;
					
					gsl_matrix_get_row(covarianceBetavector,reconstruct_init->library_collection->V,indexEbeta);
					vector2matrix(covarianceBetavector,&covarianceBetamatrix);
					gsl_vector_free(covarianceBetavector); covarianceBetavector = 0;
					tempm = gsl_matrix_submatrix(covarianceBetamatrix,0,0,vector->size,vector->size);
					gsl_matrix_memcpy(covarianceBetamatrixSHORT,&tempm.matrix);
					gsl_matrix_free(covarianceBetamatrix); covarianceBetamatrix = 0;
					perm = gsl_permutation_alloc(vector->size);
					s = 0;
					gsl_linalg_LU_decomp(covarianceBetamatrixSHORT,perm, &s);
					gsl_linalg_LU_invert(covarianceBetamatrixSHORT,perm, weightBetamatrixSHORT);
					gsl_matrix_free(covarianceBetamatrixSHORT); covarianceBetamatrixSHORT = 0;
					gsl_permutation_free(perm); perm = 0;
					
					gsl_matrix_memcpy(Wm_short,weightAlphamatrixSHORT);
					gsl_matrix_add(Wm_short,weightBetamatrixSHORT);
					gsl_matrix_scale(Wm_short,1.0/2.0);						// (Walpha + Wbeta)/2
					gsl_matrix_free(weightAlphamatrixSHORT); weightAlphamatrixSHORT = 0;
					gsl_matrix_free(weightBetamatrixSHORT); weightBetamatrixSHORT = 0;
				}
			      
				// It is not necessary to check the allocation because 'vector' size must already be > 0
				gsl_matrix *X = gsl_matrix_alloc(vector->size,2);
				gsl_matrix_set_all(X,1.0);
				gsl_matrix_set_col(X,0,Dab_short);                                              //    | x0 1 | | .  1|
														// X =| x1 1 |=|Dab 1|
														//    | .    | | .  1|
														//    | xm 1 | | .  1|

				gsl_matrix *X_transW = gsl_matrix_alloc(2,vector->size);
				if (X->size1 != Wm_short->size1)        // Because the next operation is X'W
				{
					sprintf(valERROR,"%d",__LINE__+5);
					string str(valERROR);
					message = "Wrong dimensions to compute matrix-matrix product in line " + str + " (" + __FILE__ + ")";
					EP_EXIT_ERROR(message,EPFAIL);
				}
				gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,X,Wm_short,0.0,X_transW);            // X_transW = X'W
				gsl_matrix *X_transWX = gsl_matrix_alloc(2,2);
				gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,X_transW,X,0.0,X_transWX);         // X_transWX = X'WX
				gsl_vector *X_transWY = gsl_vector_alloc(2);
				gsl_blas_dgemv(CblasNoTrans,1.0,X_transW,P_Pab,0.0,X_transWY);                  // X_transWY = X'WY
														// Y = P_Pab
			      
				perm = gsl_permutation_alloc(2);
				s=0;
				gsl_matrix *aux = gsl_matrix_alloc(2,2);
				gsl_matrix *inv = gsl_matrix_alloc(2,2);
				gsl_matrix_memcpy(aux,X_transWX);
				gsl_linalg_LU_decomp(aux, perm, &s);
				if (gsl_linalg_LU_invert(aux, perm, inv) != 0) 
				{
					sprintf(valERROR,"%d",__LINE__-2);
					string str(valERROR);
					message = "Singular matrix in line " + str + " (" + __FILE__ + ")";
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				gsl_matrix_free(aux); aux = 0;
				gsl_permutation_free(perm); perm = 0;
			      
				gsl_vector *EB = gsl_vector_alloc(2);
				gsl_blas_dgemv(CblasNoTrans,1.0,inv,X_transWY,0.0,EB);                          // |E|=((X'WX)^(-1))X'WY
														// |B|
				*calculatedEnergy = gsl_vector_get(EB,0);
				// It is calculating also the baseline, gsl_vector_get(EB,1)
				//cout<<"baseline: "<<gsl_vector_get(EB,1)<<endl;
			      
				gsl_matrix_free(X); X = 0;
				gsl_matrix_free(X_transW); X_transW = 0;
				gsl_matrix_free(X_transWX); X_transWX = 0;
				gsl_vector_free(X_transWY); X_transWY = 0;
				gsl_matrix_free(inv); inv = 0;
				gsl_vector_free(Dab); Dab = 0;
				gsl_vector_free(Pab); Pab = 0;
				gsl_vector_free(P_Pab); P_Pab = 0;
				gsl_vector_free(Dab_short); Dab_short = 0;
				gsl_matrix_free(Wm_short); Wm_short = 0;
				gsl_vector_free(EB); EB = 0;
			}
			else if (reconstruct_init->OFLib == 1)
			{
				// It is not necessary to check the allocation because 'vector' size must already be > 0 and 'reconstruct_init->pulse_length'=PulseLength(input parameter) has been checked previously
				gsl_vector *P_Pab = gsl_vector_alloc(vector->size);
				
				gsl_vector_memcpy(P_Pab,vector);      
				// It is not necessary to check the allocation because 'vector' size must already be > 0
				gsl_vector *Pab_short = gsl_vector_alloc(vector->size);
				gsl_vector_view temp;
				temp = gsl_vector_subvector(Pab,0,vector->size);
				gsl_vector_memcpy(Pab_short,&temp.vector);
				gsl_vector_sub(P_Pab,Pab_short);                                                // P-Pab
														// Pulse WITH baseline
														// Templates (PAB) WITHOUT baseline (subtracted in 'initializeReconstructionSIRENA')
				gsl_vector_free(Pab_short); Pab_short = 0;
			      
				gsl_vector *EB = gsl_vector_alloc(2);
			
				gsl_blas_dgemv(CblasNoTrans,1.0,PRCLWN,P_Pab,0.0,EB);                           // [(R'WR)^(-1)]R'W\B7D'
														// D' = P_Pab
				*calculatedEnergy = gsl_vector_get(EB,0);
			
				gsl_vector_free(P_Pab); P_Pab = 0;
				gsl_vector_free(EB); EB = 0;
			}
		}
	//}
        }
	//cout<<"calculatedEnergy="<<*calculatedEnergy<<endl;
	//std::cout << std::setprecision(9) << *calculatedEnergy << '\n';

	return EPOK;
}
/*xxxx end of SECTION B11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B12 ************************************************************
* writeFilterHDU: This function runs in RECONSTRUCTION mode and writes the optimal filter info (in the FILTER HDU) for each pulse
*                 if 'intermediate'=1.
*
* - Declare variables
* - Open intermediate FITS file
* - Create the FILTER HDU if it is the first pulse
* - Write data
* 	- OPTIMALF column
* 	- OFLENGTH column
* 	- FREQ column
* 	- OPTIMALFF column
* 	- ENERGY column
* - Modify the Primary HDU
* - Close intermediate output FITS file if it is necessary
* - Free memory
*
* Parameters:
* - reconstruct_init: Member of 'ReconstructInitSIRENA' structure to initialize the reconstruction parameters (pointer and values).
*                     In particular, this function uses 'detectFile', 'clobber' and 'pulse_length'
* - pulse_index: Index of the pulse whose info is going to be written (to know if it is the first pulse)
* - energy: Estimated energy (eV)
* - optimalfilter: Optimal filter in time domain
* - optimalfilter_f: Frequency axis of the optimal filter spectrum
* - optimalfilter_FFT: Optimal filter spectrum
* - dtcObject: FITS object for intermeadiate file name
******************************************************************************/
int writeFilterHDU(ReconstructInitSIRENA **reconstruct_init, int pulse_index, double energy, gsl_vector *optimalfilter, gsl_vector *optimalfilter_f, gsl_vector *optimalfilter_FFT, fitsfile **dtcObject)
{
	// Declare variables
	string message = "";
	int status = EPOK;

	long totalpulses = 0;

	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';

	char *tt[1];
	char *tf[1];
	char *tu[1];
	char extname[20];
	char keyname[10];
	char keyvalstr[1000];

	// Open intermediate FITS file
	if (fits_open_file(dtcObject,dtcName,READWRITE,&status))
	{
		message = "Cannot open output intermediate file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	
	// Create the FILTER HDU is if it is the first pulse
	if ( ((*reconstruct_init)->clobber == 1) && (pulse_index == 0))
	{
		strcpy(extname,"FILTER");
		if (fits_create_tbl(*dtcObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
			message = "Cannot create table " + string(extname) + " in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
	}

	strcpy(extname,"FILTER");
	if (fits_movnam_hdu(*dtcObject, ANY_HDU,extname, 0, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in output intermediate file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (fits_get_num_rows(*dtcObject,&totalpulses, &status))
	{
		message = "Cannot get number of rows in " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// It is not necessary to check the allocation because '(*reconstruct_init)->pulse_length'=PulseLength(input parameter) has been checked previously
	gsl_matrix *optimalfilter_matrix = gsl_matrix_alloc(1,(*reconstruct_init)->pulse_length);
	gsl_matrix *optimalfilter_f_matrix = gsl_matrix_alloc(1,(*reconstruct_init)->pulse_length);
	gsl_matrix *optimalfilter_FFT_matrix = gsl_matrix_alloc(1,(*reconstruct_init)->pulse_length);

	gsl_matrix_set_zero(optimalfilter_matrix);
	gsl_matrix_set_zero(optimalfilter_f_matrix);
	gsl_matrix_set_zero(optimalfilter_FFT_matrix);

	for (int i=0;i<optimalfilter->size;i++)
	{
		gsl_matrix_set(optimalfilter_matrix,0,i,gsl_vector_get(optimalfilter,i));
		gsl_matrix_set(optimalfilter_f_matrix,0,i,gsl_vector_get(optimalfilter_f,i));
		gsl_matrix_set(optimalfilter_FFT_matrix,0,i,gsl_vector_get(optimalfilter_FFT,i));
	}

	// Write data
	IOData obj;
	obj.inObject = *dtcObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"FILTER");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.type = TDOUBLE;
	obj.unit = new char [255];

	// OPTIMALF column
	obj.iniRow = totalpulses+1;
	obj.endRow = totalpulses+1;
	strcpy(obj.nameCol,"OPTIMALF");
	strcpy(obj.unit," ");
	if (writeFitsComplex (obj,optimalfilter_matrix))
	{
		message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_matrix_free(optimalfilter_matrix); optimalfilter_matrix = 0;

	// OFLENGTH column
	gsl_vector *oflength = gsl_vector_alloc(1);
	gsl_vector_set(oflength,0,optimalfilter->size);
	strcpy(obj.nameCol,"OFLENGTH");
	if (writeFitsSimple (obj,oflength))
	{
		message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) +
			    " column in FILTER";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_vector_free(oflength); oflength = 0;

	// FREQ column
	strcpy(obj.nameCol,"FREQ");
	strcpy(obj.unit,"Hz");
	if (writeFitsComplex (obj,optimalfilter_f_matrix))
	{
		message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_matrix_free(optimalfilter_f_matrix); optimalfilter_f_matrix = 0;

	// OPTIMALFF column
	strcpy(obj.nameCol,"OPTIMALFF");
	strcpy(obj.unit," ");
	if (writeFitsComplex (obj,optimalfilter_FFT_matrix))
	{
		message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_matrix_free(optimalfilter_FFT_matrix); optimalfilter_FFT_matrix = 0;

	strcpy(obj.nameTable,"PULSES");

	// ENERGY column
	gsl_vector *energygsl = gsl_vector_alloc(1);
	gsl_vector_set(energygsl,0,energy);
	gsl_vector_scale(energygsl,1.0/1e3);
	strcpy(obj.nameCol,"ENERGY");
	strcpy(obj.unit,"keV");
	if (writeFitsSimple (obj,energygsl))
	{
		message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column in PULSES";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_vector_free(energygsl); energygsl = 0;

	// Close intermediate output FITS file if it is necessary
	if (fits_close_file(*dtcObject,&status))
	{
		message = "Cannot close file " + string(dtcName);
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
        *dtcObject = 0;

	(*reconstruct_init)->clobber = 2;
	
	// Free memory
	delete [] obj.nameTable; obj.nameTable = 0;
	delete [] obj.nameCol; obj.nameCol = 0;
	delete [] obj.unit; obj.unit = 0;

	return(EPOK);
}
/*xxxx end of SECTION B12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
