/************************************************************************************************

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

   Copyright 2014:  gennoisespec has been developed by the INSTITUTO DE FISICA DE
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01 and
   ESP2013-48637-C2-1-P.

/************************************************************************************************
*                      GENNOISESPEC
 *
 *
 *
 *  File:      gennoisespec.cpp
 *  Version:   10.0.0
 *  Developer: Beatriz Cobo Martín
 *             cobo@ifca.unican.es
 *             IFCA
 *
 *  Revision History:
 *
 *  version 1.0.0: 	02/06/09    First version
 *  version 1.0.1:  02/11/09    Improvement of the error handling
 *  version 1.1.1:  02/16/09    Added EVENTCNT and FTYPE output keywords
 *  version 2.0.0:  02/18/09    Changing input (and output) keywords
 *                              Deleting some non necessary inpu parameters
 *  version 3.0.0:  02/23/09    IV input enabled
 *  version 4.0.0:  03/02/09    Using IVCAL and VVCAL
 *  version 5.0.0   03/03/09    Used the keywords libraries
 *  version 6.0.0   04/24/09	TRG_ID instead ZC_ID in getKh
 *                              Free dynamic memory
 *-------------------------------------------------------------------------------------------------------------------
 * To be according to the repository version
 *-------------------------------------------------------------------------------------------------------------------
 *  version 5.0.0   07/01/09    Changed functions names to avoid problems with new GSL functions
 *  version 6.0.0	13/11/09	If error file exists, its contests will be deleted
 *                              Different behaviour when the input FITS file type is "TESNOISE"
 *                              (it is not necessary to handle a _trg file) => Added findIntervalN function
 *  version 7.0.0   04/02/10    gsl_vector_free(EventSamples_aux) line changed
 *                  15/02/10    intervalMinBins instead eventsz to calculate the frequency step
 *  version 8.0.0   17/03/10    Factor 1/sqrt(N), in the FFT expression
 *                  20/10/10    File FITS resulting from TRIGGER is not used (not read):
 *            						getKhaux function deleted
 *             					Input parameters updated (initModule)
 *             					TESNOISE files can also have pulses
 *             						- No pulses => The whole event is used (without dividing into intervals of intervalMin size)
 *             		            finInterval and findIntervalN functions changed
 *                  10/01/11    Deleted functions to handle vectors, complex vectors and matrix
 *                              Added input keyword ASQUID
 *                              ASQUID is used from now on to assign positive polarity (instead bl)
 *                              Deleted input parameter nbins (only necessary to calculate the baseline)
 *                  01/03/11    PROCESS also includes gennoisespec version (CREATOR)
 *  version 9.0.0   24/03/11    "status" allocation and "writeLog" order changed in some cases
 *                  25/03/11    Quick Look mode no longer used
 *                  05/04/11    Pulses are not searched for if ENERGY = 0
 *                  26/05/11    New input parameters: samplesUp and nSgms
 *                              Deleted input parameter n (and n_b)
 *                              New algorithm to detect pulses: 1st derivative and threshold (findMeanSigma, medianKappaClipping and new findTstart)
 *                  25/08/11    Added input keyword PLSPOLAR
 *                              ASQUID and PLSPOLAR used to assign positive polarity to the pulses
 *                              If XRAY or IV => It looks for:
 *                              	- Single pulses (findTstart)
 *                                  - Secondary pulses (by using the pulses models library and the composed first derivative)
 *                                  - Primary pulses (findTstartPrimary)
 *                              Added 'readLib' and 'inDataIteratorLib'
 *                              'gsl_vector_view temp;' added to handle with subvectors
 *                              Code errors changed: 71xx into 710xx
 *                  27/10/11    Different input/output parameters in findTstart (tmaxDER)
 *						        Pulse model anticipated pulse and pulse maximum & next sample differ less than 5%
 *						        Composed first derivative is fixed as 0 in tstartDERgsl+tmaxDERgsl+"20"
 *							    Secondary pulse larger than 1/100 times the preceding pulse
 *					13/02/13    New input parameter scaleFactor
 *					            New parameters to handle with code hardpoints:
 *					            	stopCriteriaMKC
 *               					kappaMKC
 *               					limitCriteriaMAX
 *               					nsAftrtstart
 *               					levelPrvPulse
 *               					primaryThresholdCriteria
 *               	14/02/13   Functions 'readLib' and 'inDataIteratorLib' updated
 *                             New input parameters: b_cF, c_cF, LrsT, LbT
 *                  20/02/13   findPulses of Utils used
 *                  26/11/13   'readLib' and 'inDataIteratorLib' updated in order to not take into account the PRETRIGS
 *                             column in the library FITS file
 *                             tAftrtstart has changed from 'constant' to 'input parameter'
 *                  11/07/14   PIL/RIL/Common dependencies removed. New functions used to read input parameters
 *                  09/09/14   ADC2016 instead RAW_ADC_DATA
 *				               PXL02016 instead I0 column
 *                             Differences when reading or writing some keywords
 *                  10/09/14   50 null initial samples are added to the pulse templates read from the library ('models')
 *                  16/09/14   'intervalMinBins = floor(intervalMin*samprate);' instead 'intervalMinBins = ceil(intervalMin*samprate);'
 *                             (it produced normalizationFactor=inf in FILTER)
 *                  18/09/14   RECORDS instead ADC2016
*                              ADC instead PXL02016 column
*  version 10.0.0   13/01/15   Migrated to CFITSIO (removal of ISDC DAL)
*                              Run dependency on datatype files (xray, iv or tesnoise) deleted
 ********************************************************************************************************************/

/******************************************************************************
DESCRIPTION:

The goal of the GENNOISE task is to calculate the current noise spectral density.

 * If there are pulses in an event, the pulses are rejected and it is going to look for pulse-free intervals of a given size (intervalMin)
 * If there are no pulses in an event, the event is divided into pulse-free intervals of a given size (intervalMin)

The input data are low-pass filtered in order to look for pulses. Therefore, it is going to look for pulse-free intervals, calculate their FFT
(not filtered data) and average them.

The user must supply the following input parameters:

- inFile: Name of the input FITS file
- outFile: Name of the output FITS file
- intervalMin: Minimum length of a pulse-free interval to use (seconds)
- ntausPF: Number of tauFALLs after ending the pulse (Tend) to start the pulse-free interval
- nintervals: Number of pulse-free intervals to use
- tauFALL: Fall time of the pulses in seconds
- scaleFactor: Scale factor to apply to the fall time of the pulses in order to calculate the LPF box-car length
- samplesUp: Consecutive samples over the threshold to locate a pulse
- nSgms: Number of Sigmas to establish the threshold
- ntaus: Number of tauFALLs after starting the pulse to calculate its tend
- LrsT: Running sum length in seconds (only in notcreationlib mode)
- LbT: Baseline averaging length in seconds (only in notcreationlib mode)
- namelog: Output log file name
- verbosity: Verbosity level of the output log file

The output FITS file (_noisespec) contains three columns in two extensions, NOISE and NOISEALL:
 	- Frequency
 	- Current noise spectral density: Amount of current per unit (density) of frequency (spectral), as a function of the frequency
 	- Standard error of the mean

The NOISEALL extension contains positive and negative frequencies.

 MAP OF SECTIONS IN THIS FILE:

  - 1.- INCLUDE's
  - 2.- MAIN
  - 3.- initModule
  - 4.- inDataIterator
  - 5.- findInterval
  - 6.- findIntervalN
  - 7.- createTPSreprExten
  - 8.- writeTPSreprExten

 *******************************************************************************/

/***** SECTION 1 ************************************
*       INCLUDE's
****************************************************/
#include <gennoisespec.h>

/*xxxx end of SECTION 1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 2 ************************************
* MAIN function: This function is the main function of the GENNOISESPEC task
*
* - Read input parameters
* - Open input FITS file
* - Read and check input keywords
* - Get structure of input FITS file columns
* - Initialize variables and transform from seconds to samples
* - Pulse-free segments are divided into pulse-free intervals with intervalMinBins size
* - Every single spectrum of each pulse-free interval is stored in a row of the EventSamplesFFT array
* - Create structure to run Iteration: inDataIterator
* 	- Read columns (TIME and ADC)
* - Called iteration function: inDataIterator
* - Close input FITS file
* - Generate NOISE representation
*   - Current noise spectral density
*   - Mean
* 	- Standard error of the mean
*  - Create output FITS File: GENNOISESPEC representation file (*_noisespec.fits)
*  - Write extensions NOISE and NOISEALL (call writeTPSreprExten)
*  - Free allocate of GSL vectors
*  - Close output FITS file
*  - Free memory
*  - Finalize the task
*****************************************************/
int main (int argc, char **argv)
{
	char val[256];

	creator = "gennoisespec v.10.0.0";             //Set "CREATOR" keyword of the output FITS file
	time_t t_start = time(NULL);
	char str_stat[8];
	string message = "";
	int status=EPOK, extver=0;

	sprintf(temporalFileName,"GENNOISESPECauxfile");
	strcat(temporalFileName,".txt");
	temporalFile = fopen (temporalFileName,"w");
	if (temporalFile == NULL)
	{
		message = "Cannot open auxiliary file GENNOISESPECauxfile.txt";
		writeLog(fileRef,"Error",verbosity,message);
		EP_EXIT_ERROR(message,EPFAIL);
	}

	// Read input parameters
	if (initModule(argc, argv))
	{
	    message = "Cannot run initModule routine to get input parameters";
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	writeLog(fileRef,"Log", verbosity,"Into GENNOISESPEC task");

	// Open input FITS file
	if (fits_open_file(&infileObject, infileName,0,&status))
	{
	    message = "Cannot open file " + string(infileName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(extname,"RECORDS");
	if (fits_movnam_hdu(infileObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname) + " in " + string(infileName);
	    EP_EXIT_ERROR(message,status);
	}

	// Read and check input keywords
	strcpy(extname,"RECORDS");
	if (fits_get_num_rows(infileObject,&eventcnt, &status))
	{
	    message = "Cannot get number of rows in HDU " + string(extname);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"TRIGGSZ");
	if (fits_read_key(infileObject,TLONG,keyname, &eventsz,comment,&status))
	{
	    message = "Cannot read keyword " + string(keyname) + " in input file";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (eventsz <= 0)
	{
		message = "Legal values for TRIGGSZ (RECORDS) are integer numbers greater than 0";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	strcpy(keyname,"DELTAT");
	if (fits_read_key(infileObject,TDOUBLE,keyname, &samprate,comment,&status))
	{
	    message = "Cannot read keyword " + string(keyname) + " in input file";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (samprate <= 0)
	{
		message = "Legal values for DELTAT (RECORDS) are real numbers greater than 0";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	samprate = 1/samprate;
	ivcal=1.0;
	asquid = 1.0;
	strcpy(keyname,"MONOEN");
	if (fits_read_key(infileObject,TDOUBLE,keyname, &energy,comment,&status))
	{
	    message = "Cannot read keyword " + string(keyname) + " in input file";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (energy < 0)
	{
		message = "Legal values for MONOEN (RECORDS) are non negative real numbers";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	energy = energy*1e3;
	plspolar = 1.0;

	// Get structure of input FITS file columns
	strcpy(straux,"Time");
	if (fits_get_colnum(infileObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(infileName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"ADC");
	if (fits_get_colnum(infileObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(infileName);
	    EP_EXIT_ERROR(message,status);
	}

	message = "Open infileFits: " + string(infileName);
	writeLog(fileRef,"Error", verbosity,message);

	// Initialize variables and transform from seconds to samples
	tauFALL_b = (int) (tauFALL * samprate);
	safetyMarginTstart = safetyMarginTstart*samprate;
	Lrs = LrsT * samprate;
	Lb = LbT * samprate;

	// Pulse-free segments are divided into pulse-free intervals with intervalMinBins size
	intervalMinBins = floor(intervalMin*samprate);		//A way of casting from double to int
	if (intervalMinBins > eventsz)
	{
		message = "Illegal value in INTERVALMIN parameter. Legal values reals greater than 0 and fewer than event length";
		writeLog(fileRef,"Error", verbosity,message);
		EP_EXIT_ERROR(message,EPFAIL);
	}

	// Every single spectrum of each pulse-free interval is stored in a row of the EventSamplesFFT array
	EventSamplesFFT = gsl_matrix_alloc(nintervals,intervalMinBins);
	// Initialize to zero => Because EventSamplesFFT is allocated with maximum dimensions
	// but maybe there are not nintervals pulse-free intervals in the data
	gsl_matrix_set_zero(EventSamplesFFT);
	EventSamplesFFTMean = gsl_vector_alloc(intervalMinBins);
	gsl_vector_set_zero(EventSamplesFFTMean);
	gsl_vector *mean;
	mean = gsl_vector_alloc(intervalMinBins);
	gsl_vector_set_zero(mean);
	sigmacsdgsl = gsl_vector_alloc(intervalMinBins);
	gsl_vector_set_zero(sigmacsdgsl);

	// Create structure to run Iteration
	extern int inDataIterator(long totalrows,long offset,long firstrow,long nrows,int ncols,iteratorCol *cols,void *user_strct );

	long rows_per_loop = 0;					// 0: Use default: Optimum number of rows
	long offset = 0;						// 0: Process all the rows
	iteratorCol cols [2];					// Structure of Iteration
	int n_cols = 2; 						// Number of columns:  Time + ADC

	strcpy(extname,"RECORDS");
	if (fits_movnam_hdu(infileObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname) +" in " + string(infileName);
	    EP_EXIT_ERROR(message,status);
	}

	// Read Time column
	strcpy(straux,"TIME");
	status = fits_iter_set_by_name(&cols[0], infileObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(infileName);
	    EP_EXIT_ERROR(message,status);
	}
	
	// Read ADC column
	strcpy(straux,"ADC");
	status = fits_iter_set_by_name(&cols[1], infileObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(infileName);
	    EP_EXIT_ERROR(message,status);
	}

	// Called iteration function
	if (fits_iterate_data(n_cols,cols,offset,rows_per_loop,inDataIterator,0L,&status))
	{
	    message = "Cannot iterate data using InDataTterator";
	    EP_EXIT_ERROR(message,status);
	}

	// Close input FITS file
	if (fits_close_file(infileObject,&status))
	{
	    message = "Cannot close file " + string(infileName);
	    EP_EXIT_ERROR(message,status);
	}

	// Generate NOISE representation
	if (NumMeanSamples == 0)
	{
		message = "Pulse-free intervals not found";
		writeLog(fileRef,"Error", verbosity,message);
		EP_EXIT_ERROR(message,status);
	}
	else if	(NumMeanSamples < nintervals)
	{
		sprintf(str_stat,"%d",NumMeanSamples);
		message = "Not enough pulse-free intervals. CSD calculated with " + string(str_stat);
		writeLog(fileRef,"Log", verbosity,message);
	}

	// sqrt(sum(FFT^2)/NumMeanSamples) => sqrt(A^2) = A and sqrt(1/NumMeanSamples)=1/sqrt(Hz)
	gsl_vector_scale(EventSamplesFFTMean,(1.0/(double)NumMeanSamples));
	for (int i=0;i<EventSamplesFFTMean->size;i++)
	{
		if (gsl_vector_get(EventSamplesFFTMean,i)<0)
		{
		  message = "EventSamplesFFTMean < 0 for i=" + boost::lexical_cast<std::string>(i);
		  EP_EXIT_ERROR(message,EPFAIL); 
		}
	}
	gsl_vector_sqrtIFCA(EventSamplesFFTMean,EventSamplesFFTMean);

	// Mean
	double value_aux;
	for (int i=0;i<intervalMinBins;i++)
	{
		value_aux = 0.0;
		for (int j=0;j<NumMeanSamples;j++ )
		{
			value_aux = value_aux + gsl_matrix_get(EventSamplesFFT,j,i);
		}
		gsl_vector_set(mean,i,value_aux);
	}
	gsl_vector_scale(mean,1.0/(double)NumMeanSamples);

	// Standard error of the mean
	for (int i=0;i<intervalMinBins;i++)
	{
		value_aux = 0.0;
		for (int j=0;j<NumMeanSamples;j++)
		{
			value_aux = value_aux+pow(gsl_matrix_get(EventSamplesFFT,j,i)-gsl_vector_get(mean,i),2.0);
		}
		gsl_vector_set(sigmacsdgsl,i,value_aux);
	}
	gsl_vector_scale(sigmacsdgsl,1.0/((double)NumMeanSamples*((double)NumMeanSamples-1)));
	gsl_vector_sqrtIFCA(sigmacsdgsl,sigmacsdgsl);

	// Create output FITS File: GENNOISESPEC representation file (*_noisespec.fits)
	if(createTPSreprFile())
	{
		message = "Cannot create file " +  string(tespsName);
		EP_EXIT_ERROR(message,EPFAIL);
	}	
	if (fits_open_file(&tespsObject,tespsName,1,&status))
	{
		message = "Cannot open file " +  string(tespsName);
		EP_EXIT_ERROR(message,status);
	}
	
	// Write extensions NOISE and NOISEALL (call writeTPSreprExten)
	if(writeTPSreprExten ())
	{
		message = "Cannot write extensions in " +  string(tespsName);
		EP_EXIT_ERROR(message,EPFAIL);
	}

	// Free allocate of GSL vectors
	gsl_matrix_free(EventSamplesFFT);
	gsl_vector_free(EventSamplesFFTMean);
	gsl_vector_free(mean);
	gsl_vector_free(sigmacsdgsl);
	gsl_vector_free(startIntervalgsl);

	// Close output FITS file
	if (fits_close_file(tespsObject,&status))
	{
	    message = "Cannot close file " + string(tespsName);
	    EP_EXIT_ERROR(message,status);
	}	

	// Free memory
	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;

	fclose(temporalFile);

	// Finalize the task
	time_t t_end = time(NULL);
	sprintf(straux,"%f",(double) (t_end - t_start));
	message = "Time:" + string(straux);
	writeLog(fileRef,"Log", verbosity,message);

	writeLog(fileRef,"OK", verbosity,"Gennoisespec Module OK");

	if(fclose(fileRef))
	{
	    message = "Cannot close log file";
	    EP_EXIT_ERROR(message,EPFAIL);
	}	
	
	return EPOK;
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 ************************************************************
* initModule function: This function reads and processes the input parameters
*
* - It defines the input parameters and assigns the values to variables
* - It defines the structure for the command line options
* - It gets the parameter values from the command line or user input
* - It saves the parameter values into meaningful variables
*   (checking also if the parameter value is in the allowed range)
* - It opens the log file of the task
****************************************************************************/
int initModule(int argc, char **argv)
{
	int status = EPOK;

	// Define TRIGGER input parameters and assign values to variables
	// Parameter definition and assignation of default values
	const int npars = 14, npars1 = 15;
	inparam gennoisespecPars[npars];
	int optidx =0, par=0, fst=0, ipar;
	string message="";
	string task="gennoisespec";

	gennoisespecPars[0].name = "inFile";
	gennoisespecPars[0].description = "Enter input FITS file name";
	gennoisespecPars[0].defValStr = "a.fits";
	gennoisespecPars[0].type =  "char";
	gennoisespecPars[0].ValStr = gennoisespecPars[0].defValStr;

	gennoisespecPars[1].name = "outFile";
	gennoisespecPars[1].description = "Enter _noisespec output file name";
	gennoisespecPars[1].defValStr = "a_noisespec.fits";
	gennoisespecPars[1].type =  "char";
	gennoisespecPars[1].ValStr = gennoisespecPars[1].defValStr;

	gennoisespecPars[2].name = "intervalMin";
	gennoisespecPars[2].description = "Minimum length of a pulse-free interval (seconds)";
	gennoisespecPars[2].defValReal = 0.0025;
	gennoisespecPars[2].type = "double";
	gennoisespecPars[2].minValReal = 1.E-50;
	gennoisespecPars[2].maxValReal = +1.E+50;
	gennoisespecPars[2].ValReal = gennoisespecPars[2].defValReal;

	gennoisespecPars[3].name = "ntausPF";
	gennoisespecPars[3].description = "Number of tauFALLs after ending the pulse to start the pulse-free interval";
	gennoisespecPars[3].defValInt = 0;
	gennoisespecPars[3].type = "int";
	gennoisespecPars[3].minValInt = 0;
	gennoisespecPars[3].maxValInt = 1000;
	gennoisespecPars[3].ValInt = gennoisespecPars[3].defValInt;

	gennoisespecPars[4].name = "nintervals";
	gennoisespecPars[4].description = "Number of pulse-free intervals to use";
	gennoisespecPars[4].defValInt = 60;
	gennoisespecPars[4].type = "int";
	gennoisespecPars[4].minValInt = 1;
	gennoisespecPars[4].maxValInt = 1E4;
	gennoisespecPars[4].ValInt = gennoisespecPars[4].defValInt;

	gennoisespecPars[5].name = "tauFALL";
	gennoisespecPars[5].description = "Fall time of the pulses (seconds)";
	gennoisespecPars[5].defValReal = 3.5E-4;
	gennoisespecPars[5].type = "double";
	gennoisespecPars[5].minValReal = 1.E-50;
	gennoisespecPars[5].maxValReal = 1.E+50;
	gennoisespecPars[5].ValReal = gennoisespecPars[5].defValReal;

	gennoisespecPars[6].name = "scaleFactor";
	gennoisespecPars[6].description = "Scale factor to apply to the fall time of the pulses in order to calculate the LPF box-car length";
	gennoisespecPars[6].defValReal = 0.1;
	gennoisespecPars[6].type = "double";
	gennoisespecPars[6].minValReal = 1.E-50;
	gennoisespecPars[6].maxValReal = 1.E+50;
	gennoisespecPars[6].ValReal = gennoisespecPars[6].defValReal;

	gennoisespecPars[7].name = "samplesUp";
	gennoisespecPars[7].description = "Consecutive samples over the threshold to locate a pulse";
	gennoisespecPars[7].defValInt = 20;
	gennoisespecPars[7].type = "int";
	gennoisespecPars[7].minValInt = 1;
	gennoisespecPars[7].maxValInt = 1E4;
	gennoisespecPars[7].ValInt = gennoisespecPars[7].defValInt;

	gennoisespecPars[8].name = "nSgms";
	gennoisespecPars[8].description = "Number of Sigmas to establish the threshold";
	gennoisespecPars[8].defValInt = 3;
	gennoisespecPars[8].type = "int";
	gennoisespecPars[8].minValInt = 1;
	gennoisespecPars[8].maxValInt = 100;
	gennoisespecPars[8].ValInt = gennoisespecPars[8].defValInt;

	gennoisespecPars[9].name = "ntaus";
	gennoisespecPars[9].description = "Number of tauFALLs after starting the pulse to calculate its tend";
	gennoisespecPars[9].defValInt = 15;
	gennoisespecPars[9].type = "int";
	gennoisespecPars[9].minValInt = 1;
	gennoisespecPars[9].maxValInt = 100;
	gennoisespecPars[9].ValInt = gennoisespecPars[9].defValInt;

	gennoisespecPars[10].name = "LrsT";
	gennoisespecPars[10].description = "Running sum length (in the RS filter case) (seconds)";
	gennoisespecPars[10].defValReal = 30.E-6;
	gennoisespecPars[10].type = "double";
	gennoisespecPars[10].minValReal = 1.E-50;
	gennoisespecPars[10].maxValReal = 1.E+50;
	gennoisespecPars[10].ValReal = gennoisespecPars[10].defValReal;

	gennoisespecPars[11].name = "LbT";
	gennoisespecPars[11].description = "Baseline averaging length (in the RS filter case) (seconds)";
	gennoisespecPars[11].defValReal = 1E-3;
	gennoisespecPars[11].type = "double";
	gennoisespecPars[11].minValReal = 1.E-50;
	gennoisespecPars[11].maxValReal = 1.E+50;
	gennoisespecPars[11].ValReal = gennoisespecPars[11].defValReal;

	gennoisespecPars[12].name = "nameLog";
	gennoisespecPars[12].description = "Output log file name";
	gennoisespecPars[12].defValStr = "noise_log.txt";
	gennoisespecPars[12].type = "char";
	gennoisespecPars[12].ValStr = gennoisespecPars[12].defValStr;

	gennoisespecPars[13].name = "verbosity";
	gennoisespecPars[13].description = "Verbosity level of the output log file (in [0,3])";
	gennoisespecPars[13].defValInt = 3;
	gennoisespecPars[13].type = "int";
	gennoisespecPars[13].minValInt = 0;
	gennoisespecPars[13].maxValInt = 3;
	gennoisespecPars[13].ValInt = gennoisespecPars[13].defValInt;

	// Define structure for command line options
	static struct option long_options[npars1];
	for (optidx = 0; optidx<npars; optidx++)
	{
		long_options[optidx].name= gennoisespecPars[optidx].name.c_str();
		long_options[optidx].has_arg = 2;
		long_options[optidx].flag = 0;
		long_options[optidx].val = 0;
	}
	long_options[npars].name=0;
	long_options[npars].has_arg = 0;
	long_options[npars].flag = 0;
	long_options[npars].val = 0;

	// Get parameter values from command line or user input
	optidx = 0;
	int commandLine = 0; // initialize command line as being empty
	while ((par=getopt_long(argc, argv,"",long_options, &optidx ))!= -1) // while there are still parameters in command line
	{
		// Params have been specified in the command line
		commandLine = 1;
		switch(par)
		{
			case 0:
				// identify parameter read from command line
				for(int i=0;i<npars; i++)
				{
					if(long_options[optidx].name == gennoisespecPars[i].name.c_str())
					{
					    if (gennoisespecPars[i].type == "char") //save char value for par
					    {
	    					gennoisespecPars[i].ValStr = optarg;
					    }
					    else // check if numeric value
					    {
	    					if ((!isdigit(optarg[0]) && (optarg[0] != '-')) ||
	    							(!isdigit(optarg[0]) && (optarg[0] == '-') && (!isdigit(optarg[1]))))
	    					{
	    						message = "Invalid value for input argument " + string(gennoisespecPars[i].name);
	    						EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	    					}
	    					if (gennoisespecPars[i].type == "int")
	    					{
	    						gennoisespecPars[i].ValInt = atoi(optarg);
	    					}
	    					else
	    					{
	    						gennoisespecPars[i].ValReal= atof(optarg);
	    					}
					     }
					     break;
					} // endif
				} // endfor
				break;
		    default:
		    	message = "Invalid parameter name " + string(long_options[optidx].name);
	    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}//switch
	}//while
	// If command line is empty: ask for params interactively
	if(commandLine == 0)
	{
		if(interactivePars(gennoisespecPars,npars,task))
		{
		    message = "Error reading parameters interactively";
		    EP_PRINT_ERROR(message,EPFAIL);
		}
	}

	// Save parameter values into meaningful variables
	for(int i=0;i<npars; i++)
	{
		if(gennoisespecPars[i].name == "inFile")
		{
		    strcpy(infileName, gennoisespecPars[i].ValStr.c_str());
		}
		else if(gennoisespecPars[i].name == "outFile")
		{
		    strcpy(tespsName, gennoisespecPars[i].ValStr.c_str());
		}
		else if(gennoisespecPars[i].name == "intervalMin")
		{
		    intervalMin = gennoisespecPars[i].ValReal;
		}
		else if(gennoisespecPars[i].name == "ntausPF")
		{
		    ntausPF = gennoisespecPars[i].ValInt;
		}
		else if(gennoisespecPars[i].name == "nintervals")
		{
		    nintervals = gennoisespecPars[i].ValInt;
		}
		else if(gennoisespecPars[i].name == "tauFALL")
		{
		    tauFALL = gennoisespecPars[i].ValReal;
		}
		else if(gennoisespecPars[i].name == "scaleFactor")
		{
		    scaleFactor = gennoisespecPars[i].ValReal;
		}
		else if(gennoisespecPars[i].name == "samplesUp")
		{
		    samplesUp = gennoisespecPars[i].ValInt;
		}
		else if(gennoisespecPars[i].name == "nSgms")
		{
		    nSgms = gennoisespecPars[i].ValInt;
		}
		else if(gennoisespecPars[i].name == "ntaus")
		{
		    ntaus = gennoisespecPars[i].ValInt;
		}
		else if(gennoisespecPars[i].name == "LrsT")
		{
		    LrsT = gennoisespecPars[i].ValReal;
		}
		else if(gennoisespecPars[i].name == "LbT")
		{
		    LbT = gennoisespecPars[i].ValReal;
		}
		else if(gennoisespecPars[i].name == "nameLog")
		{
		    strcpy(nameLog,gennoisespecPars[i].ValStr.c_str());
		}
		else if(gennoisespecPars[i].name == "verbosity")
		{
		    verbosity = gennoisespecPars[i].ValInt;
		}
		
		// Check if parameter value is in allowed range
		if( gennoisespecPars[i].type == "int" &&
				(gennoisespecPars[i].ValInt < gennoisespecPars[i].minValInt ||
						gennoisespecPars[i].ValInt > gennoisespecPars[i].maxValInt))
		{
			message = "Parameter name " + string(long_options[optidx].name) + " out of range: [" +
					boost::lexical_cast<std::string>(gennoisespecPars[i].minValInt) + "," + boost::lexical_cast<std::string>(gennoisespecPars[i].maxValInt) + "]";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		else if ( gennoisespecPars[i].type == "double" &&
				(gennoisespecPars[i].ValReal < gennoisespecPars[i].minValReal ||
						gennoisespecPars[i].ValReal > gennoisespecPars[i].maxValReal))
		{
			message = "Parameter name " + string(long_options[optidx].name) + " out of range: [" +
					boost::lexical_cast<std::string>(gennoisespecPars[i].minValReal) + "," + boost::lexical_cast<std::string>(gennoisespecPars[i].maxValReal) + "]";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
	}// loop for parameters

	// Open the log file of the task
	fileRef = fopen(nameLog,"w+");	// Remove file if it already exists and open a new file to save log messages
	if (fileRef == NULL)
	{
	    message = "Cannot open Log file " + string(nameLog);
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	return (EPOK);
}
/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
* inDataIterator: This function takes the optimum number of rows to read the input FITS file
*                 and works iteratively
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each record
* 	- If there are pulses (MONOEN >0):
* 		- Assigning positive polarity (by using ASQUID and PLSPOLAR)
* 		- Low-pass filtering
*   	- Derivative after filtering
*   	- Finding the pulses: Pulses tstart are found (call findPulses)
*   - Finding the pulse-free intervals in each record read
*  		- If there are pulses => Call findInterval
*		- No pulses => The whole event is going to be used (DIVIDING into intervals of intervalMinBins size) => Call findIntervalN
*   - CSD calculus (not filtered data):
* 		- FFT calculus (each pulse-free interval)
* 		- Every single spectrum of each pulse-free interval is stored in a row of the EventSamplesFFT array
* 		- Add to mean FFT samples
*  - Free allocate of GSL vectors
****************************************************************************/
int inDataIterator(long totalrows, long offset, long firstrow, long nrows,	int ncols, iteratorCol *cols, void *user_strct)
{
	char val[256];

	string message = "";
	int status = EPOK;
	int extver=0;

	if (NumMeanSamples >= nintervals)	{return (EPOK);}

	// Declare variables
	// To read from the input FITS file
	double *time_block, *timein;	// Vector of TIME column
	double *iout_block, *ioutin;	// Vector of ADC column
	gsl_vector *timegsl_block;
	gsl_matrix *ioutgsl_block;

	// To derive
	gsl_vector *derSGN = gsl_vector_alloc(eventsz);

	gsl_vector *tstartgsl = gsl_vector_alloc(eventsz);
	gsl_vector *qualitygsl = gsl_vector_alloc(eventsz);
	gsl_vector *energygsl = gsl_vector_alloc(eventsz);
	int nPulses;
	gsl_vector *tstartDERgsl = gsl_vector_alloc(eventsz);
	gsl_vector *tmaxDERgsl = gsl_vector_alloc(eventsz);
	gsl_vector *maxDERgsl = gsl_vector_alloc(eventsz);
	gsl_vector *tendDERgsl = gsl_vector_alloc(eventsz);

	int pulseFound;	// 0->The function findTstart has not found any pulse
					// 1->The function findTstart has found at least one pulse

	// Auxiliary variables
	gsl_vector_view temp;					// In order to handle with gsl_vector_view (subvectors)

	// To calculate the FFT
	double SelectedTimeDuration;
	gsl_vector *EventSamples;
	gsl_vector *vector_aux;
	gsl_vector_complex *vector_aux1;
	vector_aux = gsl_vector_alloc(intervalMinBins);
	vector_aux1 = gsl_vector_complex_alloc(intervalMinBins);

	// Allocate input GSL vectors
	timegsl_block = gsl_vector_alloc(nrows);
	timegsl = gsl_vector_alloc(nrows); 						// GSL of input TIME column
	ioutgsl_block = gsl_matrix_alloc(nrows,eventsz);
	ioutgsl = gsl_vector_alloc(eventsz); 					// GSL of input ADC column
	gsl_vector *ioutgsl_aux = gsl_vector_alloc(eventsz);	//In case of working with ioutgsl without filtering
	gsl_vector *ioutgslNOTFIL = gsl_vector_alloc(eventsz);

	// Read iterator
	// NOTE: fits_iter_get_array because in this fits function the 1st element
	//  of the output array is the null pixel value!
	timein = (double *) fits_iter_get_array(&cols[0]);
	time_block = &timein[1];
	if (toGslVector(((void **)&time_block), &timegsl_block, nrows, 0, TDOUBLE))
	{
	    message = "Cannot run routine toGslVector for TIME vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	ioutin = (double *) fits_iter_get_array(&cols[1]);
	iout_block = &ioutin[1];
	if(toGslMatrix(((void **)&iout_block), &ioutgsl_block, eventsz, nrows, TDOUBLE,	0))
	{
	    message = "Cannot run routine toGslMatrix for IOUT ";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	
	// Pulse-free segments are divided into pulse-free intervals with intervalMinBins size
	SelectedTimeDuration = eventsz/((double)samprate);

	// Processing each record
	for (int i=0; i< nrows; i++)
	{
		if (NumMeanSamples>=nintervals)
		{
		    writeLog(fileRef,"Log", verbosity,"Enough number of pulse-free intervals to CSD calculus. No more rows read");
		    break;
		}

		sprintf(straux,"%d",ntotalrows);
		message = "-------------> Record: " + string(straux);
		sprintf(straux,"%d",eventcnt);
		message += " of " + string(straux) + " <------------------ ";
		writeLog(fileRef,"Log", verbosity,message);
		sprintf(val,"-------------> Record: %d of %d <------------------ ",ntotalrows,eventcnt);
		strcat(val,"\n");
		fputs(val,temporalFile);

		// Information has been read by blocks (with nrows per block)
		// Now, information is going to be used by rows
		gsl_vector_memcpy(timegsl,timegsl_block);
		gsl_matrix_get_row(ioutgsl,ioutgsl_block,i);
		gsl_vector_scale(ioutgsl,ivcal);		//IVCAL to change arbitrary units of voltage to non-arbitrary units of current (Amps)

		if (energy != 0)	//There are pulses in the input FITS file
		{
			// Assigning positive polarity (by using ASQUID and PLSPOLAR)
			gsl_vector_memcpy(ioutgsl_aux,ioutgsl);
			//if (asquid > 0)		gsl_vector_scale(ioutgsl_aux,-1.0);
			if (((asquid>0) && (plspolar<0)) || ((asquid<0) && (plspolar>0)))	gsl_vector_scale(ioutgsl_aux,-1.0);
			gsl_vector_memcpy(ioutgslNOTFIL,ioutgsl_aux);

			// Low-pass filtering
			status = lpf_boxcar(&ioutgsl_aux,ioutgsl_aux->size,tauFALL*scaleFactor,samprate);
			if (status == EPFAIL)
			{
			    message = "Cannot run routine lpf_boxcar for low pass filtering";
			    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
			if (status == 3)
			{
			    message = "lpf_boxcar: tauFALL*scaleFactot too small => Cut-off frequency too high => Equivalent to not filter.";
			    writeLog(fileRef,"Warning", verbosity,message);
			    status = EPOK;
			}
			if (status == 4)
			{
			    message = "lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low";
			    writeLog(fileRef,"Warning", verbosity,message);
			    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}

			// Derivative after filtering
			if (derMTHSimple (&ioutgsl_aux,&derSGN,ioutgsl_aux->size))
			{
			    message = "Cannot run routine derMTHSimple for derivative after filtering";
			    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}

			//Finding the pulses: Pulses tstart are found
			if (findPulses (ioutgslNOTFIL, ioutgsl_aux, &tstartgsl, &qualitygsl, &energygsl,
					&nPulses,
					0,
					tauFALL, scaleFactor, (int)(ntaus*tauFALL*samprate), samprate,
					samplesUp, nSgms,
					Lb, Lrs,
					library, models,
					safetyMarginTstart,
					stopCriteriaMKC,
					kappaMKC,
					levelPrvPulse,
					temporalFile,
					1))
			{
				message = "Cannot run routine findPulses";
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}

			if (nPulses != 0)	pulseFound = 1;
		}
		else
		{
			pulseFound = 0;
		}

		if (pulseFound == 1)
		{
			// Finding the pulse-free intervals in each record read
			if (findInterval(ioutgsl_aux, tstartgsl, nPulses, ntaus, ntausPF, tauFALL_b, intervalMinBins, &nIntervals, &startIntervalgsl))
			{
				message = "Cannot run routine findInterval";
			    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
		}
		else if (pulseFound == 0)
		{
			// The whole event is going to be used: It must be divided into intervals of intervalMinBins size
			if (findIntervalN(ioutgsl_aux,intervalMinBins,&nIntervals,&startIntervalgsl))
			{
				message = "Cannot run routine findIntervalN";
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
		}

		// CSD calculus (not filtered data)
		for (int i=0; i<nIntervals;i++)
		{
			sprintf(val,"startIntervalgsl %e",gsl_vector_get(startIntervalgsl,i));
			strcat(val,"\n");
			fputs(val,temporalFile);
			sprintf(val,"i %d",i);
			strcat(val,"\n");
			fputs(val,temporalFile);

			if (NumMeanSamples>=nintervals){break;}

			EventSamples = gsl_vector_alloc(intervalMinBins);
			temp = gsl_vector_subvector(ioutgsl,gsl_vector_get(startIntervalgsl,i), intervalMinBins);
			gsl_vector_memcpy(EventSamples,&temp.vector);

			// FFT calculus (EventSamplesFFT)
			if(FFT(EventSamples,vector_aux1,SelectedTimeDuration))
			{
				message = "Cannot run FFT routine for vector1";
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			gsl_vector_complex_absRIFCA(vector_aux,vector_aux1);

			// Every single spectrum is stored in a row of the EventSamplesFFT array
			gsl_matrix_set_row(EventSamplesFFT,NumMeanSamples,vector_aux);

			// Add to mean FFT samples
			gsl_vector_mul(vector_aux,vector_aux);
			gsl_vector_add(EventSamplesFFTMean,vector_aux);

			NumMeanSamples = NumMeanSamples+1;
		}

		ntotalrows++;
	}

	// Free allocate of GSL vectors
	gsl_vector_free(timegsl);
	gsl_vector_free(ioutgsl);
	gsl_vector_free(ioutgsl_aux);
	gsl_vector_free(timegsl_block);
	gsl_matrix_free(ioutgsl_block);
	gsl_vector_free(vector_aux);
	gsl_vector_complex_free(vector_aux1);
	gsl_vector_free(derSGN);
	gsl_vector_free(tstartgsl);
	gsl_vector_free(tstartDERgsl);
	gsl_vector_free(tmaxDERgsl);
	gsl_vector_free(maxDERgsl);
	gsl_vector_free(tendDERgsl);

	return (EPOK);
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* findInterval: This function finds the pulse-free intervals when the input vector has pulses.
*               The pulse-free intervals must have a minimum length (intervalMin).
*               The interval with pulse is Tstart,Tend+nPF*tau (being Tend=n*tau).
*
* Parameters:
* - invector: Input vector WITH pulses
* - startpulse: Vector with the Tstart of all the pulses of the input vector (samples)
* - npin: Number of pulses in the input vector
* - n: Number of tauFALLs after starting the pulse to ending it
* - nPF: Number of tauFALLs after ending the pulse to start the pulse-free interval
* - tau: TauFALL for each pulse (it is an input parameter to the task) (samples)
* - interval: Minimum length of the interval (samples)
* - ni: Number of pulse-free intervals in the input vector
* - startinterval: Vector with the starting time of each pulse-free interval (samples).
*
* - Declare variables
* - Processing if the input vector has more pulses
* 	- It looks for pulse-free intervals between pulses
* - Processing if there are no more pulses in the input vector
* 	-It looks for pulse-free intervals at the end of the event and the search for more pulse-free intervals is finished
*****************************************************************************/
int findInterval(gsl_vector *invector, gsl_vector *startpulse, int npin, int n, int nPF, double tau, int interval,
		int *ni, gsl_vector **startinterval)
{

	char val[256];

	// Declare variables
	long start = 0;		//Auxiliary variable, possible starting time of a pulse-free interval
	int niinc;			//Increase the number of pulse-free intervals
	*ni = 0;

	// startinterval is allocated with the maximum possible pulse-free intervals
	// Not necessary handling division by zero because of interval being always greater than zero
	*startinterval = gsl_vector_alloc(invector->size/interval+1);

	for (int i=0; i<npin; i++)
	{
		//If the input vector has pulses, it looks for pulse-free intervals between pulses
		if (gsl_vector_get(startpulse,i)-start >= 0)
		{
			niinc = (gsl_vector_get(startpulse,i)-start)/interval;

			for (int j=*ni; j<*ni+niinc; j++)
			{
				gsl_vector_set(*startinterval,j,start+(j-*ni)*interval);
			}
			start = gsl_vector_get(startpulse,i)+n*tau+nPF*tau;
			*ni = *ni+niinc;
		}
		else
			//The tend+nPF*tau of a pulse is aliased to the start of the next pulse
		{
			if (gsl_vector_get(startpulse,i)-gsl_vector_get(startpulse,i-1)<0)
			{
				start = gsl_vector_get(startpulse,i-1)+n*tau+nPF*tau;
			}
			else
			{
				start = gsl_vector_get(startpulse,i)+n*tau+nPF*tau;
			}
		}
	}

	if ((int) (invector->size-start) >=interval)
	{
		//If there are no more pulses in the event, it looks for pulse-free intervals at the end of the event
		//and the search for more pulse-free intervals is finished

		niinc = (invector->size-start)/interval;

		for (int j=*ni; j<*ni+niinc; j++)
		{
			gsl_vector_set(*startinterval,j,start+(j-*ni)*interval);
		}
		*ni = *ni+niinc;
	}

	return (EPOK);
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
* findIntervalN: This function finds the pulse-free intervals when the input vector has No pulses.
*                The pulse-free intervals must have a minimum length (intervalMin).
*
* Parameters:
* - invector: Input vector WITHOUT pulses
* - interval: Minimum length of the interval (samples)
* - ni: Number of pulse-free intervals in the input vector
* - startinterval: Vector with the starting time of each pulse-free interval (samples)
*****************************************************************************/
int findIntervalN (gsl_vector *invector, int interval, int *ni, gsl_vector **startinterval)
{
	*ni = invector->size/interval;	//Due to both numerator and denominator are integer numbers =>
									//Result is the integer part of the division

	// startinterval is allocated with the number of pulse-free intervals
	*startinterval = gsl_vector_alloc(*ni);

	for (int i=0;i<*ni;i++)
	{
		gsl_vector_set(*startinterval,i,i*interval);
	}

	return (EPOK);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* createTPSreprFile: This function creates the gennoisespec output FITS file (_noisespec.fits).
*
* - Create the NOISE representation file (if it does not exist already)
* - Create the extensions NOISE and NOISEALL
* - Write keywords
****************************************/
int createTPSreprFile ()
{
	string message = "";
	int status = EPOK, extver=0;

	// Create output FITS files: If it does not exist already
	// Create NOISE representation File
	status = fits_open_file(&tespsObject, tespsName,0,&status);
	if (status == EPOK) 
	{ // file does already exist
	    message = "Output FITS file (type 'noisespec') already exits: " + string(tespsName);
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	else {status = EPOK;}
	if (fits_create_file(&tespsObject, tespsName, &status))
	{
	    message = "Cannot create output gennoisespec file " + string(tespsName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	message = "Create gennoisespec FITS File (_noisespec): " + string(tespsName);
	writeLog(fileRef,"Log", verbosity,message);

	// Create extensions: NOISE
	strcpy(extname,"NOISE");
	if (fits_create_tbl(tespsObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
	{
	    message = "Cannot create table " + string(extname);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Create extensions: NOISEALL
	strcpy(extname,"NOISEALL");
	if (fits_create_tbl(tespsObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
	{
	    message = "Cannot create table " + string(extname);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Set PROC0 keyword
	char str_intervalMin[125];		sprintf(str_intervalMin,"%f",intervalMin);
	char str_ntausPF[125];			sprintf(str_ntausPF,"%d",ntausPF);
	char str_nintervals[125];		sprintf(str_nintervals,"%d",nintervals);
	char str_tauFALL[125];			sprintf(str_tauFALL,"%f",tauFALL);
	char str_scaleFactor[125];		sprintf(str_scaleFactor,"%f",scaleFactor);
	char str_samplesUp[125];		sprintf(str_samplesUp,"%d",samplesUp);
	char str_nSgms[125];	    		sprintf(str_nSgms,"%f",nSgms);
	char str_ntaus[125];			sprintf(str_ntaus,"%d",ntaus);
	char str_LrsT[125];			sprintf(str_LrsT,"%f",LrsT);
	char str_LbT[125];			sprintf(str_LbT,"%f",LbT);
	char str_verb[125];			sprintf(str_verb,"%d",verbosity);

	string process (string("gennoisespec")+ ' ' +
	string(infileName) 		+ ' ' + string(tespsName) 	   + ' ' +
	string(str_intervalMin) + ' ' + string(str_ntausPF) 		+ ' ' + string(str_nintervals) + ' ' +
	string(str_tauFALL)	    + ' ' + string(str_scaleFactor)     + ' ' +
	string(str_samplesUp)   + ' ' + string(str_nSgms)           + ' ' +
	string(str_ntaus) 	    + ' ' +
	string(str_LrsT)       + ' ' + string(str_LbT) + ' ' +
	string(nameLog)	        + ' ' + string(str_verb)            + ' ' +
	string("(")				+      (string) creator  +   	string(")"));

	// Create keywords in NOISE HDU
	strcpy(extname,"NOISE");
	if (fits_movnam_hdu(tespsObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname) + " in file " + string(tespsName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"EVENTCNT");
	keyvalint = intervalMinBins/2;
	if (fits_write_key(tespsObject,TINT,keyname,&keyvalint,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in file " + string(tespsName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"CREATOR");
	strcpy(keyvalstr,creator);
	if (fits_write_key(tespsObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in file " + string(tespsName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"PROC0");
	strcpy(keyvalstr,process.c_str());
	if (fits_write_key_longwarn(tespsObject,&status))
	{
        message = "Cannot write long warn in file " + string(tespsName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (fits_write_key_longstr(tespsObject,keyname,keyvalstr,comment,&status))
	{
  	    message = "Cannot write keyword " + string(keyname) + " in file " + string(tespsName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);  
	}
	
	// Create keywords in NOISEall HDU
	strcpy(extname,"NOISEALL");
	if (fits_movnam_hdu(tespsObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname) + " in file " + string(tespsName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"EVENTCNT");
	if (fits_write_key(tespsObject,TINT,keyname,&intervalMinBins,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in file " + string(tespsName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"CREATOR");
	strcpy(keyvalstr,creator);
	if (fits_write_key(tespsObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in file " + string(tespsName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}	
	strcpy(keyname,"PROC");
	strcpy(keyvalstr,process.c_str());
	if (fits_write_key_longwarn(tespsObject,&status))
	{
        message = "Cannot write long warn in file " + string(tespsName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (fits_write_key_longstr(tespsObject,keyname,keyvalstr,comment,&status))
	{
  	    message = "Cannot write keyword " + string(keyname) + " in file " + string(tespsName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);  
	}

	return EPOK;
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
* writeTPSreprExten: This function writes the noisespec output FITS file (_noisespec.fits).
*
* - Allocate GSL vectors
* - Write the data in the output FITS file (print only half of FFT to prevent aliasing)
* - Free allocate GSL vectors
*****************************************************************************/
int writeTPSreprExten ()
{
	string message = "";
	int status = EPOK;
	double SelectedTimeDuration;

	// Allocate GSL vectors
	gsl_vector *freqgsl = gsl_vector_alloc(intervalMinBins/2);
	gsl_vector *csdgsl = gsl_vector_alloc(intervalMinBins/2);
	gsl_vector *sigmacsdgslnew = gsl_vector_alloc(intervalMinBins/2);
	gsl_vector *freqALLgsl = gsl_vector_alloc(intervalMinBins);
	gsl_vector *csdALLgsl = gsl_vector_alloc(intervalMinBins);

	SelectedTimeDuration=intervalMinBins/((double)samprate);

	// Print only half of FFT to prevent aliasing
	for (int i=0; i< (intervalMinBins/2); i++)
	{
		gsl_vector_set(freqgsl,i,i/SelectedTimeDuration);
		gsl_vector_set(csdgsl,i,gsl_vector_get(EventSamplesFFTMean,i));
		gsl_vector_set(sigmacsdgslnew,i,gsl_vector_get(sigmacsdgsl,i));
		gsl_vector_set(freqALLgsl,i,i/SelectedTimeDuration);
		gsl_vector_set(csdALLgsl,i,gsl_vector_get(EventSamplesFFTMean,i));
	}
	gsl_vector_set(freqALLgsl,intervalMinBins/2,(intervalMinBins/2)/SelectedTimeDuration);
	gsl_vector_set(csdALLgsl,intervalMinBins/2,gsl_vector_get(EventSamplesFFTMean,intervalMinBins/2));
	for (int i=1; i<(intervalMinBins/2); i++)
	{
		gsl_vector_set(freqALLgsl,i+intervalMinBins/2,-(i+intervalMinBins/2-i*2)/SelectedTimeDuration);
		gsl_vector_set(csdALLgsl,i+intervalMinBins/2,gsl_vector_get(EventSamplesFFTMean,i+intervalMinBins/2));
	}

	obj.inObject = tespsObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"NOISE");
	obj.iniRow = 1;
	obj.endRow = intervalMinBins/2;
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	strcpy(obj.nameCol,"FREQ");
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit,"Hz");
	if (writeFitsSimple(obj, freqgsl))
	{
		message = "Cannot run routine writeFitsSimple for freqgsl";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	obj.inObject = tespsObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"NOISE");
	obj.iniRow = 1;
	obj.endRow = intervalMinBins/2;
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	strcpy(obj.nameCol,"CSD");
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit,"A/sqrt(Hz)");
	if (writeFitsSimple(obj, csdgsl))
	{	
		message = "Cannot run routine writeFitsSimple for csdgsl";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	obj.inObject = tespsObject;	
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"NOISE");
	obj.iniRow = 1;
	obj.endRow = intervalMinBins/2;
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	strcpy(obj.nameCol,"SIGMACSD");
	obj.type = TDOUBLE;	
	obj.unit = new char [255];
	strcpy(obj.unit,"A/sqrt(Hz)");

	if (writeFitsSimple(obj, sigmacsdgslnew))
	{
		message = "Cannot run routine writeFitsSimple for sigmacsdgsl";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	obj.inObject = tespsObject;		
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"NOISEALL");
	obj.iniRow = 1;
	obj.endRow = intervalMinBins;
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	strcpy(obj.nameCol,"FREQ");
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit,"Hz");
	if (writeFitsSimple(obj, freqALLgsl))
	{
		message = "Cannot run routine writeFitsSimple for freqALLgsl";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	
	strcpy(obj.nameCol,"CSD");
	strcpy(obj.unit,"A/sqrt(Hz)");
	if(writeFitsSimple(obj, csdALLgsl))
	{
		message = "Cannot run routine writeFitsSimple for csdALLgsl";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Free allocate of GSL vectors
	gsl_vector_free(freqgsl);
	gsl_vector_free(csdgsl);
	gsl_vector_free(sigmacsdgslnew);
	gsl_vector_free(freqALLgsl);
	gsl_vector_free(csdALLgsl);

	return (EPOK);
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
