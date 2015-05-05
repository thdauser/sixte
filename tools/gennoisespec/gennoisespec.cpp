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

   Copyright 2014:  GENNOISESPEC has been developed by the INSTITUTO DE FISICA DE
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
*                   20/01/15   Added clobbering; renaming of parameter because of tasks renaming
*                   27/01/15   The auxiliary file GENNOISESPECauxfile.txt not created
*                   10/04/15   Subtract the baseline (baseline is an input parameter)
*                              Input parameters changed
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
- intervalMinSamples: Minimum length of a pulse-free interval to use (samples) = intervalMinBins
- ntausPF: Number of tauFALLs after ending the pulse (Tend) to start the pulse-free interval
- nintervals: Number of pulse-free intervals to use
- tauFALL: Fall time of the pulses in seconds
- scaleFactor: Scale factor to apply to the fall time of the pulses in order to calculate the LPF box-car length
- samplesUp: Consecutive samples over the threshold to locate a pulse
- nSgms: Number of Sigmas to establish the threshold
- pulse_length: Pulse length (samples)
- LrsT: Running sum length in seconds (only in notcreationlib mode)
- LbT: Baseline averaging length in seconds (only in notcreationlib mode)
- baseline: Baseline (ADC units)
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
	double cutFreq = 0.;
	int boxLength = 0;

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

	// Initialize variables and transform from seconds to samples
	tauFALL_b = (int) (tauFALL * samprate);
	Lrs = LrsT * samprate;
	Lb = LbT * samprate;

	// Pulse-free segments are divided into pulse-free intervals with intervalMinBins size
	intervalMinBins = intervalMinSamples;
	if (intervalMinBins > eventsz)
	{
		message = "Illegal value in INTERVALMINSAMPLES parameter. Legal values reals greater than 0 and fewer than record length";
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

	// Check Boxlength
	cutFreq = 2 * (1/(2*pi*tauFALL*scaleFactor));
	boxLength = (int) ((1/cutFreq) * samprate);
	if (boxLength <= 1)
	{
	      message = "lpf_boxcar: tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.";
	      writeLog(fileRef,"Warning", verbosity,message);
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
		message = "Cannot create file " +  string(gnoiseName);
		EP_EXIT_ERROR(message,EPFAIL);
	}	
	if (fits_open_file(&gnoiseObject,gnoiseName,1,&status))
	{
		message = "Cannot open file " +  string(gnoiseName);
		EP_EXIT_ERROR(message,status);
	}
	
	// Write extensions NOISE and NOISEALL (call writeTPSreprExten)
	if(writeTPSreprExten ())
	{
		message = "Cannot write extensions in " +  string(gnoiseName);
		EP_EXIT_ERROR(message,EPFAIL);
	}

	// Free allocate of GSL vectors
	gsl_matrix_free(EventSamplesFFT);
	gsl_vector_free(EventSamplesFFTMean);
	gsl_vector_free(mean);
	gsl_vector_free(sigmacsdgsl);
	gsl_vector_free(startIntervalgsl);

	// Close output FITS file
	if (fits_close_file(gnoiseObject,&status))
	{
	    message = "Cannot close file " + string(gnoiseName);
	    EP_EXIT_ERROR(message,status);
	}	

	// Free memory
	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;

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

	// Define GENNOISESPEC input parameters and assign values to variables
	// Parameter definition and assignation of default values
	const int npars = 16, npars1 = 17;
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

	gennoisespecPars[2].name = "intervalMinSamples";
	gennoisespecPars[2].description = "Minimum length of a pulse-free interval (samples)";
	gennoisespecPars[2].defValInt = 512;
	gennoisespecPars[2].type = "int";
	gennoisespecPars[2].minValInt = 2;
	gennoisespecPars[2].maxValInt = 50000;
	gennoisespecPars[2].ValInt = gennoisespecPars[2].defValInt;

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

	gennoisespecPars[9].name = "pulse_length";
	gennoisespecPars[9].description = "Pulse length in samples";
	gennoisespecPars[9].defValInt = 375;
	gennoisespecPars[9].type = "int";
	gennoisespecPars[9].minValInt = 1;
	gennoisespecPars[9].maxValInt = 50000;
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

	gennoisespecPars[12].name = "baseline";
	gennoisespecPars[12].description = "Baseline (ADC units)";
	gennoisespecPars[12].defValReal = 200.0;
	gennoisespecPars[12].type = "double";
	gennoisespecPars[12].minValReal = 1.E-50;
	gennoisespecPars[12].maxValReal = 1.E+50;
	gennoisespecPars[12].ValReal = gennoisespecPars[12].defValReal;

	gennoisespecPars[13].name = "nameLog";
	gennoisespecPars[13].description = "Output log file name";
	gennoisespecPars[13].defValStr = "noise_log.txt";
	gennoisespecPars[13].type = "char";
	gennoisespecPars[13].ValStr = gennoisespecPars[13].defValStr;

	gennoisespecPars[14].name = "verbosity";
	gennoisespecPars[14].description = "Verbosity level of the output log file (in [0,3])";
	gennoisespecPars[14].defValInt = 3;
	gennoisespecPars[14].type = "int";
	gennoisespecPars[14].minValInt = 0;
	gennoisespecPars[14].maxValInt = 3;
	gennoisespecPars[14].ValInt = gennoisespecPars[14].defValInt;

	gennoisespecPars[15].name = "clobber";
	gennoisespecPars[15].description = "Re-write output files if clobber=yes";
	gennoisespecPars[15].defValStr = "no";
	gennoisespecPars[15].type = "char";
	gennoisespecPars[15].ValStr = gennoisespecPars[15].defValStr;
	
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
		    	message = "Invalid parameter name ";
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
		    strcpy(gnoiseName, gennoisespecPars[i].ValStr.c_str());
		}
		else if(gennoisespecPars[i].name == "intervalMinSamples")
		{
		    intervalMinSamples = gennoisespecPars[i].ValInt;
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
		else if(gennoisespecPars[i].name == "pulse_length")
		{
		    pulse_length = gennoisespecPars[i].ValInt;
		}
		else if(gennoisespecPars[i].name == "LrsT")
		{
		    LrsT = gennoisespecPars[i].ValReal;
		}
		else if(gennoisespecPars[i].name == "LbT")
		{
		    LbT = gennoisespecPars[i].ValReal;
		}
		else if(gennoisespecPars[i].name == "baseline")
		{
			baseline = gennoisespecPars[i].ValReal;
		}
		else if(gennoisespecPars[i].name == "nameLog")
		{
		    strcpy(nameLog,gennoisespecPars[i].ValStr.c_str());
		}
		else if(gennoisespecPars[i].name == "verbosity")
		{
		    verbosity = gennoisespecPars[i].ValInt;
		}
		else if (gennoisespecPars[i].name == "clobber")
		{
			strcpy(clobberStr, gennoisespecPars[i].ValStr.c_str());
			if(strcmp(clobberStr,"yes")==0){
			  clobber=1;
			}else{
			  clobber=0;
			}
		}
		
		// Check if parameter value is in allowed range
		if( gennoisespecPars[i].type == "int" &&
				(gennoisespecPars[i].ValInt < gennoisespecPars[i].minValInt ||
						gennoisespecPars[i].ValInt > gennoisespecPars[i].maxValInt))
		{
			message = "Parameter name " + gennoisespecPars[i].name + " out of range: [" +
					boost::lexical_cast<std::string>(gennoisespecPars[i].minValInt) + "," + boost::lexical_cast<std::string>(gennoisespecPars[i].maxValInt) + "]";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		else if ( gennoisespecPars[i].type == "double" &&
				(gennoisespecPars[i].ValReal < gennoisespecPars[i].minValReal ||
						gennoisespecPars[i].ValReal > gennoisespecPars[i].maxValReal))
		{
			message = "Parameter name " + gennoisespecPars[i].name + " out of range: [" +
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
	double threshold;
	double cutFreq;
	int boxLength;
	gsl_vector *tstartDERgsl = gsl_vector_alloc(eventsz);
	gsl_vector *tmaxDERgsl = gsl_vector_alloc(eventsz);
	gsl_vector *maxDERgsl = gsl_vector_alloc(eventsz);
	gsl_vector *index_maxDERgsl = gsl_vector_alloc(eventsz);
	gsl_vector *tendDERgsl = gsl_vector_alloc(eventsz);

	int pulseFound = 0;	// 0->The function findTstart has not found any pulse
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

		// Information has been read by blocks (with nrows per block)
		// Now, information is going to be used by rows
		gsl_vector_memcpy(timegsl,timegsl_block);
		gsl_matrix_get_row(ioutgsl,ioutgsl_block,i);
		gsl_vector_scale(ioutgsl,ivcal);		//IVCAL to change arbitrary units of voltage to non-arbitrary units of current (Amps)

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
		    status = EPOK;
		}
		if (status == 4)
		{
		    message = "lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low";
		    writeLog(fileRef,"Error", verbosity,message);
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}

		// Derivative after filtering
		if (derMTHSimple (&ioutgsl_aux,&derSGN,ioutgsl_aux->size))
		{
		    message = "Cannot run routine derMTHSimple for derivative after filtering";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}

		//Finding the pulses: Pulses tstart are found
		if (findPulsesNoise (ioutgslNOTFIL, ioutgsl_aux, &tstartgsl, &qualitygsl, &energygsl,
				&nPulses, &threshold,
				0,
				tauFALL, scaleFactor, pulse_length, samprate,
				samplesUp, nSgms,
				Lb, Lrs,
				library, models,
				stopCriteriaMKC,
				kappaMKC,
				levelPrvPulse))
		{
			message = "Cannot run routine findPulsesNoise";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		if (nPulses != 0)	pulseFound = 1;

		if (pulseFound == 1)
		{
			// Finding the pulse-free intervals in each record read
			if (findInterval(ioutgsl_aux, tstartgsl, nPulses, pulse_length, ntausPF, tauFALL_b, intervalMinBins, &nIntervals, &startIntervalgsl))
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
			if (NumMeanSamples>=nintervals){break;}

			EventSamples = gsl_vector_alloc(intervalMinBins);
			temp = gsl_vector_subvector(ioutgsl,gsl_vector_get(startIntervalgsl,i), intervalMinBins);
			gsl_vector_memcpy(EventSamples,&temp.vector);

			// Baseline estimation
			/*double baseline;
			cutFreq = 2 * (1/(2*pi*tauFALL*scaleFactor));
			boxLength = (int) ((1/cutFreq) * samprate);
			if (find_baseline(EventSamples, kappaMKC, stopCriteriaMKC, boxLength, &baseline))
			{
				message = "Cannot run find_baseline";
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}*/

			// Baseline subtraction
			gsl_vector *baselinegsl = gsl_vector_alloc(EventSamples->size);
			//gsl_vector_set_all(baselinegsl,-200.0);
			gsl_vector_set_all(baselinegsl,-1.0*baseline);
			gsl_vector_add(EventSamples,baselinegsl);
			gsl_vector_free(baselinegsl);

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
* - pulse_length: Pulse length (samples)
* - nPF: Number of tauFALLs after ending the pulse to start the pulse-free interval
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
int findInterval(gsl_vector *invector, gsl_vector *startpulse, int npin, int pulse_length, int nPF, double tau, int interval,
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
			start = gsl_vector_get(startpulse,i)+pulse_length+nPF*tau;
			*ni = *ni+niinc;
		}
		else
			//The tend+nPF*tau of a pulse is aliased to the start of the next pulse
		{
			if (gsl_vector_get(startpulse,i)-gsl_vector_get(startpulse,i-1)<0)
			{
				start = gsl_vector_get(startpulse,i-1)+pulse_length+nPF*tau;
			}
			else
			{
				start = gsl_vector_get(startpulse,i)+pulse_length+nPF*tau;
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

	if ( fileExists(string(gnoiseName)) && clobber==1)
	{
	      if (remove(gnoiseName)){
		  message = "Output FITS file ("+ string(gnoiseName)+") already exits & cannot be deleted ("+string(strerror(errno))+")";
		  EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	      }
	}
	else if(fileExists(string(gnoiseName)) && clobber==0)
	{
	      message = "Output FITS file ("+ string(gnoiseName)+") already exits: must not be overwritten";
	      EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if(!fileExists(string(gnoiseName)))
	{
	      if(fits_create_file(&gnoiseObject, gnoiseName, &status))
	      {
		  message = "Cannot create output gennoisespec file " + string(gnoiseName);
		  EP_PRINT_ERROR(message,status); return(EPFAIL);
	      }
	}
	
	message = "Create gennoisespec FITS File (_noisespec): " + string(gnoiseName);
	writeLog(fileRef,"Log", verbosity,message);

	// Create extensions: NOISE
	strcpy(extname,"NOISE");
	if (fits_create_tbl(gnoiseObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
	{
	    message = "Cannot create table " + string(extname);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Create extensions: NOISEALL
	strcpy(extname,"NOISEALL");
	if (fits_create_tbl(gnoiseObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
	{
	    message = "Cannot create table " + string(extname);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Set PROC0 keyword
	//char str_intervalMin[125];		sprintf(str_intervalMin,"%f",intervalMin);
	char str_intervalMinSamples[125];		sprintf(str_intervalMinSamples,"%d",intervalMinSamples);
	char str_ntausPF[125];			sprintf(str_ntausPF,"%d",ntausPF);
	char str_nintervals[125];		sprintf(str_nintervals,"%d",nintervals);
	char str_tauFALL[125];			sprintf(str_tauFALL,"%f",tauFALL);
	char str_scaleFactor[125];		sprintf(str_scaleFactor,"%f",scaleFactor);
	char str_samplesUp[125];		sprintf(str_samplesUp,"%d",samplesUp);
	char str_nSgms[125];	    		sprintf(str_nSgms,"%f",nSgms);
	char str_pulse_length[125];			sprintf(str_pulse_length,"%d",pulse_length);
	char str_LrsT[125];			sprintf(str_LrsT,"%f",LrsT);
	char str_LbT[125];			sprintf(str_LbT,"%f",LbT);
	char str_baseline[125];			sprintf(str_baseline,"%f",baseline);
	char str_verb[125];			sprintf(str_verb,"%d",verbosity);

	string process (string("gennoisespec")+ ' ' +
	string(infileName) 		+ ' ' + string(gnoiseName) 	   + ' ' +
	//string(str_intervalMin) + ' ' + string(str_ntausPF) 		+ ' ' + string(str_nintervals) + ' ' +
	string(str_intervalMinSamples) + ' ' + string(str_ntausPF) 		+ ' ' + string(str_nintervals) + ' ' +
	string(str_tauFALL)	    + ' ' + string(str_scaleFactor)     + ' ' +
	string(str_samplesUp)   + ' ' + string(str_nSgms)           + ' ' +
	string(str_pulse_length) 	    + ' ' +
	string(str_LrsT)       + ' ' + string(str_LbT) + ' ' +
	string(str_baseline) + ' ' +
	string(nameLog)	        + ' ' + string(str_verb)            + ' ' + string(clobberStr) + ' ' +
	string("(")				+      (string) creator  +   	string(")"));

	// Create keywords in NOISE HDU
	strcpy(extname,"NOISE");
	if (fits_movnam_hdu(gnoiseObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname) + " in file " + string(gnoiseName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"EVENTCNT");
	keyvalint = intervalMinBins/2;
	if (fits_write_key(gnoiseObject,TINT,keyname,&keyvalint,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in file " + string(gnoiseName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"CREATOR");
	strcpy(keyvalstr,creator);
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in file " + string(gnoiseName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"PROC0");
	strcpy(keyvalstr,process.c_str());
	if (fits_write_key_longwarn(gnoiseObject,&status))
	{
        message = "Cannot write long warn in file " + string(gnoiseName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (fits_write_key_longstr(gnoiseObject,keyname,keyvalstr,comment,&status))
	{
  	    message = "Cannot write keyword " + string(keyname) + " in file " + string(gnoiseName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);  
	}
	
	// Create keywords in NOISEall HDU
	strcpy(extname,"NOISEALL");
	if (fits_movnam_hdu(gnoiseObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname) + " in file " + string(gnoiseName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"EVENTCNT");
	if (fits_write_key(gnoiseObject,TINT,keyname,&intervalMinBins,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in file " + string(gnoiseName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"CREATOR");
	strcpy(keyvalstr,creator);
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in file " + string(gnoiseName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}	
	strcpy(keyname,"PROC");
	strcpy(keyvalstr,process.c_str());
	if (fits_write_key_longwarn(gnoiseObject,&status))
	{
        message = "Cannot write long warn in file " + string(gnoiseName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (fits_write_key_longstr(gnoiseObject,keyname,keyvalstr,comment,&status))
	{
  	    message = "Cannot write keyword " + string(keyname) + " in file " + string(gnoiseName);
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

	obj.inObject = gnoiseObject;
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

	obj.inObject = gnoiseObject;
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

	obj.inObject = gnoiseObject;	
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

	obj.inObject = gnoiseObject;		
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

//invector is supposed to be a pulse free interval
int find_baseline(gsl_vector *invector, double kappa, double stopCriteria, int boxLPF, double *baseline)
{
	string message = "";

	// Declare variables
	int size = invector->size; // Size of the input vector
	double mean1, sg1;
	double mean2, sg2;
	gsl_vector_view temp;
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
	//cout<<"median: "<<median<<endl;

	gsl_vector_memcpy(invectorNew,invector);

	// Iterate until no points out of the maximum excursion (kappa*sigma)
	do
	{
		temp = gsl_vector_subvector(invectorNew,0,size-boxLPF-1);
		if (findMeanSigma (&temp.vector, &mean1, &sg1))
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
			if (findMeanSigma (&temp.vector, &mean2, &sg2))
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

			break;
		}

	} while (fabs((mean2-mean1)/mean1)>(stopCriteria/100.0));	// HARDPOINT!!!!!!!!!!!!!!!!!!! (stopCriteria)

	//*baseline = mean2;
	//cout<<"baseline: "<<*baseline<<endl;

	double mean, sigma;
	findMeanSigma (invector, &mean, &sigma);
	//cout<<"mean: "<<mean<<endl;

	*baseline = median;
	//*baseline = mean;

	return(EPOK);
}


int findPulsesNoise
(
	gsl_vector *vectorin,
	gsl_vector *vectorinDER,
	gsl_vector **tstart,
	gsl_vector **quality,
	gsl_vector **energy,

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

	gsl_matrix *librarymatrix,
	gsl_matrix *modelsmatrix,

	double stopcriteriamkc,
	double kappamkc,
	double levelprvpulse)
{
	char val[256];
	char val_aux[256];

 	int status=EPOK;
	string message = "";

	const double pi = 4.0 * atan(1.0);

	// Declare variables
	int pulseFound;
	double thresholdmediankappa;	// Threshold to look for pulses in the first derivative
		// To look for single pulses during the first step
	gsl_vector *maxDERgsl = gsl_vector_alloc(vectorinDER->size);
	gsl_vector *index_maxDERgsl = gsl_vector_alloc(vectorinDER->size);
		// To look for secondary pulses during the second step
	gsl_vector *newPulsesgsl = gsl_vector_alloc(vectorinDER->size); // If a pulse is new => Look again for more pulses
	gsl_vector_set_zero(newPulsesgsl);
	gsl_vector *Lbgsl = gsl_vector_alloc(vectorinDER->size);	    // If there is no free-pulses segments longer than Lb=>
	gsl_vector_set_all(Lbgsl,lb);                                   // segments shorter than Lb will be useed and its length (< Lb)
	                                                                // must be used instead Lb in RS_filter
	gsl_vector *Bgsl;
	/*double Bprev = -999;
	gsl_vector *Bauxgsl;
	bool flagContinue = true;*/

	gsl_vector_set_zero(*quality);
	gsl_vector_set_zero(*energy);									// Estimated energy of the single pulses
																	// In order to choose the proper pulse template to calculate
																	// the adjusted derivative and to fill in the Energy column
		                                                            // in the output FITS file

	// First step to look for single pulses
	if (medianKappaClipping (vectorinDER, kappamkc, stopcriteriamkc, nsgms, (int)(pi*samplingRate*taufall*scalefactor), &thresholdmediankappa))
	{
	    message = "Cannot run medianKappaClipping looking for single pulses";
	    EP_PRINT_ERROR(message,EPFAIL);
	}
	*threshold = thresholdmediankappa;

	if (findTstart (vectorinDER, thresholdmediankappa, samplesup, 1, samplingRate, nPulses, &pulseFound, tstart, quality, &maxDERgsl, &index_maxDERgsl))
	{
	  message = "Cannot run findTstart with two rows in models";
	  EP_PRINT_ERROR(message,EPFAIL);
	}

	for (int i=0;i<*nPulses;i++)
	{
		gsl_vector_set(newPulsesgsl,i,1);
	}

	// In order to look for saturated pulses
	/*// In order to look for saturated pulses
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
	}*/

	if ((opmode == 0) && (*nPulses != 0))
	{
		if (getB(vectorin, *tstart, *nPulses, &Lbgsl, sizepulsebins, &Bgsl))
		{
		    message = "Cannot run getB routine with opmode=0 & nPulses != 0";
		    EP_PRINT_ERROR(message,EPFAIL);
		}
		double energyaux = gsl_vector_get(*energy,0);
		for (int i=0;i<*nPulses;i++)
		{
			if (i != *nPulses-1)	// Not last pulse in the record
			{
				if (getPulseHeight(vectorin, gsl_vector_get(*tstart,i), gsl_vector_get(*tstart,i+1), 0, lrs, gsl_vector_get(Lbgsl,i), gsl_vector_get(Bgsl,i), sizepulsebins, &energyaux))
				{
				    message = "Cannot run getPulseHeight routine when pulse i=" + boost::lexical_cast<std::string>(i) + " is not the last pulse";
				    EP_PRINT_ERROR(message,EPFAIL);
				}
			}
			else
			{
				if (getPulseHeight(vectorin, gsl_vector_get(*tstart,i), gsl_vector_get(*tstart,i+1), 1, lrs, gsl_vector_get(Lbgsl,i), gsl_vector_get(Bgsl,i), sizepulsebins, &energyaux))
				{
				    message = "Cannot run getPulseHeight routine when pulse i=" + boost::lexical_cast<std::string>(i) + " is the last pulse";
				    EP_PRINT_ERROR(message,EPFAIL);
				}
			}
			gsl_vector_set(*energy,i,energyaux);
		}
	}
	/*else if ((opmode == 1) && (*nPulses != 0))  // Iterative pulse searching
	{
		gsl_vector_memcpy(vectorDERComposed, vectorinDER); 	// Some parts will be overwritten
		model = gsl_vector_alloc(modelsmatrix->size2);

		int nPulsesRowAux;
		int nNewPulses;

		do
		{
			// To estimate the pulse height of each pulse
			// Sum of the Lb digitized data samples of a pulse-free interval immediately before the current pulse, B
			if (getB(vectorin, *tstart, *nPulses, &Lbgsl, sizepulsebins, &Bgsl))
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
			if (findSePulsesNoise(vectorin, vectorinDER, &vectorDERComposed,
					 thresholdmediankappa,
					 tstart, quality, energy, &maxDERgsl, &index_maxDERgsl
		             &newPulsesgsl,
					 nPulses,
					 taufall, scalefactor, sizepulsebins, samplingRate,
					 samplesup, nsgms,
					 Bgsl, lrs, Lbgsl,
					 librarymatrix, modelsmatrix, model,
					 stopcriteriamkc, kappamkc, levelprvpulse))
			{
				message = "Cannot run findPulsesNoise routine with opmode=1 & nPulses != 0";
				EP_PRINT_ERROR(message,EPFAIL);
			}

			nNewPulses = *nPulses - nPulsesRowAux;

			for (int j=0;j<*nPulses;j++)
			{
				gsl_vector_set(newPulsesgsl,j,gsl_vector_get(newPulsesgsl,j)-1.0);
			}

		} while (nNewPulses > 0);
	}*/

	// Free allocate of GSL vectors
	gsl_vector_free(maxDERgsl);
	gsl_vector_free(index_maxDERgsl);
	gsl_vector_free(newPulsesgsl);
	/*gsl_vector_free(start0);
	gsl_vector_free(end0);
	gsl_vector_free(vectorAUX);*/

	gsl_vector_free(Lbgsl);
	/*if (opmode == 1)
	{
		gsl_vector_free(model);
		gsl_vector_free(vectorDERComposed);
	}*/

	return(EPOK);
}

/*int findSePulsesNoise
(
	gsl_vector *vectorin,
	gsl_vector *vectorinDER,
	gsl_vector **vectorinDERComposed,

	double thresholdmediankappaSingle,

	gsl_vector **tstart,
	gsl_vector **quality,
	gsl_vector **energy,
	gsl_vector **maxDER,

	gsl_vector **newPulses,

	int *nPulses,

	//gsl_vector *startsaturated,
	//gsl_vector *endsaturated,
	//int nSaturated,

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

	double stopCriteriamkc,
	double kappamkc,
	double levelprvpulse)
{
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

	gsl_vector *qualitySecondary =  gsl_vector_alloc(vectorinDER->size);
	gsl_vector *maxDERSecondary =  gsl_vector_alloc(vectorinDER->size);
	gsl_vector *index_maxDERSecondary =  gsl_vector_alloc(vectorinDER->size);

	int boxLength = (int)pi*taufall*scalefactor*samplingRate;

	//int indSaturated = 0;

	double pulseheight;

	// Auxiliary variables
	gsl_vector *pulse;
	gsl_vector_view temp;					// In order to handle with gsl_vector_view (subvectors)

	// The thresholdmediankappaSingle input parameter is applied to the set of templates (from the library) in order to get
	// the value of the first sample over the threshold and the index of that sample
	if (firstSampleModelsNoise (models, thresholdmediankappaSingle, &firstSamplesgsl, &index_firstSamplesgsl))
	{
		message = "Cannot run firstSampleModelsNoise routine to get the value of the first sample over the threshold and the index of that sample";
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
			if (getPulseHeight(vectorin, gsl_vector_get(*tstart,i), gsl_vector_get(*tstart,i+1), lastOne, lrs, gsl_vector_get(lb,i), gsl_vector_get(B,i), sizepulse, &pulseheight))
			{
				message = "Cannot run getPulseHeight routine for pulse i=" + boost::lexical_cast<std::string>(i) + " when newPulses = 1";
				EP_PRINT_ERROR(message,EPFAIL);
			}
			gsl_vector_set(*energy,i,pulseheight);
			gsl_vector_set(*energy,i,-999);

			// Cut the record part corresponding to a pulse
			limInf = gsl_vector_get(*tstart,i);
			limSup =limInf+sizepulse-1;
			if (limSup>vectorin->size)	limSup = vectorin->size;
			pulse = gsl_vector_alloc(limSup-limInf);
			temp = gsl_vector_subvector(*vectorinDERComposed,limInf,limSup-limInf);
			//temp = gsl_vector_subvector(vectorinDER,limInf,limSup-limInf);
			gsl_vector_memcpy(pulse, &temp.vector);

			// Find the proper pulse model in the pulse templates library
			if (find_model1stSampleNoise(firstSampleOverThreshold, firstSamplesgsl, models, &model))
			{
				message = "Cannot run find_model1stSampleNoise routine for pulse i=" + boost::lexical_cast<std::string>(i) + " when newPulses = 1";
				EP_PRINT_ERROR(message,EPFAIL);
			}

			// Scale the template according to the maximum of the pulse
			gsl_vector_memcpy(modelScaled,model);
			//gsl_vector_scale(modelScaled,firstSampleOverThreshold/firstSampleOverThreshold_Model);

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
				limInf = gsl_vector_get(*tstart,i);
				limSup =limInf+sizepulse-1;
				if (limSup>=(vectorin->size))	limSup = vectorin->size-boxLength;

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
					limInf = gsl_vector_get(*tstart,i);
					limSup =limInf+sizepulse-1;
					if (limSup>=(vectorin->size))	limSup = vectorin->size-boxLength;

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
					limInf = gsl_vector_get(*tstart,i);
					limSup =limInf+sizepulse-1;
					if (limSup>=(vectorin->size))	limSup = vectorin->size-boxLength;

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

			//for (int j=0;j<nSaturated;j++)
			//{
			//	for (int k=gsl_vector_get(startsaturated,j);k<gsl_vector_get(endsaturated,j);k++)
			//	{
			//		gsl_vector_set(*vectorinDERComposed,k,0.0);
			//	}
			//}

			// Free allocate of GSL vectors
			gsl_vector_free(pulse);
			gsl_vector_free(modelShifted);
			gsl_vector_free(diff);
		}
	}

	gsl_vector_free(modelScaled);

	// Look for secondary pulses into 'vectorinDERComposed'
	if (medianKappaClipping (*vectorinDERComposed, kappamkc, stopCriteriamkc, nsgms, (int)(pi*samplingRate*taufall*scalefactor), &thresholdmediankappaSecondary))
	{
		message = "Cannot run medianKappaClipping routine to look for secondary pulses into 'vectorinDERComposed'";
		EP_PRINT_ERROR(message,EPFAIL);
	}

	if (findTstart (*vectorinDERComposed, thresholdmediankappaSecondary, samplesup, 1, samplingRate, &nSecondaryPulses, &istherePulse, &tstartSecondary, &qualitySecondary, &maxDERSecondary, &index_maxDERSecondary))
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
				//if ((gsl_vector_get(tstartSecondary,i) == gsl_vector_get(*tstart,j)) || (gsl_vector_get(tstartSecondary,i) == (gsl_vector_get(*tstart,j)+1)) || (gsl_vector_get(tstartDERSecondary,i) == gsl_vector_get(limSup_vector,j)))
				if (gsl_vector_get(tstartSecondary,i) == gsl_vector_get(*tstart,j))
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
				gsl_vector_set(*maxDER,ind+*nPulses,gsl_vector_get(maxDERSecondary,ind1));
				gsl_vector_set(*index_maxDER,ind+*nPulses,gsl_vector_get(maxDERSecondary,ind1));
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
		gsl_vector *tstartaux = gsl_vector_alloc(*nPulses); // 'tstartNOsmt' subvector
		temp = gsl_vector_subvector(*tstart,0,*nPulses);
		gsl_vector_memcpy(tstartaux, &temp.vector);
		gsl_permutation *perm = gsl_permutation_alloc(tstartaux->size);
		// 'gsl_sort_vector_index' indirectly sorts the elements of the vector v into ascending order, storing the resulting
		// permutation in p. The elements of p give the index of the vector element which would have been stored in that position
		// if the vector had been sorted in place. The first element of p gives the index of the least element in v, and the last
		// element of p gives the index of the greatest element in v. The vector v is not changed.
		// Example: tstartaux=(5200 6000 200 3000) tauxsorted=(200 3000 5200 6000) perm=(2 3 0 1)
		gsl_sort_vector_index(perm,tstartaux);
		gsl_vector *tstartaux1 = gsl_vector_alloc(vectorin->size);
		gsl_vector *newPulsesaux = gsl_vector_alloc(vectorin->size);
		gsl_vector *energyaux = gsl_vector_alloc(vectorin->size);
		gsl_vector *qualityaux = gsl_vector_alloc(vectorin->size);
		gsl_vector *maxDERaux = gsl_vector_alloc(vectorin->size);
		gsl_vector *index_maxDERaux = gsl_vector_alloc(vectorin->size);
		for (int i=0;i<*nPulses;i++)
		{
			gsl_vector_set(tstartaux1,i,gsl_vector_get(*tstart,gsl_permutation_get(perm,i)));
			gsl_vector_set(newPulsesaux,i,gsl_vector_get(*newPulses,gsl_permutation_get(perm,i)));
			gsl_vector_set(energyaux,i,gsl_vector_get(*energy,gsl_permutation_get(perm,i)));
			gsl_vector_set(qualityaux,i,gsl_vector_get(*quality,gsl_permutation_get(perm,i)));
			gsl_vector_set(maxDERaux,i,gsl_vector_get(*maxDER,gsl_permutation_get(perm,i)));
			gsl_vector_set(index_maxDERaux,i,gsl_vector_get(*index_maxDER,gsl_permutation_get(perm,i)));
		}
		gsl_vector_memcpy(*tstart,tstartaux1);
		gsl_vector_memcpy(*newPulses,newPulsesaux);
		gsl_vector_memcpy(*energy,energyaux);
		gsl_vector_memcpy(*quality,qualityaux);
		gsl_vector_memcpy(*maxDER,maxDERaux);
		gsl_vector_memcpy(*index_maxDER,index_maxDERaux);
		gsl_vector_free(tstartaux);
		gsl_permutation_free(perm);
		gsl_vector_free(tstartaux1);
		gsl_vector_free(newPulsesaux);
		gsl_vector_free(energyaux);
		gsl_vector_free(qualityaux);
		gsl_vector_free(maxDERaux);
		gsl_vector_free(index_maxDERaux);
	}

	gsl_vector_free(limSup_vector);
	gsl_vector_free(firstSamplesgsl);
	gsl_vector_free(index_firstSamplesgsl);


	return(EPOK);
}*/


/*int firstSampleModelsNoise (gsl_matrix *templates, double threshold, gsl_vector **firstSamples, gsl_vector **index_firstSamples)
{
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
}*/

/*int find_model1stSampleNoise(double firstSample, gsl_vector *firstSamples, gsl_matrix *models, gsl_vector **modelFound)
{
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
				if (interpolate_model(&modelAux,firstSample,model1,gsl_vector_get(firstSamples,i),model2,gsl_vector_get(firstSamples,i+1)))
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
}*/
