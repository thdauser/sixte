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

   Copyright 2014:  GENNOISESPEC has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01, 
   ESP2013-48637-C2-1-P and ESP2014-53672-C3-1-P.

/***********************************************************************
*                      GENNOISESPEC
*
*  File:       gennoisespec.cpp
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

The purpose of this package is the calculation the current noise spectral density.

MAP OF SECTIONS IN THIS FILE:

 - 1. MAIN
 - 2. initModule
 - 3. inDataIterator
 - 4. findInterval
 - 5. findIntervalN
 - 6. createTPSreprFile
 - 7. writeTPSreprExten
 - 8. find_baseline
 - 9. findPulsesNoise
 - 10. findTstartNoise
 
*******************************************************************************/

#include <gennoisespec.h>


/***** SECTION 1 ************************************
* MAIN function: This function calculates the current noise spectral density.
*                If there are pulses in a record, the pulses are rejected and it is going to look for pulse-free intervals of a given size (intervalMin).
*                If there are no pulses in a record, the event is divided into pulse-free intervals of a given size (intervalMin).
*                It is going to look for pulse-free intervals, calculate their FFT(not filtered data) and average them.
* 
* The output FITS file (_noisespec) contains three columns in two extensions, NOISE and NOISEALL:
* 	- Frequency
*  	- Current noise spectral density: Amount of current per unit (density) of frequency (spectral), as a function of the frequency
* 	- Standard error of the mean
* 	
* The user must supply the following input parameters:
*
* - inFile: Name of the input FITS file
* - outFile: Name of the output FITS file
* - intervalMinSamples: Minimum length of a pulse-free interval to use (samples) = intervalMinBins
* - ntausPF: Number of tauFALLs after ending the pulse (Tend) to start the pulse-free interval
* - nintervals: Number of pulse-free intervals to use
* - tauFALL: Fall time of the pulses in seconds
* - scaleFactor: Scale factor to apply to the fall time of the pulses in order to calculate the LPF box-car length
* - samplesUp: Consecutive samples over the threshold to locate a pulse
* - nSgms: Number of Sigmas to establish the threshold
* - pulse_length: Pulse length (samples)
* - LrsT: Running sum length in seconds (only in notcreationlib mode)
* - LbT: Baseline averaging length in seconds (only in notcreationlib mode)
* - namelog: Output log file name
* - verbosity: Verbosity level of the output log file
* 
* Steps:
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
*  - Free allocated GSL vectors
*  - Close output FITS file
*  - Free memory
*  - Finalize the task
*****************************************************/
int main (int argc, char **argv)
{
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
	
	baseline = gsl_vector_alloc(nintervals);
	gsl_vector_set_all(baseline,-999.0);
	sigma = gsl_vector_alloc(nintervals);
	gsl_vector_set_all(sigma,-999.0);

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

	noiseIntervals = gsl_matrix_alloc(nintervals,intervalMinBins);
	gsl_matrix_set_zero(noiseIntervals);

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

	long rows_per_loop = 0;		// 0: Use default: Optimum number of rows
	long offset = 0;		// 0: Process all the rows
	iteratorCol cols [2];		// Structure of Iteration
	int n_cols = 2; 		// Number of columns:  Time + ADC

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

	// Free allocated GSL vectors
	gsl_matrix_free(EventSamplesFFT);
	gsl_vector_free(EventSamplesFFTMean);
	gsl_vector_free(mean);
	gsl_vector_free(sigmacsdgsl);
	gsl_vector_free(startIntervalgsl);

	gsl_matrix_free(noiseIntervals);
	
	gsl_vector_free(baseline);
	gsl_vector_free(sigma);

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
/*xxxx end of SECTION 1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 2 ************************************************************
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
	const int npars = 15, npars1 = 16;
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
	gennoisespecPars[3].description = "Number of fall time values after the end of the pulse to start the pulse-free interval searching";
	gennoisespecPars[3].defValInt = 0;
	gennoisespecPars[3].type = "int";
	gennoisespecPars[3].minValInt = 0;
	gennoisespecPars[3].maxValInt = 1000;
	gennoisespecPars[3].ValInt = gennoisespecPars[3].defValInt;

	gennoisespecPars[4].name = "nintervals";
	gennoisespecPars[4].description = "Number of pulse-free intervals to use for the noise average";
	gennoisespecPars[4].defValInt = 60;
	gennoisespecPars[4].type = "int";
	gennoisespecPars[4].minValInt = 1;
	gennoisespecPars[4].maxValInt = 20000;
	gennoisespecPars[4].ValInt = gennoisespecPars[4].defValInt;

	gennoisespecPars[5].name = "tauFALL";
	gennoisespecPars[5].description = "Fall time of the pulses (seconds)";
	gennoisespecPars[5].defValReal = 3.5E-4;
	gennoisespecPars[5].type = "double";
	gennoisespecPars[5].minValReal = 1.E-50;
	gennoisespecPars[5].maxValReal = 1.E+50;
	gennoisespecPars[5].ValReal = gennoisespecPars[5].defValReal;

	gennoisespecPars[6].name = "scaleFactor";
	gennoisespecPars[6].description = "Scale factor to apply to the fall time of the pulses to make possible a varying cut-off frequency of the low-pass filter";
	gennoisespecPars[6].defValReal = 0.1;
	gennoisespecPars[6].type = "double";
	gennoisespecPars[6].minValReal = 1.E-50;
	gennoisespecPars[6].maxValReal = 1.E+50;
	gennoisespecPars[6].ValReal = gennoisespecPars[6].defValReal;

	gennoisespecPars[7].name = "samplesUp";
	gennoisespecPars[7].description = "Consecutive samples that the signal must cross over the threshold to trigger a pulse detection";
	gennoisespecPars[7].defValInt = 20;
	gennoisespecPars[7].type = "int";
	gennoisespecPars[7].minValInt = 1;
	gennoisespecPars[7].maxValInt = 1E4;
	gennoisespecPars[7].ValInt = gennoisespecPars[7].defValInt;

	gennoisespecPars[8].name = "nSgms";
	gennoisespecPars[8].description = "Number of quiescent-signal standard deviations to establish the threshold through the kappa-clipping algorithm";
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
	gennoisespecPars[10].description = "Running sum (RS) length for the RS-filtering for raw energy estimation, in seconds";
	gennoisespecPars[10].defValReal = 30.E-6;
	gennoisespecPars[10].type = "double";
	gennoisespecPars[10].minValReal = 1.E-50;
	gennoisespecPars[10].maxValReal = 1.E+50;
	gennoisespecPars[10].ValReal = gennoisespecPars[10].defValReal;

	gennoisespecPars[11].name = "LbT";
	gennoisespecPars[11].description = "Baseline averaging length for the RS-filtering for raw energy estimation, in seconds";
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

	gennoisespecPars[14].name = "clobber";
	gennoisespecPars[14].description = "Re-write output files if clobber=yes";
	gennoisespecPars[14].defValStr = "no";
	gennoisespecPars[14].type = "char";
	gennoisespecPars[14].ValStr = gennoisespecPars[14].defValStr;
	
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
			if(strcmp(clobberStr,"yes")==0)
			{
				clobber=1;
			}
			else
			{
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
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 ************************************************************
* inDataIterator: This function takes the optimum number of rows to read the input FITS file
*                 and works iteratively
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each record
* 	- Assigning positive polarity (by using ASQUID and PLSPOLAR)
* 	- To avoid taking into account the pulse tails at the beginning of a record as part of a pulse-free interval
* 	- Low-pass filtering
*   	- Differentiate after filtering
*   	- Finding the pulses: Pulses tstart are found (call findPulses)
* - Finding the pulse-free intervals in each record
*  	- If there are pulses => Call findInterval
*	- No pulses => The whole event is going to be used (DIVIDING into intervals of intervalMinBins size) => Call findIntervalN
* - CSD calculus (not filtered data):
* 	- Baseline subtraction
* 	- FFT calculus (each pulse-free interval)
* 	- Every single spectrum of each pulse-free interval is stored in a row of the EventSamplesFFT array
* 	- Add to mean FFT samples
* - Free allocated GSL vectors
****************************************************************************/
int inDataIterator(long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
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

	// To differentiate
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

	double baselineI;
	double mean, sg; // To handle the pulse tails at the beginning of the record
	int tail_duration;
	
	double baselineIntervalFreeOfPulses;
	double sigmaIntervalFreeOfPulses;

	// Auxiliary variables
	gsl_vector_view temp;	// In order to handle with gsl_vector_view (subvectors)

	// To calculate the FFT
	double SelectedTimeDuration;
	gsl_vector *EventSamples;
	gsl_vector *vector_aux;
	gsl_vector_complex *vector_aux1;
	vector_aux = gsl_vector_alloc(intervalMinBins);
	vector_aux1 = gsl_vector_complex_alloc(intervalMinBins);

	// Allocate input GSL vectors
	timegsl_block = gsl_vector_alloc(nrows);
	timegsl = gsl_vector_alloc(nrows); 			// GSL of input TIME column
	ioutgsl_block = gsl_matrix_alloc(nrows,eventsz);
	ioutgsl = gsl_vector_alloc(eventsz); 			// GSL of input ADC column
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
		if (((asquid>0) && (plspolar<0)) || ((asquid<0) && (plspolar>0)))	gsl_vector_scale(ioutgsl_aux,-1.0);
		gsl_vector_memcpy(ioutgslNOTFIL,ioutgsl_aux);

		// To avoid taking into account the pulse tails at the beginning of a record as part of a pulse-free interval
		tail_duration = 0;
		temp = gsl_vector_subvector(ioutgslNOTFIL,0,eventsz-(int)(pi*samprate*tauFALL*scaleFactor)-1);
		cutFreq = 2 * (1/(2*pi*tauFALL*scaleFactor));
		boxLength = (int) ((1/cutFreq) * samprate);
		if (find_baseline(&temp.vector, kappaMKC, stopCriteriaMKC, boxLength,  &mean, &sg, &baselineI))
		{
			message = "Cannot run find_baseline";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		for (int j=0;j<eventsz;j++)
		{
			if (gsl_vector_get(ioutgslNOTFIL,j) > baselineI+sg)	tail_duration = j;
			else break;
		}

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

		// Differentiate after filtering
		if (differentiate (&ioutgsl_aux,ioutgsl_aux->size))
		{
			message = "Cannot run routine differentiate for differentiating after filtering";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}

		//Finding the pulses: Pulses tstart are found
		if (findPulsesNoise (ioutgslNOTFIL, ioutgsl_aux, &tstartgsl, &qualitygsl, &energygsl,
			&nPulses, &threshold,
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

		if ((pulseFound == 1) || (tail_duration != 0))
		{
			// Finding the pulse-free intervals in each record 
			if (findInterval(tail_duration, ioutgsl_aux, tstartgsl, nPulses, pulse_length, ntausPF, tauFALL_b, intervalMinBins, &nIntervals, &startIntervalgsl))
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

			// Baseline subtraction
			findMeanSigma (EventSamples, &baselineIntervalFreeOfPulses, &sigmaIntervalFreeOfPulses);
			gsl_vector *baselinegsl = gsl_vector_alloc(EventSamples->size);
			gsl_vector_set_all(baselinegsl,-1.0*baselineIntervalFreeOfPulses);
			gsl_vector_add(EventSamples,baselinegsl);
			gsl_vector_free(baselinegsl);
			gsl_vector_set(baseline,indexBaseline,baselineIntervalFreeOfPulses);
			gsl_vector_set(sigma,indexBaseline,sigmaIntervalFreeOfPulses);
			indexBaseline++;

			gsl_matrix_set_row(noiseIntervals,NumMeanSamples,EventSamples);

			// FFT calculus (EventSamplesFFT)
			if(FFT(EventSamples,vector_aux1,SelectedTimeDuration))
			{
				message = "Cannot run FFT routine for vector1";
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			gsl_vector_complex_absIFCA(vector_aux,vector_aux1);

			// Every single spectrum is stored in a row of the EventSamplesFFT array
			gsl_matrix_set_row(EventSamplesFFT,NumMeanSamples,vector_aux);

			// Add to mean FFT samples
			gsl_vector_mul(vector_aux,vector_aux);
			gsl_vector_add(EventSamplesFFTMean,vector_aux);

			NumMeanSamples = NumMeanSamples+1;
		}

		ntotalrows++;
	}

	// Free allocated GSL vectors
	gsl_vector_free(timegsl);
	gsl_vector_free(ioutgsl);
	gsl_vector_free(ioutgsl_aux);
	gsl_vector_free(timegsl_block);
	gsl_matrix_free(ioutgsl_block);
	gsl_vector_free(vector_aux);cutFreq = 2 * (1/(2*pi*tauFALL*scaleFactor));
	boxLength = (int) ((1/cutFreq) * samprate);
	gsl_vector_complex_free(vector_aux1);
	gsl_vector_free(derSGN);
	gsl_vector_free(tstartgsl);
	gsl_vector_free(tstartDERgsl);
	gsl_vector_free(tmaxDERgsl);
	gsl_vector_free(maxDERgsl);
	gsl_vector_free(tendDERgsl);

	return (EPOK);
}
/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
* findInterval: This function finds the pulse-free intervals when the input vector has pulses.
*               The pulse-free intervals must have a minimum length (intervalMin).
*               The interval with pulse is Tstart,Tend+nPF*tau (being Tend=n*tau).
*
* Parameters:
* - tail_duration: Length of the tail of a previous pulse
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
int findInterval(int tail_duration, gsl_vector *invector, gsl_vector *startpulse, int npin, int pulse_length, int nPF, double tau, int interval,
	int *ni, gsl_vector **startinterval)
{
	// Declare variables
	long start;		// Auxiliary variable, possible starting time of a pulse-free interval
	int niinc;		// Increase the number of pulse-free intervals
	*ni = 0;

	// startinterval is allocated with the maximum possible pulse-free intervals
	// Not necessary handling division by zero because of interval being always greater than zero
	*startinterval = gsl_vector_alloc(invector->size/interval+1);

	if ((tail_duration != 0) && (tail_duration >= gsl_vector_get(startpulse,0)))	start = gsl_vector_get(startpulse,0);
	else										start = tail_duration;

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
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
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
	*ni = invector->size/interval;	// Due to both numerator and denominator are integer numbers =>
					// Result is the integer part of the division

	// startinterval is allocated with the number of pulse-free intervals
	*startinterval = gsl_vector_alloc(*ni);

	for (int i=0;i<*ni;i++)
	{
		gsl_vector_set(*startinterval,i,i*interval);
	}

	return (EPOK);
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
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
		if (remove(gnoiseName))
		{
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

	// Primary HDU
	strcpy(extname,"Primary");
	int *hdutype;
	if (fits_movabs_hdu(gnoiseObject, 1, hdutype, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strcpy(keyname,"HISTORY");
	const char * charhistory= "HISTORY Starting parameter list";
	strcpy(keyvalstr,charhistory);
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	string strhistory (string("Input File = ") + string(infileName));
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strhistory = string("Noise File = ") + string(gnoiseName);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_intervalMinSamples[125];		sprintf(str_intervalMinSamples,"%d",intervalMinSamples);
	strhistory=string("intervalMinSamples = ") + string(str_intervalMinSamples);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_ntausPF[125];			sprintf(str_ntausPF,"%d",ntausPF);
	strhistory=string("ntausPF = ") + string(str_ntausPF);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_nintervals[125];		sprintf(str_nintervals,"%d",nintervals);
	strhistory=string("nintervals = ") + string(str_nintervals);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_tauFALL[125];			sprintf(str_tauFALL,"%f",tauFALL);
	strhistory=string("tauFall = ") + string(str_tauFALL);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_scaleFactor[125];		sprintf(str_scaleFactor,"%f",scaleFactor);
	strhistory=string("scaleFactor = ") + string(str_scaleFactor);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_samplesUp[125];		sprintf(str_samplesUp,"%d",samplesUp);
	strhistory=string("samplesUp = ") + string(str_samplesUp);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_nSgms[125];	    		sprintf(str_nSgms,"%f",nSgms);
	strhistory=string("nSgms = ") + string(str_nSgms);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_pulse_length[125];			sprintf(str_pulse_length,"%d",pulse_length);
	strhistory=string("PulseLength = ") + string(str_pulse_length);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_LrsT[125];			sprintf(str_LrsT,"%f",LrsT);
	strhistory=string("LrsT = ") + string(str_LrsT);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_LbT[125];			sprintf(str_LbT,"%f",LbT);
	strhistory=string("LbT = ") + string(str_LbT);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_baseline[125];			sprintf(str_baseline,"%f",baseline);
	strhistory=string("baseline = ") + string(str_baseline);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strhistory = string("NameLog = ") + string(nameLog);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_verb[125];			sprintf(str_verb,"%d",verbosity);
	strhistory=string("verbosity = ") + string(str_verb);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_clobber[125];      sprintf(str_clobber,"%d",clobber);
	strhistory=string("clobber = ") + string(str_clobber);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	charhistory= "HISTORY Ending parameter list";
	strcpy(keyvalstr,charhistory);
	if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in noise file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	return EPOK;
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* writeTPSreprExten: This function writes the noisespec output FITS file (_noisespec.fits).
*
* - Allocate GSL vectors
* - Write the data in the output FITS file (print only half of FFT to prevent aliasing)
* - Free allocated GSL vectors
*****************************************************************************/
int writeTPSreprExten ()
{
	string message = "";
	int status = EPOK;
	int extver=0;
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
		
	strcpy(keyname,"BASELINE");
	double sumBaseline;
	gsl_vector_Sumsubvector(baseline, 0, nintervals, &sumBaseline);
	double keyvaldouble = sumBaseline/nintervals;
	if (fits_write_key(gnoiseObject,TDOUBLE,keyname,&keyvaldouble,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	
	strcpy(keyname,"NOISESTD");
	double sumSigma;
	gsl_vector_Sumsubvector(sigma, 0, nintervals, &sumSigma);
	keyvaldouble = sumSigma/nintervals;
	if (fits_write_key(gnoiseObject,TDOUBLE,keyname,&keyvaldouble,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in file " + string(gnoiseName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
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

	// Free allocated GSL vectors
	gsl_vector_free(freqgsl);
	gsl_vector_free(csdgsl);
	gsl_vector_free(sigmacsdgslnew);
	gsl_vector_free(freqALLgsl);
	gsl_vector_free(csdALLgsl);

	return (EPOK);
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
* find_baseline function: 
* 
****************************************************************************/
int find_baseline(gsl_vector *invector, double kappa, double stopCriteria, int boxLPF, double *mean, double *sigma, double *baseline)
{
	string message = "";

	// Declare variables
	int size = invector->size; // Size of the input vector
	double mean1, sg1;
	double mean2, sg2;
	gsl_vector_view temp;
	// Variables to remove input vector elements higher than the maximum excursion (kappa*sg)
	int i;							// To go through the elements of a vector
	int cnt;						// Number of points inside the excursion (mean+-excursion)
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

	*baseline = mean2;
	*mean = mean2;
	*sigma = sg2;
	//cout<<"baseline: "<<*baseline<<endl;
	//cout<<"sigma: "<<*sigma<<endl;
	
	/*double meanX, sigmaX;
	findMeanSigma (invector, &meanX, &sigmaX);
	cout<<"meanX: "<<meanX<<endl;
	cout<<"sigmaX: "<<sigmaX<<endl;*/

	//*baseline = median;
	//*baseline = mean;

	return(EPOK);
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 9 ************************************************************
* findPulsesNoise function: 
* 
****************************************************************************/
int findPulsesNoise
(
	gsl_vector *vectorin,
	gsl_vector *vectorinDER,
	gsl_vector **tstart,
	gsl_vector **quality,
	gsl_vector **energy,

	int *nPulses,
	double *threshold,

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
	string message = "";

	const double pi = 4.0 * atan(1.0);

	// Declare variables
	int pulseFound;
	double thresholdmediankappa;	// Threshold to look for pulses in the first derivative

	gsl_vector *maxDERgsl = gsl_vector_alloc(vectorinDER->size);
	gsl_vector *index_maxDERgsl = gsl_vector_alloc(vectorinDER->size);

	gsl_vector *Lbgsl = gsl_vector_alloc(vectorinDER->size);	// If there is no free-pulses segments longer than Lb=>
	gsl_vector_set_all(Lbgsl,lb);                                   // segments shorter than Lb will be useed and its length (< Lb)
	                                                                // must be used instead Lb in RS_filter
	gsl_vector *Bgsl;

	gsl_vector_set_zero(*quality);
	gsl_vector_set_zero(*energy);					// Estimated energy of the single pulses
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

	if (findTstartNoise (100,vectorinDER, thresholdmediankappa, samplesup, 1, samplingRate, nPulses, &pulseFound, tstart, quality, &maxDERgsl, &index_maxDERgsl))
	{
		message = "Cannot run findTstart with two rows in models";
		EP_PRINT_ERROR(message,EPFAIL);
	}

	if (*nPulses != 0)
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

	// Free allocated GSL vectors
	gsl_vector_free(maxDERgsl);
	gsl_vector_free(index_maxDERgsl);
	gsl_vector_free(Lbgsl);

	return(EPOK);
}
/*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 10 ************************************************************
* findTstartNoise function: 
* 
****************************************************************************/
int findTstartNoise(int maxPulsesPerRecord, gsl_vector *der, double adaptativethreshold, int nSamplesUp, int allPulsesMode, double sampling, int *numberPulses, int *thereIsPulse, gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl, gsl_vector **index_maxDERgsl)
{
	// Declare variables
	string message="";
	int szRw = der->size;		 // Size of segment of process
	*numberPulses = 0;
	*thereIsPulse = 0;
	int i = 0;			 // To go through the elements of a vector
	bool prevPulse = false;		 // false: It looks for nSamplesUp consecutive samples over the threshold
	                                 // true: It looks for nSamplesUp consecutive samples below the threshold
	int cntDown = 0;		 // To taking into account how many consecutive samples are down the threshold
	int cntUp = 0;			 // To taking into account how many consecutive samples are over the threshold
	int possibleTstart;		 // To store the first of the nSamplesUp consecutive samples over the threshold
	int possibleTend;		 // To store the first of the nSamplesUp consecutive samples below the threshold
	int maxDER_index;		 // Maximum of the derivative between tstartDER and tendDER
	gsl_vector_view temp;		 // In order to handle with gsl_vector_view (subvectors)
	bool lastPulse = true;
	bool maxFound = false;

	// Allocate GSL vectors
	gsl_vector *tstartDER = gsl_vector_alloc(maxPulsesPerRecord);	// tstarts referred to the derivative
	gsl_vector *tendDER = gsl_vector_alloc(maxPulsesPerRecord);	// tends referred to the derivative
	gsl_vector_set_zero(tendDER);
	*tstartgsl = gsl_vector_alloc(maxPulsesPerRecord);		// tstarts referred to the not filtered event
	*flagTruncated = gsl_vector_alloc(maxPulsesPerRecord);
	gsl_vector_set_zero(*flagTruncated);
	*maxDERgsl = gsl_vector_alloc(maxPulsesPerRecord);		// Maximum of the first derivative
	gsl_vector_set_all(*maxDERgsl,-1E3);
	*index_maxDERgsl = gsl_vector_alloc(maxPulsesPerRecord);	// Index of the maximum of the first derivative

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
					gsl_vector_set(*index_maxDERgsl,*numberPulses,i);
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
						gsl_vector_set(*index_maxDERgsl,*numberPulses,i);
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
						gsl_vector_set(*index_maxDERgsl,*numberPulses,i);
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
					gsl_vector_set(*index_maxDERgsl,*numberPulses-1,i);
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

	gsl_vector_memcpy(*tstartgsl,tstartDER);

	gsl_vector_free(tstartDER);
	gsl_vector_free(tendDER);

	if (*numberPulses > 0) *thereIsPulse = 1;

	return (EPOK);
}
/*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/