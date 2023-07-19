/***********************************************************************
 *   This file is part of SIXTE/SIRENA software.
 *
 *   SIXTE is free software: you can redistribute it and/or modify it
 *   under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   any later version.
 *
 *   SIXTE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 *
 *   For a copy of the GNU General Public License see
 *   <http://www.gnu.org/licenses/>.
 *
 *  Copyright 2023:  LOG has been developed by the INSTITUTO DE FISICA DE
 *  CANTABRIA (CSIC-UC) with funding under different projects:
 *  ESP2006-13608-C02-01, AYA2012-39767-C02-01, ESP2013-48637-C2-1-P,
 *  ESP2014-53672-C3-1-P, RTI2018-096686-B-C21 and PID2021-122955OB-C41.
 *
 ***********************************************************************
 *                      INTEGRASIRENA
 *
 *  File:       integraSIRENA.cpp
 *  Developers: Beatriz Cobo
 * 	            cobo@ifca.unican.es
 *              IFCA
 *              Maite Ceballos
 *              ceballos@ifca.unican.es
 *              IFCA
 *
 ***********************************************************************/
 
 /******************************************************************************
  * DESCRIPTION:
  * 
  * The purpose of this package is the integration os SIRENA in SIXTE.
  * 
  * MAP OF SECTIONS IN THIS FILE:
  * 
  * - 1. initializeReconstructionSIRENA
  * - 2. reconstructRecordSIRENA
  * - 3. newReconstructInitSIRENA
  * - 4. freeReconstructInitSIRENA
  * - 5. newPulsesCollection
  * - 6. freePulsesCollection
  * - 7. newOptimalFilterSIRENA
  * - 8. freeOptimalFilterSIRENA
  * - 9. getLibraryCollection
  * - 10. getNoiseSpec
  * - 11. IntegrafreeTesEventListSIRENA
  * - 12. checksum
  * 
  *******************************************************************************/ 
 
 #include "log.h"
 #include "scheduler.h"
 
 #include "genutils.h"
 #include "tasksSIRENA.h"
 
 #define POOLS
 const unsigned int POOL_SIZE = 200;
 const unsigned int MUL_FAC   = 2;
 const unsigned int ADDITION  = 100;
 const unsigned int MAX_SIZE  = 3201;
 int _resize_array(int size, int pulses) 
 {
     int new_size = (size < (int)MAX_SIZE) ? (size * (int)MUL_FAC) : (size + (int)ADDITION);
     return (new_size < pulses) ? pulses : new_size;
 }
 #if 0
 int _resize_array(int size, int pulses){ return (size + ADDITION < pulses) ? pulses : size + ADDITION; }
 int _resize_array(int size, int pulses){ return (size * MUL_FAC  < pulses) ? pulses : size * MUL_FAC; }
 #endif
 #define resize_array(size, pulses) _resize_array(size, pulses)
 
 /***** SECTION 1 ************************************************************
  * initializeReconstructionSIRENA: This function initializes the structure ReconstructInitSIRENA with the variables required 
  *                                 for SIRENA reconstruction. The values are taken from the input parameters.
  * 
  * - Load LibraryCollection structure if library file exists
  * - Load NoiseSpec structure
  * - Fill in the matrix tstartPulse1_i if tstartPulse1 = nameFile
  * - Fill in reconstruct_init
  * 
  * Parameters:
  * - reconstruct_init: Member of ReconstructInitSIRENA structure to initialize the reconstruction parameters (pointer and values)
  * - record_file: Filename of input data file with records
  * - fptr: FITS object with pointer to data file
  * - library_file: File name of calibration library
  * - event_file: File name of output events (with reconstructed energy)
  * - pulse_length: Pulse length
  * - scaleFactor: Detection scale factor for initial filtering
  * - samplesUp: Number of samples for threshold trespassing
  * - samplesDown: Number of samples below the threshold to look for other pulse
  * - nSgms: Number of standard deviations in the kappa-clipping process for threshold estimation
  * - detectSP: Detect secondary pulses or not
  * - opmode: Calibration run (0) or energy reconstruction run (1)
  * - detectionMode: Adjusted Derivative (AD) or Single Threshold Crossing (STC)
  * - LrsT: Running sum length for the RS raw energy estimation (seconds)
  * - LbT: Baseline averaging length for the RS raw energy estimation (seconds)
  * - noise_file: Noise file
  * - filter_domain: Filtering Domain: Time(T) or Frequency(F)
  * - filter_method: Filtering Method: F0 (deleting the zero frequency bin) or B0 (deleting the baseline) or F0B0 (deleting always the baseline)
  * - energy_method: Energy calculation Method: OPTFILT, WEIGHT, WEIGHTN, I2R, I2RFITTED or PCA
  * - filtEeV: Energy of the filters of the library to be used to calculate energy (only for OPTFILT, I2R and I2RFITTED)
  * - Ifit: Constant to apply the I2RFITTED conversion
  * - ofnoise: Noise to use with Optimal Filtering: NSD or WEIGHTM
  * - lagsornot: Lags (1) or no lags (0)
  * - nLags: Number of lags (positive odd number)
  * - Fitting35: Number of lags to analytically calculate a parabola (3) or to fit a parabola (5)
  * - ofiter: Iterate (1) or not iterate (0)
  * - oflib: Work or not with a library with optimal filters (1/0)
  * - ofinterp: Optimal Filter by using the Matched Filter or the DAB as matched filter (MF/DAB)
  *             It has been fixed in 'tesreconstruction' as 'DAB'
  * - oflength_strategy: Optimal Filter length Strategy: FREE, BASE2, BYGRADE or FIXED
  * - oflength: Optimal Filter length (taken into account if :option:`OFStrategy`=FIXED)
  * - preBuffer: Some samples added or not before the starting time of a pulse (number of added samples read from the xml file)
  * - monoenergy: Monochromatic energy of input file in eV (only for library creation)
  * - hduPRECALWN: Add or not the PRECALWN HDU in the library file (1/0) (only for library creation)
  * - hduPRCLOFWM: Add or not the PRCLOFWM HDU in the library file (1/0) (only for library creation)
  * - largeFilter: Length of the longest fixed filter
  * - interm: Write or not intermediate files (1/0)
  * - detectFile: Intermediate detections file (if intermediate=1)
  * - errorT: Additional error (in samples) added to the detected time (logically, it changes the reconstructed energies)
  * - Sum0Filt: 0-padding: Subtract the sum of the filter (1) or not (0)
  * - clobber: Overwrite or not output files if exist (1/0)
  * - maxPulsesPerRecord: Default size of the event list
  * - SaturationValue: Saturation level of the ADC curves
  * - tstartPulse1: If integer number => Sample where the first pulse starts 
  *                  or
  *                  if nameFile => File where is the tstart (in seconds) of every pulse
  * - tstartPulse2: Tstart (samples) of the second pulse
  * - tstartPulse3: Tstart (samples) of the third pulse (if 0 => PAIRS, if not 0 => TRIOS)
  * - energyPCA1: First energy (only for PCA) 
  * - energyPCA2: Second energy (only for PCA)
  * - XMLFile: File name of the XML input file with instrument definition
  * - status: Input/output status
  ******************************************************************************/
 extern "C" void initializeReconstructionSIRENA(ReconstructInitSIRENA* reconstruct_init, char* const record_file, fitsfile *fptr,
                                                char* const library_file, char* const event_file, int pulse_length, double scaleFactor, int samplesUp, int samplesDown,
                                                double nSgms, int detectSP, int opmode, char *detectionMode, double LrsT, double LbT, char* const noise_file, char* filter_domain, char* filter_method, 
                                                char* energy_method, double filtEev, double Ifit, char *ofnoise, int lagsornot, int nLags, int Fitting35, int ofiter, char oflib, char *ofinterp,
                                                char* oflength_strategy, int oflength, char preBuffer,
                                                double monoenergy, char hduPRECALWN, char hduPRCLOFWM, int largeFilter, int interm, char* const detectFile, int errorT,
                                                int Sum0Filt,
                                                char clobber, int maxPulsesPerRecord, double SaturationValue,
                                                char* const tstartPulse1, int tstartPulse2, int tstartPulse3, double energyPCA1, double energyPCA2, char * const XMLFile, int* const status)
 {  
     string message = "";

     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     // Load LibraryCollection structure if library file exists
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     int exists=0;
     if (fits_file_exists(library_file, &exists, status))
     {
         EP_EXIT_ERROR("Error checking if library file exists",*status);
     }
     
     if (opmode == 0)	// Calibration
     {
         pulse_length = pow(2,floor(log2(largeFilter)));
     }
     
     // Loading in the reconstruct_init structure values related to grading and preBuffer values from the XML file 
     if (0 != preBuffer)	reconstruct_init->preBuffer = 1;
     else			        reconstruct_init->preBuffer = 0;
     gsl_vector *pBi = gsl_vector_alloc(1);   // preBuffer values
     gsl_vector *posti = gsl_vector_alloc(1); // Filter length (including preBuffer)
                                              // post in (grading=>pre,post and pB)
                                              // filtlen in (grading=>pre,post and filtlen)
     if (reconstruct_init->preBuffer == 1)
     {
        gsl_vector_free(pBi); pBi=0;
        pBi = gsl_vector_alloc(reconstruct_init->grading->ngrades);
        gsl_matrix_get_col(pBi,reconstruct_init->grading->gradeData,2);
        reconstruct_init->preBuffer_max_value = gsl_vector_max(pBi);
        reconstruct_init->preBuffer_min_value = gsl_vector_min(pBi);
        gsl_vector_free(posti); posti=0;
        posti = gsl_vector_alloc(reconstruct_init->grading->ngrades);
        gsl_matrix_get_col(posti,reconstruct_init->grading->gradeData,1);
        reconstruct_init->post_max_value = gsl_vector_max(posti);
        reconstruct_init->post_min_value = gsl_vector_min(posti);
     }

     if (exists)
     {
         if (opmode == 1)		
         {
             if ((pulse_length < oflength) && ((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0)|| (strcmp(energy_method,"I2RDER") == 0)))
             {
                 EP_PRINT_ERROR("0-padding is going to be used",-999); // Only a warning
             }
             else if ((pulse_length > oflength) && ((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)))
             {
                 pulse_length = oflength;
             }
             largeFilter = pulse_length;
         }
         
         reconstruct_init->library_collection = getLibraryCollection(library_file, opmode, hduPRECALWN, hduPRCLOFWM, largeFilter, filter_domain, pulse_length, energy_method, ofnoise, filter_method, oflib, &ofinterp, filtEev, lagsornot, reconstruct_init->preBuffer, pBi, posti, status);
         reconstruct_init->library_collection->margin = 0.25;	// (%) Margin to be applied when several energies in the library to choose the proper filter
         if (*status)
         {
             EP_EXIT_ERROR((char*)"Error in getLibraryCollection",EPFAIL); 
         }
         
         if ((opmode == 1) && (pulse_length > reconstruct_init->library_collection->pulse_templates[0].template_duration))
         {
             if ((oflib == 0)
                 && ((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)))
             {
                 EP_EXIT_ERROR("It is not possible PulseLength>PULSE_column_length and OFLib=no",EPFAIL);
             }
             else if ((strcmp(energy_method,"WEIGHT") == 0) || (strcmp(energy_method,"WEIGHTN") == 0) || (strcmp(energy_method,"PCA") == 0))
             {
                 EP_EXIT_ERROR("Templates length in the library file must be at least as the pulse length",EPFAIL);
             }
         }
     }
     else if (!exists && opmode==1)
     {
         EP_EXIT_ERROR((char*)"Error accessing library file: it does not exists ",EPFAIL); 
     }
     
     if (reconstruct_init->preBuffer == 1)
     {
        if (pBi != NULL) {gsl_vector_free(pBi); pBi = 0;}
        if (posti != NULL) {gsl_vector_free(posti); posti = 0;}
     }
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     // Load NoiseSpec structure
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     bool baselineReadFromLibrary = true;
     if ((((((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)) && (strcmp(filter_method,"B0") == 0)) )|| (strcmp(energy_method,"WEIGHT") == 0)) && (opmode == 1) && (oflib == 1))
     {
         if (reconstruct_init->library_collection->baseline == -999.0)  baselineReadFromLibrary = false;
     }
     
     if ((opmode == 0) || 
         (((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)) && (opmode == 1) && (oflib == 0))
         || ((opmode == 1) && (strcmp(energy_method,"WEIGHT") == 0) && (oflib == 0))
         || ((opmode == 1) && (strcmp(energy_method,"WEIGHTN") == 0) && (oflib == 0))
         // If BASELINE is not in the library file, it is necessary to read its value from the noise file 
         || (((((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0))  && (strcmp(filter_method,"B0") == 0)) || (strcmp(energy_method,"WEIGHT") == 0)) && (opmode == 1) && (oflib == 1) && (baselineReadFromLibrary == false)))
     {
         exists=0;
         if(fits_file_exists(noise_file, &exists, status))
         {
             EP_EXIT_ERROR("Error checking if noise file exists",*status);
         }
         if ((opmode == 1) && (exists != 1) && (baselineReadFromLibrary == false))
         {
             EP_EXIT_ERROR("B0 chosen and BASELINE keyword not in the library file => Noise file is necessary but it does not exist",EPFAIL);
         }
         if ((exists != 1) && (baselineReadFromLibrary == true))
         {
             EP_EXIT_ERROR("The necessary noise file does not exist",EPFAIL);
         }
         reconstruct_init->noise_spectrum = getNoiseSpec(noise_file, opmode, hduPRCLOFWM, energy_method, ofnoise, filter_method, status);
         if (*status)
         {
             EP_EXIT_ERROR((char*)"Error in getNoiseSpec",EPFAIL);
         }
         if (((((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)) && (strcmp(filter_method,"B0") == 0)) || (strcmp(energy_method,"WEIGHT") == 0)) && (opmode == 1) && (oflib == 1) && (baselineReadFromLibrary == false))
         {
             reconstruct_init->library_collection->baseline = reconstruct_init->noise_spectrum->baseline;
         }
         
         if ((opmode == 1) && (strcmp(energy_method,"OPTFILT") == 0) && (strcmp(ofnoise,"WEIGHTM") == 0) 
             && (largeFilter > pow(2,reconstruct_init->noise_spectrum->weightMatrixes->size1)))
         {
             char valueAux[256];
             sprintf(valueAux,"%d",largeFilter);
             string str(valueAux);
             message = "The noise file needs to have W" + str + " (=largeFilter) in the WEIGHTMS HDU";
             str.clear();
             EP_EXIT_ERROR(message,EPFAIL);
         }
     }
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     // Fill in the matrix tstartPulse1_i because tstartPulse1 = nameFile
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     if (!isNumber(tstartPulse1))
     {
         fitsfile *tstartPulse1FileObject = NULL;
         if (fits_open_file(&tstartPulse1FileObject,tstartPulse1,READONLY,status))
         {
             message = "Cannot open tstartPulse1 file " + string(tstartPulse1);
             EP_EXIT_ERROR(message,EPFAIL);
         }
         char extname[20];
         strcpy(extname,"PIXELIMPACT");
         if (fits_movnam_hdu(tstartPulse1FileObject, ANY_HDU,extname, 0, status))
         {
             message = "Cannot move to HDU " + string(extname) + " in tstartPulse1 file " + string(tstartPulse1);
             EP_EXIT_ERROR(message,EPFAIL);
         }
         long totalpulses = 0;
         if (fits_get_num_rows(tstartPulse1FileObject,&totalpulses, status))
         {
             message = "Cannot get number of rows in " + string(tstartPulse1);
             EP_EXIT_ERROR(message,EPFAIL);
         }
         
         reconstruct_init->tstartPulse1_i = gsl_vector_alloc(totalpulses);
         
         char column_name[12];
         int time_colnum = 0;
         strcpy(column_name,"TIME");
         if (fits_get_colnum(tstartPulse1FileObject, CASEINSEN,column_name, &time_colnum, status))
         {
             message = "Cannot get column number for TIME in tstartPulse1 file";
             EP_EXIT_ERROR(message,EPFAIL);
         }
         IOData obj;
         obj.inObject = tstartPulse1FileObject;
         obj.nameTable = new char [255];
         strcpy(obj.nameTable,"PIXELIMPACT");
         obj.iniCol = 0;
         obj.nameCol = new char [255];
         obj.unit = new char [255];
         strcpy(obj.nameCol,"TIME");
         obj.type = TDOUBLE;
         obj.iniRow = 1;
         obj.endRow = totalpulses;
         if (readFitsSimple (obj,&reconstruct_init->tstartPulse1_i))
         {
             message = "Cannot read TIME column in tstartPulse1 file";
             EP_EXIT_ERROR(message,EPFAIL);
         }
         
         delete [] obj.nameTable;
         delete [] obj.nameCol;
         delete [] obj.unit;
         
         if (fits_close_file(tstartPulse1FileObject, status))
         {
             message = "Error closing tstartPulse1 file";
             EP_EXIT_ERROR(message,EPFAIL);
         }
     }
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     // Fill in reconstruct_init
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     strncpy(reconstruct_init->record_file,record_file,255);
     reconstruct_init->record_file[255]='\0';
     reconstruct_init->record_file_fptr = fptr;
     strncpy(reconstruct_init->library_file,library_file,255);
     reconstruct_init->library_file[255]='\0';
     strncpy(reconstruct_init->noise_file,noise_file,255);
     reconstruct_init->noise_file[255]='\0';
     strncpy(reconstruct_init->event_file,event_file,255);
     reconstruct_init->event_file[255]='\0';
     strncpy(reconstruct_init->detectFile,detectFile,255);
     reconstruct_init->detectFile[255]='\0';
     reconstruct_init->threshold 	= 0.0;
     reconstruct_init->pulse_length 	= pulse_length;
     reconstruct_init->scaleFactor  	= scaleFactor;
     reconstruct_init->samplesUp    	= samplesUp;
     reconstruct_init->samplesDown  	= samplesDown;
     reconstruct_init->nSgms        	= nSgms;
     reconstruct_init->detectSP      = detectSP;
     reconstruct_init->opmode	= opmode;
     strcpy(reconstruct_init->detectionMode,detectionMode);
     reconstruct_init->LrsT		= LrsT;
     reconstruct_init->LbT		= LbT;
     reconstruct_init->monoenergy 	= monoenergy;
     if (0 != hduPRECALWN)	reconstruct_init->hduPRECALWN = 1;
     else			reconstruct_init->hduPRECALWN = 0;
     if (0 != hduPRCLOFWM)	reconstruct_init->hduPRCLOFWM = 1;
     else			reconstruct_init->hduPRCLOFWM = 0;
     reconstruct_init->largeFilter 	= largeFilter;
     strcpy(reconstruct_init->FilterDomain,filter_domain);
     strcpy(reconstruct_init->FilterMethod,filter_method);
     strcpy(reconstruct_init->EnergyMethod,energy_method);
     reconstruct_init->filtEev     = filtEev;
     reconstruct_init->Ifit     = Ifit;
     strcpy(reconstruct_init->OFNoise,ofnoise);
     reconstruct_init->LagsOrNot = lagsornot;
     reconstruct_init->nLags = nLags;
     reconstruct_init->Fitting35 = Fitting35;
     reconstruct_init->OFIter = ofiter;
     if (0 != oflib)	reconstruct_init->OFLib = 1;
     else		reconstruct_init->OFLib = 0;
     strcpy(reconstruct_init->OFInterp,ofinterp);
     strcpy(reconstruct_init->OFStrategy,oflength_strategy);
     reconstruct_init->OFLength      = oflength;
     reconstruct_init->errorT = errorT;
     reconstruct_init->Sum0Filt = Sum0Filt;
     reconstruct_init->intermediate  = interm;
     reconstruct_init->SaturationValue  = SaturationValue;
     
     strncpy(reconstruct_init->tstartPulse1,tstartPulse1,255);
     reconstruct_init->tstartPulse1[255]='\0';
     
     if (tstartPulse2 != 0)
     {
         reconstruct_init->tstartPulse2 = tstartPulse2-1;	// To be consistent in the GSL indexes which start from 0
     }
     else			reconstruct_init->tstartPulse2 = tstartPulse2;
     if (tstartPulse3 != 0)
     {
         reconstruct_init->tstartPulse3 = tstartPulse3-1;	// To be consistent in the GSL indexes which start from 0
     }
     else			reconstruct_init->tstartPulse3 = tstartPulse3;
     reconstruct_init->energyPCA1 	= energyPCA1;
     reconstruct_init->energyPCA2 	= energyPCA2;
     strncpy(reconstruct_init->XMLFile,XMLFile,255);
     reconstruct_init->XMLFile[255]='\0'; 
     if (0 != clobber)	reconstruct_init->clobber = 1;
     else		        reconstruct_init->clobber = 0;
     reconstruct_init->maxPulsesPerRecord = maxPulsesPerRecord;
 }

 /*xxxx end of SECTION 1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 2 ************************************************************
  * reconstructRecordSIRENA: This function is the main wrapper function to detect, grade and calculate the energy of the pulses in the input records.
  *
  * - Inititalize PulsesCollection structure
  * - Check consistency of some input parameters
  * - If first record, read the necessary keywords and columns from the input file in order to convert from current to quasi-resistance space
  * - In case of running with threading
  * - Detect pulses in input record (runDetect()). 
  * - If in RECONSTRUCTION (:option:`opmode` = 1) and not PCA:
  *       - Filter and calculate energy of pulses (runEnergy())
  * - Fill in the pulsesAll structure
  * - Populate output event list with pulses energies, arrival time and grading
  *
  * Parameters:
  * - record: Instance of TesRecord structure that contains the input record
  * - trig_reclength: Record size (just in case threading and input files with different 'ADC' lengths but the same record size indeed)
  * - event_list: Instance of TesEventList structure that contains the information of the reconstructed pulses
  * - reconstruct_init: Member of ReconstructInitSIRENA structure to initialize the reconstruction parameters (pointer and values)
  * - lastRecord: If record being analyzed is the last one, lastRecord = 1. Otherwise it is equal to 0
  * - nRecord: Input record number
  * - pulsesAll: Member of PulsesCollection structure to successively store all the pulses used to create the library. 
  *              Re-populated after each processed record.
  * - optimalFilter: Optimal filters used in reconstruction
  * - status:Input/output status
  ******************************************************************************/
 extern "C" void reconstructRecordSIRENA(TesRecord* record, int trig_reclength, TesEventList* event_list, ReconstructInitSIRENA* reconstruct_init,  int lastRecord, int nRecord, PulsesCollection **pulsesAll, int* const status)
 {
     log_trace("reconstructRecordSIRENA: START");

     char extname[20];

     // Inititalize PulsesCollection structure
     PulsesCollection* pulsesInRecord = new PulsesCollection;
     
     PulsesCollection* pulsesAllAux = new PulsesCollection;

     // Check consistency of some input parameters
     if (record->trigger_size <= 0)
     {
         EP_EXIT_ERROR("Record size is <= 0",EPFAIL);
     }
     
     if(reconstruct_init->pulse_length > (int)record->trigger_size)
     {
         //EP_EXIT_ERROR("Warning: pulse length is larger than record size. Pulse length set to maximum value (record size)",EPFAIL);
         /*if (reconstruct_init->opmode == 0)
         {
             EP_EXIT_ERROR("largeFilter is larger than record size",EPFAIL);
         }
         else 
         {
            EP_EXIT_ERROR("Pulse length is larger than record size",EPFAIL);
         }*/
         if (reconstruct_init->opmode == 0) 
         {
             if (reconstruct_init->preBuffer == 0)
                 EP_EXIT_ERROR("largeFilter is larger than record size",EPFAIL);
             //else if (reconstruct_init->preBuffer == 1)
             //    
         }
         else 
         {
            EP_EXIT_ERROR("Pulse length is larger than record size",EPFAIL);
         }
     }
     
     // If first record, read the necessary keywords and columns from the input file in order to convert from current to quasi-resistance space
     if (nRecord == 1)
     {
         reconstruct_init->i2rdata = NULL;
         reconstruct_init->i2rdata = (I2RData*)malloc(sizeof(I2RData));
         reconstruct_init->i2rdata->I0_START = -999;    
         reconstruct_init->i2rdata->IMIN = -999;
         reconstruct_init->i2rdata->IMAX = -999;
         reconstruct_init->i2rdata->ADU_CNV = -999;
         reconstruct_init->i2rdata->I_BIAS = -999;
         reconstruct_init->i2rdata->ADU_BIAS = -999;
         reconstruct_init->i2rdata->Ifit = reconstruct_init->Ifit;
         reconstruct_init->i2rdata->V0 = -999;
         reconstruct_init->i2rdata->RL = -999;
         reconstruct_init->i2rdata->L = -999;
         //if (strcmp(reconstruct_init->EnergyMethod,"OPTFILT") != 0)
         if ((strcmp(reconstruct_init->EnergyMethod,"I2R") == 0) || (strcmp(reconstruct_init->EnergyMethod,"I2RFITTED") == 0) || (strcmp(reconstruct_init->EnergyMethod,"I2RDER") == 0))
         {
             char extname[20];
             char keyname[10];
             
             double IMIN;
             double IMAX;
             
             double ADU_CNV;
             int adu_cnv_exists = 0;
             int i_bias_exists = 0;
             int adu_bias_exists = 0;
             
             strcpy(keyname,"ADU_CNV");
             int hdunum; // Number of the current HDU (RECORDS or TESRECORDS)
             fits_get_hdu_num(reconstruct_init->record_file_fptr, &hdunum);
             strcpy(extname,"TESRECORDS");
             fits_movnam_hdu(reconstruct_init->record_file_fptr, ANY_HDU,extname, 0, status);
             if (((*status != 0) && (hdunum == 2)) || ((*status == 0) && (hdunum > 2)))
             {
                for (int i=0;i<hdunum;i++)
                {
                    fits_movabs_hdu(reconstruct_init->record_file_fptr, i+1, NULL, status);
                    fits_read_key(reconstruct_init->record_file_fptr,TDOUBLE,keyname, &ADU_CNV,NULL,status);
                    if (*status == 0)
                    {
                        adu_cnv_exists = 1;
                        reconstruct_init->i2rdata->ADU_CNV = ADU_CNV;
                        break;
                    }
                    else if ((*status != 0) && (i <= hdunum-1))
                    {
                        *status = 0;
                    }
                }
                if (adu_cnv_exists == 0)
                {
                    EP_EXIT_ERROR("ADU_CNV keyword (to be used in resistance space) is not in the input FITS file",EPFAIL);
                }
                else if (adu_cnv_exists == 1)
                {
                    double I_BIAS;

                    strcpy(keyname,"I_BIAS");
                    for (int i=0;i<hdunum;i++)
                    {
                        fits_movabs_hdu(reconstruct_init->record_file_fptr, i+1, NULL, status);
                        fits_read_key(reconstruct_init->record_file_fptr,TDOUBLE,keyname, &I_BIAS,NULL,status);
                        if (*status == 0)
                        {
                            i_bias_exists = 1;
                            reconstruct_init->i2rdata->I_BIAS = I_BIAS;
                            break;
                        }
                        else if ((*status != 0) && (i <= hdunum-1))
                        {
                            *status = 0;
                        }
                    }
                    if (i_bias_exists == 0)
                    {
                        EP_EXIT_ERROR("I_BIAS keyword (to be used in resistance space) is not in the input FITS file",EPFAIL);
                    }
                    double ADU_BIAS;

                    strcpy(keyname,"ADU_BIAS");
                    for (int i=0;i<hdunum;i++)
                    {
                        fits_movabs_hdu(reconstruct_init->record_file_fptr, i+1, NULL, status);
                        fits_read_key(reconstruct_init->record_file_fptr,TDOUBLE,keyname, &ADU_BIAS,NULL,status);
                        if (*status == 0)
                        {
                            adu_bias_exists = 1;
                            reconstruct_init->i2rdata->ADU_BIAS = ADU_BIAS;
                            break;
                        }
                        else if ((*status != 0) && (i <= hdunum-1))
                        {
                            *status = 0;
                        }
                    }
                    if (adu_bias_exists == 0)
                    {
                        EP_EXIT_ERROR("ADU_BIAS keyword (to be used in resistance space) is not in the input FITS file",EPFAIL);
                    }
                }

                if ((adu_cnv_exists == 0) || (i_bias_exists == 0) || (adu_bias_exists == 0))
                {
                    strcpy(extname,"RECORDS");
                    fits_movnam_hdu(reconstruct_init->record_file_fptr, ANY_HDU,extname, 0, status);
                    if ((*status) != 0)
                    {
                        if (((strcmp(reconstruct_init->EnergyMethod,"I2R") == 0) && ((adu_cnv_exists == 0) || (adu_bias_exists == 0) || (i_bias_exists == 0))) || (strcmp(reconstruct_init->EnergyMethod,"I2R") != 0))
                        {
                            *status = 0;

                            int hdunum; // Number of the current HDU (RECORDS or TESRECORDS)
                            int hdutype;
                            fits_get_hdu_num(reconstruct_init->record_file_fptr, &hdunum);
                            fits_get_hdu_type(reconstruct_init->record_file_fptr, &hdutype, status);

                            strcpy(extname,"ADCPARAM");
                            if (fits_movnam_hdu(reconstruct_init->record_file_fptr, ANY_HDU,extname, 0, status))
                            {
                                EP_EXIT_ERROR("ADU_CNV or I_BIAS or ADU_BIAS not in the input FITS file (to be used in resistance space). Cannot move to ADCPARAM HDU to alternatively look for IMIN and IMAX (old files with RECORDS HDU)",EPFAIL);
                            }
                            strcpy(keyname,"IMIN");
                            fits_read_key(reconstruct_init->record_file_fptr,TDOUBLE,keyname, &IMIN,NULL,status);
                            if (*status != 0)
                            {
                                EP_EXIT_ERROR("ADU_CNV or I_BIAS or ADU_BIAS not in the input FITS file (to be used in resistance space). Cannot alternatively read IMIN keyword (ADCPARAM HDU) to be used in resistance space",EPFAIL);
                            }
                            reconstruct_init->i2rdata->IMIN = IMIN;
                            strcpy(keyname,"IMAX");
                            fits_read_key(reconstruct_init->record_file_fptr,TDOUBLE,keyname, &IMAX,NULL,status);
                            reconstruct_init->i2rdata->IMAX = IMAX;
                            if (*status != 0)
                            {
                                EP_EXIT_ERROR("ADU_CNV or I_BIAS or ADU_BIAS not in the input FITS file (to be used in resistance space). Cannot alternatively read IMAX keyword (ADCPARAM HDU) to be used in resistance space",EPFAIL);
                            }

                            strcpy(extname,"TESPARAM");
                            if (fits_movnam_hdu(reconstruct_init->record_file_fptr, ANY_HDU,extname, 0, status))
                            {
                                EP_EXIT_ERROR("Cannot move to TESPARAM HDU ",EPFAIL);
                            }
                            IOData obj;
                            obj.inObject = reconstruct_init->record_file_fptr;
                            obj.nameTable = new char [255];
                            strcpy(obj.nameTable,"TESPARAM");
                            obj.iniCol = 0;
                            obj.nameCol = new char [255];
                            obj.unit = new char [255];
                            obj.type = TDOUBLE;
                            obj.iniRow = 1;
                            obj.endRow = 1;
                            gsl_vector *vector = gsl_vector_alloc(1);
                            strcpy(obj.nameCol,"I0_START");
                            if (readFitsSimple (obj,&vector))
                            {
                                EP_EXIT_ERROR("Cannot read I0_START column in records input FITS file (TESPARAM HDU)",EPFAIL);
                            }
                            reconstruct_init->i2rdata->I0_START = gsl_vector_get(vector,0);

                            if (fits_movabs_hdu(reconstruct_init->record_file_fptr, hdunum, &hdutype, status))
                            {
                                EP_EXIT_ERROR("Cannot move to RECORDS or TESRECORDS HDU in records input FITS file",EPFAIL);
                            }

                            gsl_vector_free(vector); vector = 0;
                            delete [] obj.nameTable;
                            delete [] obj.nameCol;
                            delete [] obj.unit;
                        }
                    }
                    else
                    {
                        if (((strcmp(reconstruct_init->EnergyMethod,"I2R") == 0) && ((adu_cnv_exists == 0) || (adu_bias_exists == 0) || (i_bias_exists == 0))) || (strcmp(reconstruct_init->EnergyMethod,"I2R") != 0))
                        {
                            double I0_START;

                            strcpy(keyname,"I0_START");
                            if (fits_read_key(reconstruct_init->record_file_fptr,TDOUBLE,keyname, &I0_START,NULL,status))
                                EP_EXIT_ERROR("ADU_CNV or I_BIAS or ADU_BIAS not in the input FITS file (to be used in resistance space). Cannot alternatively read I0_START keyword (TESRECORDS HDU)",EPFAIL);
                            reconstruct_init->i2rdata->I0_START = I0_START;
                            strcpy(keyname,"IMIN");
                            if (fits_read_key(reconstruct_init->record_file_fptr,TDOUBLE,keyname, &IMIN,NULL,status))
                                EP_EXIT_ERROR("ADU_CNV or I_BIAS or ADU_BIAS not in the input FITS file (to be used in resistance space). Cannot alternatively read IMIN keyword (TESRECORDS HDU)",EPFAIL);
                            reconstruct_init->i2rdata->IMIN = IMIN;
                            strcpy(keyname,"IMAX");
                            if (fits_read_key(reconstruct_init->record_file_fptr,TDOUBLE,keyname, &IMAX,NULL,status))
                                EP_EXIT_ERROR("ADU_CNV or I_BIAS or ADU_BIAS not in the input FITS file (to be used in resistance space). Cannot alternatively read IMAX keyword (TESRECORDS HDU)",EPFAIL);
                            reconstruct_init->i2rdata->IMAX = IMAX;
                        }
                    }
                }
             }
         }
         if (strcmp(reconstruct_init->EnergyMethod,"I2RDER") == 0)
         {
             *status = 0;

             char keyname[10];

             int currentHdu_num; // Number of the current HDU (RECORDS or TESRECORDS)
             int currentHdu_type;
             fits_get_hdu_num(reconstruct_init->record_file_fptr, &currentHdu_num);
             fits_get_hdu_type(reconstruct_init->record_file_fptr, &currentHdu_type, status);

             int hdunum; // How many HDU's are in the file
             fits_get_num_hdus(reconstruct_init->record_file_fptr, &hdunum,status);
             int keyword_exists = 0;
             double keywordValue;
             strcpy(keyname,"V0");
             for (int i=0;i<hdunum;i++)
             {
                 fits_movabs_hdu(reconstruct_init->record_file_fptr, i+1, NULL, status);
                 fits_read_key(reconstruct_init->record_file_fptr,TDOUBLE,keyname, &keywordValue,NULL,status);
                 if (*status == 0)
                 {
                     keyword_exists = 1;
                     reconstruct_init->i2rdata->V0 = keywordValue;
                     break;

                 }
                 else if ((*status != 0) && (i <= hdunum-1))
                 {
                    *status = 0;
                 }
             }
             if (keyword_exists == 0)
             {
                EP_EXIT_ERROR("Cannot read V0 keyword to be used in convertI2R",EPFAIL);
             }
             keyword_exists = 0;

             strcpy(keyname,"RL");
             for (int i=0;i<hdunum;i++)
             {
                fits_movabs_hdu(reconstruct_init->record_file_fptr, i+1, NULL, status);
                fits_read_key(reconstruct_init->record_file_fptr,TDOUBLE,keyname, &keywordValue,NULL,status);
                if (*status == 0)
                {
                    keyword_exists = 1;
                    reconstruct_init->i2rdata->RL = keywordValue;
                    break;
                }
                else if ((*status != 0) && (i <= hdunum-1))
                {
                    *status = 0;
                }
             }
             if (keyword_exists == 0)
             {
                EP_EXIT_ERROR("Cannot read RL keyword to be used in convertI2R",EPFAIL);
             }
             keyword_exists = 0;

             strcpy(keyname,"L");
             for (int i=0;i<hdunum;i++)
             {
                fits_movabs_hdu(reconstruct_init->record_file_fptr, i+1, NULL, status);
                fits_read_key(reconstruct_init->record_file_fptr,TDOUBLE,keyname, &keywordValue,NULL,status);
                if (*status == 0)
                {
                    keyword_exists = 1;
                    reconstruct_init->i2rdata->L = keywordValue;

                    break;
                }
                else if ((*status != 0) && (i <= hdunum-1))
                {
                    *status = 0;
                }
             }
             if (keyword_exists == 0)
             {
                EP_EXIT_ERROR("Cannot read L keyword to be used in convertI2R",EPFAIL);
             }

             if (fits_movabs_hdu(reconstruct_init->record_file_fptr, currentHdu_num, &currentHdu_type, status))
             {
                EP_EXIT_ERROR("Cannot move to RECORDS or TESRECORDS HDU ",EPFAIL);
             }
         }
     }
     int tessimOrxifusim = -999;
     // Check if input FITS file have been simulated with TESSIM or XIFUSIM
     strcpy(extname,"RECORDS");
     fits_movnam_hdu(reconstruct_init->record_file_fptr, ANY_HDU,extname, 0, status);
     if (*status != 0)
     {
         *status = 0;
         strcpy(extname,"TESRECORDS");
         if (fits_movnam_hdu(reconstruct_init->record_file_fptr, ANY_HDU,extname, 0, status))
         {
             EP_EXIT_ERROR("'TESRECORDS' HDU is not in the input FITS file",EPFAIL);
         }
         else
         {
             tessimOrxifusim = 1;
         }
     }
     else 
     {
         tessimOrxifusim = 0;
     }
     if (tessimOrxifusim == -999)
     {
         EP_EXIT_ERROR("Neither the 'RECORDS' nor 'TESRECORDS' HDUs are in the input FITS file",EPFAIL);
     }
     
     // In case of running with threading
     if (scheduler::get()->is_threading() 
         && reconstruct_init->opmode == 1
         && (strcmp(reconstruct_init->EnergyMethod, "PCA") != 0))
     {
         log_trace("reconstructRecordSIRENA:  Threading mode...");
         ReconstructInitSIRENA* rec = reconstruct_init->get_threading_object(nRecord);
         log_trace("reconstructRecordSIRENA:  Threading mode...1");
         scheduler::get()->push_detection(record, trig_reclength, nRecord, lastRecord, 
                                          *pulsesAll, &rec, &pulsesInRecord,
                                          event_list);
         freeReconstructInitSIRENA(rec);
         log_trace("reconstructRecordSIRENA:  Threading mode...2");
         return;  // The rest of 'reconstructRecordSIRENA' is not going to run: 'runDetect', 'runEnergy'...
     }
         
     // Detect pulses in record
     log_trace("Before runDetect");
     if (!scheduler::get()->is_threading())
         trig_reclength = record->trigger_size;
     runDetect(record, trig_reclength, lastRecord, nRecord, *pulsesAll, &reconstruct_init, &pulsesInRecord);
     log_trace("After runDetect");
     
     //cout<<"pulsesInRecord->ndetpulses: "<<pulsesInRecord->ndetpulses<<endl;
     //cout<<"(*pulsesAll)->ndetpulses: "<<(*pulsesAll)->ndetpulses<<endl;
     
     if ((reconstruct_init->opmode == 1) && (strcmp(reconstruct_init->EnergyMethod,"PCA") != 0))
     {
         // Filter and calculates energy
         runEnergy(record, lastRecord, nRecord, trig_reclength, &reconstruct_init, &pulsesInRecord, *pulsesAll);
     }
     log_trace("After runEnergy");
         
     // Fill in the pulsesAll structure
     if ((pulsesInRecord->ndetpulses != 0) && ((*pulsesAll)->ndetpulses == 0))
     {
         (*pulsesAll)->ndetpulses = pulsesInRecord->ndetpulses;
         if((*pulsesAll)->pulses_detected != 0 && (*pulsesAll)->size < pulsesInRecord->ndetpulses){
             delete [] (*pulsesAll)->pulses_detected; (*pulsesAll)->pulses_detected = 0;
             (*pulsesAll)->size = resize_array((*pulsesAll)->size, (*pulsesAll)->ndetpulses);
             (*pulsesAll)->pulses_detected = new PulseDetected[(*pulsesAll)->size];
         }
         
         #ifndef POOLS		
         if((*pulsesAll)->pulses_detected == 0)
         {			
             (*pulsesAll)->pulses_detected = new PulseDetected[pulsesInRecord->ndetpulses];
             (*pulsesAll)->size = pulsesInRecord->ndetpulses;
         }
         #endif
         for (int i=0;i<(*pulsesAll)->ndetpulses;i++)
         {
             (*pulsesAll)->pulses_detected[i] = pulsesInRecord->pulses_detected[i];
         }            
     }
     else
     {
         pulsesAllAux->ndetpulses = (*pulsesAll)->ndetpulses;
         
         (*pulsesAll)->ndetpulses = (*pulsesAll)->ndetpulses + pulsesInRecord->ndetpulses;
         
         if ((*pulsesAll)->pulses_detected != NULL && (*pulsesAll)->size < (*pulsesAll)->ndetpulses)
         {		
             pulsesAllAux->pulses_detected = new PulseDetected[(*pulsesAll)->ndetpulses];
             
             for (int i=0;i<pulsesAllAux->ndetpulses;i++){
                 pulsesAllAux->pulses_detected[i] = (*pulsesAll)->pulses_detected[i];
             }
             
             delete [] (*pulsesAll)->pulses_detected; (*pulsesAll)->pulses_detected = 0; 
             (*pulsesAll)->size = resize_array((*pulsesAll)->size, (*pulsesAll)->ndetpulses);								
             (*pulsesAll)->pulses_detected = new PulseDetected[(*pulsesAll)->size];			
             
             for (int i=0;i<pulsesAllAux->ndetpulses;i++){
                 (*pulsesAll)->pulses_detected[i] = pulsesAllAux->pulses_detected[i];
             }
             delete [] pulsesAllAux->pulses_detected; pulsesAllAux->pulses_detected = 0;
         }
         
         #ifndef POOLS		
         if((*pulsesAll)->pulses_detected == 0){	
             (*pulsesAll)->pulses_detected = new PulseDetected[(*pulsesAll)->ndetpulses];
             (*pulsesAll)->size = (*pulsesAll)->ndetpulses;
         }
         #endif	
         
         // Save pulses detected in current record
         for (int i=0;i<pulsesInRecord->ndetpulses;i++)
         {
             (*pulsesAll)->pulses_detected[i+pulsesAllAux->ndetpulses] = pulsesInRecord->pulses_detected[i];
         }
     }
     
     log_trace("pulsesAll: %i",(*pulsesAll)->ndetpulses);
     //cout<<"pulsesAll: "<<(*pulsesAll)->ndetpulses<<endl;
     //cout<<"pulsesInRecord: "<<pulsesInRecord->ndetpulses<<endl;
     
     // Free & Fill TesEventList structure
     event_list->index = pulsesInRecord->ndetpulses;
     event_list->energies = new double[event_list->index];
     event_list->avgs_4samplesDerivative = new double[event_list->index];
     event_list->Es_lowres = new double[event_list->index];
     event_list->phis = new double[event_list->index];
     event_list->lagsShifts = new int[event_list->index];
     event_list->bsln = new double[event_list->index];
     event_list->rmsbsln = new double[event_list->index];
     event_list->grading = new int[event_list->index];
     event_list->grades2 = new int[event_list->index];
     event_list->ph_ids = new long[event_list->index];
     event_list->ph_ids2 = new long[event_list->index];
     event_list->ph_ids3 = new long[event_list->index];
     event_list->pix_ids = new long[event_list->index];
     event_list->tends = new double[event_list->index];
     event_list->tstarts = new double[event_list->index];
     event_list->risetimes = new double[event_list->index];
     event_list->falltimes = new double[event_list->index];
     
     if (strcmp(reconstruct_init->EnergyMethod,"PCA") != 0)     // Different from PCA
     {
         for (int ip=0; ip<pulsesInRecord->ndetpulses; ip++)
         {	  			
             event_list->event_indexes[ip] = (pulsesInRecord->pulses_detected[ip].Tstart - record->time)/record->delta_t;	
             
             if (reconstruct_init->opmode == 1)
             {
                 event_list->energies[ip] = pulsesInRecord->pulses_detected[ip].energy;
             }
             else if (reconstruct_init->opmode == 0)
             {
                 event_list->energies[ip] = 999.;
             }

             event_list->avgs_4samplesDerivative[ip] = pulsesInRecord->pulses_detected[ip].avg_4samplesDerivative;
             event_list->Es_lowres[ip] = pulsesInRecord->pulses_detected[ip].E_lowres;
             event_list->phis[ip] = pulsesInRecord->pulses_detected[ip].phi;
             event_list->lagsShifts[ip] = pulsesInRecord->pulses_detected[ip].lagsShift;
             event_list->bsln[ip] = pulsesInRecord->pulses_detected[ip].bsln;
             event_list->rmsbsln[ip] = pulsesInRecord->pulses_detected[ip].rmsbsln;
             event_list->grading[ip]  = pulsesInRecord->pulses_detected[ip].grading;
             event_list->grades1[ip]  = pulsesInRecord->pulses_detected[ip].grade1;
             event_list->grades2[ip]  = pulsesInRecord->pulses_detected[ip].grade2;
             event_list->pulse_heights[ip]  = pulsesInRecord->pulses_detected[ip].pulse_height;
             event_list->pix_ids[ip]  = pulsesInRecord->pulses_detected[ip].pixid;
             event_list->ph_ids[ip]  = pulsesInRecord->pulses_detected[ip].phid;
             event_list->ph_ids2[ip]  = pulsesInRecord->pulses_detected[ip].phid2;
             event_list->ph_ids3[ip]  = pulsesInRecord->pulses_detected[ip].phid3;
             event_list->tstarts[ip]  = pulsesInRecord->pulses_detected[ip].Tstart;
             event_list->tends[ip]  = pulsesInRecord->pulses_detected[ip].Tend;
             event_list->risetimes[ip]  = pulsesInRecord->pulses_detected[ip].riseTime;
             event_list->falltimes[ip]  = pulsesInRecord->pulses_detected[ip].fallTime;
         }
         
         if (lastRecord == 1)
         {       
             double numLagsUsed_mean;
             double numLagsUsed_sigma;
             gsl_vector *numLagsUsed_vector = gsl_vector_alloc((*pulsesAll)->ndetpulses);
             //gsl_vector_set_all(numLagsUsed_vector,-999); // Debugger complains about an initialized variable in a conditional jump
                                                          // (in findMeanSigma 6 lines forward)

             for (int ip=0; ip<(*pulsesAll)->ndetpulses; ip++)
             {
                 gsl_vector_set(numLagsUsed_vector,ip,(*pulsesAll)->pulses_detected[ip].numLagsUsed);
             }
             if (findMeanSigma (numLagsUsed_vector, &numLagsUsed_mean, &numLagsUsed_sigma))
             {
                 EP_EXIT_ERROR("Cannot run findMeanSigma routine for calculating numLagsUsed statistics",EPFAIL);
             }
             if (numLagsUsed_vector != NULL) {gsl_vector_free(numLagsUsed_vector); numLagsUsed_vector = 0;}
         }
     }
     else
     {
         if (lastRecord == 1)
         {
             // Free & Fill TesEventList structure
             event_list->index = (*pulsesAll)->ndetpulses;
             //event_list->event_indexes = new double[event_list->index];
             event_list->energies = new double[event_list->index];
             event_list->avgs_4samplesDerivative = new double[event_list->index];
             event_list->Es_lowres = new double[event_list->index];
             event_list->phis = new double[event_list->index];
             event_list->lagsShifts = new int[event_list->index];
             event_list->bsln = new double[event_list->index];
             event_list->rmsbsln = new double[event_list->index];
             event_list->grading  = new int[event_list->index];
             //event_list->grades1  = new int[event_list->index];
             event_list->grades2  = new int[event_list->index];
             //event_list->pulse_heights  = new double[event_list->index];
             event_list->ph_ids   = new long[event_list->index];
             event_list->pix_ids   = new long[event_list->index];
             event_list->risetimes   = new double[event_list->index];
             event_list->falltimes   = new double[event_list->index];
             
             for (int ip=0; ip<(*pulsesAll)->ndetpulses; ip++)
             {	
                 event_list->event_indexes[ip] = ((*pulsesAll)->pulses_detected[ip].Tstart - record->time)/record->delta_t;
                 
                 if (reconstruct_init->opmode == 1)
                 {
                     event_list->energies[ip] = (*pulsesAll)->pulses_detected[ip].energy;
                 }
                 else if (reconstruct_init->opmode == 0)
                 {
                     event_list->energies[ip] = 999.;
                 }
                 
                 event_list->avgs_4samplesDerivative[ip]  = (*pulsesAll)->pulses_detected[ip].avg_4samplesDerivative;
                 event_list->Es_lowres[ip]  = (*pulsesAll)->pulses_detected[ip].E_lowres;
                 event_list->phis[ip] = (*pulsesAll)->pulses_detected[ip].phi;
                 event_list->lagsShifts[ip] = (*pulsesAll)->pulses_detected[ip].lagsShift;
                 event_list->bsln[ip] = (*pulsesAll)->pulses_detected[ip].bsln;
                 event_list->rmsbsln[ip] = (*pulsesAll)->pulses_detected[ip].rmsbsln;
                 event_list->grading[ip]  = (*pulsesAll)->pulses_detected[ip].grading;
                 event_list->grades1[ip]  = (*pulsesAll)->pulses_detected[ip].grade1;
                 event_list->grades2[ip]  = (*pulsesAll)->pulses_detected[ip].grade2;
                 event_list->pulse_heights[ip]  = (*pulsesAll)->pulses_detected[ip].pulse_height;
                 event_list->ph_ids[ip]   = (*pulsesAll)->pulses_detected[ip].phid;		    
                 event_list->pix_ids[ip]  = (*pulsesAll)->pulses_detected[ip].pixid;
                 event_list->risetimes[ip]  = (*pulsesAll)->pulses_detected[ip].riseTime;
                 event_list->falltimes[ip]  = (*pulsesAll)->pulses_detected[ip].fallTime;
             }
         }
     }
     
     delete pulsesAllAux; pulsesAllAux = 0;
     delete [] pulsesInRecord->pulses_detected; pulsesInRecord->pulses_detected = 0;
     delete pulsesInRecord; pulsesInRecord = 0;
     
     return;
 }
 /*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

 
 /***** SECTION 3 ************************************************************
  * ReconstructInitSIRENA: Constructor. It returns a pointer to an empty ReconstructInitSIRENA data structure.
  * 
  * - Initialize pointers with NULL for SIRENA
  * 
  * Parameters:
  * - status: Input/output status
  ******************************************************************************/
 extern "C" ReconstructInitSIRENA* newReconstructInitSIRENA()
 {	
     ReconstructInitSIRENA* reconstruct_init = new ReconstructInitSIRENA;
     reconstruct_init->noise_spectrum = NULL;
     
     return(reconstruct_init);
 }
 /*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 4 ************************************************************
  * freeReconstructInitSIRENA: Destructor of ReconstructInit structure
  * 
  * Parameters:
  * - reconstruct_init: Member of ReconstructInitSIRENA structure to initialize the reconstruction parameters (pointer and values)
  ******************************************************************************/
 extern "C" void freeReconstructInitSIRENA(ReconstructInitSIRENA* reconstruct_init)
 {
     if (reconstruct_init->grading != NULL)
     {
         if (reconstruct_init->grading->value != NULL)
         {
             gsl_vector_free(reconstruct_init->grading->value); reconstruct_init->grading->value = 0;
         }
         if (reconstruct_init->grading->gradeData != NULL)
         {
             gsl_matrix_free(reconstruct_init->grading->gradeData); reconstruct_init->grading->gradeData = 0;
         }
         free(reconstruct_init->grading);
         reconstruct_init->grading = 0;
     }
     if (reconstruct_init->i2rdata != NULL)
     {
         free(reconstruct_init->i2rdata); reconstruct_init->i2rdata = NULL;
     }

     delete(reconstruct_init); reconstruct_init = 0;
 }
 /*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 5 ************************************************************
  * newPulsesCollection: Constructor. It returns a pointer to an empty PulsesCollection data structure.
  * 
  * Parameters:
  * - status: Input/output status
  ******************************************************************************/
 extern "C" PulsesCollection* newPulsesCollection()
 {
     PulsesCollection* PulsesColl = new PulsesCollection;
     
     #ifndef POOLS
     // Initialize pointers with NULL for SIRENA
     PulsesColl->pulses_detected =0;
     
     // Initialize values for SIRENA
     PulsesColl->ndetpulses=0;
     PulsesColl->size = 0;
     #else
     // Initialize pointers with NULL for SIRENA
     PulsesColl->pulses_detected = new PulseDetected[POOL_SIZE];
     
     // Initialize values for SIRENA
     PulsesColl->ndetpulses=0;
     PulsesColl->size = POOL_SIZE;

     PulsesColl->pulses_detected->Tstart = -999;

     #endif
     
     return(PulsesColl);
 }
 /*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 6 ************************************************************
  * freePulsesCollection: Destructor of PulsesCollection structure.
  * 
  * Parameters:
  * - PulsesColl: Instance of PulsesCollection structure
  ******************************************************************************/
 extern "C" void freePulsesCollection(PulsesCollection* PulsesColl)
 {
     #if 0
     for(int i = 0; i < PulsesColl->ndetpulses; ++i){
         if (PulsesColl->pulses_detected[i].pulse_adc != NULL) gsl_vector_free(PulsesColl->pulses_detected[i].pulse_adc);
 }
 delete [] PulsesColl->pulses_detected;
 #endif
 delete(PulsesColl); PulsesColl = 0;
 }
 /*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 7 ************************************************************
  * newOptimalFilterSIRENA: Constructor. It returns a pointer to an empty OptimalFilterSIRENA data structure.
  * 
  * Parameters:
  * - status: Input/output status
  ******************************************************************************/
 extern "C" OptimalFilterSIRENA* newOptimalFilterSIRENA()
 {
     OptimalFilterSIRENA* OFilterColl = new OptimalFilterSIRENA;
     return(OFilterColl);
 }
 /*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 8 ************************************************************
  * freeOptimalFilterSIRENA: Destructor of OptimalFilterSIRENA structure
  * 
  * Parameters:
  * - OFilterColl: Instance of OptimalFilterSIRENA structure
  ******************************************************************************/
 extern "C" void freeOptimalFilterSIRENA(OptimalFilterSIRENA* OFilterColl)
 {
     delete(OFilterColl); OFilterColl = 0;
 }
 /*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 9 ************************************************************
  * getLibraryCollection: This funtion creates and retrieves a LibraryCollection from a file.
  * 
  * - Create LibraryCollection structure
  * - Open FITS file in READONLY mode (move to the first HDU) and get number of templates (rows)
  * - Allocate library structure
  * - Get PULSE and MF column numbers (depending on the different options)
  * - Get template duration
  * - Allocate library structure (cont.)
  * - Get matched filter duration
  * - Read different columns and populate the *LibraryCollection* structure
  * - Added new code to handle the new HDUs FIXFILTF, FIXFILTT, PRECALWN and PRCLOFWM
  * - Free allocated GSL vectors and matrices
  * 
  * Parameters:
  * - filename: File with library information
  * - opmode: Calibration run (0) or energy reconstruction run (1)
  * - hduPRECALWN: Add or not the PRECALWN HDU in the library file (1/0)
  * - hduPRCLOFWM: Add or not the PRCLOFWM HDU in the library file (1/0)
  * - largeFilter: Length of the longest fixed filter  
  * - filter_domain: Time domain ('T') or Frequency domain ('F')
  * - pulse_length: Pulse length
  * - energy_method: Energy calculation Method: OPTFILT, WEIGHT, WEIGHTN, I2R, I2RFITTED or PCA
  * - ofnoise: For optimal filtering, NSD or WEIGHTM
  * - filter_method: Filtering Method: F0 (deleting the zero frequency bin) or B0 (deleting the baseline) or F0B0 (deleting always the baseline)
  * - oflib: Work or not with a library with optimal filters (1/0)
  * - ofinterp: Optimal Filter by using the Matched Filter or the DAB as matched filter (MF/DAB) 
  * 	      It has been fixed in 'tesreconstruction' as 'DAB' (but it would be possible to work with 'MF')
  * - filtEeV: Energy of the filters of the library to be used to calculate energy
  * - lagsornot: Lags (1) or no lags (0)
  * - preBuffer: Using preBuffer (1) or not using preBuffer (0)
  * - pBi: Vector with the preBuffer values read from the XML file
  * - posti: Vector with the post values read from the XML file
  * - status: Input/output status
  ******************************************************************************/
 LibraryCollection* getLibraryCollection(const char* const filename, int opmode, int hduPRECALWN, int hduPRCLOFWM, int largeFilter, char* filter_domain, int pulse_length, char *energy_method, char *ofnoise, char *filter_method, char oflib, char **ofinterp, double filtEev, int lagsornot, int preBuffer, gsl_vector *pBi, gsl_vector *posti, int* const status)
 {
     string message = "";

     // Create LibraryCollection structure
     LibraryCollection* library_collection = new LibraryCollection;
     
     // Open FITS file in READONLY mode
     fitsfile* fptr = NULL;
     if (fits_open_file(&fptr, filename, READONLY, status))
     {
         EP_PRINT_ERROR("Error opening library file",*status);
         return(library_collection);
     }
     
     // Move to the first HDU
     int extver = 1;
     char HDUname[12];
     strcpy(HDUname,"LIBRARY");
     if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
     {
         EP_PRINT_ERROR("Error moving to HDU LIBRARY in library file",*status);
         return(library_collection);
     }
     
     library_collection->baseline = -999.0;
     if ((strcmp(filter_method,"B0") == 0) || (strcmp(energy_method,"WEIGHT") == 0))
     {
         if (fits_read_key(fptr,TDOUBLE,"BASELINE", &library_collection->baseline,NULL,status))
         {
             EP_PRINT_ERROR("Cannot read BASELINE keyword from HDU LIBRARY in library file => Check the noise file",-999);
             *status = 0;
         }
     }
     
     // Get number of templates (rows)
     long ntemplates;
     if (fits_get_num_rows(fptr,&ntemplates, status))
     {
         EP_PRINT_ERROR("Cannot get number of rows in library file",*status);
         return(library_collection);
     }
     if (ntemplates == 0)	
     {
         EP_PRINT_ERROR("The library has no rows",EPFAIL); 
         *status = EPFAIL; return(library_collection);
     }
     
     if ((opmode == 1) && 
         (((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0))))
     {
         if (ntemplates == 1)
         {
             //if (strcmp(*ofinterp,"DAB") == 0)  *ofinterp = "MF";
             if (strcmp(*ofinterp,"DAB") == 0)  *ofinterp = (char *)("MF");
         }
         else
         {
             if (filtEev != 0)
             {
                 //if (strcmp(*ofinterp,"DAB") == 0)  *ofinterp = "MF";
                 if (strcmp(*ofinterp,"DAB") == 0)  *ofinterp = (char *)("MF");
                 
                 EP_PRINT_ERROR("The library has several rows, but only the row related to filtEev is going to be used in reconstruction",-999); // Only a warning
             }
             else if ((filtEev == 0) && (lagsornot == 1))
             {
                 EP_PRINT_ERROR("filtEev=0 (filters interpolation) and LagsOrNot=1 is not developed yet => Please, change your choice to LagsOrNot=0",EPFAIL); 
                 *status = EPFAIL; return(library_collection);
             }
         }
     }
     
     // Allocate library structure
     // It is not necessary to check the allocation because 'ntemplates' has been checked previously
     library_collection->ntemplates                = (int)ntemplates;
     library_collection->energies           	      = gsl_vector_alloc(ntemplates);
     library_collection->pulse_heights      	      = gsl_vector_alloc(ntemplates);
     library_collection->maxDERs      	      = gsl_vector_alloc(ntemplates);
     library_collection->samp1DERs      	      = gsl_vector_alloc(ntemplates);
     library_collection->pulse_templates           = new PulseTemplate[ntemplates];
     library_collection->pulse_templatesMaxLengthFixedFilter = new PulseTemplate[ntemplates];
     library_collection->pulse_templates_filder    = new PulseTemplate[ntemplates];
     library_collection->pulse_templates_B0        = new PulseTemplate[ntemplates];
     library_collection->matched_filters           = new MatchedFilter[ntemplates];
     library_collection->matched_filters_B0        = new MatchedFilter[ntemplates];
     library_collection->optimal_filters           = new OptimalFilterSIRENA[ntemplates];
     library_collection->optimal_filtersab         = new OptimalFilterSIRENA[ntemplates];
     library_collection->optimal_filtersFREQ       = new OptimalFilterSIRENA[ntemplates];
     library_collection->optimal_filtersTIME       = new OptimalFilterSIRENA[ntemplates];
     library_collection->optimal_filtersabFREQ     = new OptimalFilterSIRENA[ntemplates];
     library_collection->optimal_filtersabTIME     = new OptimalFilterSIRENA[ntemplates];
     
     // Get PULSE and MF column numbers (depending the different options)
     char column_name[12];
     int template_colnum = 0;
     int mfilter_colnum = 0;
     
     strcpy(column_name,"PULSE");
     if (fits_get_colnum(fptr, CASEINSEN,column_name, &template_colnum, status))
     {
         EP_PRINT_ERROR("Cannot get column number for PULSE in library file",*status);
         *status=EPFAIL; return(library_collection);
     }

     if ((opmode == 0) ||
     (((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)) && (oflib == 0) && (strcmp(*ofinterp,"MF") == 0) && (opmode == 1)))
     {
         strcpy(column_name,"MF");
         if (fits_get_colnum(fptr, CASEINSEN,column_name, &mfilter_colnum, status))
         {
             EP_PRINT_ERROR("Cannot get column number for MF in library file",*status);
             *status=EPFAIL; return(library_collection);
         }
     }
     
     // Get template duration
     int naxis = 0;
     long naxes = 0;
     if (fits_read_tdim(fptr, template_colnum, 1, &naxis, &naxes, status))
     {
         EP_PRINT_ERROR("Cannot read dim of column PULSE",*status);
         *status=EPFAIL; return(library_collection);
     }
     int template_duration = naxes;
     if (template_duration == 0)	
     {
         EP_PRINT_ERROR("PULSE column vectors length is 0",EPFAIL); 
         *status = EPFAIL; return(library_collection);
     }
     
     if (preBuffer == 0)
     {
        if ((opmode == 0) && (pulse_length != template_duration))
        {
            EP_PRINT_ERROR("It is not possible the PulseLength provided because it does not match with the PulseLength (PULSE_column_length) of the existing library",EPFAIL);
            *status=EPFAIL; return(library_collection);
        }
     }
     else   // preBuffer = 1
     {
        if ((opmode == 0) && (gsl_vector_max(posti) != template_duration))
        {
            EP_PRINT_ERROR("It is not possible the post and pB values provided in the XML file because it does not match with the template duration of the existing library",EPFAIL);
            *status=EPFAIL; return(library_collection);
        }
     }
     
     int template_durationPLSMXLFF = -999;
     int ncols;
     if (fits_get_num_cols(fptr,&ncols, status))
     {
         EP_PRINT_ERROR("Cannot get number of cols in library file",*status);
         return(library_collection);
     }
     
     if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
     {
         int template_colnumPLSMXLFF = 0;
         strcpy(column_name,"PLSMXLFF");
         if (fits_get_colnum(fptr, CASEINSEN,column_name, &template_colnumPLSMXLFF, status))
         {
             EP_PRINT_ERROR("Cannot get column number for PLSMXLFF in library file",*status);
             *status=EPFAIL; return(library_collection);
         }
         if (fits_read_tdim(fptr, template_colnumPLSMXLFF, 1, &naxis, &naxes, status))
         {
             EP_PRINT_ERROR("Cannot read dim of column PLSMXLFF",*status);
             *status=EPFAIL; return(library_collection);
         }
         template_durationPLSMXLFF = naxes;
         if (template_durationPLSMXLFF == 0)	
         {
             EP_PRINT_ERROR("PLSMXLFF column vectors length is 0",EPFAIL); 
             *status=EPFAIL; return(library_collection);
         }
         
         if ((opmode == 0) && (largeFilter != template_durationPLSMXLFF))
         {
             EP_PRINT_ERROR("It is not possible the largeFilter provided because it does not match with the largeFilter (PLSMXLFF_column_length) of the existing library",EPFAIL);
             *status=EPFAIL; return(library_collection);
         }
     }
     else template_durationPLSMXLFF = 999;
     
     // Allocate library structure (cont.)
     // It is not necessary to check the allocation because 'ntemplates' has been checked previously
     if ((((strcmp(energy_method,"WEIGHT") == 0) || (strcmp(energy_method,"WEIGHTN") == 0)) && (opmode == 1)) || ((opmode == 0) && (hduPRECALWN == 1)))
     {
         library_collection->V = gsl_matrix_alloc(ntemplates,template_duration*template_duration);
         library_collection->W = gsl_matrix_alloc(ntemplates,template_duration*template_duration);
         if (ntemplates > 1)
         {
             library_collection->WAB = gsl_matrix_alloc(ntemplates-1,template_duration*template_duration);
             library_collection->T = gsl_matrix_alloc(ntemplates-1,template_duration);
             library_collection->t = gsl_vector_alloc(ntemplates-1);
             library_collection->X = gsl_matrix_alloc(ntemplates-1,template_duration*template_duration);
             library_collection->Y = gsl_matrix_alloc(ntemplates-1,template_duration);
             library_collection->Z = gsl_matrix_alloc(ntemplates-1,template_duration);
             library_collection->r = gsl_vector_alloc(ntemplates-1);
         }
         else
         {
             library_collection->WAB = gsl_matrix_alloc(1,template_duration*template_duration);
             library_collection->T = gsl_matrix_alloc(1,template_duration);
             library_collection->t = gsl_vector_alloc(1);
             library_collection->X = gsl_matrix_alloc(1,template_duration*template_duration);
             library_collection->Y = gsl_matrix_alloc(1,template_duration);
             library_collection->Z = gsl_matrix_alloc(1,template_duration);
             library_collection->r = gsl_vector_alloc(1);
         }
     }
     if (ntemplates > 1)
     {
         library_collection->PAB = gsl_matrix_alloc(ntemplates-1,template_duration);
         library_collection->PABMXLFF = gsl_matrix_alloc(ntemplates-1,template_durationPLSMXLFF);
         library_collection->DAB = gsl_matrix_alloc(ntemplates-1,template_duration);
     }
     else
     {
         library_collection->PAB = gsl_matrix_alloc(1,template_duration);
         library_collection->PABMXLFF = gsl_matrix_alloc(1,template_durationPLSMXLFF);
         library_collection->DAB = gsl_matrix_alloc(1,template_duration);
     }
     
     // Get matched filter duration
     int mfilter_duration = -999;
     if ((opmode == 0) ||
     (((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0)|| (strcmp(energy_method,"I2RDER") == 0)) && (oflib == 0) && (strcmp(*ofinterp,"MF") == 0)))
     {
         if (fits_read_tdim(fptr, mfilter_colnum, 1, &naxis, &naxes, status))
         {
             EP_PRINT_ERROR("Cannot read dim of column MF",*status);
             *status=EPFAIL; return(library_collection);
         }
         mfilter_duration = naxes;
         if (mfilter_duration == 0) 	
         {
             EP_PRINT_ERROR("MF column vectors length is 0",EPFAIL); 
             return(library_collection);
         }
         
         //Check that TEMPLATE & MFILTER have the same duration
         if(template_duration != mfilter_duration)
         {
             EP_PRINT_ERROR("Error: template and matched filter have different durations",EPFAIL);
             *status=EPFAIL; return(library_collection);
         }
     }
     
     // Read different columns and populate the *LibraryCollection* structure
     // Read ENERGY column
     IOData obj;
     obj.inObject = fptr;
     obj.nameTable = new char [255];
     strcpy(obj.nameTable,"LIBRARY");
     obj.iniCol = 0;
     obj.nameCol = new char [255];
     obj.unit = new char [255];
     strcpy(obj.nameCol,"ENERGY");
     obj.type = TDOUBLE;
     obj.iniRow = 1;
     obj.endRow = ntemplates;
     if (readFitsSimple (obj,&library_collection->energies))
     {
         EP_PRINT_ERROR("Cannot read ENERGY column in library FITS file",*status);
         *status=EPFAIL; return(library_collection);
     }
     if ((opmode == 1) && (oflib == 1))
     {
         int filtEevIsAEnergy = 0;
         if (filtEev != 0)
         {
             for (int i=0;i<ntemplates;i++)
             {
                 if (filtEev == gsl_vector_get(library_collection->energies,i))
                 {
                     filtEevIsAEnergy = 1;
                     break;
                 }
             }
         }
         if ((filtEev !=0) && (filtEevIsAEnergy == 0) )
         {
             EP_EXIT_ERROR("filtEv is not one of the energies in the library",EPFAIL);    
         }
     }
     
     // Read PHEIGHT column
     strcpy(obj.nameCol,"PHEIGHT");
     if (readFitsSimple (obj,&library_collection->pulse_heights))
     {
         EP_PRINT_ERROR("Cannot read PHEIGHT column in library FITS file",*status);
         *status=EPFAIL; return(library_collection);
     }
     
     // It is not necessary to check the allocation because 'ntemplates' and 'template_duration' have been checked previously
     gsl_matrix *matrixAux_PULSEMaxLengthFixedFilter = gsl_matrix_alloc(ntemplates,template_durationPLSMXLFF);
     if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
     {
         strcpy(obj.nameCol,"PLSMXLFF");
         if (readFitsComplex (obj,&matrixAux_PULSEMaxLengthFixedFilter))
         {
             EP_PRINT_ERROR("Cannot read PLSMXLFF column in library FITS file",*status);
             *status=EPFAIL; return(library_collection);
         }
     }
     
     // It is not necessary to check the allocation because 'ntemplates' and 'template_duration' have been checked previously
     gsl_matrix *matrixAux_PULSE = gsl_matrix_alloc(ntemplates,template_duration);
     strcpy(obj.nameCol,"PULSE");
     if (readFitsComplex (obj,&matrixAux_PULSE))
     {
         EP_PRINT_ERROR("Cannot read PULSE column in library FITS file",*status);
         *status=EPFAIL; return(library_collection);
     }
     
     // It is not necessary to check the allocation because 'ntemplates' and 'template_duration' have been checked previously
     gsl_matrix *matrixAux_PULSEB0 = gsl_matrix_alloc(ntemplates,template_duration);
     strcpy(obj.nameCol,"PULSEB0");
     if (readFitsComplex (obj,&matrixAux_PULSEB0))
     {
         EP_PRINT_ERROR("Cannot read PULSEB0 column in library FITS file",*status);
         *status=EPFAIL; return(library_collection);
     }
     
     gsl_matrix *matrixAux_MF = NULL;
     gsl_matrix *matrixAux_MFB0 = NULL;
     if ((opmode == 0) ||
     (((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)) && (oflib == 0) && (strcmp(*ofinterp,"MF") == 0) && (opmode == 1)))
     {
         if ((opmode == 0) || ((strcmp(filter_method,"F0") == 0) || (strcmp(filter_method,"F0B0") == 0)))
         {
             // It is not necessary to check the allocation because 'ntemplates' and 'mfilter_duration' have been checked previously
             matrixAux_MF = gsl_matrix_alloc(ntemplates,mfilter_duration);
             strcpy(obj.nameCol,"MF");
             if (readFitsComplex (obj,&matrixAux_MF))
             {
                 EP_PRINT_ERROR("Cannot read MF column in library FITS file",*status);
                 *status=EPFAIL; return(library_collection);
             }
         }
         if ((opmode == 0) || (strcmp(filter_method,"B0") == 0))
         {
             // It is not necessary to check the allocation because 'ntemplates' and 'mfilter_duration' have been checked previously
             matrixAux_MFB0 = gsl_matrix_alloc(ntemplates,mfilter_duration);
             strcpy(obj.nameCol,"MFB0");
             if (readFitsComplex (obj,&matrixAux_MFB0))
             {
                 EP_PRINT_ERROR("Cannot read MFB0 column in library FITS file",*status);
                 *status=EPFAIL; return(library_collection);
             }  
         }  
     }
     
     gsl_matrix *matrixAux_WAB = NULL;
     gsl_vector *vectorAux_WAB = NULL;
     gsl_matrix *matrixAux_V = NULL;
     gsl_vector *vectorAux_V = NULL;
     gsl_matrix *matrixAux_W = NULL;
     gsl_vector *vectorAux_W = NULL;
     gsl_matrix *matrixAux_T = NULL;
     gsl_vector *vectorAux_T = NULL;
     gsl_vector *vectorAux_t = NULL;
     gsl_matrix *matrixAux_X = NULL;
     gsl_vector *vectorAux_X = NULL;
     gsl_matrix *matrixAux_Y = NULL;
     gsl_vector *vectorAux_Y = NULL;
     gsl_matrix *matrixAux_Z = NULL;
     gsl_vector *vectorAux_Z = NULL;
     gsl_vector *vectorAux_r = NULL;
     gsl_matrix *matrixAux_PAB = NULL;
     gsl_matrix *matrixAux_PABMaxLengthFixedFilter = NULL;
     gsl_vector *vectorAux_PAB = NULL;
     gsl_vector *vectorAux_PABMaxLengthFixedFilter = NULL;
     gsl_matrix *matrixAux_DAB = NULL;
     gsl_vector *vectorAux_DAB = NULL;
     
     int dim;
     if ((((strcmp(energy_method,"WEIGHT") == 0) || (strcmp(energy_method,"WEIGHTN") == 0)) && (opmode == 1)) || ((opmode == 0) && (hduPRECALWN == 1)))
     {
         // It is not necessary to check the allocation because 'ntemplates' and 'template_duration' have been checked previously
         matrixAux_V = gsl_matrix_alloc(ntemplates,template_duration*template_duration);
         vectorAux_V = gsl_vector_alloc(template_duration*template_duration);
         strcpy(obj.nameCol,"COVARM");
         if (readFitsComplex (obj,&matrixAux_V))
         {
             EP_PRINT_ERROR("Cannot read COVARM column in library FITS file",*status);
             *status=EPFAIL; return(library_collection);
         }
         
         // It is not necessary to check the allocation because 'ntemplates' and 'ofilter_durationab' have been checked previously
         matrixAux_W = gsl_matrix_alloc(ntemplates,template_duration*template_duration);
         vectorAux_W = gsl_vector_alloc(template_duration*template_duration);
         strcpy(obj.nameCol,"WEIGHTM");
         if (readFitsComplex (obj,&matrixAux_W))
         {
             EP_PRINT_ERROR("Cannot read WEIGHTM column in library FITS file",*status);
             *status=EPFAIL; return(library_collection);
         }
         
         if ((opmode == 1) || ((opmode == 0) && (ntemplates >1)))
         {
             if (ntemplates == 1)
             {
                 obj.endRow = 1;
                 dim = 1;
             }
             else
             {
                 obj.endRow = ntemplates-1;
                 dim = ntemplates-1;
             }
             
             if (((opmode == 1) && (strcmp(energy_method,"WEIGHT") == 0)) || ((opmode == 0) && (ntemplates >1)))
             {
                 // It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
                 matrixAux_T = gsl_matrix_alloc(dim,template_duration);
                 vectorAux_T = gsl_vector_alloc(template_duration);
                 strcpy(obj.nameCol,"TV");
                 if (readFitsComplex (obj,&matrixAux_T))
                 {
                     EP_PRINT_ERROR("Cannot read TV column in library FITS file",*status);
                     *status=EPFAIL; return(library_collection);
                 }
                 
                 // It is not necessary to check the allocation because dim > 0
                 vectorAux_t = gsl_vector_alloc(dim);
                 strcpy(obj.nameCol,"tE");
                 if (readFitsSimple (obj,&vectorAux_t))
                 {
                     EP_PRINT_ERROR("Cannot read tE column in library FITS file",*status);
                     *status=EPFAIL; return(library_collection);
                 }
                 
                 // It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
                 matrixAux_X = gsl_matrix_alloc(dim,template_duration*template_duration);
                 vectorAux_X = gsl_vector_alloc(template_duration*template_duration);
                 strcpy(obj.nameCol,"XM");
                 if (readFitsComplex (obj,&matrixAux_X))
                 {
                     EP_PRINT_ERROR("Cannot read XM column in library FITS file",*status);
                     *status=EPFAIL; return(library_collection);
                 }
                 
                 // It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
                 matrixAux_Y = gsl_matrix_alloc(dim,template_duration);
                 vectorAux_Y = gsl_vector_alloc(template_duration);
                 strcpy(obj.nameCol,"YV");
                 if (readFitsComplex (obj,&matrixAux_Y))
                 {
                     EP_PRINT_ERROR("Cannot read YV column in library FITS file",*status);
                     *status=EPFAIL; return(library_collection);
                 }
                 
                 // It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
                 matrixAux_Z = gsl_matrix_alloc(dim,template_duration);
                 vectorAux_Z = gsl_vector_alloc(template_duration);
                 strcpy(obj.nameCol,"ZV");
                 if (readFitsComplex (obj,&matrixAux_Z))
                 {
                     EP_PRINT_ERROR("Cannot read ZV column in library FITS file",*status);
                     *status=EPFAIL; return(library_collection);
                 }
                 
                 // It is not necessary to check the allocation because dim > 0
                 vectorAux_r = gsl_vector_alloc(dim);
                 strcpy(obj.nameCol,"rE");
                 if (readFitsSimple (obj,&vectorAux_r))
                 {
                     EP_PRINT_ERROR("Cannot read rE column in library FITS file",*status);
                     *status=EPFAIL; return(library_collection);
                 }
             }
             
             if (((opmode == 1) && (strcmp(energy_method,"WEIGHTN") == 0)) || ((opmode == 0) && (ntemplates >1)))
             {
                 // It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
                 matrixAux_WAB = gsl_matrix_alloc(dim,template_duration*template_duration);
                 vectorAux_WAB = gsl_vector_alloc(template_duration*template_duration);
                 strcpy(obj.nameCol,"WAB");
                 if (readFitsComplex (obj,&matrixAux_WAB))
                 {
                     EP_PRINT_ERROR("Cannot read WAB column in library FITS file",*status);
                     *status=EPFAIL; return(library_collection);
                 }
             }
         }
     }
     
     if (((opmode == 1) && ((strcmp(energy_method,"WEIGHTN") == 0) || (((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)) && (strcmp(*ofinterp,"DAB") == 0) && (strcmp(ofnoise,"NSD") == 0)))) || ((opmode == 0) && (ntemplates >1)))
     {
         if (ntemplates == 1)
         {
             obj.endRow = 1;
             dim = 1;
         }
         else
         {
             obj.endRow = ntemplates-1;
             dim = ntemplates-1;
         }
         // It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
         matrixAux_PAB = gsl_matrix_alloc(dim,template_duration);
         vectorAux_PAB = gsl_vector_alloc(template_duration);
         strcpy(obj.nameCol,"PAB");
         if (readFitsComplex (obj,&matrixAux_PAB))
         {
             EP_PRINT_ERROR("Cannot read PAB column in library FITS file",*status);
             *status=EPFAIL; return(library_collection);
         }
         
         matrixAux_PABMaxLengthFixedFilter = gsl_matrix_alloc(dim,template_durationPLSMXLFF);
         vectorAux_PABMaxLengthFixedFilter = gsl_vector_alloc(template_durationPLSMXLFF);
         if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
         {
             strcpy(obj.nameCol,"PABMXLFF");
             if (readFitsComplex (obj,&matrixAux_PABMaxLengthFixedFilter))
             {
                 EP_PRINT_ERROR("Cannot read PABMXLFF column in library FITS file",*status);
                 *status=EPFAIL; return(library_collection);
             }
         }
         
         // It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
         matrixAux_DAB = gsl_matrix_alloc(dim,template_duration);
         vectorAux_DAB = gsl_vector_alloc(template_duration);
         strcpy(obj.nameCol,"DAB");
         if (readFitsComplex (obj,&matrixAux_DAB))
         {
             EP_PRINT_ERROR("Cannot read DAB column in library FITS file",*status);
             *status=EPFAIL; return(library_collection);
         }
     }
     
     for (int it = 0 ; it < ntemplates ; it++)
     {
         // It is not necessary to check the allocation because 'template_duration' has been checked previously
         library_collection->pulse_templates[it].ptemplate    		= gsl_vector_alloc(template_duration);
         library_collection->pulse_templates_filder[it].ptemplate   = gsl_vector_alloc(template_duration);
         library_collection->pulse_templates_B0[it].ptemplate 		= gsl_vector_alloc(template_duration);
         library_collection->matched_filters[it].mfilter      		= gsl_vector_alloc(template_duration);
         library_collection->matched_filters_B0[it].mfilter   		= gsl_vector_alloc(template_duration);

         library_collection->pulse_templatesMaxLengthFixedFilter[it].template_duration = template_durationPLSMXLFF;
         library_collection->pulse_templates[it].template_duration		= template_duration;
         library_collection->pulse_templates_filder[it].template_duration    	= -1;
         library_collection->pulse_templates_B0[it].template_duration 		= template_duration;
         library_collection->matched_filters[it].mfilter_duration     		= template_duration;
         library_collection->matched_filters_B0[it].mfilter_duration  		= template_duration;
         
         library_collection->pulse_templatesMaxLengthFixedFilter[it].energy = gsl_vector_get(library_collection->energies,it);
         library_collection->pulse_templates[it].energy    	= gsl_vector_get(library_collection->energies,it);
         library_collection->pulse_templates_filder[it].energy	= gsl_vector_get(library_collection->energies,it);
         library_collection->pulse_templates_B0[it].energy 	= gsl_vector_get(library_collection->energies,it);
         library_collection->matched_filters[it].energy    	= gsl_vector_get(library_collection->energies,it);
         library_collection->matched_filters_B0[it].energy 	= gsl_vector_get(library_collection->energies,it);
         
         if (((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19)) && ((opmode == 0) ||
         ((opmode == 1) && ((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) ||
         (strcmp(energy_method,"I2RDER") == 0)))))
         {
             library_collection->pulse_templatesMaxLengthFixedFilter[it].ptemplate    = gsl_vector_alloc(template_durationPLSMXLFF);
             gsl_matrix_get_row(library_collection->pulse_templatesMaxLengthFixedFilter[it].ptemplate,matrixAux_PULSEMaxLengthFixedFilter,it);
         }
         
         gsl_matrix_get_row(library_collection->pulse_templates[it].ptemplate,matrixAux_PULSE,it);
         
         gsl_matrix_get_row(library_collection->pulse_templates_B0[it].ptemplate,matrixAux_PULSEB0,it);
         
         if ((opmode == 0) || 
         (((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)) && (oflib == 0) && (strcmp(*ofinterp,"MF") == 0)))
         {
             if ((opmode == 0) || ((strcmp(filter_method,"F0") == 0) || (strcmp(filter_method,"F0B0") == 0)))
             {
                 gsl_matrix_get_row(library_collection->matched_filters[it].mfilter,matrixAux_MF,it);
             }
             if ((opmode == 0) || (strcmp(filter_method,"B0") == 0))
             {
                 gsl_matrix_get_row(library_collection->matched_filters_B0[it].mfilter,matrixAux_MFB0,it);
             }
         }
         
         if (((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)) && (oflib == 0) && (strcmp(*ofinterp,"DAB") == 0))
         {
             if ((opmode == 1)  && (it < ntemplates-1))
             {
                 gsl_matrix_get_row(vectorAux_DAB,matrixAux_DAB,it);
                 gsl_vector_memcpy(library_collection->matched_filters[it].mfilter,vectorAux_DAB);
             }
         }
         
         if (((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)) && 
         (strcmp(*ofinterp,"DAB") == 0) && (opmode == 1))
         {
             if ((opmode == 1)  && (it < ntemplates-1))
             {
                 gsl_matrix_get_row(vectorAux_PAB,matrixAux_PAB,it);
                 gsl_matrix_set_row(library_collection->PAB,it,vectorAux_PAB);
                 
                 if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                 {
                     gsl_matrix_get_row(vectorAux_PABMaxLengthFixedFilter,matrixAux_PABMaxLengthFixedFilter,it);
                     gsl_matrix_set_row(library_collection->PABMXLFF,it,vectorAux_PABMaxLengthFixedFilter);
                 }
                 
                 gsl_matrix_get_row(vectorAux_DAB,matrixAux_DAB,it);
                 gsl_matrix_set_row(library_collection->DAB,it,vectorAux_DAB);
             }
         }
         
         if (((opmode == 0) && (hduPRECALWN == 1)) || ((strcmp(energy_method,"WEIGHT") == 0) || (strcmp(energy_method,"WEIGHTN") == 0)))
         {
             gsl_matrix_get_row(vectorAux_V,matrixAux_V,it);
             gsl_matrix_set_row(library_collection->V,it,vectorAux_V);
             
             gsl_matrix_get_row(vectorAux_W,matrixAux_W,it);
             gsl_matrix_set_row(library_collection->W,it,vectorAux_W);
         }
         
         if (((opmode == 1) && (it < ntemplates-1) && (strcmp(energy_method,"WEIGHT") == 0)) ||
             ((opmode == 0) && (hduPRECALWN == 1) && (ntemplates > 1) && (it < ntemplates-1)))
         {
             gsl_matrix_get_row(vectorAux_T,matrixAux_T,it);
             gsl_matrix_set_row(library_collection->T,it,vectorAux_T);
             
             gsl_vector_set(library_collection->t,it,gsl_vector_get(vectorAux_t,it));
             
             gsl_matrix_get_row(vectorAux_X,matrixAux_X,it);
             gsl_matrix_set_row(library_collection->X,it,vectorAux_X);
             
             gsl_matrix_get_row(vectorAux_Y,matrixAux_Y,it);
             gsl_matrix_set_row(library_collection->Y,it,vectorAux_Y);
             
             gsl_matrix_get_row(vectorAux_Z,matrixAux_Z,it);
             gsl_matrix_set_row(library_collection->Z,it,vectorAux_Z);
             
             gsl_vector_set(library_collection->r,it,gsl_vector_get(vectorAux_r,it));
         }	
         
         if (((opmode == 1)  && (it < ntemplates-1) && (strcmp(energy_method,"WEIGHTN") == 0)) ||
             ((opmode == 0) && (hduPRECALWN == 1) && (ntemplates > 1) && (it < ntemplates-1)))
         {
             gsl_matrix_get_row(vectorAux_WAB,matrixAux_WAB,it);
             gsl_matrix_set_row(library_collection->WAB,it,vectorAux_WAB);
         }
         
         if (((opmode == 1)  && (it < ntemplates-1) &&(strcmp(energy_method,"WEIGHTN") == 0)) ||
             ((opmode == 0) && (ntemplates > 1) && (it < ntemplates-1)))
         {
             
             
             gsl_matrix_get_row(vectorAux_PAB,matrixAux_PAB,it);
             gsl_matrix_set_row(library_collection->PAB,it,vectorAux_PAB);
             
             if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
             {
                 gsl_matrix_get_row(vectorAux_PABMaxLengthFixedFilter,matrixAux_PABMaxLengthFixedFilter,it);
                 gsl_matrix_set_row(library_collection->PABMXLFF,it,vectorAux_PABMaxLengthFixedFilter);
             }
             
             gsl_matrix_get_row(vectorAux_DAB,matrixAux_DAB,it);
             gsl_matrix_set_row(library_collection->DAB,it,vectorAux_DAB);
         }
     }
     
     // Free allocated GSL vectors and matrices
     gsl_matrix_free(matrixAux_PULSE); matrixAux_PULSE = 0;
     gsl_matrix_free(matrixAux_PULSEMaxLengthFixedFilter); matrixAux_PULSEMaxLengthFixedFilter = 0;
     gsl_matrix_free(matrixAux_PULSEB0); matrixAux_PULSEB0 = 0;
     if (matrixAux_MF != NULL) {gsl_matrix_free(matrixAux_MF); matrixAux_MF = 0;}
     if (matrixAux_MFB0 != NULL) {gsl_matrix_free(matrixAux_MFB0); matrixAux_MFB0 = 0;}
     if (matrixAux_V != NULL) {gsl_matrix_free(matrixAux_V); matrixAux_V = 0;}
     if (vectorAux_V != NULL) {gsl_vector_free(vectorAux_V); vectorAux_V = 0;}
     if (matrixAux_W != NULL) {gsl_matrix_free(matrixAux_W); matrixAux_W = 0;}
     if (vectorAux_W != NULL) {gsl_vector_free(vectorAux_W); vectorAux_W = 0;}
     if (matrixAux_WAB != NULL) {gsl_matrix_free(matrixAux_WAB); matrixAux_WAB = 0;}
     if (vectorAux_WAB != NULL) {gsl_vector_free(vectorAux_WAB); vectorAux_WAB = 0;}
     if (matrixAux_T != NULL) {gsl_matrix_free(matrixAux_T); matrixAux_T = 0;}
     if (vectorAux_T != NULL) {gsl_vector_free(vectorAux_T); vectorAux_T = 0;}
     if (vectorAux_t != NULL) {gsl_vector_free(vectorAux_t); vectorAux_t = 0;}
     if (matrixAux_X != NULL) {gsl_matrix_free(matrixAux_X); matrixAux_X = 0;}
     if (vectorAux_X != NULL) {gsl_vector_free(vectorAux_X); vectorAux_X = 0;}
     if (matrixAux_Y != NULL) {gsl_matrix_free(matrixAux_Y); matrixAux_Y = 0;}
     if (vectorAux_Y != NULL) {gsl_vector_free(vectorAux_Y); vectorAux_Y = 0;}
     if (matrixAux_Z != NULL) {gsl_matrix_free(matrixAux_Z); matrixAux_Z = 0;}
     if (vectorAux_Z != NULL) {gsl_vector_free(vectorAux_Z); vectorAux_Z = 0;}
     if (vectorAux_r != NULL) {gsl_vector_free(vectorAux_r); vectorAux_r = 0;}
     if (matrixAux_PAB != NULL) {gsl_matrix_free(matrixAux_PAB); matrixAux_PAB = 0;}
     if (vectorAux_PAB != NULL) {gsl_vector_free(vectorAux_PAB); vectorAux_PAB = 0;}
     if (matrixAux_PABMaxLengthFixedFilter != NULL) {gsl_matrix_free(matrixAux_PABMaxLengthFixedFilter); matrixAux_PABMaxLengthFixedFilter = 0;}
     if (vectorAux_PABMaxLengthFixedFilter != NULL) {gsl_vector_free(vectorAux_PABMaxLengthFixedFilter); vectorAux_PABMaxLengthFixedFilter = 0;}
     if (matrixAux_DAB != NULL) {gsl_matrix_free(matrixAux_DAB); matrixAux_DAB = 0;}
     if (vectorAux_DAB != NULL) {gsl_vector_free(vectorAux_DAB); vectorAux_DAB = 0;}
     
     // Added new code to handle the new HDUs FIXFILTF, FIXFILTT, PRECALWN and PRCLOFWM
     //if (opmode == 0)
     if ((opmode == 0) && (preBuffer == 0)) 
     {
         // FIXFILTF HDU
         strcpy(HDUname,"FIXFILTF");
         if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
         {
             EP_PRINT_ERROR("Error moving to HDU FIXFILTF in library file",*status);
             return(library_collection);
         }
         
         // Get number of fixed optimal filters (columnss)
         int nOFs;
         if (fits_get_num_cols(fptr,&nOFs, status))
         {
             EP_PRINT_ERROR("Cannot get number of rows in library file",*status);
             return(library_collection);
         }
         if (ntemplates == 1)	nOFs = nOFs-1;		// -1 because the ENERGYcolumn
         else 			        nOFs = (nOFs-1)/2;	// -1 because the ENERGYcolumn and /2 because the AB column
         
         if (nOFs == 0)	
         {
             EP_PRINT_ERROR("The library has no fixed optimal filters",EPFAIL); 
             *status = EPFAIL; return(library_collection);
         }
         library_collection->nfixedfilters = nOFs;
         
         int lengthALL_F = 0;
         int lengthALL_T = 0;
         int lengthALL_PRCLWN = 0;
         int lengthALL_PRCLOFWM = 0;
         for (int i=0;i<nOFs;i++)
         {
             if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
             {
                 if (i==0)
                 {
                     lengthALL_F = largeFilter*2;
                     lengthALL_T = largeFilter;
                 }
                 else
                 {
                     lengthALL_F = lengthALL_F + pow(2,floor(log2(pulse_length))-i+1)*2;
                     lengthALL_T = lengthALL_T + pow(2,floor(log2(pulse_length))-i+1);
                 }
             }
             else
             {
                 lengthALL_F = lengthALL_F + pow(2,floor(log2(pulse_length))-i)*2;
                 lengthALL_T = lengthALL_T + pow(2,floor(log2(pulse_length))-i);
             }
         }
         if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))            lengthALL_PRCLWN = lengthALL_F-largeFilter*2;
         else              					                                            lengthALL_PRCLWN = lengthALL_F;
         
         lengthALL_PRCLOFWM = lengthALL_F;
         
         gsl_matrix *matrixALL_OFFx = gsl_matrix_alloc(ntemplates,lengthALL_F);
         gsl_matrix *matrixALL_OFTx = gsl_matrix_alloc(ntemplates,lengthALL_T);
         gsl_matrix *matrixALLab_OFFx = gsl_matrix_alloc(ntemplates,lengthALL_F);
         gsl_matrix *matrixALLab_OFTx = gsl_matrix_alloc(ntemplates,lengthALL_T);
         gsl_matrix *matrixALL_PRCLWNx = gsl_matrix_alloc(ntemplates,lengthALL_PRCLWN);
         gsl_matrix *matrixALL_PRCLOFWMx = gsl_matrix_alloc(ntemplates,lengthALL_PRCLOFWM);
         
         char str_length[125];
         
         gsl_matrix *matrixAux_OFFx = NULL;
         gsl_matrix *matrixAuxab_OFFx = NULL;
         int index = 0;
         strcpy(obj.nameTable,"FIXFILTF");
         obj.iniRow = 1;
         obj.endRow = ntemplates;
         for (int i=0;i<nOFs;i++)
         {
             if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
             {
                 if (i==0)
                 {
                     snprintf(str_length,125,"%d",largeFilter);
                     matrixAux_OFFx = gsl_matrix_alloc(ntemplates,largeFilter*2);
                 }
                 else
                 {	
                     snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(pulse_length))-i+1)));	//-1 because of the largeFilter-length filter
                     matrixAux_OFFx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i+1)*2);
                 }
             }
             else
             {
                 snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(pulse_length))-i)));	
                 matrixAux_OFFx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i)*2);
             }
             strcpy(obj.nameCol,(string("F")+string(str_length)).c_str());
             if (readFitsComplex (obj,&matrixAux_OFFx))
             {
                 message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                 EP_PRINT_ERROR(message,*status);
                 *status=EPFAIL; return(library_collection);
             }
             for (int j=0;j<(int)(matrixAux_OFFx->size1);j++)
             {
                 for (int k=0;k<(int)(matrixAux_OFFx->size2);k++)
                 {
                     gsl_matrix_set(matrixALL_OFFx,j,k+index,gsl_matrix_get(matrixAux_OFFx,j,k));
                 }
             }
             
             if (ntemplates > 1)
             {
                 strcpy(obj.nameCol,(string("ABF")+string(str_length)).c_str());
                 if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                 {
                     if (i==0)	matrixAuxab_OFFx = gsl_matrix_alloc(ntemplates,largeFilter*2);
                     else 		matrixAuxab_OFFx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i+1)*2);
                 }
                 else	matrixAuxab_OFFx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i)*2);
                 if (readFitsComplex (obj,&matrixAuxab_OFFx))
                 {
                     message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                     EP_PRINT_ERROR(message,*status);
                     *status=EPFAIL; return(library_collection);
                 }
                 for (int j=0;j<(int)(matrixAuxab_OFFx->size1);j++)
                 {
                     for (int k=0;k<(int)(matrixAuxab_OFFx->size2);k++)
                     {
                         gsl_matrix_set(matrixALLab_OFFx,j,k+index,gsl_matrix_get(matrixAuxab_OFFx,j,k));
                     }
                 }
             }
             
             if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
             {
                 if (i==0) 	index = index + largeFilter*2;
                 else 		index = index + pow(2,floor(log2(pulse_length))-i+1)*2;
             }
             else	index = index + pow(2,floor(log2(pulse_length))-i)*2;
             
             gsl_matrix_free(matrixAux_OFFx); matrixAux_OFFx = 0;
             gsl_matrix_free(matrixAuxab_OFFx); matrixAuxab_OFFx = 0;
         }
         
         // FIXFILTT HDU
         strcpy(HDUname,"FIXFILTT");
         if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
         {
             EP_PRINT_ERROR("Error moving to HDU FIXFILTT in library file",*status);
             return(library_collection);
         }
         
         gsl_matrix *matrixAux_OFTx = NULL;
         gsl_matrix *matrixAuxab_OFTx = NULL;
         index = 0;
         strcpy(obj.nameTable,"FIXFILTT");
         for (int i=0;i<nOFs;i++)
         {
             if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
             {
                 if (i==0)
                 {
                     snprintf(str_length,125,"%d",largeFilter);
                     matrixAux_OFTx = gsl_matrix_alloc(ntemplates,largeFilter);
                 }
                 else
                 {
                     snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(pulse_length))-i+1)));
                     matrixAux_OFTx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i+1));
                 }
             }
             else
             {
                 snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(pulse_length))-i)));
                 matrixAux_OFTx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i));
             }
             strcpy(obj.nameCol,(string("T")+string(str_length)).c_str());
             if (readFitsComplex (obj,&matrixAux_OFTx))
             {
                 message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                 EP_PRINT_ERROR(message,*status);
                 *status=EPFAIL; return(library_collection);
             }
             for (int j=0;j<(int)(matrixAux_OFTx->size1);j++)
             {
                 for (int k=0;k<(int)(matrixAux_OFTx->size2);k++)
                 {
                     gsl_matrix_set(matrixALL_OFTx,j,k+index,gsl_matrix_get(matrixAux_OFTx,j,k));
                 }
             }
             
             if (ntemplates > 1)
             {
                 strcpy(obj.nameCol,(string("ABT")+string(str_length)).c_str());
                 if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                 {
                     if (i==0)	matrixAuxab_OFTx = gsl_matrix_alloc(ntemplates,largeFilter);
                     else 		matrixAuxab_OFTx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i+1));
                 }
                 else	matrixAuxab_OFTx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i));
                 if (readFitsComplex (obj,&matrixAuxab_OFTx))
                 {
                     message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                     EP_PRINT_ERROR(message,*status);
                     *status=EPFAIL; return(library_collection);
                 }
                 for (int j=0;j<(int)(matrixAuxab_OFTx->size1);j++)
                 {
                     for (int k=0;k<(int)(matrixAuxab_OFTx->size2);k++)
                     {
                         gsl_matrix_set(matrixALLab_OFTx,j,k+index,gsl_matrix_get(matrixAuxab_OFTx,j,k));
                     }
                 }
             }
             
             if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
             {
                 if (i==0)	index = index + largeFilter;
                 else		index = index + pow(2,floor(log2(pulse_length))-i+1);
             }
             else 	index = index + pow(2,floor(log2(pulse_length))-i);
             
             gsl_matrix_free(matrixAux_OFTx); matrixAux_OFTx = 0;
             gsl_matrix_free(matrixAuxab_OFTx); matrixAuxab_OFTx = 0;
         }
         for (int it=0;it<ntemplates;it++)
         {
             library_collection->optimal_filtersFREQ[it].energy		= gsl_vector_get(library_collection->energies,it);
             library_collection->optimal_filtersFREQ[it].ofilter_duration	= lengthALL_F;
             library_collection->optimal_filtersFREQ[it].ofilter    		= gsl_vector_alloc(lengthALL_F);
             
             gsl_matrix_get_row(library_collection->optimal_filtersFREQ[it].ofilter,matrixALL_OFFx,it);
             
             library_collection->optimal_filtersTIME[it].energy		= gsl_vector_get(library_collection->energies,it);
             library_collection->optimal_filtersTIME[it].ofilter_duration 	= lengthALL_T;
             library_collection->optimal_filtersTIME[it].ofilter    		= gsl_vector_alloc(lengthALL_T);
             
             gsl_matrix_get_row(library_collection->optimal_filtersTIME[it].ofilter,matrixALL_OFTx,it);
             
             if (it < ntemplates-1)
             {
                 library_collection->optimal_filtersabFREQ[it].energy		= gsl_vector_get(library_collection->energies,it);
                 library_collection->optimal_filtersabFREQ[it].ofilter_duration  = lengthALL_F;
                 library_collection->optimal_filtersabFREQ[it].ofilter    	= gsl_vector_alloc(lengthALL_F);
                 
                 gsl_matrix_get_row(library_collection->optimal_filtersabFREQ[it].ofilter,matrixALLab_OFFx,it);
                 
                 library_collection->optimal_filtersabTIME[it].energy		= gsl_vector_get(library_collection->energies,it);
                 library_collection->optimal_filtersabTIME[it].ofilter_duration	= lengthALL_F;
                 library_collection->optimal_filtersabTIME[it].ofilter    	= gsl_vector_alloc(lengthALL_T);
                 
                 gsl_matrix_get_row(library_collection->optimal_filtersabTIME[it].ofilter,matrixALLab_OFTx,it);
             }
         }
         
         strcpy(HDUname,"PRECALWN");
         if (hduPRECALWN == 1)
         {
             if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
             {
                 if (*status != 0)
                 {
                     EP_PRINT_ERROR("The input hduPRECALWN parameter does not match the existing library",1);
                     return(library_collection);
                 }
             }
         }
         if ((ntemplates > 1) && (hduPRECALWN == 1)) 
         {
             // PRECALWN HDU
             if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
             {
                 EP_PRINT_ERROR("Error moving to HDU PRECALWN in library file",*status);
                 return(library_collection);
             }
             
             if (fits_get_num_cols(fptr,&nOFs, status))
             {
                 EP_PRINT_ERROR("Cannot get number of rows in library file",*status);
                 return(library_collection);
             }
             nOFs = nOFs-1;		// -1 because the ENERGY column 
             gsl_matrix *matrixAux_PRCLWNx = NULL;
             index = 0;
             strcpy(obj.nameTable,"PRECALWN");
             for (int i=0;i<nOFs;i++)
             {
                 snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(pulse_length))-i)));
                 strcpy(obj.nameCol,(string("PCL")+string(str_length)).c_str());
                 matrixAux_PRCLWNx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i)*2);
                 if (readFitsComplex (obj,&matrixAux_PRCLWNx))
                 {
                     message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                     EP_PRINT_ERROR(message,*status);
                     *status=EPFAIL; return(library_collection);
                 }
                 for (int j=0;j<(int)(matrixAux_PRCLWNx->size1);j++)
                 {
                     for (int k=0;k<(int)(matrixAux_PRCLWNx->size2);k++)
                     {
                         gsl_matrix_set(matrixALL_PRCLWNx,j,k+index,gsl_matrix_get(matrixAux_PRCLWNx,j,k));
                     }
                 }
                 
                 index = index + pow(2,floor(log2(pulse_length))-i)*2;
                 
                 gsl_matrix_free(matrixAux_PRCLWNx); matrixAux_PRCLWNx = 0;
             }
             
             gsl_vector *vectorAux_PRCLWNx = gsl_vector_alloc(lengthALL_PRCLWN);
             library_collection->PRECALWN = gsl_matrix_alloc(ntemplates,lengthALL_PRCLWN);
             for (int it=0;it<ntemplates;it++)
             {
                 if (it < ntemplates-1)
                 {
                     gsl_matrix_get_row(vectorAux_PRCLWNx,matrixALL_PRCLWNx,it);
                     gsl_matrix_set_row(library_collection->PRECALWN,it,vectorAux_PRCLWNx);
                 }
             }
             gsl_vector_free(vectorAux_PRCLWNx); vectorAux_PRCLWNx = 0;
         }
         
         gsl_matrix_free(matrixALL_OFFx); matrixALL_OFFx = 0;
         gsl_matrix_free(matrixALL_OFTx); matrixALL_OFTx = 0;
         gsl_matrix_free(matrixALLab_OFFx); matrixALLab_OFFx = 0;
         gsl_matrix_free(matrixALLab_OFTx); matrixALLab_OFTx = 0;
         gsl_matrix_free(matrixALL_PRCLWNx); matrixALL_PRCLWNx = 0;
         
         // PRCLOFWM HDU
         strcpy(HDUname,"PRCLOFWM");
         if (hduPRCLOFWM == 1)
         {
             if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))	
             {
                 if (*status != 0)
                 {
                     EP_PRINT_ERROR("The input hduPRCLOFWM parameter does not match the existing library",1);
                     return(library_collection);
                 }
             }
         }
         if (hduPRCLOFWM == 1)
         {
             if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
             {
                 EP_PRINT_ERROR("Error moving to HDU PRCLOFWM in library file",*status);
                 return(library_collection);
             }
             
             if (fits_get_num_cols(fptr,&nOFs, status))
             {
                 EP_PRINT_ERROR("Cannot get number of rows in library file",*status);
                 return(library_collection);
             }
             nOFs = nOFs-1;		// -1 because the ENERGY column 
             gsl_matrix *matrixAux_PRCLOFWMx = NULL;
             index = 0;
             strcpy(obj.nameTable,"PRCLOFWM");
             for (int i=0;i<nOFs;i++)
             {
                 if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                 {
                     if (i == 0)
                     {
                         snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(largeFilter))-i)));
                         matrixAux_PRCLOFWMx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(largeFilter))-i)*2);
                     }
                     else
                     {
                         snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(pulse_length))-i+1)));
                         matrixAux_PRCLOFWMx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i+1)*2);
                     }
                 }
                 else
                 {
                     snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(pulse_length))-i)));
                     matrixAux_PRCLOFWMx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i)*2);
                 }
                 strcpy(obj.nameCol,(string("OFW")+string(str_length)).c_str());
                 if (readFitsComplex (obj,&matrixAux_PRCLOFWMx))
                 {
                     message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                     EP_PRINT_ERROR(message,*status);
                     *status=EPFAIL; return(library_collection);
                 }
                 
                 for (int j=0;j<(int)(matrixAux_PRCLOFWMx->size1);j++)
                 {
                     for (int k=0;k<(int)(matrixAux_PRCLOFWMx->size2);k++)
                     {
                         gsl_matrix_set(matrixALL_PRCLOFWMx,j,k+index,gsl_matrix_get(matrixAux_PRCLOFWMx,j,k));
                     }
                 }
                 
                 if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                 {
                     if (i==0)	index = index + pow(2,floor(log2(largeFilter))-i)*2;
                     else		index = index + pow(2,floor(log2(pulse_length))-i+1)*2;
                 }
                 else 	index = index + pow(2,floor(log2(pulse_length))-i)*2;
                 
                 gsl_matrix_free(matrixAux_PRCLOFWMx); matrixAux_PRCLOFWMx = 0;
             }
             
             gsl_vector *vectorAux_PRCLOFWMx = gsl_vector_alloc(lengthALL_PRCLOFWM);
             library_collection->PRCLOFWM = gsl_matrix_alloc(ntemplates,lengthALL_PRCLOFWM);
             for (int it=0;it<ntemplates;it++)
             {
                 gsl_matrix_get_row(vectorAux_PRCLOFWMx,matrixALL_PRCLOFWMx,it);
                 gsl_matrix_set_row(library_collection->PRCLOFWM,it,vectorAux_PRCLOFWMx);
             }
             gsl_vector_free(vectorAux_PRCLOFWMx); vectorAux_PRCLOFWMx = 0;
         }
         gsl_matrix_free(matrixALL_PRCLOFWMx); matrixALL_PRCLOFWMx = 0;
     }
     else if ((opmode == 0) && (preBuffer == 1))
     {
         // FIXFILTF HDU
         strcpy(HDUname,"FIXFILTF");
         if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
         {
             EP_PRINT_ERROR("Error moving to HDU FIXFILTF in library file",*status);
             return(library_collection);
         }
         
         library_collection->nfixedfilters = posti->size;

         int lengthALL_F = 0;
         int lengthALL_T = 0;         

         for (int i=0;i<(int)(posti->size);i++)
         {
             lengthALL_T = lengthALL_T + gsl_vector_get(posti,i);
         }
         
         lengthALL_F = lengthALL_T*2;
         
         gsl_matrix *matrixALL_OFFx = gsl_matrix_alloc(ntemplates,lengthALL_F);
         gsl_matrix *matrixALL_OFTx = gsl_matrix_alloc(ntemplates,lengthALL_T);
         gsl_matrix *matrixALLab_OFFx = gsl_matrix_alloc(ntemplates,lengthALL_F);
         gsl_matrix *matrixALLab_OFTx = gsl_matrix_alloc(ntemplates,lengthALL_T);
         
         char str_length[125];
         
         gsl_matrix *matrixAux_OFFx = NULL;
         gsl_matrix *matrixAuxab_OFFx = NULL;
         int index = 0;
         strcpy(obj.nameTable,"FIXFILTF");
         obj.iniRow = 1;
         obj.endRow = ntemplates;
         for (int i=0;i<(int)(posti->size);i++)
         {
             snprintf(str_length,125,"%d",(int) gsl_vector_get(posti,i));	
             matrixAux_OFFx = gsl_matrix_alloc(ntemplates,gsl_vector_get(posti,i)*2);
             
             strcpy(obj.nameCol,(string("F")+string(str_length)).c_str());
             if (readFitsComplex (obj,&matrixAux_OFFx))
             {
                 message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                 EP_PRINT_ERROR(message,*status);
                 *status=EPFAIL; return(library_collection);
             }
             for (int j=0;j<(int)(matrixAux_OFFx->size1);j++)
             {
                 for (int k=0;k<(int)(matrixAux_OFFx->size2);k++)
                 {
                     gsl_matrix_set(matrixALL_OFFx,j,k+index,gsl_matrix_get(matrixAux_OFFx,j,k));
                 }
             }
            
             index = index + gsl_vector_get(posti,i)*2;
             
             if (ntemplates > 1)
             {
                 EP_PRINT_ERROR("Code not prepared for ntemplates>1 and preBuffer=yes",*status);
                 return(library_collection);
             }
             
             gsl_matrix_free(matrixAux_OFFx); matrixAux_OFFx = 0;
             gsl_matrix_free(matrixAuxab_OFFx); matrixAuxab_OFFx = 0;
         }
         
         // FIXFILTT HDU
         strcpy(HDUname,"FIXFILTT");
         if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
         {
             EP_PRINT_ERROR("Error moving to HDU FIXFILTT in library file",*status);
             return(library_collection);
         }
         
         gsl_matrix *matrixAux_OFTx = NULL;
         gsl_matrix *matrixAuxab_OFTx = NULL;
         index = 0;
         strcpy(obj.nameTable,"FIXFILTT");
         for (int i=0;i<(int)(posti->size);i++)
         {
             snprintf(str_length,125,"%d",(int) gsl_vector_get(posti,i));
             matrixAux_OFTx = gsl_matrix_alloc(ntemplates,gsl_vector_get(posti,i));
             
             strcpy(obj.nameCol,(string("T")+string(str_length)).c_str());
             if (readFitsComplex (obj,&matrixAux_OFTx))
             {
                 message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                 EP_PRINT_ERROR(message,*status);
                 *status=EPFAIL; return(library_collection);
             }
             for (int j=0;j<(int)(matrixAux_OFTx->size1);j++)
             {
                 for (int k=0;k<(int)(matrixAux_OFTx->size2);k++)
                 {
                     gsl_matrix_set(matrixALL_OFTx,j,k+index,gsl_matrix_get(matrixAux_OFTx,j,k));
                 }
             }
             
             if (ntemplates > 1)
             {
                 EP_PRINT_ERROR("Code not prepared for ntemplates>1 and preBuffer=yes",*status);
                 return(library_collection);
             }
             
             index = index + gsl_vector_get(posti,i);
             
             gsl_matrix_free(matrixAux_OFTx); matrixAux_OFTx = 0;
             gsl_matrix_free(matrixAuxab_OFTx); matrixAux_OFTx = 0;
         }
         for (int it=0;it<ntemplates;it++)
         {
             library_collection->optimal_filtersFREQ[it].energy		= gsl_vector_get(library_collection->energies,it);
             library_collection->optimal_filtersFREQ[it].ofilter_duration	= lengthALL_F;
             library_collection->optimal_filtersFREQ[it].ofilter    		= gsl_vector_alloc(lengthALL_F);
             
             gsl_matrix_get_row(library_collection->optimal_filtersFREQ[it].ofilter,matrixALL_OFFx,it);
             
             library_collection->optimal_filtersTIME[it].energy		= gsl_vector_get(library_collection->energies,it);
             library_collection->optimal_filtersTIME[it].ofilter_duration 	= lengthALL_T;
             library_collection->optimal_filtersTIME[it].ofilter    		= gsl_vector_alloc(lengthALL_T);
             
             gsl_matrix_get_row(library_collection->optimal_filtersTIME[it].ofilter,matrixALL_OFTx,it);
             
             if (it < ntemplates-1)
             {
                 library_collection->optimal_filtersabFREQ[it].energy		= gsl_vector_get(library_collection->energies,it);
                 library_collection->optimal_filtersabFREQ[it].ofilter_duration  = lengthALL_F;
                 library_collection->optimal_filtersabFREQ[it].ofilter    	= gsl_vector_alloc(lengthALL_F);
                 
                 gsl_matrix_get_row(library_collection->optimal_filtersabFREQ[it].ofilter,matrixALLab_OFFx,it);
                 
                 library_collection->optimal_filtersabTIME[it].energy		= gsl_vector_get(library_collection->energies,it);
                 library_collection->optimal_filtersabTIME[it].ofilter_duration	= lengthALL_F;
                 library_collection->optimal_filtersabTIME[it].ofilter    	= gsl_vector_alloc(lengthALL_T);
                 
                 gsl_matrix_get_row(library_collection->optimal_filtersabTIME[it].ofilter,matrixALLab_OFTx,it);
             }
         }
         
         gsl_matrix_free(matrixALL_OFFx); matrixALL_OFFx = 0;
         gsl_matrix_free(matrixALL_OFTx); matrixALL_OFTx = 0;
         gsl_matrix_free(matrixALLab_OFFx); matrixALLab_OFFx = 0;
         gsl_matrix_free(matrixALLab_OFTx); matrixALLab_OFTx = 0;
     }
     
     if (opmode == 1)
     { 
         char str_length[125];
         obj.iniRow = 1;
         obj.endRow = ntemplates;
         int index = 0;
         
         if (strcmp(filter_domain,"F") == 0)
         {
             // FIXFILTF HDU
             strcpy(HDUname,"FIXFILTF");
             if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
             {
                 EP_PRINT_ERROR("Error moving to HDU FIXFILTF in library file",*status);
                 return(library_collection);
             }
             
             // Get number of fixed optimal filters (columnss)
             int nOFs;
             if (fits_get_num_cols(fptr,&nOFs, status))
             {
                 EP_PRINT_ERROR("Cannot get number of columns in library file",*status);
                 return(library_collection);
             }
             
             if (preBuffer == 0)
             {
                if (ntemplates == 1)
                {
                    int nOFs_aux;
                    nOFs_aux = nOFs-1;		// -1 because the ENERGYcolumn
                    if ((nOFs_aux == floor(log2(template_duration))) || (nOFs_aux == 1 + floor(log2(template_duration))))		nOFs = nOFs-1;		// -1 because the ENERGYcolumn
                    else                                                                            							nOFs = (nOFs-1)/2;	// /2 because the AB column
                }
                else 								                                                                            nOFs = (nOFs-1)/2;	// /2 because the AB column
             }
             else // preBuffer =1
             {
                 if (ntemplates == 1)   nOFs = nOFs-1;		// -1 because the ENERGYcolumn
                 else                   nOFs = (nOFs-1)/2;	// /2 because the AB column             
             }
             
             if (nOFs == 0)	
             {
                 EP_PRINT_ERROR("The library has no fixed optimal filters",EPFAIL); 
                 *status=EPFAIL; return(library_collection);
             }
             library_collection->nfixedfilters = nOFs;
             
             int lengthALL_F = 0;
             int lengthALL_T = 0;
             if (preBuffer == 0)
             {
                 for (int i=0;i<nOFs;i++)
                {
                    if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                    {
                        if (i==0)	lengthALL_F = lengthALL_F + template_durationPLSMXLFF*2;
                        else		lengthALL_F = lengthALL_F + pow(2,floor(log2(template_duration))-i+1)*2;
                    }
                    else	lengthALL_F = lengthALL_F + pow(2,floor(log2(template_duration))-i)*2;
                }
             }
             else  // preBuffer = 1
             {
                if (nOFs != (int)(posti->size))
                {
                    EP_PRINT_ERROR("The number of optimal filters in the library does not match the grading info in the XML file",EPFAIL); 
                    *status=EPFAIL; return(library_collection);
                }
                for (int i=0;i<(int)(posti->size);i++)
                {
                    lengthALL_T = lengthALL_T + gsl_vector_get(posti,i) + gsl_vector_get(pBi,i);
                }
                
                lengthALL_F = lengthALL_T*2;  
             }
                  
             strcpy(obj.nameTable,"FIXFILTF");
             
             if (strcmp(*ofinterp,"MF") == 0)
             {
                 gsl_matrix *matrixALL_OFFx = gsl_matrix_alloc(ntemplates,lengthALL_F);
                 gsl_matrix *matrixAux_OFFx = NULL;
                 
                 if (preBuffer == 0)
                 {
                    for (int i=0;i<nOFs;i++)
                    {
                        if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                        {  
                            if (i==0)
                            {
                                snprintf(str_length,125,"%d",template_durationPLSMXLFF);
                                matrixAux_OFFx = gsl_matrix_alloc(ntemplates,template_durationPLSMXLFF*2);
                            }
                            else
                            {
                                snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i+1)));
                                matrixAux_OFFx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i+1)*2);
                            }
                        }
                        else
                        {
                            snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i)));
                            matrixAux_OFFx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i)*2);
                        }
                        strcpy(obj.nameCol,(string("F")+string(str_length)).c_str());
                        if (readFitsComplex (obj,&matrixAux_OFFx))
                        {
                            message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                            EP_PRINT_ERROR(message,*status);
                            *status=EPFAIL; return(library_collection);
                        }
                        for (int j=0;j<(int)(matrixAux_OFFx->size1);j++)
                        {
                            for (int k=0;k<(int)(matrixAux_OFFx->size2);k++)
                            {
                                gsl_matrix_set(matrixALL_OFFx,j,k+index,gsl_matrix_get(matrixAux_OFFx,j,k));
                            }
                        }
                        
                        if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                        {
                            if (i==0)	index = index + template_durationPLSMXLFF*2;
                            else 		index = index + pow(2,floor(log2(template_duration))-i+1)*2;
                        }
                        else    index = index + pow(2,floor(log2(template_duration))-i)*2;
                        
                        gsl_matrix_free(matrixAux_OFFx); matrixAux_OFFx = 0;
                    }
                 }
                 else // preBuffer=1
                 {
                    for (int i=0;i<nOFs;i++)
                    {
                        snprintf(str_length,125,"%d",(int) gsl_vector_get(posti,i));
                        matrixAux_OFFx = gsl_matrix_alloc(ntemplates,gsl_vector_get(posti,i)*2);
                        
                        strcpy(obj.nameCol,(string("F")+string(str_length)).c_str());
                        if (readFitsComplex (obj,&matrixAux_OFFx))
                        {
                            message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                            EP_PRINT_ERROR(message,*status);
                            *status=EPFAIL; return(library_collection);
                        }
                        for (int j=0;j<(int)(matrixAux_OFFx->size1);j++)
                        {
                            for (int k=0;k<(int)(matrixAux_OFFx->size2);k++)
                            {
                                gsl_matrix_set(matrixALL_OFFx,j,k+index,gsl_matrix_get(matrixAux_OFFx,j,k));
                            }
                        }
                        
                        index = index + gsl_vector_get(posti,i)*2;
                        
                        gsl_matrix_free(matrixAux_OFFx); matrixAux_OFFx = 0;
                    }
                 }
                 
                 for (int it=0;it<ntemplates;it++)
                 {
                     library_collection->optimal_filters[it].energy 			= gsl_vector_get(library_collection->energies,it);
                     library_collection->optimal_filters[it].ofilter_duration 	= lengthALL_F;
                     library_collection->optimal_filters[it].ofilter    		= gsl_vector_alloc(lengthALL_F);
                     
                     gsl_matrix_get_row(library_collection->optimal_filters[it].ofilter,matrixALL_OFFx,it);
                 }
                 
                 gsl_matrix_free(matrixALL_OFFx); matrixALL_OFFx = 0;
             }
             else if (strcmp(*ofinterp,"DAB") == 0)
             {
                 gsl_matrix *matrixALLab_OFFx = gsl_matrix_alloc(ntemplates,lengthALL_F);
                 gsl_matrix *matrixAuxab_OFFx = NULL;
                 
                 if (preBuffer == 0)
                 {
                    for (int i=0;i<nOFs;i++)
                    {
                        if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                        {
                            if (i==0)
                            {
                                snprintf(str_length,125,"%d",template_durationPLSMXLFF);
                                matrixAuxab_OFFx = gsl_matrix_alloc(ntemplates,template_durationPLSMXLFF*2);
                            }
                            else
                            {
                                snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i+1)));
                                matrixAuxab_OFFx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i+1)*2);
                            }
                        }
                        else
                        {
                            snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i)));
                            matrixAuxab_OFFx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i)*2);
                        }
                        strcpy(obj.nameCol,(string("ABF")+string(str_length)).c_str());
                        if (readFitsComplex (obj,&matrixAuxab_OFFx))
                        {
                            message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                            EP_PRINT_ERROR(message,*status);
                            *status=EPFAIL; return(library_collection);
                        }
                        for (int j=0;j<(int)(matrixAuxab_OFFx->size1);j++)
                        {
                            for (int k=0;k<(int)(matrixAuxab_OFFx->size2);k++)
                            {
                                gsl_matrix_set(matrixALLab_OFFx,j,k+index,gsl_matrix_get(matrixAuxab_OFFx,j,k));
                            }
                        }
                        
                        if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                        {
                            if (i==0)	index = index + template_durationPLSMXLFF*2;
                            else		index = index + pow(2,floor(log2(template_duration))-i+1)*2;
                        }
                        else    index = index + pow(2,floor(log2(template_duration))-i)*2;
                        
                        gsl_matrix_free(matrixAuxab_OFFx); matrixAuxab_OFFx = 0;
                    }
                 }
                 else //preBuffer = 1
                 {
                     for (int i=0;i<nOFs;i++)
                     {
                         snprintf(str_length,125,"%d",(int) gsl_vector_get(posti,i));
                         matrixAuxab_OFFx = gsl_matrix_alloc(ntemplates,gsl_vector_get(posti,i)*2);
                         
                         strcpy(obj.nameCol,(string("ABF")+string(str_length)).c_str());
                         if (readFitsComplex (obj,&matrixAuxab_OFFx))
                         {
                             message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                             EP_PRINT_ERROR(message,*status);
                             *status=EPFAIL; return(library_collection);
                         }
                         for (int j=0;j<(int)(matrixAuxab_OFFx->size1);j++)
                         {
                             for (int k=0;k<(int)(matrixAuxab_OFFx->size2);k++)
                             {
                                 gsl_matrix_set(matrixALLab_OFFx,j,k+index,gsl_matrix_get(matrixAuxab_OFFx,j,k));
                             }
                         }
                         
                         index = index + gsl_vector_get(posti,i)*2;
                         
                         gsl_matrix_free(matrixAuxab_OFFx); matrixAuxab_OFFx = 0;
                     }
                 }
                 
                 for (int it=0;it<ntemplates;it++)
                 {
                     library_collection->optimal_filters[it].energy			= gsl_vector_get(library_collection->energies,it);
                     library_collection->optimal_filters[it].ofilter_duration 	= lengthALL_F;
                     library_collection->optimal_filters[it].ofilter    		= gsl_vector_alloc(lengthALL_F);
                     
                     gsl_matrix_get_row(library_collection->optimal_filters[it].ofilter,matrixALLab_OFFx,it);
                 }
                 
                 gsl_matrix_free(matrixALLab_OFFx); matrixALLab_OFFx = 0;
             }
         }
         else if (strcmp(filter_domain,"T") == 0)
         {
             // FIXFILTT HDU
             strcpy(HDUname,"FIXFILTT");
             if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
             {
                 EP_PRINT_ERROR("Error moving to HDU FIXFILTT in library file",*status);
                 return(library_collection);
             }
             
             // Get number of fixed optimal filters (columnss)
             int nOFs;
             if (fits_get_num_cols(fptr,&nOFs, status))
             {
                 EP_PRINT_ERROR("Cannot get number of columns in library file",*status);
                 return(library_collection);
             }
             
             if (preBuffer == 0)
             {
                if (ntemplates == 1)	
                {
                    int nOFs_aux;
                    nOFs_aux = nOFs-1;		// -1 because the ENERGYcolumn
                    if ((nOFs_aux == floor(log2(template_duration))) || (nOFs_aux == 1 + floor(log2(template_duration))))		nOFs = nOFs-1;		// -1 because the ENERGYcolumn
                    else 								nOFs = (nOFs-1)/2;	// /2 because the AB column
                }
                else	  nOFs = (nOFs-1)/2;	// /2 because the AB column
             }
             else // preBuffer = 1
             {
                 if (ntemplates == 1)   nOFs = nOFs-1;		// -1 because the ENERGYcolumn
                 else                   nOFs = (nOFs-1)/2;	// /2 because the AB column
             }
             
             if (nOFs == 0)	
             {
                 EP_PRINT_ERROR("The library has no fixed optimal filters",EPFAIL); 
                 *status=EPFAIL; return(library_collection);
             }
             library_collection->nfixedfilters = nOFs;
             
             int lengthALL_T = 0;
             if (preBuffer == 0)
             {
                for (int i=0;i<nOFs;i++)
                {
                    if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                    {
                        if (i==0)	lengthALL_T = lengthALL_T + template_durationPLSMXLFF;
                        else  		lengthALL_T = lengthALL_T + pow(2,floor(log2(template_duration))-i+1);
                    }
                    else    
                    {
                        lengthALL_T = lengthALL_T + pow(2,floor(log2(template_duration))-i);
                    }
                }
             }
             else // preBuffer =1
             {
                if (nOFs != (int)(posti->size))
                {
                    EP_PRINT_ERROR("The number of optimal filters in the library does not match the grading info in the XML file",EPFAIL); 
                    *status=EPFAIL; return(library_collection);
                }
                for (int i=0;i<nOFs;i++)
                {
                    lengthALL_T = lengthALL_T + gsl_vector_get(posti,i) + gsl_vector_get(pBi,i);
                }
             }
             
             strcpy(obj.nameTable,"FIXFILTT");
             
             if (strcmp(*ofinterp,"MF") == 0)
             {
                 gsl_matrix *matrixALL_OFTx = gsl_matrix_alloc(ntemplates,lengthALL_T);
                 gsl_matrix *matrixAux_OFTx = NULL;
                 
                 if (preBuffer == 0)
                 {
                    for (int i=0;i<nOFs;i++)
                    {
                        if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                        {
                            if (i==0)
                            {
                                snprintf(str_length,125,"%d",template_durationPLSMXLFF);
                                matrixAux_OFTx = gsl_matrix_alloc(ntemplates,template_durationPLSMXLFF);
                            }
                            else
                            {
                                snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i+1)));
                                matrixAux_OFTx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i+1));
                            }
                        }
                        else
                        {
                            snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i)));
                            matrixAux_OFTx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i));
                        }
                        strcpy(obj.nameCol,(string("T")+string(str_length)).c_str());
                        if (readFitsComplex (obj,&matrixAux_OFTx))
                        {
                            message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                            EP_PRINT_ERROR(message,*status);
                            *status=EPFAIL; return(library_collection);
                        }
                        
                        for (int j=0;j<(int)(matrixAux_OFTx->size1);j++)
                        {
                            for (int k=0;k<(int)(matrixAux_OFTx->size2);k++)
                            {
                                gsl_matrix_set(matrixALL_OFTx,j,k+index,gsl_matrix_get(matrixAux_OFTx,j,k));
                            }
                        }
                        
                        if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                        {
                            if (i==0)	index = index + template_durationPLSMXLFF;
                            else		index = index + pow(2,floor(log2(template_duration))-i+1);
                        }
                        else    index = index + pow(2,floor(log2(template_duration))-i);
                        
                        gsl_matrix_free(matrixAux_OFTx); matrixAux_OFTx = 0;
                    }
                 }
                 else // preBuffer =1
                 {
                    for (int i=0;i<nOFs;i++)
                    {
                        snprintf(str_length,125,"%d",(int) gsl_vector_get(posti,i));
                        matrixAux_OFTx = gsl_matrix_alloc(ntemplates,gsl_vector_get(posti,i));
                        
                        strcpy(obj.nameCol,(string("T")+string(str_length)).c_str());
                        if (readFitsComplex (obj,&matrixAux_OFTx))
                        {
                            message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                            EP_PRINT_ERROR(message,*status);
                            *status=EPFAIL; return(library_collection);
                        }
                        
                        for (int j=0;j<(int)(matrixAux_OFTx->size1);j++)
                        {
                            for (int k=0;k<(int)(matrixAux_OFTx->size2);k++)
                            {
                                gsl_matrix_set(matrixALL_OFTx,j,k+index,gsl_matrix_get(matrixAux_OFTx,j,k));
                            }
                        }
                        
                        index = index + gsl_vector_get(posti,i);
                        
                        gsl_matrix_free(matrixAux_OFTx); matrixAux_OFTx = 0;
                    }
                 }
                 
                 for (int it=0;it<ntemplates;it++)
                 {
                     library_collection->optimal_filters[it].energy			  = gsl_vector_get(library_collection->energies,it);
                     library_collection->optimal_filters[it].ofilter_duration = lengthALL_T;
                     library_collection->optimal_filters[it].ofilter    	  = gsl_vector_alloc(lengthALL_T);
                     
                     gsl_matrix_get_row(library_collection->optimal_filters[it].ofilter,matrixALL_OFTx,it);
                 }
                 
                 gsl_matrix_free(matrixALL_OFTx); matrixALL_OFTx = 0;
             }
             else if (strcmp(*ofinterp,"DAB") == 0)
             {
                 gsl_matrix *matrixALLab_OFTx = gsl_matrix_alloc(ntemplates,lengthALL_T);
                 gsl_matrix *matrixAuxab_OFTx = NULL;
                 
                 if (preBuffer == 0)
                 {
                    for (int i=0;i<nOFs;i++)
                    {
                        if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                        {
                            if (i==0)
                            {
                                snprintf(str_length,125,"%d",template_durationPLSMXLFF);
                                matrixAuxab_OFTx = gsl_matrix_alloc(ntemplates,template_durationPLSMXLFF);
                            }
                            else
                            {
                                snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i+1)));
                                matrixAuxab_OFTx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i+1));
                            }
                        }
                        else
                        {
                            snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i)));
                            matrixAuxab_OFTx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i));
                        }
                        strcpy(obj.nameCol,(string("ABT")+string(str_length)).c_str());
                        if (readFitsComplex (obj,&matrixAuxab_OFTx))
                        {
                            message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                            EP_PRINT_ERROR(message,*status);
                            *status=EPFAIL; return(library_collection);
                        }
                        for (int j=0;j<(int)(matrixAuxab_OFTx->size1);j++)
                        {
                            for (int k=0;k<(int)(matrixAuxab_OFTx->size2);k++)
                            {
                                gsl_matrix_set(matrixALLab_OFTx,j,k+index,gsl_matrix_get(matrixAuxab_OFTx,j,k));
                            }
                        }
                        
                        if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                        {
                            if (i==0)	index = index + template_durationPLSMXLFF;
                            else		index = index + pow(2,floor(log2(template_duration))-i+1);
                        }
                        else    index = index + pow(2,floor(log2(template_duration))-i);
                        
                        gsl_matrix_free(matrixAuxab_OFTx); matrixAuxab_OFTx = 0;
                    }
                 }
                 else // preBuffer = 1
                 {
                    for (int i=0;i<nOFs;i++)
                    {
                        snprintf(str_length,125,"%d",(int) gsl_vector_get(posti,i));
                        matrixAuxab_OFTx = gsl_matrix_alloc(ntemplates,gsl_vector_get(posti,i)+gsl_vector_get(pBi,i));
                        
                        strcpy(obj.nameCol,(string("ABT")+string(str_length)).c_str());
                        if (readFitsComplex (obj,&matrixAuxab_OFTx))
                        {
                            message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                            EP_PRINT_ERROR(message,*status);
                            *status=EPFAIL; return(library_collection);
                        }
                        for (int j=0;j<(int)(matrixAuxab_OFTx->size1);j++)
                        {
                            for (int k=0;k<(int)(matrixAuxab_OFTx->size2);k++)
                            {
                                gsl_matrix_set(matrixALLab_OFTx,j,k+index,gsl_matrix_get(matrixAuxab_OFTx,j,k));
                            }
                        }
                        
                        index = index + gsl_vector_get(posti,i);
                        
                        gsl_matrix_free(matrixAuxab_OFTx); matrixAuxab_OFTx = 0;
                    }
                 }
                 
                 for (int it=0;it<ntemplates;it++)
                 {
                     library_collection->optimal_filters[it].energy			= gsl_vector_get(library_collection->energies,it);
                     library_collection->optimal_filters[it].ofilter_duration 	= lengthALL_T;
                     library_collection->optimal_filters[it].ofilter    		= gsl_vector_alloc(lengthALL_T);
                     
                     gsl_matrix_get_row(library_collection->optimal_filters[it].ofilter,matrixALLab_OFTx,it);
                 }
                 
                 gsl_matrix_free(matrixALLab_OFTx); matrixALLab_OFTx = 0;
             }	
         }  
     }
     
     if ((opmode == 1) && ((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)) && (strcmp(ofnoise,"WEIGHTM") == 0) && (oflib == 1))
     {
         char str_length[125];
         obj.iniRow = 1;
         obj.endRow = ntemplates;
         strcpy(obj.nameTable,"PRCLOFWM");
         int index = 0;
         
         // PRCLOFWM HDU
         strcpy(HDUname,"PRCLOFWM");
         if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
         {
             * status = 1;
             EP_PRINT_ERROR("Error moving to HDU PRCLOFWM in library file, probably PRCLOFWM HDU does not exist",*status);
             return(library_collection);
         }
         
         // Get number of columns
         int nOFs;
         if (fits_get_num_cols(fptr,&nOFs, status))
         {
             EP_PRINT_ERROR("Cannot get number of columns in library file",*status);
             return(library_collection);
         }
         nOFs = nOFs-1;	// -1 because the ENERGYcolumn
         if (nOFs == 0)	
         {
             EP_PRINT_ERROR("The library has no columns (PRCLOFWM)",EPFAIL); 
             *status=EPFAIL; return(library_collection);
         }
         library_collection->nfixedfilters = nOFs;
         
         int lengthALL_PRCLOFWM = 0;
         if (preBuffer == 0)
         {
            for (int i=0;i<nOFs;i++)
            {
                if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                {
                    if (i==0)	lengthALL_PRCLOFWM = lengthALL_PRCLOFWM + template_durationPLSMXLFF*2;
                    else  		lengthALL_PRCLOFWM = lengthALL_PRCLOFWM + pow(2,floor(log2(template_duration))-i+1)*2;
                }
                else    lengthALL_PRCLOFWM = lengthALL_PRCLOFWM + pow(2,floor(log2(template_duration))-i)*2;
            }
         }
         else
         {
            for (int i=0;i<(int)(posti->size);i++)
            {
                lengthALL_PRCLOFWM = lengthALL_PRCLOFWM + gsl_vector_get(posti,i);
            }
            lengthALL_PRCLOFWM = lengthALL_PRCLOFWM*2;
         }

         gsl_matrix *matrixALL_PRCLOFWMx = gsl_matrix_alloc(ntemplates,lengthALL_PRCLOFWM);
         gsl_matrix *matrixAux_PRCLOFWMx = NULL;
         for (int i=0;i<nOFs;i++)
         {
             if (preBuffer == 0)
             {
                if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                {
                    if (i==0)
                    {
                        snprintf(str_length,125,"%d",template_durationPLSMXLFF);
                        matrixAux_PRCLOFWMx = gsl_matrix_alloc(ntemplates,template_durationPLSMXLFF*2);
                    }
                    else
                    {
                        snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i+1)));
                        matrixAux_PRCLOFWMx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i+1)*2);
                    }
                }
                else
                {
                    snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i)));
                    matrixAux_PRCLOFWMx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i)*2);
                }
             }
             else if (preBuffer == 1)
             {
                 snprintf(str_length,125,"%d",(int) (gsl_vector_get(posti,i)));
                 matrixAux_PRCLOFWMx = gsl_matrix_alloc(ntemplates,gsl_vector_get(posti,i)*2);
             }
             
             strcpy(obj.nameCol,(string("OFW")+string(str_length)).c_str());
             if (readFitsComplex (obj,&matrixAux_PRCLOFWMx))
             {
                 message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                 EP_PRINT_ERROR(message,*status);
                 *status=EPFAIL; return(library_collection);
             }
             
             for (int j=0;j<(int)(matrixAux_PRCLOFWMx->size1);j++)
             {
                 for (int k=0;k<(int)(matrixAux_PRCLOFWMx->size2);k++)
                 {
                     gsl_matrix_set(matrixALL_PRCLOFWMx,j,k+index,gsl_matrix_get(matrixAux_PRCLOFWMx,j,k));
                 }
             }
             
             if (preBuffer == 0)
             {
                if ((ncols == 7) || (ncols == 9) || (ncols == 10) || (ncols == 19))
                {
                    if (i==0)	index = index + template_durationPLSMXLFF*2;
                    else		index = index + pow(2,floor(log2(template_duration))-i+1)*2;
                }
                else    index = index + pow(2,floor(log2(template_duration))-i)*2;
             }
             else
             {
                 index = index + gsl_vector_get(posti,i)*2;
             }
             
             gsl_matrix_free(matrixAux_PRCLOFWMx); matrixAux_PRCLOFWMx = 0;
         }
         
         library_collection->PRCLOFWM = gsl_matrix_alloc(ntemplates, lengthALL_PRCLOFWM);
         gsl_matrix_memcpy(library_collection->PRCLOFWM,matrixALL_PRCLOFWMx);
         gsl_matrix_free(matrixALL_PRCLOFWMx); matrixALL_PRCLOFWMx = 0;
     }
     
     if ((opmode == 1) && (strcmp(energy_method,"WEIGHTN") == 0) && (ntemplates > 1))
     {
         char str_length[125];
         obj.iniRow = 1;
         obj.endRow = ntemplates;
         int index = 0;
         
         // PRECALWN HDU
         strcpy(HDUname,"PRECALWN");
         if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
         {
             EP_PRINT_ERROR("Error moving to HDU PRECALWN in library file",*status);
             return(library_collection);
         }
         
         // Get number of columns
         int nOFs;
         if (fits_get_num_cols(fptr,&nOFs, status))
         {
             EP_PRINT_ERROR("Cannot get number of columns in library file",*status);
             return(library_collection);
         }
         nOFs = nOFs-1;	// -1 because the ENERGYcolumn
         if (nOFs == 0)	
         {
             EP_PRINT_ERROR("The library has no fixed columns (PRECALWN)",EPFAIL); 
             *status=EPFAIL; return(library_collection);
         }
         library_collection->nfixedfilters = nOFs;
         
         int lengthALL_PRCLWN = 0;
         for (int i=0;i<nOFs;i++)
         {
             lengthALL_PRCLWN = lengthALL_PRCLWN + pow(2,floor(log2(template_duration))-i)*2;
         }
         
         strcpy(obj.nameTable,"PRECALWN");
         
         gsl_matrix *matrixALL_PRCLWNx = gsl_matrix_alloc(ntemplates,lengthALL_PRCLWN);
         gsl_matrix *matrixAux_PRCLWNx = NULL;
         for (int i=0;i<nOFs;i++)
         {
             snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i)));
             strcpy(obj.nameCol,(string("PCL")+string(str_length)).c_str());
             matrixAux_PRCLWNx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i)*2);
             if (readFitsComplex (obj,&matrixAux_PRCLWNx))
             {
                 message = "Cannot read " + string(obj.nameCol) + " column in library FITS file";
                 EP_PRINT_ERROR(message,*status);
                 *status=EPFAIL; return(library_collection);
             }
             
             for (int j=0;j<(int)(matrixAux_PRCLWNx->size1);j++)
             {
                 for (int k=0;k<(int)(matrixAux_PRCLWNx->size2);k++)
                 {
                     gsl_matrix_set(matrixALL_PRCLWNx,j,k+index,gsl_matrix_get(matrixAux_PRCLWNx,j,k));
                 }
             }
             
             index = index + pow(2,floor(log2(template_duration))-i)*2;
             
             gsl_matrix_free(matrixAux_PRCLWNx); matrixAux_PRCLWNx = 0;
         }
         
         library_collection->PRECALWN = gsl_matrix_alloc(ntemplates, lengthALL_PRCLWN);
         gsl_matrix_memcpy(library_collection->PRECALWN,matrixALL_PRCLWNx);
         gsl_matrix_free(matrixALL_PRCLWNx); matrixALL_PRCLWNx = 0;
     }
     
     delete [] obj.nameTable;
     delete [] obj.nameCol;
     delete [] obj.unit;
     
     if (fits_close_file(fptr, status))
     {
         EP_PRINT_ERROR("Error closing library file",*status);
         return(library_collection);
     }
     
     return(library_collection);
 }
 /*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 10 ************************************************************
  * getNoiseSpec: This function creates and retrieves a NoiseSpec from a file.
  * 
  * - Create *NoiseSpec* structure
  * - Open FITS file, move to the NOISE, NOISEALL and WEIGHTMS HDUs and get necessary keywords
  * - Allocate *NoiseSpec* structure
  * - Get noise spectrum (CSD), and noise frequencies (FREQ) column numbers
  * - Read column CSD and save it into the structure
  * - Read column FREQ and save it into the structure
  * - Read columns Wx with the noise weight matrix from noise intervals and save them into the structure
  * - Return noise spectrum
  *    
  * Parameters:
  * - filename: File name with noise
  * - opmode: Calibration run (0) or energy reconstruction run (1)
  * - hduPRCLOFWM: Add or not the PRCLOFWM HDU in the library file (1/0) 
  * - energy_method: Energy calculation Method: OPTFILT, WEIGHT, WEIGHTN, I2R, I2RFITTED or PCA
  * - ofnoise: Noise to use with Optimal Filtering: NSD or WEIGHTM
  * - filter_method: Filtering Method: F0 (deleting the zero frequency bin) or B0 (deleting the baseline) or F0B0 (deleting always the baseline)
  * - status: Input/output status
  ******************************************************************************/
 NoiseSpec* getNoiseSpec(const char* const filename, int opmode, int hduPRCLOFWM, char *energy_method, char *ofnoise, char *filter_method, int* const status)
 {
     string message = "";

     // Create NoiseSpec structure
     NoiseSpec* noise_spectrum = new NoiseSpec;
     
     // Open FITS file in READONLY mode
     fitsfile* fptr = NULL;
     if(fits_open_file(&fptr, filename, READONLY, status))
     {
         EP_PRINT_ERROR("Error opening noise file",*status);
         *status=EPFAIL; return(noise_spectrum);
     }
     
     // Move to the NOISE HDU
     int extver = 1;
     char HDUname[12];
     char keyname[10];
     strcpy(HDUname,"NOISE");
     if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
     {
         EP_PRINT_ERROR("Error moving to HDU NOISE in noise file",*status);
         *status=EPFAIL;return(noise_spectrum);
     }
     strcpy(keyname,"NOISESTD");
     if (fits_read_key(fptr,TDOUBLE,keyname, &noise_spectrum->noiseStd,NULL,status))
     {
         EP_PRINT_ERROR("Cannot read NOISESTD keyword",*status);
         *status=EPFAIL;return(noise_spectrum);
     }
     
     if ((opmode == 0) ||
         ((opmode == 1) && (strcmp(filter_method,"B0") == 0) && ((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0)|| (strcmp(energy_method,"I2RDER") == 0)))
         || ((opmode == 1) && (strcmp(energy_method,"WEIGHT") == 0)))
     {
         strcpy(keyname,"BASELINE");
         
         if (fits_read_key(fptr,TDOUBLE,keyname, &noise_spectrum->baseline,NULL,status))
         {
             EP_PRINT_ERROR("Cannot read BASELINE keyword",*status);
             *status=EPFAIL;return(noise_spectrum);
         }
     }
     
     if (strcmp(energy_method,"WEIGHT") != 0)
     {
         // Move to the NOISEALL hdu
         strcpy(HDUname,"NOISEALL");
         if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
         {
             EP_PRINT_ERROR("Error moving to HDU NOISEALL in noise file",*status);
             *status=EPFAIL;return(noise_spectrum);
         }
         
         // Get number of rows
         long noise_duration;
         if (fits_get_num_rows(fptr,&noise_duration, status))
         {
             EP_PRINT_ERROR("Cannot get number of rows in noise file",*status);
             *status=EPFAIL;return(noise_spectrum);
         }
         if (noise_duration == 0)	
         {	
             EP_PRINT_ERROR("Number of rows in noise file is 0",EPFAIL); 
             *status=EPFAIL; return(noise_spectrum);
         }
         noise_spectrum->noise_duration = noise_duration/2;
         
         // Allocate NoiseSpec structure
         // It is not necessary to check the allocation because 'noise_duration' has been checked previously
         noise_spectrum->noisespec  = gsl_vector_alloc(noise_duration);
         noise_spectrum->noisefreqs = gsl_vector_alloc(noise_duration);
         
         // Get noise spectrum (CSD), and noise frequencies (FREQ) column numbers
         char column_name[12];
         int CSD_colnum = 0;
         int FREQ_colnum = 0;
         strcpy(column_name,"CSD");
         if (fits_get_colnum(fptr, CASEINSEN,column_name, &CSD_colnum, status))
         {
             EP_PRINT_ERROR("Cannot get column number for CSD in noise file",*status);
             *status=EPFAIL; return(noise_spectrum);
         }
         strcpy(column_name,"FREQ");
         if (fits_get_colnum(fptr, CASEINSEN,column_name, &FREQ_colnum, status))
         {
             EP_PRINT_ERROR("Cannot get column number for FREQ in noise file", *status);
             *status=EPFAIL; return(noise_spectrum);
         }	
         
         // Read column CSD and save it into the structure
         IOData obj;
         obj.inObject = fptr;
         obj.nameTable = new char [255];
         strcpy(obj.nameTable,"NOISEALL");
         obj.iniCol = 0;
         obj.nameCol = new char [255];
         obj.unit = new char [255];
         strcpy(obj.nameCol,"CSD");
         obj.type = TDOUBLE;
         obj.iniRow = 1;
         obj.endRow = noise_duration;
         if (readFitsSimple (obj,&noise_spectrum->noisespec))
         {
             EP_PRINT_ERROR("Cannot read CSD colum in noise file",*status);
             *status=EPFAIL; return(noise_spectrum);
         }
         
         // Read column FREQ and save it into the structure
         strcpy(obj.nameCol,"FREQ");
         if (readFitsSimple (obj,&noise_spectrum->noisefreqs))
         {
             EP_PRINT_ERROR("Cannot read FREQ colum in noise file",*status);
             *status=EPFAIL; return(noise_spectrum);
         }
         
         if (((opmode == 0) && (hduPRCLOFWM == 1))
             || ((opmode == 1) && ((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RFITTED") == 0) || (strcmp(energy_method,"I2RDER") == 0)) && (strcmp(ofnoise,"WEIGHTM") == 0)))
         {
             // Move to the WEIGHTMS hdu
             strcpy(HDUname,"WEIGHTMS");
             if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
             {
                 EP_PRINT_ERROR("Error moving to HDU WEIGHTMS in noise file",*status);
                 *status=EPFAIL;return(noise_spectrum);
             }
             
             // Get number of columns
             int noiseW_numcols;
             if (fits_get_num_cols(fptr,&noiseW_numcols, status))
             {
                 EP_PRINT_ERROR("Cannot get number of columns in noise file (WEIGHTMS)",*status);
                 *status=EPFAIL;return(noise_spectrum);
             }
             noise_spectrum->weightMatrixes = gsl_matrix_alloc(noiseW_numcols,pow(2,noiseW_numcols)*pow(2,noiseW_numcols));
             gsl_matrix_set_all(noise_spectrum->weightMatrixes,-999.0);
             gsl_vector *weightpoints = gsl_vector_alloc(noiseW_numcols);
             for (int i=0;i<(int)(weightpoints->size);i++)	gsl_vector_set(weightpoints,i,pow(2,noiseW_numcols-i));
             
             strcpy(obj.nameTable,"WEIGHTMS");
             gsl_matrix *weightMatrixi;
             char str_length[125];
             for (int i=0;i<(int)(weightpoints->size);i++)
             {
                 snprintf(str_length,125,"%d",(int) gsl_vector_get(weightpoints,i));
                 strcpy(obj.nameCol,(string("W")+string(str_length)).c_str());
                 obj.type = TDOUBLE;
                 obj.iniRow = 1;
                 obj.endRow = 1;
                 weightMatrixi = gsl_matrix_alloc(1,gsl_vector_get(weightpoints,i)*gsl_vector_get(weightpoints,i));
                 if (readFitsComplex (obj,&weightMatrixi))
                 {
                     message = "Cannot read " + string(obj.nameCol) + " column in noise FITS file";
                     EP_PRINT_ERROR(message,*status);
                     *status=EPFAIL; return(noise_spectrum);
                 }
                 
                 if (i == 0) 	
                 {
                     gsl_vector *vectoraux = gsl_vector_alloc(weightMatrixi->size2);
                     gsl_matrix_get_row(vectoraux,weightMatrixi,0);
                     gsl_matrix_set_row(noise_spectrum->weightMatrixes,0,vectoraux);
                     gsl_vector_free(vectoraux); vectoraux = 0;
                 }
                 else
                 {
                     for (int j=0;j<gsl_vector_get(weightpoints,i)*gsl_vector_get(weightpoints,i);j++)
                     {
                         gsl_matrix_set(noise_spectrum->weightMatrixes,i,j,gsl_matrix_get(weightMatrixi,0,j));
                     }
                 }
                 
                 gsl_matrix_free(weightMatrixi); weightMatrixi = 0;
             }
         }
         
         delete [] obj.nameTable;
         delete [] obj.nameCol;
         delete [] obj.unit;
     }
     
     if (fits_close_file(fptr, status))
     {
         EP_PRINT_ERROR("Error closing noise file",*status);
         return(noise_spectrum);
         
     }
     fptr = 0;
     // Return noise spectrum
     return(noise_spectrum);
 }
 /*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


 /***** SECTION 11 ************************************************************
  * IntegrafreeTesEventListSIRENA: This function frees the structure in the input parameter.
  *
  ******************************************************************************/
 extern "C" void IntegrafreeTesEventListSIRENA(TesEventList* event_list)
 {
    if (event_list->energies != NULL) 	delete [] event_list->energies;
    if (event_list->avgs_4samplesDerivative != NULL) 	delete [] event_list->avgs_4samplesDerivative;
    if (event_list->Es_lowres != NULL) 	delete [] event_list->Es_lowres;
    if (event_list->phis != NULL) 	        delete [] event_list->phis;
    if (event_list->lagsShifts != NULL) 	delete [] event_list->lagsShifts;
    if (event_list->bsln != NULL) 	        delete [] event_list->bsln;
    if (event_list->rmsbsln != NULL) 	delete [] event_list->rmsbsln;
    if (event_list->grading != NULL) 	delete [] event_list->grading;
    if (event_list->grades2 != NULL) 	delete [] event_list->grades2;
    if (event_list->ph_ids != NULL) 	delete [] event_list->ph_ids;
    if (event_list->ph_ids2 != NULL) 	delete [] event_list->ph_ids2;
    if (event_list->ph_ids3 != NULL) 	delete [] event_list->ph_ids3;
    if (event_list->pix_ids != NULL) 	delete [] event_list->pix_ids;
    if (event_list->tstarts != NULL) 	delete [] event_list->tstarts;
    if (event_list->tends != NULL) 	        delete [] event_list->tends;
    if (event_list->risetimes != NULL)       delete [] event_list->risetimes;
    if (event_list->falltimes != NULL)       delete [] event_list->falltimes;
}
/*xxxx end of SECTION 11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 12 ************************************************************
* checksum: Calculate the checksum
*
* Parameters:
* - buffer:
* - len:
* - seed:
******************************************************************************/
extern "C" unsigned checksum(void *buffer, size_t len, unsigned int seed)
{
      unsigned char *buf = (unsigned char *)buffer;
      size_t i;

      for (i = 0; i < len; ++i)
            seed += (unsigned int)(*buf++);
      return seed;
}

/*xxxx end of SECTION 12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

 
 // It waits until all threads finish and it builds the 'event_list' by using the results
 void th_end(ReconstructInitSIRENA* reconstruct_init,
             PulsesCollection** pulsesAll)
 {
     //log_trace("Ending the reconstruction...");
     if (!scheduler::get()->is_threading()) { 
         delete scheduler::get();
         return;
     }
     
     if(strcmp(reconstruct_init->EnergyMethod,"PCA") != 0){
         scheduler::get()->set_is_running_energy(true);
     }
     
     scheduler::get()->finish_reconstruction_v2(pulsesAll);
     //scheduler::get()->finish_reconstruction(reconstruct_init, 
     //                                        pulsesAll, optimalFilter);
 }
 
 // It returns the current 'event_list'
 int th_get_event_list(TesEventList** test_event, TesRecord** record)
 {
     //log_trace("Getting event list...");
     if(!scheduler::get()->has_records()) {
         delete scheduler::get();
         *test_event = 0;
         *record = 0;
         return 0;
     }
     
     scheduler::get()->get_test_event(test_event, record);
     
     return 1;
 }
 
 // It returns 'true' if THREADING mode
 int is_threading(){
     return scheduler::get()->is_threading();
 }
 
 /* structs constructors and destructors */
 ReconstructInitSIRENA::ReconstructInitSIRENA():
 library_collection(0),
 threshold(0.0f),
 record_file_fptr(0),
 pulse_length(0),
 scaleFactor(0.0f),
 samplesUp(0),
 samplesDown(0),
 nSgms(0.0f),
 detectSP(0),
 monoenergy(0.0f),
 hduPRECALWN(0),
 hduPRCLOFWM(0),
 largeFilter(0),
 LrsT(0.0f),
 LbT(0.0f),
 opmode(0),
 noise_spectrum(0),
 filtEev(0.0f),
 Ifit(0.0f),
 LagsOrNot(0),
 nLags(0),
 Fitting35(0),
 OFIter(0),
 OFLib(0),
 OFLength(0),
 preBuffer(0),
 intermediate(0),
 errorT(0),
 Sum0Filt(0),
 clobber(0),
 maxPulsesPerRecord(0),
 SaturationValue(0.0f),
 tstartPulse2(0),
 tstartPulse3(0),
 energyPCA1(0.0f),
 energyPCA2(0.0f),
 grading(0),
 i2rdata(0)
 {
     
 }
 
 ReconstructInitSIRENA::ReconstructInitSIRENA(const ReconstructInitSIRENA& other):
 library_collection(0),
 threshold(other.threshold),
 record_file_fptr(0),
 pulse_length(other.pulse_length),
 scaleFactor(other.scaleFactor),
 samplesUp(other.samplesUp),
 samplesDown(other.samplesDown),
 nSgms(other.nSgms),
 detectSP(other.detectSP),
 monoenergy(other.monoenergy),
 hduPRECALWN(other.hduPRECALWN),
 hduPRCLOFWM(other.hduPRCLOFWM),
 largeFilter(other.largeFilter),
 LrsT(other.LrsT),
 LbT(other.LbT),
 opmode(other.opmode),
 noise_spectrum(0),
 filtEev(other.filtEev),
 Ifit(other.Ifit),
 LagsOrNot(other.LagsOrNot),
 nLags(other.nLags),
 Fitting35(other.Fitting35),
 OFIter(other.OFIter),
 OFLib(other.OFLib),
 OFLength(other.OFLength),
 preBuffer(other.preBuffer),
 intermediate(other.intermediate),
 errorT(other.errorT),
 Sum0Filt(other.Sum0Filt),
 clobber(other.clobber),
 maxPulsesPerRecord(other.maxPulsesPerRecord),
 SaturationValue(other.SaturationValue),
 tstartPulse2(other.tstartPulse2),
 tstartPulse3(other.tstartPulse3),
 energyPCA1(other.energyPCA1),
 energyPCA2(other.energyPCA2),
 grading(0),
 i2rdata(0)
 {
     strcpy(sirenaVersion, other.sirenaVersion);

     if(other.library_collection){
         library_collection = new LibraryCollection();
         *library_collection = (*other.library_collection);
     }

     strcpy(library_file,other.library_file);
     strcpy(record_file,other.record_file);

     //record_file_fptr
     // Here we copy the ptr of the fits file, this is NOT thread safe,
     // even allowing reentrant here we should open the file again,
     // and even by doing that the writing through multiple threads
     // won't be safe
     record_file_fptr = other.record_file_fptr;

     strcpy(noise_file,other.noise_file);
     strcpy(event_file,other.event_file);

     strcpy(detectionMode, other.detectionMode);

     if(other.noise_spectrum){
         noise_spectrum = new NoiseSpec();
         *noise_spectrum = (*other.noise_spectrum);
     }

     strcpy(FilterDomain, other.FilterDomain);
     strcpy(FilterMethod, other.FilterMethod);
     strcpy(EnergyMethod, other.EnergyMethod);

     strcpy(OFNoise, other.OFNoise);

     strcpy(OFInterp, other.OFInterp);
     strcpy(OFStrategy, other.OFStrategy);

     strcpy(detectFile, other.detectFile);

     strcpy(XMLFile, other.XMLFile);

     if(other.grading){
         grading = new Grading();
         *grading = (*other.grading);
     }

     if(other.i2rdata){
         i2rdata = new I2RData();
         *i2rdata = (*other.i2rdata);
     }
 }

 ReconstructInitSIRENA&
 ReconstructInitSIRENA::operator=(const ReconstructInitSIRENA& other)
 {
     strcpy(sirenaVersion, other.sirenaVersion);

     if (this != &other){
         if(library_collection) {
             delete library_collection; library_collection = 0;
         }
         if(other.library_collection){
             library_collection = new LibraryCollection();
             *library_collection = (*other.library_collection);
         }
         
         threshold = other.threshold;
         strcpy(library_file,other.library_file);
         strcpy(record_file,other.record_file);
         
         //record_file_fptr
         // Here we copy the ptr of the fits file, this is NOT thread safe,
         // even allowing reentrant here we should open the file again,
         // and even by doing that the writing through multiple threads
         // won't be safe
         record_file_fptr = other.record_file_fptr;
         
         strcpy(noise_file,other.noise_file);
         strcpy(event_file,other.event_file);
         
         pulse_length = other.pulse_length;
         scaleFactor = other.scaleFactor;
         samplesUp = other.samplesUp;
         samplesDown = other.samplesDown;
         nSgms = other.nSgms;
         detectSP = other.detectSP;
         monoenergy = other.monoenergy;
         hduPRECALWN = other.hduPRECALWN;
         hduPRCLOFWM = other.hduPRCLOFWM;
         largeFilter = other.largeFilter;
         LrsT = other.LrsT;
         LbT = other.LbT;
         opmode = other.opmode;
         strcpy(detectionMode, other.detectionMode);
         
         //NoiseSpec
         if(noise_spectrum) {
             delete noise_spectrum; noise_spectrum = 0;
         }
         if(other.noise_spectrum){
             noise_spectrum = new NoiseSpec();
             *noise_spectrum = (*other.noise_spectrum);
         }
         
         strcpy(FilterDomain, other.FilterDomain);
         strcpy(FilterMethod, other.FilterMethod);
         strcpy(EnergyMethod, other.EnergyMethod);
         
         filtEev = other.filtEev;
         
         Ifit = other.Ifit;
         
         strcpy(OFNoise, other.OFNoise);
         
         LagsOrNot = other.LagsOrNot;
         nLags = other.nLags;
         Fitting35 = other.Fitting35;
         OFIter = other.OFIter;
         OFLib = other.OFLib;
         
         strcpy(OFInterp, other.OFInterp);
         strcpy(OFStrategy, other.OFStrategy);
         
         OFLength = other.OFLength;
         preBuffer = other.preBuffer;
         intermediate = other.intermediate;
         strcpy(detectFile, other.detectFile);

         errorT = other.errorT;
         Sum0Filt = other.Sum0Filt;

         clobber = other.clobber;
         maxPulsesPerRecord = other.maxPulsesPerRecord;
         SaturationValue = other.SaturationValue;
         strcpy(tstartPulse1,other.tstartPulse1);
         tstartPulse2 = other.tstartPulse2;
         tstartPulse3 = other.tstartPulse3;
         energyPCA1 = other.energyPCA1;
         energyPCA2 = other.energyPCA2;

         strcpy(XMLFile, other.XMLFile);

         //Grading
         if(grading) {
             delete grading; grading = 0;
         }
         if(other.grading){
             grading = new Grading();
             *grading = (*other.grading);
         }

         //I2RData
         if(i2rdata) {
             delete i2rdata; i2rdata = 0;
         }
         if(other.i2rdata){
             i2rdata = new I2RData();
             *i2rdata = (*other.i2rdata);
         }

     }
     return *this;
 }
 
 ReconstructInitSIRENA::~ReconstructInitSIRENA()
 {
     if(library_collection) {
         delete library_collection; library_collection = 0;
     }
     if(noise_spectrum) {
         delete noise_spectrum; noise_spectrum = 0;
     }
     if(grading) {
         //delete grading; grading = 0;

         //gsl_vector_free(grading->value); grading->value = 0; //Error
         //gsl_matrix_free(grading->gradeData); grading->gradeData = 0;

         //delete [] grading; // Aumentan en 2 los errores en valgrind

          //delete [] grading->value; grading->value = 0; // Aumentan en 2 los errores en valgrind
          //delete [] grading->gradeData; grading->gradeData = 0;
          //delete grading; grading = 0;
     }
     if(i2rdata) {
         delete i2rdata; i2rdata = 0;
     }
 }
 
 /* This method copies the data from the object to a new object except
  *   for the library and the record file */
 ReconstructInitSIRENA* ReconstructInitSIRENA::get_threading_object(int n_record)
 {
     //log_trace("Getting the reconstruction structure from record %i", n_record);
     ReconstructInitSIRENA* ret = new ReconstructInitSIRENA;
     if(this->library_collection){
         ret->library_collection = this->library_collection;
     }
     
     ret->threshold = this->threshold;
     strcpy(ret->library_file, this->library_file);
     strcpy(ret->record_file, this->record_file);
     
     //record_file_fptr
     // Here we copy the ptr of the fits file, this is NOT thread safe,
     // even allowing reentrant here we should open the file again,
     // and even by doing that the writing through multiple threads
     // won't be safe
     ret->record_file_fptr = this->record_file_fptr;
     
     strcpy(ret->noise_file, this->noise_file);
     strcpy(ret->event_file, this->event_file);
     
     ret->pulse_length = this->pulse_length;
     ret->scaleFactor = this->scaleFactor;
     ret->samplesUp = this->samplesUp;
     ret->samplesDown = this->samplesDown;
     ret->nSgms = this->nSgms;
     ret->detectSP = this->detectSP;
     ret->monoenergy = this->monoenergy;
     ret->hduPRECALWN = this->hduPRECALWN;
     ret->hduPRCLOFWM = this->hduPRCLOFWM;
     ret->largeFilter = this->largeFilter;
     ret->LrsT = this->LrsT;
     ret->LbT = this->LbT;
     ret->opmode = this->opmode;
     strcpy(ret->detectionMode, this->detectionMode);
     
     //NoiseSpec
     if(this->noise_spectrum){
         ret->noise_spectrum = new NoiseSpec();
         *ret->noise_spectrum = (*this->noise_spectrum);
     }
     
     strcpy(ret->FilterDomain, this->FilterDomain);
     strcpy(ret->FilterMethod, this->FilterMethod);
     strcpy(ret->EnergyMethod, this->EnergyMethod);
     
     ret->filtEev = this->filtEev;
     
     ret->Ifit = this->Ifit;
     
     strcpy(ret->OFNoise, this->OFNoise);
     
     ret->LagsOrNot = this->LagsOrNot;
     ret->nLags = this->nLags;
     ret->Fitting35 = this->Fitting35;
     ret->OFIter = this->OFIter;
     ret->OFLib = this->OFLib;
     
     strcpy(ret->OFInterp, this->OFInterp);
     strcpy(ret->OFStrategy, this->OFStrategy);
     
     ret->OFLength = this->OFLength;
     ret->preBuffer = this->preBuffer;
     ret->intermediate = this->intermediate;
     strcpy(ret->detectFile, this->detectFile);
     //sprintf(ret->detectFile, "%s_%i", ret->detectFile, n_record);
     strcat(ret->detectFile,"_");
     strcat(ret->detectFile,to_string(n_record).c_str());

     ret->errorT = this->errorT;
     ret->Sum0Filt = this->Sum0Filt;

     ret->clobber = this->clobber;
     ret->maxPulsesPerRecord = this->maxPulsesPerRecord;
     ret->SaturationValue = this->SaturationValue;
     strcpy(ret->tstartPulse1, this->tstartPulse1);
     ret->tstartPulse2 = this->tstartPulse2;
     ret->tstartPulse3 = this->tstartPulse3;
     ret->energyPCA1 = this->energyPCA1;
     ret->energyPCA2 = this->energyPCA2;

     strcpy(ret->XMLFile, this->XMLFile);

     //Grading
     if(this->grading){
         ret->grading = new Grading();
         *ret->grading = (*this->grading);
     }
     //Grading
     if(this->i2rdata){
         ret->i2rdata = new I2RData();
         *ret->i2rdata = (*this->i2rdata);
     }
     return ret;
 }
 
 // LibraryCollection
 LibraryCollection::LibraryCollection():
 ntemplates(0),
 nfixedfilters(0),
 energies(0),
 pulse_heights(0),
 pulse_templatesMaxLengthFixedFilter(0),
 pulse_templates(0),
 pulse_templates_filder(0),
 maxDERs(0),
 samp1DERs(0),
 pulse_templates_B0(0),
 matched_filters(0),
 matched_filters_B0(0),
 optimal_filters(0),
 optimal_filtersFREQ(0),
 optimal_filtersTIME(0),
 V(0),
 W(0),
 WAB(0),
 T(0),
 t(0),
 X(0),
 Y(0),
 Z(0),
 r(0),
 PAB(0),
 PABMXLFF(0),
 DAB(0),
 optimal_filtersab(0),
 optimal_filtersabTIME(0),
 optimal_filtersabFREQ(0),
 PRECALWN(0),
 PRCLOFWM(0)
 {
     
 }
 
 LibraryCollection::LibraryCollection(const LibraryCollection& other):
 ntemplates(other.ntemplates),
 nfixedfilters(other.nfixedfilters),
 energies(0),
 pulse_heights(0),
 pulse_templatesMaxLengthFixedFilter(0),
 pulse_templates(0),
 pulse_templates_filder(0),
 maxDERs(0),
 samp1DERs(0),
 pulse_templates_B0(0),
 matched_filters(0),
 matched_filters_B0(0),
 optimal_filters(0),
 optimal_filtersFREQ(0),
 optimal_filtersTIME(0),
 V(0),
 W(0),
 WAB(0),
 T(0),
 t(0),
 X(0),
 Y(0),
 Z(0),
 r(0),
 PAB(0),
 PABMXLFF(0),
 DAB(0),
 optimal_filtersab(0),
 optimal_filtersabTIME(0),
 optimal_filtersabFREQ(0),
 PRECALWN(0),
 PRCLOFWM(0)
 {
     if(other.energies){
         energies = gsl_vector_alloc(other.energies->size);
         gsl_vector_memcpy(energies, other.energies);
     }
     if(other.pulse_heights){
         pulse_heights = gsl_vector_alloc(other.pulse_heights->size);
         gsl_vector_memcpy(pulse_heights, other.pulse_heights);
     }
     if(other.pulse_templatesMaxLengthFixedFilter){
         pulse_templatesMaxLengthFixedFilter = new PulseTemplate[ntemplates];
         for (int i = 0; i < ntemplates; ++i){
             pulse_templatesMaxLengthFixedFilter[i] = 
             other.pulse_templatesMaxLengthFixedFilter[i];
         }
     }
     if(pulse_templates){
         pulse_templates = new PulseTemplate[ntemplates];
         for (int i = 0; i < ntemplates; ++i){
             pulse_templates[i] = other.pulse_templates[i];
         }
     }
     if (other.pulse_templates_filder){
         pulse_templates_filder = new PulseTemplate[ntemplates];
         for (int i = 0; i < ntemplates; ++i){
             pulse_templates_filder[i] = other.pulse_templates_filder[i];
         }
     } 
     if (other.maxDERs){
         maxDERs = gsl_vector_alloc(other.maxDERs->size);
         gsl_vector_memcpy(maxDERs, other.maxDERs);
     }
     if(other.samp1DERs){
         samp1DERs = gsl_vector_alloc(other.samp1DERs->size);
         gsl_vector_memcpy(samp1DERs, other.samp1DERs);
     }
     if(other.pulse_templates_B0){
         pulse_templates_B0 = new PulseTemplate[ntemplates];
         for (int i = 0; i < ntemplates; ++i){
             pulse_templates_B0[i] = other.pulse_templates_B0[i];
         }
     }
     if(other.matched_filters){
         matched_filters = new MatchedFilter[ntemplates];
         for (int i = 0; i < ntemplates; ++i){
             matched_filters[i] = other.matched_filters[i];
         }
     }
     if(other.matched_filters_B0){
         matched_filters_B0 = new MatchedFilter[ntemplates];
         for (int i = 0; i < ntemplates; ++i){
             matched_filters_B0[i] = other.matched_filters_B0[i];
         }
     }
     if(other.optimal_filters){
         optimal_filters = new OptimalFilterSIRENA[ntemplates];
         for (int i = 0; i < ntemplates; ++i){
             optimal_filters[i] = other.optimal_filters[i];
         }
     }
     if(other.optimal_filtersFREQ){
         optimal_filtersFREQ = new OptimalFilterSIRENA[ntemplates];
         for (int i = 0; i < ntemplates; ++i){
             optimal_filtersFREQ[i] = other.optimal_filtersFREQ[i];
         }
     }
     if (other.optimal_filtersTIME){
         optimal_filtersTIME = new OptimalFilterSIRENA[ntemplates];
         for (int i = 0; i < ntemplates; ++i){
             optimal_filtersTIME[i] = other.optimal_filtersTIME[i];
         }
     }
     
     if(other.V){
         V = gsl_matrix_alloc(other.V->size1,
                              other.V->size2);
         gsl_matrix_memcpy(V, other.V);
     }
     if (other.W){
         W = gsl_matrix_alloc(other.W->size1,
                              other.W->size2);
         gsl_matrix_memcpy(W, other.W);    
     }
     
     if(other.WAB){
         WAB = gsl_matrix_alloc(other.WAB->size1,
                                other.WAB->size2);
         gsl_matrix_memcpy(WAB, other.WAB);
     }
     
     if(other.T){
         T = gsl_matrix_alloc(other.T->size1,
                              other.T->size2);
         gsl_matrix_memcpy(T, other.T);
     }
     
     if(other.t){
         t = gsl_vector_alloc(other.t->size);
         gsl_vector_memcpy(t, other.t);
     }
     
     if(other.X){
         X = gsl_matrix_alloc(other.X->size1,
                              other.X->size2);
         gsl_matrix_memcpy(X, other.X);
     }
     
     if(other.Y){
         Y = gsl_matrix_alloc(other.Y->size1,
                              other.Y->size2);
         gsl_matrix_memcpy(Y, other.Y);
     }
     
     if(other.Z){
         Z = gsl_matrix_alloc(other.Z->size1,
                              other.Z->size2);
         gsl_matrix_memcpy(Z, other.Z);
     }
     
     if(other.r){
         r = gsl_vector_alloc(other.r->size);
         gsl_vector_memcpy(r, other.r);
     }
     
     if(other.PAB){
         PAB = gsl_matrix_alloc(other.PAB->size1,
                                other.PAB->size2);
         gsl_matrix_memcpy(PAB, other.PAB);
     }
     
     if(other.PABMXLFF){
         PABMXLFF = gsl_matrix_alloc(other.PABMXLFF->size1,
                                     other.PABMXLFF->size2);
         gsl_matrix_memcpy(PABMXLFF, other.PABMXLFF);
     }
     
     if(other.DAB){
         DAB = gsl_matrix_alloc(other.DAB->size1,
                                other.DAB->size2);
         gsl_matrix_memcpy(DAB, other.DAB);
     }
     
     if(other.optimal_filtersab){
         optimal_filtersab = new OptimalFilterSIRENA[ntemplates];
         for (int i = 0; i < ntemplates; ++i){
             optimal_filtersab[i] = other.optimal_filtersab[i];
         }
     }
     if(other.optimal_filtersabTIME){
         optimal_filtersabTIME = new OptimalFilterSIRENA[ntemplates];
         for (int i = 0; i < ntemplates; ++i){
             optimal_filtersabTIME[i] = other.optimal_filtersabTIME[i];
         }
     }
     if(other.optimal_filtersabFREQ){
         optimal_filtersabFREQ = new OptimalFilterSIRENA[ntemplates];
         for (int i = 0; i < ntemplates; ++i){
             optimal_filtersabFREQ[i] = other.optimal_filtersabFREQ[i];
         }
     }
     
     if(other.PRECALWN){
         PRECALWN = gsl_matrix_alloc(other.PRECALWN->size1,
                                     other.PRECALWN->size2);
         gsl_matrix_memcpy(PRECALWN, other.PRECALWN);
     }
     
     if(other.PRCLOFWM){
         PRCLOFWM = gsl_matrix_alloc(other.PRCLOFWM->size1,
                                     other.PRCLOFWM->size2);
         gsl_matrix_memcpy(PRCLOFWM, other.PRCLOFWM);
     }
 }
 
 LibraryCollection& LibraryCollection::operator=(const LibraryCollection& other)
 {
     if (this != &other){
         ntemplates = other.ntemplates;
         nfixedfilters = other.nfixedfilters;
         
         if(energies) {
             gsl_vector_free(energies); energies = 0;
         }
         if(other.energies){
             energies = gsl_vector_alloc(other.energies->size);
             gsl_vector_memcpy(energies, other.energies);
         }
         if(pulse_heights) { 
             gsl_vector_free(pulse_heights); pulse_heights = 0;
         }
         if(other.pulse_heights){
             pulse_heights = gsl_vector_alloc(other.pulse_heights->size);
             gsl_vector_memcpy(pulse_heights, other.pulse_heights);
         }
         if(ntemplates > 0 && pulse_templatesMaxLengthFixedFilter) {
             delete [] pulse_templatesMaxLengthFixedFilter;
             pulse_templatesMaxLengthFixedFilter = 0;
         }
         if(other.pulse_templatesMaxLengthFixedFilter){
             pulse_templatesMaxLengthFixedFilter = new PulseTemplate[ntemplates];
             for (int i = 0; i < ntemplates; ++i){
                 pulse_templatesMaxLengthFixedFilter[i] = 
                 other.pulse_templatesMaxLengthFixedFilter[i];
             }
         }
         if(ntemplates > 0 && pulse_templates) {
             delete [] pulse_templates; pulse_templates = 0;
         }
         if(pulse_templates){
             pulse_templates = new PulseTemplate[ntemplates];
             for (int i = 0; i < ntemplates; ++i){
                 pulse_templates[i] = other.pulse_templates[i];
             }
         }
         if(ntemplates > 0 && pulse_templates_filder) {
             delete [] pulse_templates_filder; pulse_templates_filder = 0;
         }
         if (other.pulse_templates_filder){
             pulse_templates_filder = new PulseTemplate[ntemplates];
             for (int i = 0; i < ntemplates; ++i){
                 pulse_templates_filder[i] = other.pulse_templates_filder[i];
             }
         }
         
         if(maxDERs) {
             gsl_vector_free(maxDERs); maxDERs = 0;
         }
         if (other.maxDERs){
             maxDERs = gsl_vector_alloc(other.maxDERs->size);
             gsl_vector_memcpy(maxDERs, other.maxDERs);
         }
         
         if(samp1DERs) {
             gsl_vector_free(samp1DERs); samp1DERs = 0;
         }
         if(other.samp1DERs){
             samp1DERs = gsl_vector_alloc(other.samp1DERs->size);
             gsl_vector_memcpy(samp1DERs, other.samp1DERs);
         }
         if(ntemplates > 0 && pulse_templates_B0){
             delete [] pulse_templates_B0; pulse_templates_B0 = 0;
         }
         if(other.pulse_templates_B0){
             pulse_templates_B0 = new PulseTemplate[ntemplates];
             for (int i = 0; i < ntemplates; ++i){
                 pulse_templates_B0[i] = other.pulse_templates_B0[i];
             }
         }
         if(ntemplates > 0 && matched_filters) {
             delete [] matched_filters; matched_filters = 0;
         }
         if(other.matched_filters){
             matched_filters = new MatchedFilter[ntemplates];
             for (int i = 0; i < ntemplates; ++i){
                 matched_filters[i] = other.matched_filters[i];
             }
         }
         if(ntemplates > 0 && matched_filters_B0){ 
             delete [] matched_filters_B0; matched_filters_B0 = 0;
         }
         if(other.matched_filters_B0){
             matched_filters_B0 = new MatchedFilter[ntemplates];
             for (int i = 0; i < ntemplates; ++i){
                 matched_filters_B0[i] = other.matched_filters_B0[i];
             }
         }
         if(ntemplates > 0 && optimal_filters) {
             delete [] optimal_filters; optimal_filters = 0;
         }
         if(other.optimal_filters){
             optimal_filters = new OptimalFilterSIRENA[ntemplates];
             for (int i = 0; i < ntemplates; ++i){
                 optimal_filters[i] = other.optimal_filters[i];
             }
         }
         if(ntemplates > 0 && optimal_filtersFREQ) {
             delete [] optimal_filtersFREQ; optimal_filtersFREQ = 0;
         }
         if(other.optimal_filtersFREQ){
             optimal_filtersFREQ = new OptimalFilterSIRENA[ntemplates];
             for (int i = 0; i < ntemplates; ++i){
                 optimal_filtersFREQ[i] = other.optimal_filtersFREQ[i];
             }
         }
         if(ntemplates > 0 && optimal_filtersTIME) {
             delete [] optimal_filtersTIME; optimal_filtersTIME = 0;
         }
         if (other.optimal_filtersTIME){
             optimal_filtersTIME = new OptimalFilterSIRENA[ntemplates];
             for (int i = 0; i < ntemplates; ++i){
                 optimal_filtersTIME[i] = other.optimal_filtersTIME[i];
             }
         }
         
         if(V) {
             gsl_matrix_free(V); V = 0;
         }
         if(other.V){
             V = gsl_matrix_alloc(other.V->size1,
                                  other.V->size2);
             gsl_matrix_memcpy(V, other.V);
         }
         
         if(W) {
             gsl_matrix_free(W); W = 0;
         }
         if (other.W){
             W = gsl_matrix_alloc(other.W->size1,
                                  other.W->size2);
             gsl_matrix_memcpy(W, other.W);    
         }
         
         if(WAB) {
             gsl_matrix_free(WAB); WAB = 0;
         }
         if(other.WAB){
             WAB = gsl_matrix_alloc(other.WAB->size1,
                                    other.WAB->size2);
             gsl_matrix_memcpy(WAB, other.WAB);
         }
         
         if(T) {
             gsl_matrix_free(T); T = 0;
         }
         if(other.T){
             T = gsl_matrix_alloc(other.T->size1,
                                  other.T->size2);
             gsl_matrix_memcpy(T, other.T);
         }
         
         if(t) {
             gsl_vector_free(t); t = 0;
         }
         if(other.t){
             t = gsl_vector_alloc(other.t->size);
             gsl_vector_memcpy(t, other.t);
         }
         
         if(X) {
             gsl_matrix_free(X); X = 0;
         }
         if(other.X){
             X = gsl_matrix_alloc(other.X->size1,
                                  other.X->size2);
             gsl_matrix_memcpy(X, other.X);
         }
         
         if(Y) {
             gsl_matrix_free(Y); Y = 0;
         }
         if(other.Y){
             Y = gsl_matrix_alloc(other.Y->size1,
                                  other.Y->size2);
             gsl_matrix_memcpy(Y, other.Y);
         }
         
         if(Z) {
             gsl_matrix_free(Z); Z = 0;
         }
         if(other.Z){
             Z = gsl_matrix_alloc(other.Z->size1,
                                  other.Z->size2);
             gsl_matrix_memcpy(Z, other.Z);
         }
         
         if(r) {
             gsl_vector_free(r); r = 0;
         }
         if(other.r){
             r = gsl_vector_alloc(other.r->size);
             gsl_vector_memcpy(r, other.r);
         }
         
         if(PAB) {
             gsl_matrix_free(PAB); PAB = 0;
         }
         if(other.PAB){
             PAB = gsl_matrix_alloc(other.PAB->size1,
                                    other.PAB->size2);
             gsl_matrix_memcpy(PAB, other.PAB);
         }
         
         if(PABMXLFF) {
             gsl_matrix_free(PABMXLFF); PABMXLFF = 0;
         }
         if(other.PABMXLFF){
             PABMXLFF = gsl_matrix_alloc(other.PABMXLFF->size1,
                                         other.PABMXLFF->size2);
             gsl_matrix_memcpy(PABMXLFF, other.PABMXLFF);
         }
         
         if(DAB) {
             gsl_matrix_free(DAB); DAB = 0;
         }
         if(other.DAB){
             DAB = gsl_matrix_alloc(other.DAB->size1,
                                    other.DAB->size2);
             gsl_matrix_memcpy(DAB, other.DAB);
         }
         
         if(ntemplates > 0 && optimal_filtersab) {
             delete [] optimal_filtersab; optimal_filtersab = 0;
         }
         if(other.optimal_filtersab){
             optimal_filtersab = new OptimalFilterSIRENA[ntemplates];
             for (int i = 0; i < ntemplates; ++i){
                 optimal_filtersab[i] = other.optimal_filtersab[i];
             }
         }
         if(ntemplates > 0 && optimal_filtersabTIME) {
             delete [] optimal_filtersabTIME; optimal_filtersabTIME = 0;
         }
         if(other.optimal_filtersabTIME){
             optimal_filtersabTIME = new OptimalFilterSIRENA[ntemplates];
             for (int i = 0; i < ntemplates; ++i){
                 optimal_filtersabTIME[i] = other.optimal_filtersabTIME[i];
             }
         }
         if(ntemplates > 0 && optimal_filtersabFREQ) {
             delete [] optimal_filtersabFREQ; optimal_filtersabFREQ = 0;
         }
         if(other.optimal_filtersabFREQ){
             optimal_filtersabFREQ = new OptimalFilterSIRENA[ntemplates];
             for (int i = 0; i < ntemplates; ++i){
                 optimal_filtersabFREQ[i] = other.optimal_filtersabFREQ[i];
             }
         }
         
         if(PRECALWN) {
             gsl_matrix_free(PRECALWN); PRECALWN = 0;
         }
         if(other.PRECALWN){
             PRECALWN = gsl_matrix_alloc(other.PRECALWN->size1,
                                         other.PRECALWN->size2);
             gsl_matrix_memcpy(PRECALWN, other.PRECALWN);
         }
         
         if(PRCLOFWM) {
             gsl_matrix_free(PRCLOFWM); PRCLOFWM = 0;
         }
         if(other.PRCLOFWM){
             PRCLOFWM = gsl_matrix_alloc(other.PRCLOFWM->size1,
                                         other.PRCLOFWM->size2);
             gsl_matrix_memcpy(PRCLOFWM, other.PRCLOFWM);
         }
     }
     return *this;
 }
 
 LibraryCollection::~LibraryCollection()
 {
     if(energies) {
         gsl_vector_free(energies); energies = 0;
     }
     if(pulse_heights) {
         gsl_vector_free(pulse_heights); pulse_heights = 0;
     }
     
     if(ntemplates > 0 && pulse_templatesMaxLengthFixedFilter){
         delete [] pulse_templatesMaxLengthFixedFilter;
         pulse_templatesMaxLengthFixedFilter = 0;
     }
     if(ntemplates > 0 && pulse_templates){
         delete [] pulse_templates;
         pulse_templates = 0;
     }
     if(ntemplates > 0 && pulse_templates_filder){
         delete [] pulse_templates_filder;
         pulse_templates_filder = 0;
     }
     
     if(maxDERs) {
         gsl_vector_free(maxDERs); maxDERs = 0;
     }
     if(samp1DERs) {
         gsl_vector_free(samp1DERs); samp1DERs = 0;
     }
     
     if(ntemplates > 0 && pulse_templates_B0){
         delete [] pulse_templates_B0;
         pulse_templates_B0 = 0;
     }
     if(ntemplates > 0 && matched_filters){
         delete [] matched_filters;
         matched_filters = 0;
     }
     if(ntemplates > 0 && matched_filters_B0){
         delete [] matched_filters_B0;
         matched_filters_B0 = 0;
     }
     if(ntemplates > 0 && optimal_filters){
         delete [] optimal_filters;
         optimal_filters = 0;
     }
     if(ntemplates > 0 && optimal_filtersFREQ){
         delete [] optimal_filtersFREQ;
         optimal_filtersFREQ = 0;
     }
     if(ntemplates > 0 && optimal_filtersTIME){
         delete [] optimal_filtersTIME;
         optimal_filtersTIME = 0;
     }
     
     if(V) {
         gsl_matrix_free(V); V = 0;
     }
     if(W) {
         gsl_matrix_free(W); W = 0;
     }
     if(WAB) { 
         gsl_matrix_free(WAB); WAB = 0;
     }
     if(T) {
         gsl_matrix_free(T); T = 0;
     }
     if(t) {
         gsl_vector_free(t); t = 0;
     }
     if(X) {
         gsl_matrix_free(X); X = 0;
     }
     if(Y) {
         gsl_matrix_free(Y); Y = 0;
     }
     if(Z) {
         gsl_matrix_free(Z); Z = 0;
     }
     if(r) {
         gsl_vector_free(r); r = 0;
     }
     if(PAB) {
         gsl_matrix_free(PAB); PAB = 0;
     }
     if(PABMXLFF) {
         gsl_matrix_free(PABMXLFF); PABMXLFF = 0;
     }
     if(DAB) {
         gsl_matrix_free(DAB); DAB = 0;
     }
     
     if(ntemplates > 0 && optimal_filtersab){
         delete [] optimal_filtersab;
     }
     if(ntemplates > 0 && optimal_filtersabFREQ){
         delete [] optimal_filtersabFREQ;
     }
     if(ntemplates > 0 && optimal_filtersabTIME){
         delete [] optimal_filtersabTIME;
     }
     
     if(PRECALWN) { 
         gsl_matrix_free(PRECALWN); PRECALWN = 0;
     }
     if(PRCLOFWM) {
         gsl_matrix_free(PRCLOFWM); PRCLOFWM = 0;
     }
 }
 
 // PulseDetected
 PulseDetected::PulseDetected():
 pulse_duration(0),
 grade1(0),
 grade2(0),
 grade2_1(0),
 pixid(0),
 phid(0),
 phid2(0),
 phid3(0),
 pulse_adc(0),
 pulse_adc_preBuffer(0),
 Tstart(0.0f),
 TstartSamples(0.0f),
 Tend(0.0f),
 riseTime(0.0f),
 fallTime(0.0f),
 pulse_height(0.0f),
 maxDER(0.0f),
 samp1DER(0.0f),
 energy(0.0f),
 grading(0),
 avg_4samplesDerivative(0.0f),
 E_lowres(0.0f),
 phi(0.0f),
 lagsShift(0),
 quality(0.0f),
 numLagsUsed(0),
 bsln(0.0f),
 rmsbsln(0.0f)
 {

 }
 
 PulseDetected::PulseDetected(const PulseDetected& other):
 pulse_duration(other.pulse_duration),
 grade1(other.grade1),
 grade2(other.grade2),
 grade2_1(other.grade2_1),
 pixid(other.pixid),
 phid(other.phid),
 phid2(other.phid2),
 phid3(other.phid3),
 pulse_adc(0),
 pulse_adc_preBuffer(0),
 Tstart(other.Tstart),
 TstartSamples(other.TstartSamples),
 Tend(other.Tend),
 riseTime(other.riseTime),
 fallTime(other.fallTime),
 pulse_height(other.pulse_height),
 maxDER(other.maxDER),
 samp1DER(other.samp1DER),
 energy(other.energy),
 grading(other.grading),
 avg_4samplesDerivative(other.avg_4samplesDerivative),
 E_lowres(other.E_lowres),
 phi(other.phi),
 lagsShift(other.lagsShift),
 quality(other.quality),
 numLagsUsed(other.numLagsUsed),
 bsln(other.bsln),
 rmsbsln(other.rmsbsln)
 {
     if(other.pulse_adc){
         pulse_adc = gsl_vector_alloc(other.pulse_adc->size);
         gsl_vector_memcpy(pulse_adc, other.pulse_adc);
     }
     if(other.pulse_adc_preBuffer){
         pulse_adc_preBuffer = gsl_vector_alloc(other.pulse_adc_preBuffer->size);
         gsl_vector_memcpy(pulse_adc_preBuffer, other.pulse_adc_preBuffer);
     }
 }
 
 PulseDetected& PulseDetected::operator=(const PulseDetected& other)
 {
     if(this != &other){
         
         pulse_duration = other.pulse_duration;
         grade1 = other.grade1;
         grade2 = other.grade2;
         grade2_1 = other.grade2_1;
         pixid = other.pixid;
         phid = other.phid;
         phid2 = other.phid2;
         phid3 = other.phid3;
         if(pulse_adc) {
             gsl_vector_free(pulse_adc); pulse_adc = 0;
         }
         if(other.pulse_adc){
             pulse_adc = gsl_vector_alloc(other.pulse_adc->size);
             gsl_vector_memcpy(pulse_adc, other.pulse_adc);
         }
         if(pulse_adc_preBuffer) {
             gsl_vector_free(pulse_adc_preBuffer); pulse_adc_preBuffer = 0;
         }
         if(other.pulse_adc_preBuffer){
             pulse_adc_preBuffer = gsl_vector_alloc(other.pulse_adc_preBuffer->size);
             gsl_vector_memcpy(pulse_adc_preBuffer, other.pulse_adc_preBuffer);
         }
         Tstart = other.Tstart;
         TstartSamples = other.TstartSamples;
         Tend = other.Tend;
         riseTime = other.riseTime;
         fallTime = other.fallTime;
         pulse_height = other.pulse_height;
         maxDER = other.maxDER;
         samp1DER = other.samp1DER;
         energy = other.energy;
         grading = other.grading;
         avg_4samplesDerivative = other.avg_4samplesDerivative;
         E_lowres = other.E_lowres;
         phi = other.phi;
         lagsShift = other.lagsShift;
         quality = other.quality;
         numLagsUsed = other.numLagsUsed;
         bsln = other.bsln;
         rmsbsln = other.rmsbsln;
     }
     return *this;
 }
 
 PulseDetected::~PulseDetected()
 {
     if(pulse_adc) {
         gsl_vector_free(pulse_adc); pulse_adc = 0;
     }
     if(pulse_adc_preBuffer) {
         gsl_vector_free(pulse_adc_preBuffer); pulse_adc_preBuffer = 0;
     }
 }
 
 MatrixStruct::MatrixStruct():
 matrixRows(0), 
 matrixColumns(0), 
 matrixBody(0)
 {
     
 }
 
 MatrixStruct::MatrixStruct(const MatrixStruct& other)
 :matrixRows(other.matrixRows),
 matrixColumns(other.matrixColumns),
 matrixBody(0)
 {
     if(other.matrixBody){
         matrixBody = gsl_matrix_alloc(matrixRows, matrixColumns);
         gsl_matrix_memcpy(matrixBody, other.matrixBody);
     }
 }
 
 MatrixStruct& MatrixStruct::operator=(const MatrixStruct& other)
 {
     if(this != &other){
         matrixRows = other.matrixRows;
         matrixColumns = other.matrixColumns;
         if (matrixBody) {
             gsl_matrix_free(matrixBody); matrixBody = 0;
         }
         if(other.matrixBody){
             matrixBody = gsl_matrix_alloc(matrixRows, matrixColumns);
             gsl_matrix_memcpy(matrixBody, other.matrixBody);
         }
     }
     return *this;
 } 
 
 MatrixStruct::~MatrixStruct()
 {
     if(matrixBody) {
         gsl_matrix_free(matrixBody); matrixBody = 0;
     }
 }
 
 PulseTemplate::PulseTemplate():
 template_duration(0), 
 ptemplate(0), 
 energy(0.0f), 
 pulse_height(0.0f)
 {
     
 }
 
 PulseTemplate::PulseTemplate(const PulseTemplate& other):
 template_duration(other.template_duration),
 ptemplate(0),
 energy(other.energy),
 pulse_height(other.pulse_height)
 {
     if (other.ptemplate){
         ptemplate = gsl_vector_alloc(other.ptemplate->size);
         gsl_vector_memcpy(ptemplate, other.ptemplate);
     }
 }
 
 PulseTemplate& PulseTemplate::operator=(const PulseTemplate& other)
 {
     if(this != &other){
         template_duration = other.template_duration;
         if (ptemplate) {
             gsl_vector_free(ptemplate); ptemplate = 0;
         }
         if (other.ptemplate){
             ptemplate = gsl_vector_alloc(other.ptemplate->size);
             gsl_vector_memcpy(ptemplate, other.ptemplate);
         }
         energy = other.energy;
         pulse_height = other.pulse_height;
     }
     return *this;
 }
 
 PulseTemplate::~PulseTemplate()
 {
     if(ptemplate) {
         gsl_vector_free(ptemplate); ptemplate = 0;
     }
 }
 
 MatchedFilter::MatchedFilter():
 mfilter_duration(0),
 mfilter(0),
 energy(0.0f), 
 pulse_height(0.0f)
 {
     
 }
 
 MatchedFilter::MatchedFilter(const MatchedFilter& other):
 mfilter_duration(other.mfilter_duration),
 mfilter(0),
 energy(other.energy),
 pulse_height(other.pulse_height)
 {
     if(other.mfilter){
         mfilter = gsl_vector_alloc(other.mfilter->size);
         gsl_vector_memcpy(mfilter, other.mfilter);
     }
 }
 
 MatchedFilter& MatchedFilter::operator=(const MatchedFilter& other)
 {
     if(this != &other){
         mfilter_duration = other.mfilter_duration;
         if (mfilter) {
             gsl_vector_free(mfilter); mfilter = 0;
         }
         if(other.mfilter){
             mfilter = gsl_vector_alloc(other.mfilter->size);
             gsl_vector_memcpy(mfilter, other.mfilter);
         }
         energy = other.energy;
         pulse_height = other.pulse_height;
     }
     return *this;
 }
 
 MatchedFilter::~MatchedFilter()
 {
     if(mfilter) {
         gsl_vector_free(mfilter); mfilter = 0;
     }
 }
 
 OptimalFilterSIRENA::OptimalFilterSIRENA():
 ofilter_duration(0),
 ofilter(0), 
 energy(0.0f)
 {
     
 }
 
 OptimalFilterSIRENA::OptimalFilterSIRENA(const OptimalFilterSIRENA& other):
 ofilter_duration(other.ofilter_duration),
 ofilter(0),
 energy(other.energy)
 {
     if(other.ofilter){
         ofilter = gsl_vector_alloc(other.ofilter->size);
         gsl_vector_memcpy(ofilter, other.ofilter);
     }
 }
 
 OptimalFilterSIRENA&
 OptimalFilterSIRENA::operator=(const OptimalFilterSIRENA& other)
 {
     if(this != &other){
         ofilter_duration = other.ofilter_duration;
         if (ofilter) {
             gsl_vector_free(ofilter); ofilter = 0;
         }
         if(other.ofilter){
             ofilter = gsl_vector_alloc(other.ofilter->size);
             gsl_vector_memcpy(ofilter, other.ofilter);
         }
         energy = other.energy;
     }
     return *this;
 }
 
 OptimalFilterSIRENA::~OptimalFilterSIRENA()
 {
     if(ofilter) {
         gsl_vector_free(ofilter); ofilter = 0;
     }
 }
 
 NoiseSpec::NoiseSpec():
 noiseStd(0.0f), 
 baseline(0.0f), 
 noise_duration(0), 
 noisespec(0),
 noisefreqs(0), 
 weightMatrixes(0)
 {
     
 }  
 
 NoiseSpec::NoiseSpec(const NoiseSpec& other):
 noiseStd(other.noiseStd), 
 baseline(other.baseline), 
 noise_duration(other.noise_duration), 
 noisespec(0),
 noisefreqs(0), 
 weightMatrixes(0)
 {
     if(other.noisespec){
         noisespec = gsl_vector_alloc(other.noisespec->size);
         gsl_vector_memcpy(noisespec, other.noisespec);
     }
     if(other.noisefreqs){
         noisefreqs = gsl_vector_alloc(other.noisefreqs->size);
         gsl_vector_memcpy(noisefreqs, other.noisefreqs);
     }
     if(other.weightMatrixes){
         weightMatrixes = gsl_matrix_alloc(other.weightMatrixes->size1,
                                           other.weightMatrixes->size2);
         gsl_matrix_memcpy(weightMatrixes, other.weightMatrixes);
     }
 }
 
 NoiseSpec& NoiseSpec::operator=(const NoiseSpec& other)
 {
     if(this != &other){
         noiseStd = other.noiseStd;
         baseline = other.baseline;
         noise_duration = other.noise_duration;
         if (noisespec) {
             gsl_vector_free(noisespec); noisespec = 0;
         }
         if(other.noisespec){
             noisespec = gsl_vector_alloc(other.noisespec->size);
             gsl_vector_memcpy(noisespec, other.noisespec);
         }
         if (noisefreqs) {
             gsl_vector_free(noisefreqs); noisefreqs = 0;
         }
         if(other.noisefreqs){
             noisefreqs = gsl_vector_alloc(other.noisefreqs->size);
             gsl_vector_memcpy(noisefreqs, other.noisefreqs);
         }
         if (weightMatrixes) {
             gsl_matrix_free(weightMatrixes); weightMatrixes = 0;
         }
         if(other.weightMatrixes){
             weightMatrixes = gsl_matrix_alloc(other.weightMatrixes->size1,
                                               other.weightMatrixes->size2);
             gsl_matrix_memcpy(weightMatrixes, other.weightMatrixes);
         }
     }
     return *this;
 }
 NoiseSpec::~NoiseSpec()
 {
     if(noisespec) {
         gsl_vector_free(noisespec); noisespec = 0;
     }
     if(noisefreqs){
         gsl_vector_free(noisefreqs); noisefreqs = 0;
     }
     if(weightMatrixes) {
         gsl_matrix_free(weightMatrixes); weightMatrixes = 0;
     }
 }
 
 Grading::Grading():
 ngrades(0), 
 value(0), 
 gradeData(0)
 {
     
 }
 
 Grading::Grading(const Grading& other):
 ngrades(other.ngrades), 
 value(0), 
 gradeData(0)
 {
     if(other.value){
         value = gsl_vector_alloc(other.value->size);
         gsl_vector_memcpy(value, other.value);
     }
     if(other.gradeData){
         gradeData = gsl_matrix_alloc(other.gradeData->size1,
                                      other.gradeData->size2);
         gsl_matrix_memcpy(gradeData, other.gradeData);
     }
 }
 
 Grading& Grading::operator=(const Grading& other)
 {
     if(this != &other){
         ngrades = other.ngrades;
         if (value) {
             gsl_vector_free(value); value = 0;
         }
         if(other.value){
             value = gsl_vector_alloc(other.value->size);
             gsl_vector_memcpy(value, other.value);
         }
         if (gradeData) {
             gsl_matrix_free(gradeData); gradeData = 0;
         }
         if(other.gradeData){
             gradeData = gsl_matrix_alloc(other.gradeData->size1,
                                          other.gradeData->size2);
             gsl_matrix_memcpy(gradeData, other.gradeData);
         }
     }
     return *this;
 }
 Grading::~Grading()
 {
     if(value) {
         gsl_vector_free(value); value = 0;
     }
     if(gradeData) {
         gsl_matrix_free(gradeData); gradeData = 0;
     }
 }
 
 I2RData::I2RData():
 I0_START(0.0f),
 IMIN(0.0f), 
 IMAX(0.0f)
 {
     
 }
 
 I2RData::I2RData(const I2RData& other):
 I0_START(other.I0_START), 
 IMIN(other.IMIN),
 IMAX(other.IMAX)
 {
     
     
 }
 
 I2RData& I2RData::operator=(const I2RData& other)
 {
     if(this != &other){
         I0_START = other.I0_START;
         IMIN = other.IMIN;
         IMAX = other.IMAX;
     }
     return *this;
 }
 
 I2RData::~I2RData()
 {
 }
 
 PulsesCollection::PulsesCollection():
 ndetpulses(0),
 size(0),
 pulses_detected(0)
 {

 }
 
 PulsesCollection::PulsesCollection(const PulsesCollection& other):
 ndetpulses(other.ndetpulses),
 size(other.size),
 pulses_detected(0)
 {
     if(other.pulses_detected){
         pulses_detected = new PulseDetected[size];
         for (int i = 0; i < ndetpulses; ++i){
             pulses_detected[i] = other.pulses_detected[i];
         }
     }
 }
 
 PulsesCollection& PulsesCollection::operator=(const PulsesCollection& other)
 {
     if(this != &other){
         ndetpulses = other.ndetpulses;
         size = other.size;
         if (pulses_detected) {
             delete [] pulses_detected; pulses_detected = 0;
         }
         if(other.pulses_detected){
             pulses_detected = new PulseDetected[size];
             for (int i = 0; i < ndetpulses; ++i){
                 pulses_detected[i] = other.pulses_detected[i];
             }
         }
     }
     return *this;
 }
 
 PulsesCollection::~PulsesCollection()
 {
     if(ndetpulses > 0 && pulses_detected) {
         delete [] pulses_detected; pulses_detected = 0;
     }
 }
 
