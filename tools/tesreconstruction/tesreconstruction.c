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

   Copyright 2015 Philippe Peille, IRAP
   Copyright 2014:  TASKSSIRENA has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01, 
   ESP2013-48637-C2-1-P, ESP2014-53672-C3-1-P and RTI2018-096686-B-C21.

***********************************************************************
*                      TESRECONSTRUCTION
*
*  File:       tesreconstruction.c
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*              Philippe Peille, IRAP
*                                                                     
***********************************************************************/

#include "tesreconstruction.h"

/***** SECTION 1 ************************************
* MAIN function: This function is mainly a wrapper to pass a data file to the SIRENA tasks in order to reconstruct the energy of the incoming X-ray photons after their detection.
* 
* It can run the SIRENA tasks or the Philippe Peille's tasks depending on the 'Rcmethod' selected.
* 	
* The user must supply the following input parameters (.par file).
* 
* Common parameters:
* 
* - Rcmethod: Reconstruction method (PP or SIRENA)
* - SIRENA: If Rcmethod starts with '@' it provides a file text containing several record input FITS files
* - RecordFile: Record FITS file
* - TesEventFile: Output event list file
* - EventListSize: Default size of the event list
* - clobber:Overwrite or not output files if exist (1/0)
* - history: write program parameters into output file
* 
* PP parameters:
* 
* - SaturationValue: Saturation level of the ADC curves
* - OptimalFilterFile: Optimal filters file
* - PulseTemplateFile: Pulse template file
* - Threshold: Threshold level
* - Calfac: Calibration factor (should be read from the xml file)
* - NormalExclusion: Minimal distance before using OFs after a misreconstruction
* - DerivateExclusion: Minimal distance before reconstructing any event after a misreconstruction
*
* SIRENA parameters:
* 
* - LibraryFile: File with calibration library
* - scaleFactor: Detection scale factor for initial filtering
* - samplesUp: Number of consecutive samples up for threshold trespassing (only used in calibration run, and in production run with STC detection mode)
* - samplesDown: Number of consecutive samples below the threshold to look for other pulse (only used in production run with STC detection mode)
* - nSgms: Number of quiescent-signal standard deviations to establish the threshold through the kappa-clipping algorithm
* - detectSP: Detect secondary pulses (1) or not (0)
* - LrsT: Running sum length for the RS raw energy estimation (seconds) (only for library creation)
* - LbT: Baseline averaging length (seconds)
* - monoenergy: Monochromatic energy of the pulses in the input FITS file in eV (only for library creation)
* - hduPRECALWN: Add or not the PRECALWN HDU in the library file (1/0) (only for library creation)
* - hduPRCLOFWM: Add or not the PRCLOFWM HDU in the library file (1/0) (only for library creation)
* - largeFilter: Length of the longest fixed filter (only for library creation)
* - opmode: Calibration run (0) or energy reconstruction run (1)
* - detectionMode: Adjusted Derivative (AD) or Single Threshold Crossing (STC)
* - NoiseFile: Noise FITS file with noise spectrum
* - FilterDomain: Filtering Domain: Time (T) or Frequency (F)
* - FilterMethod: Filtering Method: F0 (deleting the zero frequency bin) or B0 (deleting the baseline)
* - EnergyMethod: Energy calculation Method: OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL, I2RFITTED, I2RDER or PCA
* - filtEeV: Energy of the filters of the library to be used to calculate energy (only for OPTFILT, I2R, I2RFITTED and I2RDER)
* - Ifit: Constant to apply the I2RFITTED conversion
* - OFNoise: Noise to use with Optimal Filtering: NSD or WEIGHTM
* - LagsOrNot: Lags or no lags (1/0)
* - nLags: Number of lags (positive odd number)
* - Fitting35: Number of lags to analytically calculate a parabola (3) or to fit a parabola (5)
* - OFIter: Iterate or not iterate (1/0) 
* - OFLib: Work or not with a library (1/0)
* - OFStrategy: Optimal Filter length Strategy: FREE, BYGRADE or FIXED
* - OFLength: Optimal Filter length (taken into account if OFStrategy=FIXED)
* - OFLengthNotPadded: Filter length not padded with 0s (only necessary when reconstructing with 0-padding)
* - preBuffer: Some samples added before the starting time of a pulse (number of samples added read from the xml file)
*              SIRENA's format XML file (grading=>pre,post and pB) or new format XML file (grading=>pre,post and filtlen)
*                                      pre=494, post=8192, pB=1000                          pre=494, post=7192, filtlen=8192
*                                                                                             preBuffer=filtlen-post
* - intermediate: Write or not intermediate files (1/0)
* - detectFile: Intermediate detections file (if intermediate=1)
* - errorT: Additional error (in samples) added to the detected time (Logically, it changes the reconstructed energies )
* - Sum0Filt: 0-padding: Subtract the sum of the filter (1) or not (0)
* - tstartPulse1: Integer number: Sample where the first pulse starts or nameFile: File where the tstart (in seconds) of every pulse is
* - tstartPulse2: Tstart (samples) of the second pulse
* - tstartPulse3: Tstart (samples) of the third pulse (if 0 => PAIRS, if not 0 => TRIOS)
* - energyPCA1: First energy (only for PCA)
* - energyPCA2: Second energy (only for PCA)
* - XMLFile: XML input FITS file with instrument definition
* 
* Steps:
* 
* - Register HEATOOL
* - Reading all programm parameters by using PIL
* - Read XML info
* - Sixt standard keywords structure
* - Open output FITS file
* - Initialize PP data structures needed for pulse filtering
* - Initialize SIRENA data structures needed for pulse filtering
* - Read the grading data from the XML file and store it in 'reconstruct_init_sirena->grading'
* - Obtain the 'trig_reclength' and the sampling rate:
*   - If Rcmethod starts with '@' => List of record input FITS files. For every FITS file:
*       - Open FITS file
*       - Check if input FITS file have been simulated with TESSIM or XIFUSIM
*       - Get the sampling rate from the HISTORY keyword from the input FITS file
*       - If it is a xifusim simulated file
*           - Obtain 'trig_reclength' from the HISTORY block
*   - If Rcemethod doesn't start with '@' => Single record input FITS file
*       - Open FITS file
*       - Check if input FITS file have been simulated with TESSIM or XIFUSIM
*       - Get the sampling rate from the HISTORY keyword from the input FITS file
*       - If it is a xifusim simulated file
*           - Obtain 'trig_reclength' from the HISTORY block
* - Build up TesEventList to recover the results of the reconstruction
* - Reconstruct the input record FITS file:
*   - If Rcmethod starts with '@' => List of record input FITS files. For every FITS file:
*       - Open record file
*       - Initialize: initializeReconstruction or initializeReconstructionSIRENA
*       - Build up TesRecord to read the file
*       - Iterate of records and do the reconstruction
*           - Reconstruct: reconstructRecord or reconstructRecordSIRENA
*           - Save events to the event_list
*           - Copy trigger keywords to event file
*           - Close file
*   - If Rcemethod doesn't start with '@' => Single record input FITS file
*       - Open record file
*       - Initialize: initializeReconstruction or initializeReconstructionSIRENA
*       - Build up TesRecord to read the file
*       - Iterate of records and do the reconstruction
*           - Reconstruct: reconstructRecord or reconstructRecordSIRENA
*           - Save events to the event_list
*           - Copy trigger keywords to event file
*           - Close file
* - Save GTI extension to event file
* - Free memory
*****************************************************/
int tesreconstruction_main() {
  time_t ttstart = time(0);
  
  // Containing all programm parameters read by PIL.
  struct Parameters par;
  par.hduPRCLOFWM = 0;  // Debugger complains about an initialized variable (only the boolean type)
  par.hduPRECALWN = 0;  // Debugger complains about an initialized variable (only the boolean type)
  par.preBuffer = 0;    // Debugger complains about an initialized variable (only the boolean type)
  par.OFLib = 1;        // Debugger complains about an initialized variable (only the boolean type)
  
  // Error status.
  int status=EXIT_SUCCESS;

  // Register HEATOOL:
  set_toolname("tesreconstruction");
  set_toolversion("0.05");
  
  do { // Beginning of the ERROR handling loop (will at
       // most be run once).
    headas_chat(3, "initialize ...\n");

    // Get program parameters.
    status=getpar(&par);
    CHECK_STATUS_BREAK(status);
    
    if ((strcmp(par.Rcmethod,"SIRENA") == 0) && (strcmp(par.EnergyMethod,"I2RFITTED") == 0) && (par.Ifit == 0.0))
    {
        SIXT_ERROR("Ifit value must be provided");
        return(EXIT_FAILURE);
    }

    if (strcmp(par.Rcmethod,"SIRENA") && (par.OFLengthNotPadded < par.OFLength))
        printf("%s","Running 0-padding: 0 < OFLengthNotPadded < OFLength\n");
        
    double sf = -999.;
    double sampling_rate = -999.0;
    AdvDet *det = newAdvDet(&status);
    CHECK_STATUS_BREAK(status);
    
    if ((par.opmode ==1) || ((par.opmode == 0) && (par.preBuffer == 1)))
    {
        // Read XML info
        //--------------
        det = loadAdvDet(par.XMLFile, &status);
        CHECK_STATUS_BREAK(status);
        
        sf = det->SampleFreq;
    }
    
    if ((par.preBuffer == 1) && (par.opmode == 0))
    {
        printf("%s","Attention: preBuffer used => Parameters of library filters read from XML file\n");
    }
    
    int trig_reclength = -999;
    
    char firstchar = par.RecordFile[0];
    char firstchar2[2] = {firstchar , '\0'};

    // Check if input file header is complete to work with xifusim/tessim simulated files
    // -------------------------------------------------------------------------------
    fitsfile* fptr = NULL;
    //fitsfile* libptr = NULL;
    int numfits;
    int hdunum; // Number of HDUs (RECORDS-file or TESRECORDS-file)
    int tessimOrxifusim = -999;     // 0: tessim, 1: xifusim
    int numberkeywords;
    char *headerPrimary = NULL;
    char *sample_rate_pointer = NULL;
    char *trig_reclength_pointer = NULL;
    char each_character_after_srate[125];
    char each_character_after_treclength[125];
    char characters_after_srate[125];
    char characters_after_treclength[125];
    char extname[20];
    int extver = 0;
    
    if (strcmp(firstchar2,"@") == 0)
    {
            FILE *filetxt = fopen(strndup(par.RecordFile+1, strlen(par.RecordFile)-1), "r");
            if (filetxt == NULL)    
            {
                    printf("%s","File given in RecordFile does not exist\n");
                    status = 104;
            }
            CHECK_STATUS_BREAK(status);
            
            char filefits[256];
            numfits = 0;
            while(fscanf(filetxt,"%s",filefits)!=EOF)
            {
                numfits++;
            }
            fclose(filetxt);
            
            filetxt = fopen(strndup(par.RecordFile+1, strlen(par.RecordFile)-1), "r");
        
            for (int j=0;j<numfits;j++)   // For every FITS file
            {
                    fgets(filefits, 256, filetxt);
                    strtok(filefits, "\n");     // To delete '/n' from filefits (if not, 'fits_open_file' can not open the file)
            
                    fits_open_file(&fptr, filefits, READONLY, &status);
                    if (status != 0)    printf("%s","FITS file read from ASCII file does not exist\n");
                    CHECK_STATUS_BREAK(status);  
                    
                    fits_get_num_hdus(fptr, &hdunum,&status);
                    
                    // Check if input FITS file have been simulated with TESSIM or XIFUSIM
                    strcpy(extname,"RECORDS");
                    fits_movnam_hdu(fptr, ANY_HDU,extname, extver, &status);
                    if (status != 0)
                    {
                        status = 0;
                        strcpy(extname,"TESRECORDS");
                        fits_movnam_hdu(fptr, ANY_HDU,extname, extver, &status);
                        if (status != 0)
                        {
                            SIXT_ERROR("Cannot move to TESRECORDS HDU in input FITS file");
                            return(EXIT_FAILURE);
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
                        SIXT_ERROR("Neither the 'RECORDS' nor 'TESRECORDS' HDUs are in the input FITS file");
                        return(EXIT_FAILURE);
                    }
                    
                    // Get the sampling rate from the HISTORY keyword from the input FITS file
                    // Move to "Primary" HDU
                    fits_movabs_hdu(fptr, 1, NULL, &status); 
                    CHECK_STATUS_BREAK(status);
                    // and read full Primary HDU and store it in 'headerPrimary'
                    fits_hdr2str(fptr, 0, NULL, 0,&headerPrimary, &numberkeywords, &status); 
                    CHECK_STATUS_BREAK(status);
                    
                    // Pointer to where the text "sample_rate=" is in HISTORY block          
                    sample_rate_pointer = strstr (headerPrimary,"sample_rate=");    
                    if(!sample_rate_pointer)
                    {
                        // read it from xml file
                        sampling_rate = sf;
                    }
                    else
                    {
                        // Pointer to the next character to "sample_rate=" (12 characters)   
                        sample_rate_pointer = sample_rate_pointer + 12; 
                        snprintf(each_character_after_srate,125,"%c",*sample_rate_pointer);
                        snprintf(characters_after_srate,125,"%c",*sample_rate_pointer);
                        
                        while (*sample_rate_pointer != ' ')
                        {
                            sample_rate_pointer = sample_rate_pointer + 1;
                            snprintf(each_character_after_srate,125,"%c",*sample_rate_pointer);
                            strcat(characters_after_srate,each_character_after_srate); 
                        }
                        
                        sampling_rate = atof(characters_after_srate);
                        
                        if (((sampling_rate != -999.0) && (sf != -999.0)) && (sampling_rate != sf))
                        {
                            SIXT_ERROR("Sampling rate from input FITS file and from XML file do not match");
                            return(EXIT_FAILURE);
                        }
                    }
                    
                    if (tessimOrxifusim == 1) //xifusim simulated file (with TESRECORDS)
                    {    
                        /*// Move to "Primary" HDU 
                        fits_movabs_hdu(fptr, 1, NULL, &status); 
                        CHECK_STATUS_BREAK(status);
                        // and read full Primary HDU and store it in 'headerPrimary'
                        fits_hdr2str(fptr, 0, NULL, 0,&headerPrimary, &numberkeywords, &status); 
                        CHECK_STATUS_BREAK(status);*/
                    
                        // Pointer to where the text "trig_reclength=" is in HISTORY block
                        trig_reclength = -999.0;
                        trig_reclength_pointer = strstr (headerPrimary,"trig_reclength=");    
                        if(!trig_reclength_pointer)
                        {
                                //SIXT_ERROR("'trig_reclength' is not in the input FITS file (necessary if SIRENA is going tu run in THREADING mode)");
                                //return(EXIT_FAILURE);
                                printf("%s","Attention: 'trig_reclength' is not in the input FITS file (necessary if SIRENA is going to run in THREADING mode)\n");
                        }
                        else
                        {
                                // Pointer to the next character to "trig_reclength=" (15 characters)   
                                trig_reclength_pointer = trig_reclength_pointer + 15; 
                                snprintf(each_character_after_treclength,125,"%c",*trig_reclength_pointer);
                                
                                snprintf(characters_after_treclength,125,"%c",*trig_reclength_pointer);
                                
                                while (*trig_reclength_pointer != ' ')
                                {
                                    trig_reclength_pointer = trig_reclength_pointer + 1;
                                    snprintf(each_character_after_treclength,125,"%c",*trig_reclength_pointer);
                                    strcat(characters_after_treclength,each_character_after_treclength); 
                                }
                                
                                trig_reclength = atoi(characters_after_treclength);
                        }
                    }
                    
                    fits_close_file(fptr,&status);
                    CHECK_STATUS_BREAK(status);
            }
            
            fclose(filetxt);
    }
    else
    {
            numfits = 1;
            fits_open_file(&fptr, par.RecordFile, READONLY, &status);
            if (status != 0)    printf("%s","File given in RecordFile does not exist\n");
            
            fits_get_num_hdus(fptr, &hdunum,&status);
            
            // Check if input FITS file have been simulated with TESSIM or XIFUSIM
            strcpy(extname,"RECORDS");
            fits_movnam_hdu(fptr, ANY_HDU,extname, extver, &status);
            if (status != 0)
            {
                status = 0;
                strcpy(extname,"TESRECORDS");
                fits_movnam_hdu(fptr, ANY_HDU,extname, extver, &status);
                if (status != 0)
                {
                    SIXT_ERROR("Cannot move to TESRECORDS HDU in input FITS file");
                    return(EXIT_FAILURE);
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
                SIXT_ERROR("Neither the 'RECORDS' nor 'TESRECORDS' HDUs are in the input FITS file");
                return(EXIT_FAILURE);
            }

            if (par.opmode == 1)
            {
                status = checkXmls(&par);
            }

            // Get the sampling rate from the HISTORY keyword from the input FITS file
            // Move to "Primary" HDU
            fits_movabs_hdu(fptr, 1, NULL, &status); 
            CHECK_STATUS_BREAK(status);
            // and read full Primary HDU and store it in 'headerPrimary'
            if (fits_hdr2str(fptr, 0, NULL, 0,&headerPrimary, &numberkeywords, &status))
            {
                free(headerPrimary);
                CHECK_STATUS_BREAK(status);
            }
            
            // Pointer to where the text "sample_rate=" is in HISTORY block          
            sample_rate_pointer = strstr (headerPrimary,"sample_rate=");    
            if(!sample_rate_pointer)
            {
                // read it from xml file
                sampling_rate = sf;
            }
            else
            {
                // Pointer to the next character to "sample_rate=" (12 characters)
                sample_rate_pointer = sample_rate_pointer + 12; 
                snprintf(each_character_after_srate,125,"%c",*sample_rate_pointer);
                snprintf(characters_after_srate,125,"%c",*sample_rate_pointer);
                
                while (*sample_rate_pointer != ' ')
                {
                    sample_rate_pointer = sample_rate_pointer + 1;
                    snprintf(each_character_after_srate,125,"%c",*sample_rate_pointer);
                    strcat(characters_after_srate,each_character_after_srate); 
                }
                
                sampling_rate = atof(characters_after_srate);
                
                if (((sampling_rate != -999.0) && (sf != -999.0)) && (sampling_rate != sf))
                {
                    SIXT_ERROR("Sampling rate from input FITS file and from XML file do not match");
                    return(EXIT_FAILURE);
                }
            }
            
            if (tessimOrxifusim == 1) //xifusim simulated file (with TESRECORDS)
            {    
               /* // Move to "Primary" HDU
                fits_movabs_hdu(fptr, 1, NULL, &status); 
                CHECK_STATUS_BREAK(status);
                // and read full Primary HDU and store it in 'headerPrimary'
                fits_hdr2str(fptr, 0, NULL, 0,&headerPrimary, &numberkeywords, &status); 
                CHECK_STATUS_BREAK(status);*/
            
                // Pointer to where the text "trig_reclength=" is in HISTORY block
                trig_reclength = -999.0;
                trig_reclength_pointer = strstr (headerPrimary,"trig_reclength=");
                if(!trig_reclength_pointer)
                {
                    //SIXT_ERROR("'trig_reclength' is not in the input FITS file (necessary if SIRENA is going tu run in THREADING mode)");
                    //return(EXIT_FAILURE);
                    printf("%s","Attention: 'trig_reclength' is not in the input FITS file (necessary if SIRENA is going to run in THREADING mode)\n");
                }
                else
                {
                    // Pointer to the next character to "trig_reclength=" (15 characters)   
                    trig_reclength_pointer = trig_reclength_pointer + 15; 
                    snprintf(each_character_after_treclength,125,"%c",*trig_reclength_pointer);
                    
                    snprintf(characters_after_treclength,125,"%c",*trig_reclength_pointer);
                    
                    while (*trig_reclength_pointer != ' ')
                    {
                        trig_reclength_pointer = trig_reclength_pointer + 1;
                        snprintf(each_character_after_treclength,125,"%c",*trig_reclength_pointer);
                        strcat(characters_after_treclength,each_character_after_treclength); 
                    }
                    
                    trig_reclength = atoi(characters_after_treclength);
                }
            }
            fits_close_file(fptr,&status);
            CHECK_STATUS_BREAK(status);
    }

    if (sample_rate_pointer != NULL)
    {
        sample_rate_pointer = NULL;
        free(sample_rate_pointer);
    }
    if (trig_reclength_pointer != NULL)
    {
        trig_reclength_pointer = NULL;
        free(trig_reclength_pointer);
    }
    if (headerPrimary != NULL)
    {
        free(headerPrimary);
    }
    
    // Sixt standard keywords structure
    //----------------------------------
    SixtStdKeywords* keywords = newSixtStdKeywords(&status);
    CHECK_STATUS_BREAK(status);
    
    //Open outfile 
    //------------
    TesEventFile * outfile = opennewTesEventFileSIRENA(par.TesEventFile,
                                                 keywords,
                                                 SIRENA_VERSION,
                                                 par.clobber,
                                                 &status);
    CHECK_STATUS_BREAK(status);
    
    // Initialize PP data structures needed for pulse filtering
    //---------------------------------------------------------
    ReconstructInit* reconstruct_init = newReconstructInit(&status);
    CHECK_STATUS_BREAK(status);
    
    // Initialize SIRENA data structures needed for pulse filtering
    //-------------------------------------------------------------
    ReconstructInitSIRENA* reconstruct_init_sirena = newReconstructInitSIRENA(&status);
    CHECK_STATUS_BREAK(status);
    PulsesCollection* pulsesAll = newPulsesCollection(&status);
    CHECK_STATUS_BREAK(status);  
    OptimalFilterSIRENA* optimalFilter = newOptimalFilterSIRENA(&status);
    CHECK_STATUS_BREAK(status);// define a second structure for calibration
    
    if ((par.opmode ==1) || ((par.opmode == 0) && (par.preBuffer == 1)))
    {
        // Read the grading data from the XML file and store it in 'reconstruct_init_sirena->grading'
        reconstruct_init_sirena->grading = NULL;
        reconstruct_init_sirena->grading = (Grading*)malloc(sizeof(Grading));
        
        reconstruct_init_sirena->grading->ngrades = 0;
        reconstruct_init_sirena->grading->value  = NULL;
        reconstruct_init_sirena->grading->gradeData = NULL;
        
        if ((det->nrecons == 0) && (det->npix == 0))
        {
            SIXT_ERROR("The provided XMLFile does not have the grading info");
            return(EXIT_FAILURE);
        }
        else if ((det->nrecons == 0) && (det->npix != 0))
        {
            if (det->pix->grades == NULL)
            {
                SIXT_ERROR("The provided XMLFile does not have the grading info");
                return(EXIT_FAILURE);
            }
            reconstruct_init_sirena->grading->ngrades=det->pix->ngrades;
            reconstruct_init_sirena->grading->gradeData = gsl_matrix_alloc(det->pix->ngrades,3);
            for (int i=0;i<det->pix->ngrades;i++)
            {
                if ((int) (det->pix->grades[i].gradelim_post) != 8)
                {
                    if ((int) (det->pix->grades[i].grade_preBuffer) >= (int) (det->pix->grades[i].gradelim_post))
                    {
                        printf("%s %d %s %d %s","preBuffer=",(int) (det->pix->grades[i].grade_preBuffer)," for filter length ",(int) (det->pix->grades[i].gradelim_post),"\n");
                        SIXT_ERROR("preBuffer values provided in the XML file should be lower than corresponding filter lengths");
                        return(EXIT_FAILURE);
                    }
                }
                gsl_matrix_set(reconstruct_init_sirena->grading->gradeData,i,0,(int) (det->pix->grades[i].gradelim_pre));
                gsl_matrix_set(reconstruct_init_sirena->grading->gradeData,i,1,(int) (det->pix->grades[i].gradelim_post));
                gsl_matrix_set(reconstruct_init_sirena->grading->gradeData,i,2,(int) (det->pix->grades[i].grade_preBuffer));
                if (((int) (det->pix->grades[i].gradelim_post) == 8) && ((int) (det->pix->grades[i].grade_preBuffer) != 0))
                {
                    printf("%s","Attention: preBuffer!=0 for low resolution (filter length 8) but preBuffer=0 is going to be used\n");
                    gsl_matrix_set(reconstruct_init_sirena->grading->gradeData,i,2,0);
                }
                /*printf("%s %f %s","0pre=",gsl_matrix_get(reconstruct_init_sirena->grading->gradeData,i,0),"\n");
                printf("%s %f %s","post=",gsl_matrix_get(reconstruct_init_sirena->grading->gradeData,i,1),"\n");
                printf("%s %f %s","pB=",gsl_matrix_get(reconstruct_init_sirena->grading->gradeData,i,2),"\n");*/
            }
        }
        else if(((det->nrecons != 0) && (det->npix == 0)) || ((det->nrecons == 1) && (det->npix == 1)))
        {
            if (det->recons->grades == NULL)
            {
                SIXT_ERROR("The provided XMLFile does not have the grading info");
                return(EXIT_FAILURE);
            }
            reconstruct_init_sirena->grading->ngrades=det->recons->ngrades;
            reconstruct_init_sirena->grading->gradeData = gsl_matrix_alloc(det->recons->ngrades,3);
            for (int i=0;i<det->recons->ngrades;i++)
            {
                if ((int) (det->recons->grades[i].gradelim_post) != 8)
                {
                    if ((int) (det->recons->grades[i].grade_preBuffer) >= (int) (det->recons->grades[i].gradelim_post))
                    {
                        printf("%s %d %s %d %s","preBuffer=",(int) (det->recons->grades[i].grade_preBuffer)," for filter length ",(int) (det->recons->grades[i].gradelim_post),"\n");
                        SIXT_ERROR("preBuffer values provided in the XML file should be lower than corresponding filter lengths");
                        return(EXIT_FAILURE);
                    }
                }
                gsl_matrix_set(reconstruct_init_sirena->grading->gradeData,i,0,(int) (det->recons->grades[i].gradelim_pre));
                gsl_matrix_set(reconstruct_init_sirena->grading->gradeData,i,1,(int) (det->recons->grades[i].gradelim_post));
                gsl_matrix_set(reconstruct_init_sirena->grading->gradeData,i,2,(int) (det->recons->grades[i].grade_preBuffer));
                if (((int) (det->recons->grades[i].gradelim_post) == 8) && ((int) (det->recons->grades[i].grade_preBuffer) != 0))
                {
                    printf("%s","Attention: preBuffer!=0 for low resolution (filter length 8) but preBuffer=0 is going to be used\n");
                    gsl_matrix_set(reconstruct_init_sirena->grading->gradeData,i,2,0);
                }
                /*printf("%s %f %s","1pre=",gsl_matrix_get(reconstruct_init_sirena->grading->gradeData,i,0),"\n");
                printf("%s %f %s","post=",gsl_matrix_get(reconstruct_init_sirena->grading->gradeData,i,1),"\n");
                printf("%s %f %s","pB=",gsl_matrix_get(reconstruct_init_sirena->grading->gradeData,i,2),"\n");*/
            }
        }
        
        int OFlengthvsposti = 0;
        if ((par.preBuffer == 1) && (par.opmode == 1))
        {
            for (int i=0;i<reconstruct_init_sirena->grading->ngrades;i++) 
            {
                if (par.OFLength == gsl_matrix_get(reconstruct_init_sirena->grading->gradeData,i,1))
                {
                    OFlengthvsposti = 1;
                    break;
                }
            }
            if (OFlengthvsposti == 0)
            {
                SIXT_ERROR("The grading/preBuffer info of the XML file does not match the OFLength input parameter");
                return(EXIT_FAILURE);
            }
        }
    }
    destroyAdvDet(&det);

    // Build up TesEventList to recover the results of the reconstruction
    // SIRENA method
    TesEventList* event_list = newTesEventListSIRENA(&status);
    allocateTesEventListTrigger(event_list,par.EventListSize,&status);
    //allocateWholeTesEventListTrigger(event_list,par.EventListSize,&status);
    CHECK_STATUS_BREAK(status);
    // PP method
    TesEventList* event_listPP = newTesEventList(&status);
    allocateTesEventListTrigger(event_listPP,par.EventListSize,&status);
    CHECK_STATUS_BREAK(status);
            
    TesTriggerFile* record_file;
    TesRecord* record;
    int lastRecord = 0, nrecord = 0, startRecordGroup = 0,nrecord_filei = 0;    //last record required for SIRENA library creation
    
    if (strcmp(firstchar2,"@") == 0)
    {
            FILE *filetxt = fopen(strndup(par.RecordFile+1, strlen(par.RecordFile)-1), "r");
            if (status != 0)    printf("%s","FITS file read from ASCII file does not exist\n");
            CHECK_STATUS_BREAK(status);  
            
            char filefits[256];
        
            for (int j=0;j<numfits;j++)   // For every FITS file
            {
                    fgets(filefits, 256, filetxt);
                    strtok(filefits, "\n");     // To delete '/n' from filefits (if not, 'fits_open_file' can not open the file)
                    //printf("%s %s %s","FITS file i: ",filefits,"\n");
            
                    // Open record file
                    // ----------------
                    record_file = openexistingTesTriggerFile(filefits,keywords,&status);
                    CHECK_STATUS_BREAK(status);
                    
                    if(!strcmp(par.Rcmethod,"PP"))
                    {
                            initializeReconstruction(reconstruct_init,par.OptimalFilterFile,par.OFLengthNotPadded,
                                                    par.PulseTemplateFile,par.Threshold,par.Calfac,par.NormalExclusion,
                                                    par.DerivateExclusion,par.SaturationValue,&status);
                    }
                    else
                    {
                            initializeReconstructionSIRENA(reconstruct_init_sirena, par.RecordFile, record_file->fptr, 
                                                    par.LibraryFile, par.TesEventFile, par.OFLengthNotPadded, par.scaleFactor, par.samplesUp,
                                                    par.samplesDown, par.nSgms, par.detectSP, par.opmode, par.detectionMode, par.LrsT, 
                                                    par.LbT, par.NoiseFile, par.FilterDomain, par.FilterMethod, par.EnergyMethod, 
                                                    par.filtEev, par.Ifit, par.OFNoise, par.LagsOrNot, par.nLags, par.Fitting35, par.OFIter, 
                                                    par.OFLib, par.OFInterp, par.OFStrategy, par.OFLength, par.preBuffer,par.monoenergy, 
                                                    par.hduPRECALWN, par.hduPRCLOFWM, par.largeFilter, par.intermediate, par.detectFile, 
                                                    par.errorT, par.Sum0Filt, par.clobber, par.EventListSize, par.SaturationValue, par.tstartPulse1, 
                                                    par.tstartPulse2, par.tstartPulse3, par.energyPCA1, par.energyPCA2, par.XMLFile, &status);
                                            
                    }  
                    CHECK_STATUS_BREAK(status);
                    
                    // Build up TesRecord to read the file
                    record = newTesRecord(&status);
                    if (record_file->delta_t == -999) record_file->delta_t = 1./sampling_rate;
                    if (trig_reclength == -999) trig_reclength = record_file->trigger_size;
                    allocateTesRecord(record,trig_reclength,record_file->delta_t,0,&status);
                    //allocateTesRecord(record,record_file->trigger_size,record_file->delta_t,0,&status);
                    CHECK_STATUS_BREAK(status);
                    
                    // Iterate of records and do the reconstruction
                    //int lastRecord = 0, nrecord = 0;    //last record required for SIRENA library creation
                    nrecord_filei = 0;
                    while(getNextRecord(record_file,record,&lastRecord,&startRecordGroup,&status))
                    {
                            if(!strcmp(par.Rcmethod,"PP"))
                            {
                                    reconstructRecord(record,event_listPP,reconstruct_init,0,&status);
                            }
                            else
                            {
                                    nrecord = nrecord + 1;
                                    nrecord_filei = nrecord_filei + 1;
                                    if ((nrecord_filei == record_file->nrows) && (j == numfits-1)) lastRecord=1;  // lastRecord of all the FITS files
                                   
                                    if ((strcmp(par.EnergyMethod,"I2R") == 0) || (strcmp(par.EnergyMethod,"I2RFITTED") == 0)
                                        || (strcmp(par.EnergyMethod,"I2RDER") == 0))
                                    {
                                            strcpy(reconstruct_init_sirena->EnergyMethod,par.EnergyMethod);
                                    }
                                
                                    //printf("%s %d %s","**TESRECONSTRUCTION nrecord = ",nrecord,"\n");
                                    reconstructRecordSIRENA(record,trig_reclength,event_list,reconstruct_init_sirena,
                                                            lastRecord, startRecordGroup, &pulsesAll, &optimalFilter, &status);
                            }
                            CHECK_STATUS_BREAK(status);

                            if ((strcmp(par.EnergyMethod,"PCA") != 0) || ((strcmp(par.EnergyMethod,"PCA") == 0) && lastRecord == 1))
                            {
                                    // In THREADING mode, saveEventListToFileSIRENA is not called until finishing with calculus 
                                    // (ordering is necessary previously)  
                                    if(!is_threading()){    
                                            saveEventListToFileSIRENA(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
                                            CHECK_STATUS_BREAK(status);
                                            //Reinitialize event list
                                            event_list->index=0;//////////////////////!!!!!!!!!!!!!!!!!!!!!OJO!!!!!!!!!!!!!!!!!!! Igual no hay que inicializarlo a 0
                                    }
                                    //else
                                    //        printf("%s","Not prepared to run in THREADING mode with a input ASCII file (with several FITS files)\n");
                            }
                    } // while getNextRecord
                    if(is_threading()) 
                    {
                            th_end(&reconstruct_init_sirena, &pulsesAll, &optimalFilter);
                            int i = 1;
                            int aux = 1;
                            while((aux = th_get_event_list(&event_list, &record)) == 1)
                            {
                                    saveEventListToFileSIRENA(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
                                    CHECK_STATUS_BREAK(status);
                                    ++i;
                            }
                    }
                    
                    if ((!strcmp(par.Rcmethod,"SIRENA")) && (pulsesAll->ndetpulses == 0)) 
                            printf("%s %s %s","WARNING: no pulses have been detected in the current FITS file: ", filefits,"\n");
                    
                    if (numfits == 0)   // Only one time
                    {
                            // Copy trigger keywords to event file
                            copyTriggerKeywords(record_file->fptr,outfile->fptr,&status);
                            CHECK_STATUS_BREAK(status);
                            
                            // Messages providing info of some columns
                            char keywordvalue[9];
                            char comment[MAXMSG];
                            
                            fits_movnam_hdu(outfile->fptr, ANY_HDU,"EVENTS", 0, &status);
                            CHECK_STATUS_BREAK(status);
                            
                            fits_read_key(outfile->fptr, TSTRING, "TTYPE1", &keywordvalue, NULL, &status);
                            strcpy(comment, "Starting time");
                            fits_update_key(outfile->fptr, TSTRING, "TTYPE1", keywordvalue, comment, &status);
                            
                            fits_read_key(outfile->fptr, TSTRING, "TTYPE2", &keywordvalue, NULL, &status);
                            strcpy(comment, "Reconstructed-uncalibrated energy");
                            fits_update_key(outfile->fptr, TSTRING, "TTYPE2", keywordvalue, comment, &status);      
                            
                            fits_read_key(outfile->fptr, TSTRING, "TTYPE3", &keywordvalue, NULL, &status);
                            strcpy(comment, "Average first 4 samples (derivative)");
                            fits_update_key(outfile->fptr, TSTRING, "TTYPE3", keywordvalue, comment, &status);      
                            
                            fits_read_key(outfile->fptr, TSTRING, "TTYPE4", &keywordvalue, NULL, &status);
                            strcpy(comment, "Optimal filter length");
                            fits_update_key(outfile->fptr, TSTRING, "TTYPE4", keywordvalue, comment, &status);      
                            
                            fits_read_key(outfile->fptr, TSTRING, "TTYPE5", &keywordvalue, NULL, &status);
                            strcpy(comment, "Starting time-starting time previous event");
                            fits_update_key(outfile->fptr, TSTRING, "TTYPE5", keywordvalue, comment, &status);
                    }
                    
                    freeTesTriggerFile(&record_file,&status);   // The record_file (every FITS file) is closed
                    
                    CHECK_STATUS_BREAK(status);
            
            }   // for every FITS file
            
            fclose(filetxt);
    }
    else
    {
            // Open record file
            // ----------------
            record_file = openexistingTesTriggerFile(par.RecordFile,keywords,&status);
            CHECK_STATUS_BREAK(status);

            if(!strcmp(par.Rcmethod,"PP")){
                initializeReconstruction(reconstruct_init,par.OptimalFilterFile,par.OFLengthNotPadded,
                        par.PulseTemplateFile,par.Threshold,par.Calfac,par.NormalExclusion,
                        par.DerivateExclusion,par.SaturationValue,&status);
            }else{
                initializeReconstructionSIRENA(reconstruct_init_sirena, par.RecordFile, record_file->fptr, 
                        par.LibraryFile, par.TesEventFile, par.OFLengthNotPadded, par.scaleFactor, par.samplesUp,
                        par.samplesDown, par.nSgms, par.detectSP, par.opmode, par.detectionMode, par.LrsT, 
                        par.LbT, par.NoiseFile, par.FilterDomain, par.FilterMethod, par.EnergyMethod, 
                        par.filtEev, par.Ifit, par.OFNoise, par.LagsOrNot, par.nLags, par.Fitting35, par.OFIter, 
                        par.OFLib, par.OFInterp, par.OFStrategy, par.OFLength, par.preBuffer, par.monoenergy, 
                        par.hduPRECALWN, par.hduPRCLOFWM, par.largeFilter, par.intermediate, par.detectFile, 
                        par.errorT, par.Sum0Filt, par.clobber, par.EventListSize, par.SaturationValue, par.tstartPulse1, 
                        par.tstartPulse2, par.tstartPulse3, par.energyPCA1, par.energyPCA2, par.XMLFile, &status);
            }
            CHECK_STATUS_BREAK(status);
            
            // Build up TesRecord to read the file
            record = newTesRecord(&status);
            if ((record_file->delta_t == -999) && (sampling_rate == -999))
            {
                SIXT_ERROR("Cannot read or get the sampling rate neither from the input FITS file nor the XML file. Please, include the DELTAT keyword (inverse of sampling rate) in the input FITS file before running TESRECONSTRUCTION again");
                return(EXIT_FAILURE);
            }
            if (record_file->delta_t == -999) record_file->delta_t = 1./sampling_rate;
            if (trig_reclength == -999) trig_reclength = record_file->trigger_size;
            //allocateTesRecord(record,record_file->trigger_size,record_file->delta_t,0,&status);
            allocateTesRecord(record,trig_reclength,record_file->delta_t,0,&status);
            CHECK_STATUS_BREAK(status);
                
            // Iterate of records and do the reconstruction
            lastRecord = 0, nrecord = 0;    //last record required for SIRENA library creation
            while(getNextRecord(record_file,record,&lastRecord,&startRecordGroup,&status))
            {
                    if(!strcmp(par.Rcmethod,"PP"))
                    {
                            reconstructRecord(record,event_listPP,reconstruct_init,0,&status);
                    }
                    else
                    {
                            nrecord = nrecord + 1;
                            //if(nrecord == record_file->nrows) lastRecord=1;
                            /*if(nrecord < 5887) 
                            {
                              continue;
                            }
                            else if(nrecord > 5920)
                            {
                              status=1;
                              CHECK_STATUS_BREAK(status);
                            }*/
                            /*if(nrecord > 1)
                            {
                            	status=1;
                                CHECK_STATUS_BREAK(status);
                            }*/
                            if ((strcmp(par.EnergyMethod,"I2R") == 0) || (strcmp(par.EnergyMethod,"I2RFITTED") == 0)
                                || (strcmp(par.EnergyMethod,"I2RDER") == 0))
                            {
                                strcpy(reconstruct_init_sirena->EnergyMethod,par.EnergyMethod);
                            }
                        
                            //printf("%s %d %s","**TESRECONSTRUCTION nrecord = ",nrecord,"\n");
                            reconstructRecordSIRENA(record,trig_reclength, event_list,reconstruct_init_sirena,
                                                    lastRecord, startRecordGroup, &pulsesAll, &optimalFilter, &status);
                    }
                    CHECK_STATUS_BREAK(status);

                    if ((strcmp(par.EnergyMethod,"PCA") != 0) || ((strcmp(par.EnergyMethod,"PCA") == 0) && lastRecord == 1))
                    {
                            // In THREADING mode, saveEventListToFileSIRENA is not called until finishing with calculus 
                            // (ordering is necessary previously)  
                            if(!is_threading()){    
                                    saveEventListToFileSIRENA(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
                                    CHECK_STATUS_BREAK(status);
                                    //Reinitialize event list
                                    event_list->index=0;
                            }

                            IntegrafreeTesEventListSIRENA(event_list);
                            //if (NULL!=event_list->event_indexes) free(event_list->event_indexes);
                            //if (NULL!=event_list->grades1) free(event_list->grades1);
                            //if (NULL!=event_list->pulse_heights) free(event_list->pulse_heights);
                    }
            }
            
            if(is_threading()) 
            {
                    //printf("%s","**Threading...waiting \n");
                    th_end(&reconstruct_init_sirena, &pulsesAll, &optimalFilter);
                    int i = 1;
                    int aux = 1;
                    while((aux = th_get_event_list(&event_list, &record)) == 1)
                    {
                            saveEventListToFileSIRENA(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
                            CHECK_STATUS_BREAK(status);
                            ++i;
                    }
                    IntegrafreeTesEventListSIRENA(event_list);
                    //if (NULL!=event_list->event_indexes) free(event_list->event_indexes);
                    //if (NULL!=event_list->grades1) free(event_list->grades1);
                    //if (NULL!=event_list->pulse_heights) free(event_list->pulse_heights);
            }
            
            if ((!strcmp(par.Rcmethod,"SIRENA")) && (pulsesAll->ndetpulses == 0)) 
            printf("%s","WARNING: no pulses have been detected\n");
            
            // Copy trigger keywords to event file
            copyTriggerKeywords(record_file->fptr,outfile->fptr,&status);
            CHECK_STATUS_BREAK(status);
            
            // Messages providing info of some columns
            char keywordvalue[9];
            char comment[MAXMSG];
            
            fits_movnam_hdu(outfile->fptr, ANY_HDU,"EVENTS", 0, &status);
            CHECK_STATUS_BREAK(status);
            
            fits_read_key(outfile->fptr, TSTRING, "TTYPE1", &keywordvalue, NULL, &status);
            strcpy(comment, "Starting time");
            fits_update_key(outfile->fptr, TSTRING, "TTYPE1", keywordvalue, comment, &status);
            
            fits_read_key(outfile->fptr, TSTRING, "TTYPE2", &keywordvalue, NULL, &status);
            strcpy(comment, "Reconstructed-uncalibrated energy");
            fits_update_key(outfile->fptr, TSTRING, "TTYPE2", keywordvalue, comment, &status);      
            
            fits_read_key(outfile->fptr, TSTRING, "TTYPE3", &keywordvalue, NULL, &status);
            strcpy(comment, "Average first 4 samples (derivative)");
            fits_update_key(outfile->fptr, TSTRING, "TTYPE3", keywordvalue, comment, &status);      
            
            fits_read_key(outfile->fptr, TSTRING, "TTYPE4", &keywordvalue, NULL, &status);
            strcpy(comment, "Optimal filter length");
            fits_update_key(outfile->fptr, TSTRING, "TTYPE4", keywordvalue, comment, &status);      
            
            fits_read_key(outfile->fptr, TSTRING, "TTYPE5", &keywordvalue, NULL, &status);
            strcpy(comment, "Starting time-starting time previous event");
            fits_update_key(outfile->fptr, TSTRING, "TTYPE5", keywordvalue, comment, &status);
            
            freeTesTriggerFile(&record_file,&status);
    }
    
    // Save GTI extension to event file
    GTI* gti=getGTIFromFileOrContinuous("none",keywords->tstart, keywords->tstop,keywords->mjdref, &status);
    saveGTIExt(outfile->fptr, "STDGTI", gti, &status);    
    CHECK_STATUS_BREAK(status);
    
    //Free memory
    freeReconstructInit(reconstruct_init);
    if (reconstruct_init_sirena->grading != NULL)
    {
        gsl_vector_free(reconstruct_init_sirena->grading->value); reconstruct_init_sirena->grading->value = 0;
        gsl_matrix_free(reconstruct_init_sirena->grading->gradeData); reconstruct_init_sirena->grading->gradeData = 0;
    }
    free(reconstruct_init_sirena->grading);
    reconstruct_init_sirena->grading = 0;
    freeReconstructInitSIRENA(reconstruct_init_sirena);
    freePulsesCollection(pulsesAll);
    freeOptimalFilterSIRENA(optimalFilter);
    freeTesEventFile(outfile,&status);
    freeTesEventListSIRENA(event_list);
    freeTesEventList(event_listPP);
    freeTesRecord(&record);
    freeSixtStdKeywords(keywords);
    CHECK_STATUS_BREAK(status);
 
  } while(0); // END of the error handling loop.
  
  if (EXIT_SUCCESS==status) 
  {
	headas_chat(3, "finished successfully!\n\n");
        time_t ttcurrent = time(0);
        printf("Elapsed time: %f\n", ((float)(ttcurrent - ttstart)));
	return(EXIT_SUCCESS);
  } 
  else 
  {
        time_t ttcurrent = time(0);
        printf("Elapsed time: %f\n", ((float)(ttcurrent - ttstart)));
	return(status);
  }
}
/*xxxx end of SECTION 1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 2 ************************************************************
* getpar function: This function gets the input parameter from the command line or their default values from the tesreconstruction .par file
*
* Parameters:
* - par: Structure containing the input parameters
******************************************************************************/
int getpar(struct Parameters* const par)

{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  status=ape_trad_query_string("Rcmethod", &sbuffer);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the reconstruction method");
	  return(status);
  }
  strcpy(par->Rcmethod, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("RecordFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the optimal filter file");
    return(status);
  }
  strcpy(par->RecordFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("TesEventFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event file");
    return(status);
  }
  strcpy(par->TesEventFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("EventListSize", &par->EventListSize);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the EventListSize parameter");
	  return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the clobber parameter");
	  return(status);
  }

  status=ape_trad_query_bool("history", &par->history);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the history parameter");
    return(status);
  }

  status=ape_trad_query_double("SaturationValue", &par->SaturationValue);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the SaturationValue parameter");
    return(status);
  }

  if(strcmp(par->Rcmethod,"PP")==0){
	// PP  reconstruction method
	status=ape_trad_query_string("OptimalFilterFile", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the optimal filter file");
		return(status);
	}
	strcpy(par->OptimalFilterFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("PulseTemplateFile", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the pulse template file");
		return(status);
	}
	strcpy(par->PulseTemplateFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_double("Threshold", &par->Threshold);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the Threshold parameter");
		return(status);
	}

	status=ape_trad_query_double("Calfac", &par->Calfac);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the Calfac parameter");
		return(status);
	}

	status=ape_trad_query_int("NormalExclusion", &par->NormalExclusion);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the NormalExclusion parameter");
		return(status);
	}

	status=ape_trad_query_int("DerivateExclusion", &par->DerivateExclusion);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the DerivateExclusion parameter");
		return(status);
	}
  }else if(strcmp(par->Rcmethod,"SIRENA")==0){
	
      // SIRENA parameters
      status=ape_trad_query_string("LibraryFile", &sbuffer);
      strcpy(par->LibraryFile, sbuffer);
      free(sbuffer);
      
      status=ape_trad_query_double("scaleFactor", &par->scaleFactor);
      
      status=ape_trad_query_int("samplesUp", &par->samplesUp);
      
      status=ape_trad_query_int("samplesDown", &par->samplesDown);
      
      status=ape_trad_query_double("nSgms", &par->nSgms);

      status=ape_trad_query_string("detectionMode", &sbuffer);
      strcpy(par->detectionMode, sbuffer);
      free(sbuffer);
      
      status=ape_trad_query_int("detectSP", &par->detectSP);
      
      status=ape_trad_query_int("opmode", &par->opmode);
      
      status=ape_trad_query_string("detectionMode", &sbuffer);
      strcpy(par->detectionMode, sbuffer);
      free(sbuffer);
      
      status=ape_trad_query_double("LrsT", &par->LrsT);
      
      status=ape_trad_query_double("LbT", &par->LbT);
      
      status=ape_trad_query_int("intermediate", &par->intermediate);
      
      status=ape_trad_query_string("detectFile", &sbuffer);
      strcpy(par->detectFile, sbuffer);
      free(sbuffer);
      
      status=ape_trad_query_double("monoenergy", &par->monoenergy);
      
      status=ape_trad_query_bool("hduPRECALWN", &par->hduPRECALWN);
      status=ape_trad_query_bool("hduPRCLOFWM", &par->hduPRCLOFWM);
      
      status=ape_trad_query_int("largeFilter", &par->largeFilter);
      
      status=ape_trad_query_string("NoiseFile", &sbuffer);
      strcpy(par->NoiseFile, sbuffer);
      free(sbuffer);
      
      status=ape_trad_query_string("FilterDomain", &sbuffer);
      strcpy(par->FilterDomain, sbuffer);
      free(sbuffer);
      
      status=ape_trad_query_string("FilterMethod", &sbuffer);
      strcpy(par->FilterMethod, sbuffer);
      free(sbuffer);
      
      status=ape_trad_query_string("EnergyMethod", &sbuffer);
      strcpy(par->EnergyMethod, sbuffer);
      free(sbuffer);

      status=ape_trad_query_int("OFLengthNotPadded", &par->OFLengthNotPadded);
      if (EXIT_SUCCESS!=status) {
        SIXT_ERROR("failed reading the OFLengthNotPadded parameter");
        return(status);
      }
      assert(par->OFLengthNotPadded > 0);
      //MyAssert((par->OFLengthNotPadded > 0) && (par->OFLengthNotPadded <= par->OFLength), "0-padding: 0 < OFLengthNotPadded <= OFLength");

      MyAssert((par->opmode == 0) || (par->opmode == 1), "opmode must be 0 or 1");
      
      status=ape_trad_query_double("filtEev", &par->filtEev);
      
      status=ape_trad_query_double("Ifit", &par->Ifit);
      
      status=ape_trad_query_string("OFNoise", &sbuffer);
      strcpy(par->OFNoise, sbuffer);
      free(sbuffer);
      
      status=ape_trad_query_int("LagsOrNot", &par->LagsOrNot);
      status=ape_trad_query_int("nLags", &par->nLags);
      status=ape_trad_query_int("Fitting35", &par->Fitting35);
      
      status=ape_trad_query_int("OFIter", &par->OFIter);
      
      status=ape_trad_query_bool("OFLib", &par->OFLib);
      
      strcpy(par->OFInterp, "DAB");
      
      status=ape_trad_query_string("OFStrategy", &sbuffer);
      strcpy(par->OFStrategy, sbuffer);
      free(sbuffer);

      if ((strcmp(par->EnergyMethod,"OPTFILT") == 0) || (strcmp(par->EnergyMethod,"I2R") == 0) || (strcmp(par->EnergyMethod,"I2RFITTED") == 0) || (strcmp(par->EnergyMethod,"I2RDER") == 0))
      {
        if (strcmp(par->OFStrategy,"FREE") == 0)  par->OFLib = 0;
        else if (strcmp(par->OFStrategy,"FIXED") == 0)    par->OFLib = 1;
        else if (strcmp(par->OFStrategy,"BYGRADE") == 0)  par->OFLib = 1;
      }

      status=ape_trad_query_int("OFLength", &par->OFLength);
      
      status=ape_trad_query_bool("preBuffer", &par->preBuffer);
      
      status=ape_trad_query_int("errorT", &par->errorT);

      status=ape_trad_query_int("Sum0Filt", &par->Sum0Filt);
      
      //status=ape_trad_query_int("tstartPulse1", &par->tstartPulse1);
      status=ape_trad_query_string("tstartPulse1", &sbuffer);
      strcpy(par->tstartPulse1, sbuffer);
      free(sbuffer);
      
      status=ape_trad_query_int("tstartPulse2", &par->tstartPulse2);
      
      status=ape_trad_query_int("tstartPulse3", &par->tstartPulse3);
      
      status=ape_trad_query_double("energyPCA1", &par->energyPCA1);
      
      status=ape_trad_query_double("energyPCA2", &par->energyPCA2);
      
      status=ape_trad_query_string("XMLFile", &sbuffer);
      strcpy(par->XMLFile, sbuffer);
      free(sbuffer);
      
      if (EXIT_SUCCESS!=status) {
          SIXT_ERROR("failed reading some SIRENA parameter");
          return(status);
      }
      
      MyAssert((par->opmode == 0) || (par->opmode == 1), "opmode must be 0 or 1");
      int isNumber = 1;
      for (int i = 0; i < strlen(par->tstartPulse1); i++) 
      {
          if (isdigit(par->tstartPulse1[i]) == 0)    
          {
              isNumber = 0;
              break;
          }
      }
      if ((isNumber == 0) && (par->opmode == 0))
      {
          SIXT_ERROR("tstartPulse1 can not be a file if CALIBRATION mode");
          return(EXIT_FAILURE);
      }
      if ((isNumber == 0) && (strcmp(par->FilterDomain,"F") == 0))    // It is only implemented tstartPulse1 as a file for time domain
      {
          SIXT_ERROR("It is not possible to work in FREQUENCY domain if tstartPulse1 is a file => Change FilterDomain to TIME domain (T) ");
          return(EXIT_FAILURE);
      }
      
      MyAssert((par->intermediate == 0) || (par->intermediate == 1), "intermediate must be 0 or 1");
      
      if (par->opmode == 0) MyAssert(par->monoenergy > 0, "monoenergy must be greater than 0");
      
      MyAssert((strcmp(par->FilterDomain,"T") == 0) || (strcmp(par->FilterDomain,"F") == 0), "FilterDomain must be T or F");
      
      MyAssert((strcmp(par->FilterMethod,"F0") == 0) || (strcmp(par->FilterMethod,"B0") == 0),"FilterMethod must be F0 or B0");
      
      MyAssert((strcmp(par->EnergyMethod,"OPTFILT") == 0) || (strcmp(par->EnergyMethod,"WEIGHT") == 0) || (strcmp(par->EnergyMethod,"WEIGHTN") == 0) ||
      (strcmp(par->EnergyMethod,"I2R") == 0) ||	(strcmp(par->EnergyMethod,"I2RFITTED") == 0) ||	(strcmp(par->EnergyMethod,"I2RDER") == 0)
      || (strcmp(par->EnergyMethod,"PCA") == 0), "EnergyMethod must be OPTFILT, WEIGHT, WEIGHTN, I2R, I2RFITTED, I2RDER or PCA");
      
      MyAssert((strcmp(par->OFNoise,"NSD") == 0) || (strcmp(par->OFNoise,"WEIGHTM") == 0), "OFNoise must be NSD or WEIGHTM");
      
      MyAssert((strcmp(par->detectionMode,"AD") == 0) || (strcmp(par->detectionMode,"STC") == 0), "detectionMode must be AD or STC");
      
      MyAssert((par->LagsOrNot ==0) || (par->LagsOrNot ==1), "LagsOrNot must me 0 or 1");
      if ((par->nLags)%2 == 0)
      {
          SIXT_ERROR("parameter error: nLags must be odd");
          return(EXIT_FAILURE);
      }
      MyAssert((par->Fitting35 ==3) || (par->Fitting35 ==5), "Fitting35 must me 3 or 5");
      if ((par->Fitting35 ==3) && (par->nLags<3))
      {
          SIXT_ERROR("parameter error: nLags must be at least 3");
          return(EXIT_FAILURE);
      }
      if ((par->Fitting35 ==5) && (par->nLags<5))
      {
          SIXT_ERROR("parameter error: nLags must be at least 5");
          return(EXIT_FAILURE);
      }
      
      MyAssert((par->Sum0Filt ==0) || (par->Sum0Filt ==1), "Sum0Filt must be 0 or 1");
      
      if ((strcmp(par->EnergyMethod,"WEIGHT") == 0) && (par->LagsOrNot == 1))
      {
          SIXT_ERROR("parameter error: EnergyMethod=WEIGHT and Lags not implemented yet");
          return(EXIT_FAILURE);
      }
      
      MyAssert((par->OFIter ==0) || (par->OFIter ==1), "OFIter must be 0 or 1");
      
      // It was in order to not ask for the noise file if OFLib=1
      /*if ((par->OFLib == 1) && (strcmp(par->FilterMethod,"F0") != 0))
       {                                                       *
       SIXT_ERROR("parameter error: If OFLib=yes => FilterMethod must be F0");
       return(EXIT_FAILURE);
  }*/
      
      if ((par->OFLengthNotPadded < par->OFLength) && (strcmp(par->FilterDomain,"F") == 0))
      {
          SIXT_ERROR("Code is not prepared to run 0-padding in Frequency domain");
          // To run 0-padding in Frequency domain the steps should be:
          //1. Take the 8192-samples-length filter in Time domain
          //2. Cut the 0-padding length first samples (the first 4096 samples, or the first 2048 samples...) => 0-padding filter
          //3. FFT of the 0-padding filter
          //4. FFT of the 0-padding pulse (pulse cut according the 0-padding)
          //5. Scalar product in Frequency domain
          return(EXIT_FAILURE);
      }
      
      if ((strcmp(par->EnergyMethod,"WEIGHT") == 0) && (par->OFLib == 1))
      {
          SIXT_ERROR("parameter error: EnergyMethod=WEIGHT => OFLib should be 'no'");
          return(EXIT_FAILURE);
      }
      
      if ((strcmp(par->EnergyMethod,"OPTFILT") == 0) && (strcmp(par->OFNoise,"WEIGHTM") == 0) && (par->OFLib == 0))
      {
          SIXT_ERROR("parameter error: EnergyMethod=OPTFILT && OFNoise=WEIGHTM => OFLib should be 'yes'");
          return(EXIT_FAILURE);
      }
      
      MyAssert((strcmp(par->OFStrategy,"FREE") == 0) || (strcmp(par->OFStrategy,"BYGRADE") == 0) || (strcmp(par->OFStrategy,"FIXED") == 0), 
               "OFStrategy must be FREE, BYGRADE or FIXED");
      
      MyAssert(par->OFLength > 0, "OFLength must be greater than 0");
      
      MyAssert(par->energyPCA1 > 0, "energyPCA1 must be greater than 0");
      MyAssert(par->energyPCA2 > 0, "energyPCA2 must be greater than 0");
      
      MyAssert(par->LbT > 0, "LbT must be greater than 0");
	
  } else {
	SIXT_ERROR("failed reading the Rcmethod parameter");
	return(EXIT_FAILURE);
  }
  return(status);
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 ************************************************************
* MyAssert function: This function displays an error message if the condition in 'expr' is true
*
* Parameters:
* - expr: Condition to be true in order to display the error message
* - msg: Message to be displayed
******************************************************************************/
void MyAssert(int expr, char* msg)
{
    if (expr == 0)
    {
        printf("%s %s %s"," Assertion failure: ",msg,"\n");
        abort();
    }
}
/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
* checkxmls function: Check if the XML file used to build the library is the same to be used to recconstruct (by checking the checksums)
*
* Parameters:
* - par: Structure containing the input parameters
******************************************************************************/
int checkXmls(struct Parameters* const par)
{
    // Error status.
    int status=EXIT_SUCCESS;

    fitsfile* libptr = NULL;
    int numberkeywords;
    char *libheaderPrimary = NULL;
    char *xml_pointer = NULL;
    char *XMLFile_pointer = NULL;
    char *HISTORY_pointer = NULL;
    char *space_pointer = NULL;
    //char *slash_pointer = NULL;
    char libXMLfile[1024] = "x";
    char libXMLfile1[1024] = "x";
    char libXMLfile2[1024] = "x";
    char wholelibXMLfile[1024] = "x";
    char charAux[1024] = "x";
    char reconsXMLfile[1024] = "x";
    int lengthstr = 0;

    // Move to "Primary" HDU of the library file
    fits_open_file(&libptr, par->LibraryFile, READONLY, &status);
    if (status != 0)
    {
        SIXT_ERROR("File given in LibraryFile does not exist");
        return(EXIT_FAILURE);
    }

    if (fits_movabs_hdu(libptr, 1, NULL, &status))
    {
        return(EXIT_FAILURE);
    }
    // and read full Primary HDU and store it in 'libheaderPrimary'
    if (fits_hdr2str(libptr, 0, NULL, 0,&libheaderPrimary, &numberkeywords, &status))
    {
        free(libheaderPrimary);
        return(EXIT_FAILURE);
    }

    xml_pointer = strstr(libheaderPrimary,".xml");
    //printf("%s %s","xml_pointer0: ","\n");
    //printf(xml_pointer);
    if(!xml_pointer)
    {
        SIXT_ERROR("XML file info not included in Primary HDU in library file");
        return(EXIT_FAILURE);
    }
    XMLFile_pointer = strstr(libheaderPrimary,"XMLFile = ");
    //printf("%s %s","XMLFile_pointer0: ","\n");
    //printf(XMLFile_pointer);
    if(!XMLFile_pointer)
    {
        SIXT_ERROR("XML file info not included in Primary HDU in library file");
        return(EXIT_FAILURE);
    }
    strncpy(wholelibXMLfile,XMLFile_pointer+10,xml_pointer-XMLFile_pointer-10+4);  // 10 -> "XMLFile = ", 4 -> ".xml"
    //strcpy(wholelibXMLfile,"GSHISTORY P53 FCHISTORY P133 _rl8192_pB450.xml");
    //strcpy(wholelibXMLfile,"GSFCHISTORY P133 _rl8192_pB450.xml");
    //strcpy(wholelibXMLfile,"GSFC_rl8192_pB450.xml");
    //printf("%s %s","wholelibXMLfile: ","\n");
    //printf("%s %s",wholelibXMLfile,"\n");

    int HISTORYnum = 0;
    int k;
    int lengthTOTAL = (int)strlen(wholelibXMLfile);
    int spaces[lengthTOTAL];
    //printf("%s %d %s","lengthTOTAL=: ",lengthTOTAL,"\n");
    k=0;
    do
    {
        //printf("%s %c %s","wholelibXMLfile(k)",wholelibXMLfile[k],"\n");
        if(wholelibXMLfile[k]==' ')
        {
            if (HISTORYnum % 2 == 0)
            {  // Par
                spaces[HISTORYnum] = k-7;
            }
            else
            {
                spaces[HISTORYnum] = k;
            }

            HISTORYnum++;
            //printf("%s %d %s","k_' '=: ",k,"\n");
        }
        k++;
    }while(k<=lengthTOTAL);

    //for (int i=0;i<HISTORYnum;i++)
    //    printf("%s %d %s","white spaces: ",spaces[i],"\n");

    if (HISTORYnum != 0) HISTORYnum = HISTORYnum/2;

    //printf("%s %d %s","Numero de HISTORY: ",HISTORYnum,"\n");

    if (HISTORYnum == 0)
    {
        strcpy(libXMLfile,wholelibXMLfile);
    }
    else
    {
        for (int i=0;i<=HISTORYnum;i++)
        {
            if (i == 0) // First
            {
                subString (wholelibXMLfile, 0, spaces[i], libXMLfile);
            }
            else
            {
                if (i == HISTORYnum) // Last
                {
                    subString (wholelibXMLfile, spaces[2*i-1]+1, lengthTOTAL-spaces[2*i-1]-1, libXMLfile2);
                    strcat(libXMLfile,libXMLfile2);
                    memset(libXMLfile2,0,1024);
                }
                else
                {
                    subString (wholelibXMLfile, spaces[2*i-1]+1, spaces[2*i]-spaces[2*i-1]-1, libXMLfile2);
                    strcat(libXMLfile,libXMLfile2);
                    memset(libXMLfile2,0,1024);
                }
            }
        }
    }

    strcpy(reconsXMLfile,par->XMLFile);

    if (libheaderPrimary != NULL)
    {
        free(libheaderPrimary);
    }

    fits_close_file(libptr,&status);

    FILE *fp_libXMLfile, *fp_reconsXMLfile;
    size_t len_libXMLfile, len_reconsXMLfile;
    char buf_libXMLfile[4096], buf_reconsXMLfile[4096];
    unsigned checksum_libXMLfile, checksum_reconsXMLfile;

    if (NULL == (fp_libXMLfile = fopen(libXMLfile, "rb")))
    {
          printf("Unable to open XML from library, %s, for reading it and calculate its checksum.\n", libXMLfile);
          return -1;
    }
    len_libXMLfile = fread(buf_libXMLfile, sizeof(char), sizeof(buf_libXMLfile), fp_libXMLfile);
    //printf("%d bytes read\n", len_libXMLfile);
    checksum_libXMLfile = checksum(buf_libXMLfile, len_libXMLfile, 0);


    if (NULL == (fp_reconsXMLfile = fopen(reconsXMLfile, "rb")))
    {
          printf("Unable to open provided XML, %s, for reading it and calculate its checksum\n", reconsXMLfile);
          return -1;
    }
    len_reconsXMLfile = fread(buf_reconsXMLfile, sizeof(char), sizeof(buf_reconsXMLfile), fp_reconsXMLfile);
    //printf("%d bytes read\n", len_reconsXMLfile);
    checksum_reconsXMLfile = checksum(buf_reconsXMLfile, len_reconsXMLfile, 0);

    if (checksum_libXMLfile == checksum_reconsXMLfile) status = 0;
    else status = 1;

    if (status != 0)
    {
        SIXT_ERROR("XML file from library FITS file and from input parameter do not match (checksum)");
        printf("The checksum of XML from library, %s, is %#x\n", libXMLfile, checksum_libXMLfile);
        printf("The checksum of provided XML, %s, is %#x\n", reconsXMLfile, checksum_reconsXMLfile);
        return(EXIT_FAILURE);
    }

    return(status);
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* subString: Extract some elements from an array of characters
*
* Parameters:
* - input: Array of characters from which extract some elements
* - offset: Offset
* - len: Length (number of elements to extract)
* - dest: Array of characters where write the elements extracted
******************************************************************************/
char* subString (const char* input, int offset, int len, char* dest)
{
  int input_len = strlen (input);

  if (offset + len > input_len)
  {
     return NULL;
  }

  strncpy (dest, input + offset, len);
  return dest;
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
