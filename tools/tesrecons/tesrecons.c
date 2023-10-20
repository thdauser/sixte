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
*                      TESRECONS
*
*  File:       tesLIB.c
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*              Philippe Peille, IRAP
*                                                                     
***********************************************************************/

#include "tesrecons.h"

/***** SECTION 1 ************************************
* MAIN function: This function is mainly a wrapper to pass a data file to the SIRENA tasks in order to reconstruct the energies.
*
* The user must supply the following input parameters (.par file).
* 
* Parameters:
* 
* - RecordFile: Record FITS file
*   - If RecordFile starts with '@' it provides a file text containing several record input FITS files
* - TesEventFile: Output event list file
* - LibraryFile: Library FITS file to be created
* - XMLFile: XML input FITS file with instrument definition
* - preBuffer: Some samples added before the starting time of a pulse (number of samples added read from the XML file)
*              SIRENA's format XML file (grading=>pre,post and pB) or new format XML file (grading=>pre,post and filtlen)
*                                      pre=494, post=8192, pB=1000                          pre=494, post=7192, filtlen=8192
*                                                                                             preBuffer=filtlen-post
* - EventListSize: Default size of the event list
* - clobber:Overwrite or not output files if exist (1/0)
* - history: write program parameters into output file
* - scaleFactor: Detection scale factor for initial filtering
* - samplesUp: Number of consecutive samples up for threshold trespassing (only used with STC detection mode)
* - samplesDown: Number of consecutive samples below the threshold to look for other pulse (only used with STC detection mode)
* - nSgms: Number of quiescent-signal standard deviations to establish the threshold through the kappa-clipping algorithm
* - LbT: Baseline averaging length (seconds)


* - FilterDomain: Filtering Domain: Time (T) or Frequency (F)
******* - FilterMethod: Filtering Method: F0 (deleting the zero frequency bin) or B0 (deleting the baseline) or F0B0 (deleting always the baseline)
* - FilterMethod: Filtering Method: F0 (deleting the zero frequency bin) or B0 (deleting the baseline)
* - EnergyMethod: Energy calculation Method: OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL, I2RFITTED, I2RDER or PCA
* - Ifit: Constant to apply the I2RFITTED conversion

* - intermediate: Write or not intermediate files (1/0)
* - detectFile: Intermediate detections file (if intermediate=1)
* 
* Steps:
* 
* - Register HEATOOL
* - Reading all programm parameters by using PIL
* - Read XML info
* - getSamplingrate_trigreclength => Obtain the 'trig_reclength' and the sampling rate
* - Sixt standard keywords structure
* - Open output FITS file
* - Initialize data structures needed for pulse filtering
* - Read the grading data from the XML file and store it in 'reconstruct_init_sirena->grading'
* - Build up TesEventList
* - Call SIRENA to build reconstruct the energies
* - Save GTI extension to event file
* - Free memory
*****************************************************/
int tesrecons_main() {
  printf("Running TESRECONS\n");
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
  set_toolname("tesrecons");
  //set_toolversion("0.05");
  
  do { // Beginning of the ERROR handling loop (will at
       // most be run once).
    headas_chat(3, "initialize ...\n");

    // Get program parameters
    status=getpar_tesrecons(&par);
    //CHECK_STATUS_BREAK(status);
    if (status != EXIT_SUCCESS)
    {
    	//printf("getpar_tesrecons error\n");
    	return(status);
    }
    par.opmode = 1; // Reconstructing the energies

    if ((strcmp(par.EnergyMethod,"I2RFITTED") == 0) && (par.Ifit == 0.0))
    {
        SIXT_ERROR("Ifit value must be provided");
        return(EXIT_FAILURE);
    }
        
    AdvDet *det = newAdvDet(&status);
    CHECK_STATUS_BREAK(status);
    double sampling_rate_XML = -999.; // Sampling rate from XML
    // Read XML info
    det = loadAdvDet(par.XMLFile, &status);
    CHECK_STATUS_BREAK(status);
    sampling_rate_XML = det->SampleFreq;
    if (par.preBuffer == 1)    printf("%s","Attention: preBuffer used => Parameters of library filters read from XML file\n");

    // Obtain the 'trig_reclength' and the sampling rate
    double sampling_rate = -999.0;
    int trig_reclength = -999;
    sampling_rate = sampling_rate_XML;
    int numfits = -999;
    status = getSamplingrate_trigreclength (par.RecordFile,par,&sampling_rate,&trig_reclength, &numfits);
    if (status != EXIT_SUCCESS)
    {
        SIXT_ERROR("Error in 'getSamplingrate_trigreclength' function");
        return(EXIT_FAILURE);
    }

    // Sixt standard keywords structure
    SixtStdKeywords* keywords = newSixtStdKeywords(&status);
    CHECK_STATUS_BREAK(status);
    
    //Open outfile 
    TesEventFile * outfile = opennewTesEventFileSIRENA(par.TesEventFile,
                                                 keywords,
                                                 SIRENA_VERSION,
                                                 par.clobber,
                                                 &status);
    CHECK_STATUS_BREAK(status);
    
    // Initialize data structures needed
    ReconstructInitSIRENA* reconstruct_init_sirena = newReconstructInitSIRENA();
    CHECK_STATUS_BREAK(status);
    PulsesCollection* pulsesAll = newPulsesCollection();
    CHECK_STATUS_BREAK(status);  
    OptimalFilterSIRENA* optimalFilter = newOptimalFilterSIRENA();
    CHECK_STATUS_BREAK(status);// define a second structure for calibration
    
    // Read the grading data from the XML file and store it in 'reconstruct_init_sirena->grading'
    status = fillReconstructInitSIRENAGrading (par, det, &reconstruct_init_sirena);
    destroyAdvDet(&det);

    // Build up TesEventList
    TesEventList* event_list = newTesEventListSIRENA(&status);
    allocateTesEventListTrigger(event_list,par.EventListSize,&status);
    CHECK_STATUS_BREAK(status);
            
    // Call SIRENA to reconstruct
    status = callSIRENA(par.RecordFile, keywords, reconstruct_init_sirena, par, sampling_rate, &trig_reclength, pulsesAll, outfile);

    // Save GTI extension to event file
    GTI* gti=getGTIFromFileOrContinuous("none",keywords->tstart, keywords->tstop,keywords->mjdref, &status);
    saveGTIExt(outfile->fptr, "STDGTI", gti, &status);    
    CHECK_STATUS_BREAK(status);
    
    //Free memory
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
* getpar_tesrecons function: This function gets the input parameter from the command line or their default values from the tesrecons.par file
*
* Parameters:
* - par: Structure containing the input parameters
******************************************************************************/
int getpar_tesrecons(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

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

  status=ape_trad_query_string("LibraryFile", &sbuffer);
  strcpy(par->LibraryFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file");
    return(status);
  }
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("preBuffer", &par->preBuffer);

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
    //printf("failed reading the history parameter\n");
    return(status);
  }

  status=ape_trad_query_double("scaleFactor", &par->scaleFactor);
  status=ape_trad_query_int("samplesUp", &par->samplesUp);
  status=ape_trad_query_int("samplesDown", &par->samplesDown);
  status=ape_trad_query_double("nSgms", &par->nSgms);

  status=ape_trad_query_string("detectionMode", &sbuffer);
  strcpy(par->detectionMode, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("detectSP", &par->detectSP);
      
  status=ape_trad_query_double("LbT", &par->LbT);

  status=ape_trad_query_int("intermediate", &par->intermediate);
  status=ape_trad_query_string("detectFile", &sbuffer);
  strcpy(par->detectFile, sbuffer);
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

  status=ape_trad_query_int("OFLengthNotPadded", &par->OFLengthNotPadded);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the OFLengthNotPadded parameter");
    return(status);
  }
  assert(par->OFLengthNotPadded > 0);

  status=ape_trad_query_int("errorT", &par->errorT);
  status=ape_trad_query_int("Sum0Filt", &par->Sum0Filt);
  status=ape_trad_query_string("tstartPulse1", &sbuffer);
  strcpy(par->tstartPulse1, sbuffer);
  free(sbuffer);
  int isNumber = 1;
  for (int i = 0; i < (int)(strlen(par->tstartPulse1)); i++)
  {
      if (isdigit(par->tstartPulse1[i]) == 0)
      {
          isNumber = 0;
          break;
      }
  }
  status=ape_trad_query_int("tstartPulse2", &par->tstartPulse2);
  status=ape_trad_query_int("tstartPulse3", &par->tstartPulse3);
  status=ape_trad_query_double("energyPCA1", &par->energyPCA1);
  status=ape_trad_query_double("energyPCA2", &par->energyPCA2);
      
  if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("failed reading some TESRECONS parameter");
      return(status);
  }

  MyAssert(par->LbT > 0, "LbT must be greater than 0");

  MyAssert((strcmp(par->FilterDomain,"T") == 0) || (strcmp(par->FilterDomain,"F") == 0), "FilterDomain must be T or F");

  //MyAssert((strcmp(par->FilterMethod,"F0") == 0) || (strcmp(par->FilterMethod,"B0") == 0) || (strcmp(par->FilterMethod,"F0B0") == 0),"FilterMethod must be F0 or B0 or F0B0");
  MyAssert((strcmp(par->FilterMethod,"F0") == 0) || (strcmp(par->FilterMethod,"B0") == 0),"FilterMethod must be F0 or B0");

  MyAssert((strcmp(par->EnergyMethod,"OPTFILT") == 0) || (strcmp(par->EnergyMethod,"WEIGHT") == 0) || (strcmp(par->EnergyMethod,"WEIGHTN") == 0) ||
  (strcmp(par->EnergyMethod,"I2R") == 0) ||	(strcmp(par->EnergyMethod,"I2RFITTED") == 0) || (strcmp(par->EnergyMethod,"PCA") == 0), "EnergyMethod must be OPTFILT, WEIGHT, WEIGHTN, I2R, I2RFITTED or PCA");

  if ((isNumber == 0) && (strcmp(par->FilterDomain,"F") == 0))    // It is only implemented tstartPulse1 as a file for time domain
  {
      SIXT_ERROR("It is not possible to work in FREQUENCY domain if tstartPulse1 is a file => Change FilterDomain to TIME domain (T) ");
      return(EXIT_FAILURE);
  }

  MyAssert((par->intermediate == 0) || (par->intermediate == 1), "intermediate must be 0 or 1");

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

  return(status);
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
