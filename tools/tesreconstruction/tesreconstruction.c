/*
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


   Copyright 2015 Philippe Peille, IRAP
*/

#include "tesreconstruction.h"

////////////////////////////////////
/** Main procedure. */
int tesreconstruction_main() {
  // Containing all programm parameters read by PIL.
  struct Parameters par;

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
    
    // Sixt standard keywords structure
    SixtStdKeywords* keywords = newSixtStdKeywords(&status);
    CHECK_STATUS_BREAK(status);

    // Open record file
    TesTriggerFile* record_file = openexistingTesTriggerFile(par.RecordFile,keywords,&status);
    CHECK_STATUS_BREAK(status);

    //Open outfile
    TesEventFile * outfile = opennewTesEventFile(par.TesEventFile,
    		keywords,
    		par.clobber,
    		&status);
    CHECK_STATUS_BREAK(status);

    // Initialize PP data structures needed for pulse filtering
    ReconstructInit* reconstruct_init = newReconstructInit(&status);
    CHECK_STATUS_BREAK(status);
    
    // Initialize SIRENA data structures needed for pulse filtering
    ReconstructInitSIRENA* reconstruct_init_sirena = newReconstructInitSIRENA(&status);
    CHECK_STATUS_BREAK(status);
    PulsesCollection* pulsesAll = newPulsesCollection(&status);
    CHECK_STATUS_BREAK(status);  
    OptimalFilterSIRENA* optimalFilter = newOptimalFilterSIRENA(&status);
    CHECK_STATUS_BREAK(status);// define a second structure for calibration
    
    if(!strcmp(par.Rcmethod,"PP")){
	  initializeReconstruction(reconstruct_init,par.OptimalFilterFile,par.PulseLength,
    		par.PulseTemplateFile,par.Threshold,par.Calfac,par.NormalExclusion,
    		par.DerivateExclusion,par.SaturationValue,&status);
    }else{
	  initializeReconstructionSIRENA(reconstruct_init_sirena, par.RecordFile, record_file->fptr, par.LibraryFile, par.TesEventFile, 
		par.PulseLength, par.scaleFactor, par.samplesUp, par.samplesDown, par.nSgms, par.detectSP, par.mode, par.detectionMode, par.LrsT, par.LbT, par.NoiseFile, 
		par.FilterDomain, par.FilterMethod, par.EnergyMethod, par.filtEev, par.OFNoise, par.LagsOrNot, par.OFIter, par.OFLib, par.OFInterp, par.OFStrategy, par.OFLength,
		par.monoenergy, par.hduPRECALWN, par.hduPRCLOFWM, par.largeFilter, par.intermediate, par.detectFile, par.filterFile, par.clobber, par.EventListSize, par.SaturationValue,
		par.tstartPulse1, par.tstartPulse2, par.tstartPulse3, par.energyPCA1, par.energyPCA2, par.XMLFile, &status);
	  
	  // Read the grading data from the XML file and store it in 'reconstruct_init_sirena->grading'
	  reconstruct_init_sirena->grading = NULL;
	  reconstruct_init_sirena->grading = (Grading*)malloc(sizeof(Grading));
	  
	  reconstruct_init_sirena->grading->ngrades = 0;
	  reconstruct_init_sirena->grading->value  = NULL;
	  reconstruct_init_sirena->grading->gradeData = NULL;
	  
	  AdvDet *det = newAdvDet(&status);
	  CHECK_STATUS_BREAK(status);
	  det = loadAdvDet(par.XMLFile, &status);
	  CHECK_STATUS_BREAK(status);
	  if (det->pix->grades == NULL)
	  {
		SIXT_ERROR("The provided XMLFile does not have the grading info");
		return(EXIT_FAILURE);
  	  }
	  
	  reconstruct_init_sirena->grading->ngrades=det->pix->ngrades;
	  reconstruct_init_sirena->grading->value = gsl_vector_alloc(det->pix->ngrades);
	  reconstruct_init_sirena->grading->gradeData = gsl_matrix_alloc(det->pix->ngrades,2);
	  for (int i=0;i<det->pix->ngrades;i++)
	  {
	      gsl_vector_set(reconstruct_init_sirena->grading->value,i,det->pix->grades[i].value);
	      gsl_matrix_set(reconstruct_init_sirena->grading->gradeData,i,0,det->pix->grades[i].gradelim_pre);
	      gsl_matrix_set(reconstruct_init_sirena->grading->gradeData,i,1,det->pix->grades[i].gradelim_post);
	  }
	  destroyAdvDet(&det);
    }  
    CHECK_STATUS_BREAK(status);

    // Build up TesRecord to read the file
    TesRecord* record = newTesRecord(&status);
    allocateTesRecord(record,record_file->trigger_size,record_file->delta_t,0,&status);
    CHECK_STATUS_BREAK(status);

    //Build up TesEventList to recover the results of the reconstruction
    TesEventList* event_list = newTesEventList(&status);
    allocateTesEventListTrigger(event_list,par.EventListSize,&status);
    CHECK_STATUS_BREAK(status);

    // Iterate of records and do the reconstruction
    int lastRecord = 0, nrecord = 0; //last record required for SIRENA library creation
    while(getNextRecord(record_file,record,&status))
    {
      if(!strcmp(par.Rcmethod,"PP"))
      {
	    /*nrecord = nrecord + 1;
	    if(nrecord == record_file->nrows) lastRecord=1;
	    printf("%s %d %s","**TESRECONSTRUCTION nrecord = ",nrecord,"\n");*/
	    
	    reconstructRecord(record,event_list,reconstruct_init,0,&status);
      }
      else
      {
	    nrecord = nrecord + 1;
	    if(nrecord == record_file->nrows) lastRecord=1;
	    /*if(nrecord < 349) 
            {
	      continue;
	    }
            else if(nrecord > 349)
            {
	      status=1;
	      CHECK_STATUS_BREAK(status);
	    }*/
	    /*if(nrecord > 1)
	    {
	    	status=1;
	        CHECK_STATUS_BREAK(status);
	    }*/
	    if ((strcmp(par.EnergyMethod,"I2R") == 0) || (strcmp(par.EnergyMethod,"I2RALL") == 0) || (strcmp(par.EnergyMethod,"I2RNOL") == 0) || (strcmp(par.EnergyMethod,"I2RFITTED") == 0))
	    {
	    	strcpy(reconstruct_init_sirena->EnergyMethod,par.EnergyMethod);
	    }
	
	    //printf("%s %d %s","**TESRECONSTRUCTION nrecord = ",nrecord,"\n");
	    reconstructRecordSIRENA(record,event_list,reconstruct_init_sirena,
				    lastRecord, nrecord, &pulsesAll, &optimalFilter, &status);
      }
      CHECK_STATUS_BREAK(status);

      if ((strcmp(par.EnergyMethod,"PCA") != 0) || ((strcmp(par.EnergyMethod,"PCA") == 0) && lastRecord == 1))
      {
	  saveEventListToFile(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
	  CHECK_STATUS_BREAK(status);
	  
	  //Reinitialize event list
	  event_list->index=0;
      }
    }
    
    if ((!strcmp(par.Rcmethod,"SIRENA")) && (pulsesAll->ndetpulses == 0))  printf("%s","WARNING: no pulses have been detected\n");
    
    // Copy trigger keywords to event file
    copyTriggerKeywords(record_file->fptr,outfile->fptr,&status);
    CHECK_STATUS_BREAK(status);
    // Save GTI extension to event file
    GTI* gti=getGTIFromFileOrContinuous("none",keywords->tstart, keywords->tstop,keywords->mjdref, &status);
    saveGTIExt(outfile->fptr, "STDGTI", gti, &status);    
    CHECK_STATUS_BREAK(status);

    //Free memory
    freeReconstructInit(reconstruct_init);
    freeReconstructInitSIRENA(reconstruct_init_sirena);
    freePulsesCollection(pulsesAll);
    freeOptimalFilterSIRENA(optimalFilter);
    freeTesTriggerFile(&record_file,&status);
    freeTesEventFile(outfile,&status);
    freeTesEventList(event_list);
    freeTesRecord(&record);
    freeSixtStdKeywords(keywords);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.
  
  if (EXIT_SUCCESS==status) 
  {
	headas_chat(3, "finished successfully!\n\n");
	return(EXIT_SUCCESS);
  } 
  else 
  {
	return(status);
  }
}

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

  status=ape_trad_query_int("PulseLength", &par->PulseLength);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the PulseLength parameter");
	  return(status);
  }
  assert(&par->PulseLength > 0);

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

	/*status=ape_trad_query_double("SaturationValue", &par->SaturationValue);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the SaturationValue parameter");
		return(status);
	}*/
	
  }else if(strcmp(par->Rcmethod,"SIRENA")==0){
	
	// SIRENA parameters
	status=ape_trad_query_string("LibraryFile", &sbuffer);
	strcpy(par->LibraryFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_double("scaleFactor", &par->scaleFactor);
    
	status=ape_trad_query_double("samplesUp", &par->samplesUp);
        
        status=ape_trad_query_double("samplesDown", &par->samplesDown);
  
	status=ape_trad_query_double("nSgms", &par->nSgms);
        
        status=ape_trad_query_int("detectSP", &par->detectSP);
  
	status=ape_trad_query_int("mode", &par->mode);
        
        status=ape_trad_query_string("detectionMode", &sbuffer);
	strcpy(par->detectionMode, sbuffer);
	free(sbuffer);
  
	status=ape_trad_query_double("LrsT", &par->LrsT);

	status=ape_trad_query_double("LbT", &par->LbT);

	status=ape_trad_query_int("intermediate", &par->intermediate);

	status=ape_trad_query_string("detectFile", &sbuffer);
	strcpy(par->detectFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("filterFile", &sbuffer);
	strcpy(par->filterFile, sbuffer);
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
        
        status=ape_trad_query_double("filtEev", &par->filtEev);

	status=ape_trad_query_string("OFNoise", &sbuffer);
	strcpy(par->OFNoise, sbuffer);
	free(sbuffer);
	
	status=ape_trad_query_int("LagsOrNot", &par->LagsOrNot);

	status=ape_trad_query_int("OFIter", &par->OFIter);

	status=ape_trad_query_bool("OFLib", &par->OFLib);
	
	//status=ape_trad_query_string("OFInterp", &sbuffer);
        //strcpy(par->OFInterp, sbuffer);
        //free(sbuffer);
	strcpy(par->OFInterp, "DAB");
	//strcpy(par->OFInterp, "MF");
	
	status=ape_trad_query_string("OFStrategy", &sbuffer);
	strcpy(par->OFStrategy, sbuffer);
	free(sbuffer);
	
	status=ape_trad_query_int("OFLength", &par->OFLength);

	status=ape_trad_query_int("tstartPulse1", &par->tstartPulse1);
	
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
	
	//assert((par->mode ==0) || (par->mode ==1));
	/*int mode_0_1;
	mode_0_1 = ((par->mode ==0) || (par->mode ==1));
	assert(mode_0_1);*/
	MyAssert((par->mode == 0) || (par->mode == 1), "mode must be 0 or 1");
	  
	//assert((par->intermediate == 0) || (par->intermediate == 1));
	MyAssert((par->intermediate == 0) || (par->intermediate == 1), "intermediate must be 0 or 1");
	
	//assert(&par->monoenergy > 0);
	MyAssert(&par->monoenergy > 0, "monoenergy must be greater than 0");
	
	//assert((strcmp(par->FilterDomain,"T") == 0) || (strcmp(par->FilterDomain,"F") == 0));
	MyAssert((strcmp(par->FilterDomain,"T") == 0) || (strcmp(par->FilterDomain,"F") == 0), "FilterDomain must be T or F");
	
	//assert((strcmp(par->FilterMethod,"F0") == 0) || (strcmp(par->FilterMethod,"B0") == 0));
	MyAssert((strcmp(par->FilterMethod,"F0") == 0) || (strcmp(par->FilterMethod,"B0") == 0),"FilterMethod must be F0 or B0");
	
	MyAssert((strcmp(par->EnergyMethod,"OPTFILT") == 0) || (strcmp(par->EnergyMethod,"WEIGHT") == 0) || (strcmp(par->EnergyMethod,"WEIGHTN") == 0) ||
		(strcmp(par->EnergyMethod,"I2R") == 0) || (strcmp(par->EnergyMethod,"I2RALL") == 0) || (strcmp(par->EnergyMethod,"I2RNOL") == 0) || 
		(strcmp(par->EnergyMethod,"I2RFITTED") == 0) || (strcmp(par->EnergyMethod,"PCA") == 0), "EnergyMethod must be OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL, I2RFITTED or PCA");
	
	MyAssert((strcmp(par->OFNoise,"NSD") == 0) || (strcmp(par->OFNoise,"WEIGHTM") == 0), "OFNoise must be NSD or WEIGHTM");
	
	//assert((par->LagsOrNot ==0) || (par->LagsOrNot ==1));
	MyAssert((par->LagsOrNot ==0) || (par->LagsOrNot ==1), "LagsOrNot must me 0 or 1");

	if (((strcmp(par->EnergyMethod,"WEIGHT") == 0) || (strcmp(par->EnergyMethod,"WEIGHTN") == 0)) && (par->LagsOrNot == 1))
	{
		SIXT_ERROR("parameter error: EnergyMethod=WEIGHT/WEIGHTN and Lags not implemented yet");
		return(EXIT_FAILURE);
	}
	
	//assert((par->OFIter ==0) || (par->OFIter ==1));
	MyAssert((par->OFIter ==0) || (par->OFIter ==1), "OFIter must be 0 or 1");
	
	if ((par->OFLib == 1) && (strcmp(par->FilterMethod,"F0") != 0))
	{
		SIXT_ERROR("parameter error: If OFLib=yes => FilterMethod must be F0");
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
	
	//assert((strcmp(par->OFStrategy,"FREE") == 0) || (strcmp(par->OFStrategy,"BASE2") == 0) || (strcmp(par->OFStrategy,"BYGRADE") == 0) || (strcmp(par->OFStrategy,"FIXED") == 0));
	MyAssert((strcmp(par->OFStrategy,"FREE") == 0) || (strcmp(par->OFStrategy,"BASE2") == 0) || (strcmp(par->OFStrategy,"BYGRADE") == 0) || (strcmp(par->OFStrategy,"FIXED") == 0), 
		 "OFStrategy must be FREE, BASE2, BYGRADE or FIXED");
	
	//assert(&par->OFLength > 0);
	MyAssert(&par->OFLength > 0, "OFLength must be greater than 0");
	
	/*if (((strcmp(par->EnergyMethod,"I2R") == 0) || (strcmp(par->EnergyMethod,"I2RALL") == 0) || (strcmp(par->EnergyMethod,"I2RNOL") == 0)) && (par->tstartPulse1 == 0))
	{
		printf("%s %d %s","Error",status,"\n");
		SIXT_ERROR("parameter error: EnergyMethod=I2R/I2RALL/I2RNOL and tstartPulse1=0 (If I2R/I2RALL/I2RNOL, tstartPulsex must be always provided)");
		return(EXIT_FAILURE);
	}*/	
	
	//assert(&par->energyPCA1 > 0);
	MyAssert(&par->energyPCA1 > 0, "energyPCA1 must be greater than 0");
	//assert(&par->energyPCA2 > 0);
	MyAssert(&par->energyPCA2 > 0, "energyPCA2 must be greater than 0");
	
  } else {
	SIXT_ERROR("failed reading the Rcmethod parameter");
	return(EXIT_FAILURE);
  }
  return(status);
}

void MyAssert(int expr, char* msg)
{
    if (expr == 0)
    {
        printf("%s %s %s"," Assertion failure: ",msg,"\n");
        abort();
    }
}
