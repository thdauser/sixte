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
    PulsesCollection* pulsesAll2 = newPulsesCollection(&status);
    CHECK_STATUS_BREAK(status);    
    
    if(!strcmp(par.Rcmethod,"PP")){
    
	  initializeReconstruction(reconstruct_init,par.OptimalFilterFile,par.PulseLength,
    		par.PulseTemplateFile,par.Threshold,par.Calfac,par.NormalExclusion,
    		par.DerivateExclusion,par.SaturationValue,&status);
	  
    }else{
	  initializeReconstructionSIRENA(reconstruct_init_sirena, par.RecordFile, par.LibraryFile,
		par.tauFall, par.PulseLength, par.scaleFactor, par.samplesUp, par.nSgms, par.crtLib,
		par.mode, par.LrsT, par.LbT, par.NoiseFile, par.FilterDomain, par.FilterMethod,
		par.calibLQ, par.b_cF,par.b_cF, par.monoenergy, par.intermediate, par.detectFile,
		par.filterFile, par.RecordFileCalib2, par.monoenergy2, par.clobber, &status);
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
    while(getNextRecord(record_file,record,&status)){

	if(!strcmp(par.Rcmethod,"PP")){
	      reconstructRecord(record,event_list,reconstruct_init,0,&status);
	}else{
	    nrecord = nrecord + 1;
	    if(nrecord == record_file->nrows) lastRecord=1;
	    reconstructRecordSIRENA(record,event_list,reconstruct_init_sirena,
				    lastRecord, nrecord, &pulsesAll, &optimalFilter, &status);
	}	
 	CHECK_STATUS_BREAK(status);

	saveEventListToFile(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
 	CHECK_STATUS_BREAK(status);
	//Reinitialize event list
	event_list->index=0;
	
    }
    
    // If SIRENA processing is in calibration (b,c calibration factors calculations)
    if(!strcmp(par.Rcmethod,"SIRENA") && par.mode==0 && par.crtLib==0){
        // Open second record file (calibration)
        TesTriggerFile* record_file2 = openexistingTesTriggerFile(par.RecordFileCalib2,keywords,&status);
        CHECK_STATUS_BREAK(status);

	// Build up TesRecord to read the 2nd file
	TesRecord* record2 = newTesRecord(&status);
	allocateTesRecord(record2,record_file2->trigger_size,record_file2->delta_t,0,&status);
	CHECK_STATUS_BREAK(status);
	
	// Iterate of records and do the reconstruction
	lastRecord = 0, nrecord = 0; 
	while(getNextRecord(record_file2,record2,&status)){
	    nrecord = nrecord + 1;
	    if(nrecord == record_file2->nrows) lastRecord=1;
	    reconstructRecordSIRENA(record2, event_list, reconstruct_init_sirena,
				    lastRecord, nrecord, &pulsesAll2, &optimalFilter, &status);
	    saveEventListToFile(outfile,event_list,record2->time,record_file2->delta_t,
				record2->pixid,&status);
	    //Reinitialize event list
	    event_list->index=0;
	}
	CHECK_STATUS_BREAK(status);
	double b_cF = 0.0, c_cF = 0.0;
	runEnergyCalib(reconstruct_init_sirena, pulsesAll, pulsesAll2, &b_cF, &c_cF);
	
	// Copy calib keywords to event file
	writeCalibKeywords(outfile->fptr,b_cF, c_cF,&status);
	CHECK_STATUS_BREAK(status);
	
	// free calibration memory
	freeTesTriggerFile(&record_file2,&status);
	freeTesRecord(&record2);
    }
    
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
    freePulsesCollection(pulsesAll2);
    freeOptimalFilterSIRENA(optimalFilter);
    freeTesTriggerFile(&record_file,&status);
    freeTesEventFile(outfile,&status);
    freeTesEventList(event_list);
    freeTesRecord(&record);
    freeSixtStdKeywords(keywords);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
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

	status=ape_trad_query_double("SaturationValue", &par->SaturationValue);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the SaturationValue parameter");
		return(status);
	}
	
  }else if(strcmp(par->Rcmethod,"SIRENA")==0){
	
	// SIRENA parameters
	status=ape_trad_query_string("LibraryFile", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the library file");
		return(status);
	}
	strcpy(par->LibraryFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_double("scaleFactor", &par->scaleFactor);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the scaleFactor parameter");
		return(status);
	}
  
	status=ape_trad_query_double("tauFall", &par->tauFall);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the tauFall parameter");
		return(status);
	}
  
	status=ape_trad_query_double("samplesUp", &par->samplesUp);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the samplesUp parameter");
		return(status);
	}
  
	status=ape_trad_query_double("nSgms", &par->nSgms);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the nSgms parameter");
		return(status);
	}
  
	status=ape_trad_query_int("mode", &par->mode);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the mode parameter");
		return(status);
	}
	
	status=ape_trad_query_int("crtLib", &par->crtLib);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the crtLib parameter");
		return(status);
	}

	if(par->mode==1 && par->crtLib==1){
		SIXT_ERROR("parameter error: mode=1 (production) and crtLib=1 (library creation - calibration)");
		return(status);
	}
  
	status=ape_trad_query_double("LrsT", &par->LrsT);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the LrsT parameter");
		return(status);
	}

	status=ape_trad_query_double("LbT", &par->LbT);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the LbT parameter");
		return(status);
	}
	
	if(par->mode==0){
		status=ape_trad_query_double("monoenergy", &par->monoenergy);
		if (EXIT_SUCCESS!=status) {
		      SIXT_ERROR("failed reading the monoenergy parameter");
		      return(status);
		}
	}
  
	if(par->crtLib==0){
	 
		status=ape_trad_query_string("NoiseFile", &sbuffer);
		if (EXIT_SUCCESS!=status) {
		      SIXT_ERROR("failed reading the name of the noise file");
		      return(status);
		}
		strcpy(par->NoiseFile, sbuffer);
		free(sbuffer);
		
		status=ape_trad_query_string("FilterDomain", &sbuffer);
		if (EXIT_SUCCESS!=status) {
		  SIXT_ERROR("failed reading the Filter domain");
		  return(status);
		}
		strcpy(par->FilterDomain, sbuffer);
		free(sbuffer);
		
		status=ape_trad_query_string("FilterMethod", &sbuffer);
		if (EXIT_SUCCESS!=status) {
		  SIXT_ERROR("failed reading the Filter method");
		  return(status);
		}
		strcpy(par->FilterMethod, sbuffer);
		free(sbuffer);
  
		status=ape_trad_query_int("calibLQ", &par->calibLQ);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the calibLQ parameter");
			return(status);
		}
	  
		status=ape_trad_query_double("b_cF", &par->b_cF);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the b_cF parameter");
			return(status);
		}
	  
		status=ape_trad_query_double("c_cF", &par->c_cF);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the c_cF parameter");
			return(status);
		}

		status=ape_trad_query_int("intermediate", &par->intermediate);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the intermediate parameter");
			return(status);
		}
		if(par->intermediate){
			status=ape_trad_query_string("detectFile", &sbuffer);
			if (EXIT_SUCCESS!=status) {
			    SIXT_ERROR("failed reading the name of the detections output intermediate file");
			    return(status);
			}
			strcpy(par->detectFile, sbuffer);
			free(sbuffer);
			
			status=ape_trad_query_string("filterFile", &sbuffer);
			if (EXIT_SUCCESS!=status) {
			    SIXT_ERROR("failed reading the name of the filters output intermediate file");
			    return(status);
			}
			strcpy(par->filterFile, sbuffer);
			free(sbuffer); 
		} // if intermediate
		if(par->mode == 0){
			status=ape_trad_query_string("RecordFileCalib2", &sbuffer);
			if (EXIT_SUCCESS!=status) {
			    SIXT_ERROR("failed reading the name of second calibration file");
			    return(status);
			}
			strcpy(par->RecordFileCalib2, sbuffer);
			free(sbuffer);
			
			status=ape_trad_query_double("monoenergy2", &par->monoenergy2);
			if (EXIT_SUCCESS!=status) {
			    SIXT_ERROR("failed reading the monoenergy2 parameter");
			  return(status);
			}
		} // if calibration
		
	} // if crtLib

  } else {
	SIXT_ERROR("failed reading the Rcmethod parameter");
	return(EXIT_FAILURE);
  }
  return(status);
}
