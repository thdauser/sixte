#include "analyse_xms_events.h"


int analyse_xms_events_main() {
  struct Parameters par;
  EventListFile* elf=NULL;
  PatternFile* plf=NULL;

  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("analyse_xms_events");
  set_toolversion("0.05");


  do { // ERROR handling loop

    // Event grade counters.
    long nphotons=0, ngrade0=0, ngrade1=0, ngrade2=0, ngrade3=0;

    // Read parameters by PIL:
    status=analyse_xms_events_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Set the input event file.
    elf=openEventListFile(par.EventList, READWRITE, &status);
    CHECK_STATUS_BREAK(status);

    // Determine mission keywords.
    char comment[MAXMSG], telescop[MAXMSG]={""}, instrume[MAXMSG]={""};
    fits_read_key(elf->fptr, TSTRING, "TELESCOP", &telescop, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(elf->fptr, TSTRING, "INSTRUME", &instrume, comment, &status);
    CHECK_STATUS_BREAK(status);

    // Create and open a new pattern file.
    plf=openNewPatternFile(par.PatternList, 
			   telescop, instrume, "Normal",
			   par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Upate the TLMIN and TLMAX keywords.
    char keystr[MAXMSG];
    long value;
    sprintf(keystr, "TLMIN%d", elf->cpi);
    fits_read_key(elf->fptr, TLONG, keystr, &value, comment, &status);
    CHECK_STATUS_BREAK(status);
    sprintf(keystr, "TLMIN%d", plf->cpi);
    fits_update_key(plf->fptr, TLONG, keystr, &value, "", &status);
    CHECK_STATUS_BREAK(status);
    sprintf(keystr, "TLMAX%d", elf->cpi);
    fits_read_key(elf->fptr, TLONG, keystr, &value, comment, &status);
    CHECK_STATUS_BREAK(status);
    sprintf(keystr, "TLMAX%d", plf->cpi);
    fits_update_key(plf->fptr, TLONG, keystr, &value, "", &status);
    CHECK_STATUS_BREAK(status);

    // Loop over all events in the event file.
    long row;
    for (row=0; row<elf->nrows; row++) {
      
      // Read a new event from the file.
      Event event;
      getEventFromFile(elf, row+1, &event, &status);
      CHECK_STATUS_BREAK(status);

      // Count the total number of photon events in the event list.
      nphotons++;

      // Check the events before and after the current one 
      // within the specified time spans.
      int nbefore_short=0, nbefore_long=0, nbefore_veryshort=0;
      int nafter_short=0, nafter_long=0, nafter_veryshort=0;

      // Former events:
      long row2;
      for (row2=row-1; row2>=0; row2--) {
	Event event2; // Buffer.
	getEventFromFile(elf, row2+1, &event2, &status);
	CHECK_STATUS_BREAK(status);
	
	if (event.time-event2.time > par.PostTrigger*par.TimeUnit) break;
	if ((event.rawx==event2.rawx)&&(event.rawy==event2.rawy)) {
	  nbefore_long++;
	  if (event.time-event2.time < par.PreTrigger*par.TimeUnit) {
	    nbefore_short++;
	  }	
	  if (event.time-event2.time < par.PileupTime) {
	    nbefore_veryshort++;
	  }	
	}
	// Avoid too many unnecessary loop runs.
	if ((nbefore_short>0) || (nbefore_long>0)) break;
      }
      CHECK_STATUS_BREAK(status);

      // Subsequent events:
      for (row2=row+1; row2<elf->nrows; row2++) {
	Event event2; // Buffer.
	getEventFromFile(elf, row2+1, &event2, &status);
	CHECK_STATUS_BREAK(status);

	if (event2.time-event.time > par.PostTrigger*par.TimeUnit) break;
	if ((event.rawx==event2.rawx)&&(event.rawy==event2.rawy)) {
	  nafter_long++;
	  if (event2.time-event.time < par.PreTrigger*par.TimeUnit) {
	    nafter_short++;
	  }
	  if (event2.time-event.time < par.PileupTime) {
	    nafter_veryshort++;
	  }
	}
	// Avoid too many unnecessary loop runs.
	if ((nafter_short>0) || (nafter_long>0)) break;
      }
      CHECK_STATUS_BREAK(status);

      Pattern* pattern=getPattern(&status);
      CHECK_STATUS_BREAK(status);

      // Copy the pattern data.
      pattern->time=event.time;
      pattern->frame=event.frame;
      pattern->signal=event.signal;
      pattern->pi=event.pi;
      pattern->rawx=event.rawx;
      pattern->rawy=event.rawy;
      pattern->ra=0.;
      pattern->dec=0.;
      pattern->npixels=1;
      pattern->pileup=0;
      
      int ii;
      for (ii=0; (ii<NEVENTPHOTONS)&&(ii<NPATTERNPHOTONS); ii++){
	pattern->ph_id[ii] =event.ph_id[ii];
	pattern->src_id[ii]=event.src_id[ii];
	if ((ii>0)&&(pattern->ph_id[ii]!=0)) {
	  pattern->pileup=1;
	}
      }
    
      pattern->signals[0]=0.;
      pattern->signals[1]=0.;
      pattern->signals[2]=0.;
      pattern->signals[3]=0.;
      pattern->signals[4]=event.signal;
      pattern->signals[5]=0.;
      pattern->signals[6]=0.;
      pattern->signals[7]=0.;
      pattern->signals[8]=0.;

      // Determine the event grade.
      if ((nbefore_veryshort>0)||(nafter_veryshort>0)) {
	pattern->type = 3;
      } else if ((nbefore_short>0)||(nafter_short>0)) {
	pattern->type = 2;
      } else if ((nbefore_short==0) && (nafter_long==0)) {
	pattern->type = 0;
      } else {
	pattern->type = 1;
      } 

      switch (pattern->type) {
      case 0: ngrade0++; break;
      case 1: ngrade1++; break;
      case 2: ngrade2++; break;
      case 3: ngrade3++; break;
      }
      
      // Write the data to the output file.
      addPattern2File(plf, pattern, &status);	  
      CHECK_STATUS_BREAK(status);

      // Release memory.
      freePattern(&pattern);

    } // End of loop over all events in the event file
    CHECK_STATUS_BREAK(status);

    // Write header keywords.
    fits_write_key(plf->fptr, TLONG, "NPHOTONS", &nphotons, "", &status);
    fits_write_key(plf->fptr, TLONG, "NGRADE0", &ngrade0, "", &status);
    fits_write_key(plf->fptr, TLONG, "NGRADE1", &ngrade1, "", &status);
    fits_write_key(plf->fptr, TLONG, "NGRADE2", &ngrade2, "", &status);
    fits_write_key(plf->fptr, TLONG, "NGRADE3", &ngrade3, "", &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // End of error handling loop

  // --- Clean Up ---

  // Close the files.
  freeEventListFile(&elf, &status);
  destroyPatternFile(&plf, &status);

  return(status);
}


int analyse_xms_events_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  status=ape_trad_query_file_name("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the event list!\n", status);
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);


  status=ape_trad_query_string("PatternList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the pattern list!\n", status);
    return(status);
  } 
  strcpy(par->PatternList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_double("TimeUnit", &par->TimeUnit);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the time unit!\n", status);
    return(status);
  }

  status=ape_trad_query_int("PreTrigger", &par->PreTrigger);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the pre-trigger!\n", status);
    return(status);
  }

  status=ape_trad_query_int("PostTrigger", &par->PostTrigger);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the post-trigger!\n", status);
    return(status);
  }

  status=ape_trad_query_double("PileupTime", &par->PileupTime);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the pile-up time!\n", status);
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}





