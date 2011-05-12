#include "analyse_xms_events.h"


int analyse_xms_events_main() {
  struct Parameters par;
  EventListFile* elf=NULL;
  GenPatternFile* plf=NULL;

  int status = EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("analyse_xms_events");
  set_toolversion("0.02");

  do { // ERROR handling loop

    // Event grade counters.
    long nphotons=0, ngrade0=0, ngrade1=0, ngrade2=0, ngrade3=0;

    // Read parameters by PIL:
    status = analyse_xms_events_getpar(&par);
    if (EXIT_SUCCESS!=status) break;

    // Set the input event file.
    elf=openEventListFile(par.EventList, READWRITE, &status);
    if (EXIT_SUCCESS!=status) break;

    // Create and open the output event file.
    // Filename of the template file.
    char template[MAXMSG];
    strcpy(template, par.data_path);
    strcat(template, "/templates/patternlist.tpl");
    // Open a new pattern file from the specified template.
    plf=openNewGenPatternFile(par.PatternList, template, elf->mjdref, &status);
    if (EXIT_SUCCESS!=status) break;


    // Loop over all events in the event file.
    long row;
    for (row=0; row<elf->nrows; row++) {
      
      // Read a new event from the file.
      Event event;
      getEventFromFile(elf, row+1, &event, &status);
      CHECK_STATUS_BREAK(status);

      // Count the total number of photon events in the event list.
      nphotons++;

      // Neglect the inner pixel(s).
      //if ((event.rawx>=14)&&(event.rawx<=16)&&
      //(event.rawy>=14)&&(event.rawy<=16)) continue;

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

      GenPattern pattern = {
	.pat_type= 0,
	.pileup  = 0,
	.event   = event
      };

      // Determine the event grade.
      if ((nbefore_veryshort>0) || (nafter_veryshort>0)) {
	pattern.pat_type = 3;
      } else if ((nbefore_short>0)||(nafter_short>0)) {
	pattern.pat_type = 2;
      } else if ((nbefore_short==0) && (nafter_long==0)) {
	pattern.pat_type = 0;
      } else {
	pattern.pat_type = 1;
      } 

      switch (pattern.pat_type) {
      case 0: ngrade0++; break;
      case 1: ngrade1++; break;
      case 2: ngrade2++; break;
      case 3: ngrade3++; break;
      }
      
      // Write the data to the output file.
      addGenPattern2File(plf, &pattern, &status);	  
      CHECK_STATUS_BREAK(status);

    } // End of loop over all events in the event file
    CHECK_STATUS_BREAK(status);

    // Write header keywords.
    fits_write_key(plf->eventlistfile->fptr, TLONG, "NPHOTONS", &nphotons, "", &status);
    fits_write_key(plf->eventlistfile->fptr, TLONG, "NGRADE0", &ngrade0, "", &status);
    fits_write_key(plf->eventlistfile->fptr, TLONG, "NGRADE1", &ngrade1, "", &status);
    fits_write_key(plf->eventlistfile->fptr, TLONG, "NGRADE2", &ngrade2, "", &status);
    fits_write_key(plf->eventlistfile->fptr, TLONG, "NGRADE3", &ngrade3, "", &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // End of error handling loop

  // --- Clean Up ---

  // Close the files.
  freeEventListFile(&elf, &status);
  destroyGenPatternFile(&plf, &status);

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

  // Get the name of the FITS directory containing the data
  // required for the simulation from the environment variable.
  if (NULL!=(sbuffer=getenv("SIXT_DATA_PATH"))) {
    strcpy(par->data_path, sbuffer);
    // Note: the char* pointer returned by getenv should not
    // be modified nor free'd.
  } else {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error reading the environment variable 'SIXT_DATA_PATH'!\n", 
		   status);
    return(status);
  }

  return(status);
}





