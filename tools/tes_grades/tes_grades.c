#include "tes_grades.h"


int tes_grades_main() {
  struct Parameters par;
  EventFile* elf=NULL;

  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("tes_grades");
  set_toolversion("0.07");


  do { // ERROR handling loop

    // Grade counters.
    long ngrade0=0, ngrade1=0, ngrade2=0, ngrade3=0;

    // Read parameters by PIL:
    status=tes_grades_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Open the event file.
    elf=openEventFile(par.EventList, READWRITE, &status);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "analyse events ...\n");

    // Loop over all entries in the event file.
    long row1;
    for (row1=0; row1<elf->nrows; row1++) {
      
      // Read the time of the event from the file.
      Event ev1;
      getEventFromFile(elf, row1+1, &ev1, &status);
      CHECK_STATUS_BREAK(status);

      // Check the events before and after the current one 
      // within the specified time spans.
      int nbefore_short=0, nbefore_long=0, nbefore_veryshort=0;
      int nafter_short=0, nafter_long=0, nafter_veryshort=0;

      // Former events:
      long row2;
      for (row2=row1-1; row2>=0; row2--) {
	Event ev2; // Buffer.
	getEventFromFile(elf, row2+1, &ev2, &status);
	CHECK_STATUS_BREAK(status);
	
	if (ev1.time-ev2.time > par.PostTrigger*par.TimeUnit) break;
	if ((ev1.rawx==ev2.rawx)&&(ev1.rawy==ev2.rawy)) {
	  nbefore_long++;
	  if (ev1.time-ev2.time < par.PreTrigger*par.TimeUnit) {
	    nbefore_short++;
	  }	
	  if (ev1.time-ev2.time < par.PileupTime) {
	    nbefore_veryshort++;
	  }	
	}
	// Avoid too many unnecessary loop runs.
	if ((nbefore_short>0) || (nbefore_long>0)) break;
      }
      CHECK_STATUS_BREAK(status);

      // Subsequent events:
      for (row2=row1+1; row2<elf->nrows; row2++) {
	Event ev2; // Buffer.
	getEventFromFile(elf, row2+1, &ev2, &status);
	CHECK_STATUS_BREAK(status);

	if (ev2.time-ev1.time > par.PostTrigger*par.TimeUnit) break;
	if ((ev1.rawx==ev2.rawx)&&(ev1.rawy==ev2.rawy)) {
	  nafter_long++;
	  if (ev2.time-ev1.time < par.PreTrigger*par.TimeUnit) {
	    nafter_short++;
	  }
	  if (ev2.time-ev1.time < par.PileupTime) {
	    nafter_veryshort++;
	  }
	}
	// Avoid too many unnecessary loop runs.
	if ((nafter_short>0) || (nafter_long>0)) break;
      }
      CHECK_STATUS_BREAK(status);

      // Determine the event grade.
      if ((nbefore_veryshort>0)||(nafter_veryshort>0)) {
	ev1.type=3;
      } else if ((nbefore_short>0)||(nafter_short>0)) {
	ev1.type=2;
      } else if ((nbefore_short==0)&&(nafter_long==0)) {
	ev1.type=0;
      } else {
	ev1.type=1;
      } 

      switch (ev1.type) {
      case 0: ngrade0++; break;
      case 1: ngrade1++; break;
      case 2: ngrade2++; break;
      case 3: ngrade3++; break;
      }
      
      // Update the event information in the file.
      fits_write_col(elf->fptr, TINT, elf->ctype, row1+1, 
		     1, 1, &ev1.type, &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);
    // End of loop over all events in the event file

    // Write header keywords.
    fits_update_key(elf->fptr, TLONG, "NGRADE0", &ngrade0, "", &status);
    fits_update_key(elf->fptr, TLONG, "NGRADE1", &ngrade1, "", &status);
    fits_update_key(elf->fptr, TLONG, "NGRADE2", &ngrade2, "", &status);
    fits_update_key(elf->fptr, TLONG, "NGRADE3", &ngrade3, "", &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // End of error handling loop

  // --- Clean Up ---

  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  freeEventFile(&elf, &status);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int tes_grades_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  status=ape_trad_query_file_name("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event list");
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);


  status=ape_trad_query_double("TimeUnit", &par->TimeUnit);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the time unit");
    return(status);
  }

  status=ape_trad_query_int("PreTrigger", &par->PreTrigger);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the pre-trigger");
    return(status);
  }

  status=ape_trad_query_int("PostTrigger", &par->PostTrigger);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the post-trigger");
    return(status);
  }

  status=ape_trad_query_double("PileupTime", &par->PileupTime);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the pile-up time");
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}

