#include "tes_grades.h"


int tes_grades_main() {
  struct Parameters par;
  PatternFile* plf=NULL;

  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("tes_grades");
  set_toolversion("0.06");


  do { // ERROR handling loop

    // Grade counters.
    long ngrade0=0, ngrade1=0, ngrade2=0, ngrade3=0;

    // Read parameters by PIL:
    status=tes_grades_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Set the input pattern file.
    plf=openPatternFile(par.PatternList, READWRITE, &status);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "analyse events ...\n");

    // Loop over all entries in the pattern file.
    long row1;
    for (row1=0; row1<plf->nrows; row1++) {
      
      // Read the time of the pattern from the file.
      Pattern pat1;
      getPatternFromFile(plf, row1+1, &pat1, &status);
      CHECK_STATUS_BREAK(status);

      // Check the events before and after the current one 
      // within the specified time spans.
      int nbefore_short=0, nbefore_long=0, nbefore_veryshort=0;
      int nafter_short=0, nafter_long=0, nafter_veryshort=0;

      // Former events:
      long row2;
      for (row2=row1-1; row2>=0; row2--) {
	Pattern pat2; // Buffer.
	getPatternFromFile(plf, row2+1, &pat2, &status);
	CHECK_STATUS_BREAK(status);
	
	if (pat1.time-pat2.time > par.PostTrigger*par.TimeUnit) break;
	if ((pat1.rawx==pat2.rawx)&&(pat1.rawy==pat2.rawy)) {
	  nbefore_long++;
	  if (pat1.time-pat2.time < par.PreTrigger*par.TimeUnit) {
	    nbefore_short++;
	  }	
	  if (pat1.time-pat2.time < par.PileupTime) {
	    nbefore_veryshort++;
	  }	
	}
	// Avoid too many unnecessary loop runs.
	if ((nbefore_short>0) || (nbefore_long>0)) break;
      }
      CHECK_STATUS_BREAK(status);

      // Subsequent events:
      for (row2=row1+1; row2<plf->nrows; row2++) {
	Pattern pat2; // Buffer.
	getPatternFromFile(plf, row2+1, &pat2, &status);
	CHECK_STATUS_BREAK(status);

	if (pat2.time-pat1.time > par.PostTrigger*par.TimeUnit) break;
	if ((pat1.rawx==pat2.rawx)&&(pat1.rawy==pat2.rawy)) {
	  nafter_long++;
	  if (pat2.time-pat1.time < par.PreTrigger*par.TimeUnit) {
	    nafter_short++;
	  }
	  if (pat2.time-pat1.time < par.PileupTime) {
	    nafter_veryshort++;
	  }
	}
	// Avoid too many unnecessary loop runs.
	if ((nafter_short>0) || (nafter_long>0)) break;
      }
      CHECK_STATUS_BREAK(status);

      // Determine the event grade.
      if ((nbefore_veryshort>0)||(nafter_veryshort>0)) {
	pat1.type=3;
      } else if ((nbefore_short>0)||(nafter_short>0)) {
	pat1.type=2;
      } else if ((nbefore_short==0)&&(nafter_long==0)) {
	pat1.type=0;
      } else {
	pat1.type=1;
      } 

      switch (pat1.type) {
      case 0: ngrade0++; break;
      case 1: ngrade1++; break;
      case 2: ngrade2++; break;
      case 3: ngrade3++; break;
      }
      
      // Update the pattern information in the file.
      fits_write_col(plf->fptr, TINT, plf->ctype, row1+1, 
		     1, 1, &pat1.type, &status);
      CHECK_STATUS_BREAK(status);

    }
    CHECK_STATUS_BREAK(status);
    // End of loop over all events in the event file

    // Write header keywords.
    fits_update_key(plf->fptr, TLONG, "NGRADE0", &ngrade0, "", &status);
    fits_update_key(plf->fptr, TLONG, "NGRADE1", &ngrade1, "", &status);
    fits_update_key(plf->fptr, TLONG, "NGRADE2", &ngrade2, "", &status);
    fits_update_key(plf->fptr, TLONG, "NGRADE3", &ngrade3, "", &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // End of error handling loop

  // --- Clean Up ---

  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  destroyPatternFile(&plf, &status);

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int tes_grades_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  status=ape_trad_query_file_name("PatternList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the pattern list");
    return(status);
  } 
  strcpy(par->PatternList, sbuffer);
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

