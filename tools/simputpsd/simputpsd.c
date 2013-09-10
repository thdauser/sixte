#include "simputpsd.h"


int simputpsd_main() 
{
  // Program parameters.
  struct Parameters par;

  // Input ASCII file containing the PSD.
  FILE* asciipsd=NULL;

  // Output SimputPSD.
  SimputPSD* simputpsd=NULL;

  // SimputCtlg the PSD should be attached to.
  SimputCtlg* cat=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("simputpsd");
  set_toolversion("0.01");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=simputpsd_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Open the ASCII file with the PSD.
    asciipsd=fopen(par.PSDFile,"r");
    CHECK_NULL_BREAK(asciipsd, status, "could not open input PSD file");

    // ---- END of Initialization ----


    // ---- Main Part ----

    // Determine the number of rows.
    long nlines=0;
    char c=0;
    while(!feof(asciipsd)) {
      c=fgetc(asciipsd);
      if ('\n'==c) {
	nlines++;
      }
    }
    // Check if the last line has been empty.
    if('\n'==c) {
      nlines--;
    }

    // Allocate memory.
    simputpsd=newSimputPSD(&status);
    CHECK_STATUS_BREAK(status);
    simputpsd->nentries=nlines;
    simputpsd->frequency=(float*)malloc(nlines*sizeof(float));
    CHECK_NULL_BREAK(simputpsd->frequency, status, "memory allocation failed");
    simputpsd->power=(float*)malloc(nlines*sizeof(float));
    CHECK_NULL_BREAK(simputpsd->power, status, "memory allocation failed");

    // Reset the file pointer, read the data and store them in
    // the SimputPSD data structure.
    rewind(asciipsd);
    long ii;
    for (ii=0; ii<nlines; ii++) {
      fscanf(asciipsd, "%f %f\n",
	     &(simputpsd->frequency[ii]), &(simputpsd->power[ii]));
    }

    // Store the PSD in the SIMPUT file.
    saveSimputPSD(simputpsd, par.Simput, "TIMING", 1, &status);
    CHECK_STATUS_BREAK(status);


    // Open the SimputCtlg.
    cat=openSimputCtlg(par.Simput, READWRITE, 32, 32, 32, 32, &status);
    CHECK_STATUS_BREAK(status);

    // Set the timing reference in the source catalog.
    char* timeref=(char*)malloc(32*sizeof(char));
    CHECK_NULL_BREAK(timeref, status, "memory allocation failed");
    strcpy(timeref, "[TIMING,1]");
    fits_write_col(cat->fptr, TSTRING, cat->ctiming, 1, 1, 1, 
    		   &timeref, &status);
    CHECK_STATUS_BREAK(status);

    // ---- END of Main Part ----

  } while(0); // END of error handling loop.

  // Close open files.
  if (NULL!=asciipsd) {
    fclose(asciipsd);
    asciipsd=NULL;
  }

  // Release memory.
  freeSimputPSD(&simputpsd);
  freeSimputCtlg(&cat, &status);

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int simputpsd_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("Simput", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the name of the SIMPUT catalog failed");
    return(status);
  } 
  strcpy(par->Simput, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("PSDFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the name of the PSD file failed");
    return(status);
  }
  strcpy(par->PSDFile, sbuffer);
  free(sbuffer);

  return(status);
}


