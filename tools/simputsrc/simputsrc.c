#include "simputsrc.h"


int simputsrc_main() 
{
  // Program parameters.
  struct Parameters par;

  // Simput data structure.
  SimputCatalog* cat=NULL;
  SimputSource* src=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("simputsrc");
  set_toolversion("0.01");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----
    
    // Read the parameters using PIL.
    status=simputsrc_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // ---- END of Initialization ----


    // ---- Main Part ----

    // Create a new SIMPUT catalog.
    remove(par.Simput);
    cat=openSimputCatalog(par.Simput, READWRITE, &status);
    CHECK_STATUS_BREAK(status);

    // Insert a point-like source.
    src=getSimputSourceV(1, "", 0., 0., 0., 1., 
			 par.Emin, par.Emax, par.Flux, 
			 par.Spectrum, "", "", &status);
    CHECK_STATUS_BREAK(status);
    appendSimputSource(cat, src, &status);
    CHECK_STATUS_BREAK(status);

    // ---- END of Main Part ----

  } while(0); // END of error handling loop.


  // Release memory.
  freeSimputSource(&src);
  freeSimputCatalog(&cat, &status);

  if (status==EXIT_SUCCESS) headas_chat(0, "finished successfully!\n\n");
  return(status);
}



int simputsrc_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_float("Flux", &par->Flux);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the Flux parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Emin", &par->Emin);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the Emin parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Emax", &par->Emax);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the Emax parameter failed");
    return(status);
  }

  status=ape_trad_query_file_name("Spectrum", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the name of the spectrum file failed");
    return(status);
  } 
  strcpy(par->Spectrum, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Simput", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the name of the output SIMPUT catalog file failed");
    return(status);
  } 
  strcpy(par->Simput, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }

  return(status);
}


