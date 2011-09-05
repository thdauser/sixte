#include "simputspec.h"


int simputspec_main() 
{
  // Program parameters.
  struct Parameters par;

  // Temporary file for ISIS interaction.
  FILE* tmpfile=NULL;

  // Simput data structures (used as buffers).
  SimputMissionIndepSpec* spec=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("simputspec");
  set_toolversion("0.01");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----
    
    // Read the parameters using PIL.
    status=simputspec_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // ---- END of Initialization ----


    // ---- Main Part ----

    // Run the command.
    system("isis");
    // TODO Process return value.

    // ---- END of Main Part ----

  } while(0); // END of error handling loop.

  if (status==EXIT_SUCCESS) headas_chat(0, "finished successfully!\n\n");
  return(status);
}



int simputspec_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_float("plIndex", &par->plIndex);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the plIndex parameter failed");
    return(status);
  }

  status=ape_trad_query_float("plNorm", &par->plNorm);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the plNorm parameter failed");
    return(status);
  }

  status=ape_trad_query_float("bbkT", &par->bbkT);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the bbkT parameter failed");
    return(status);
  }

  status=ape_trad_query_float("bbNorm", &par->bbNorm);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the bbNorm parameter failed");
    return(status);
  }

  status=ape_trad_query_float("flSigma", &par->flSigma);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the flSigma parameter failed");
    return(status);
  }

  status=ape_trad_query_float("flNorm", &par->flNorm);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the flNorm parameter failed");
    return(status);
  }

  status=ape_trad_query_float("rflSpin", &par->rflSpin);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the rflSpin parameter failed");
    return(status);
  }

  status=ape_trad_query_float("rflNorm", &par->rflNorm);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the rflNorm parameter failed");
    return(status);
  }

  status=ape_trad_query_file_name("Outfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the name of the output file failed");
    return(status);
  } 
  strcpy(par->Outfile, sbuffer);
  free(sbuffer);


  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }

  return(status);
}


