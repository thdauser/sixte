#include "attitudecatalog.h"


AttitudeCatalog* get_AttitudeCatalog(const char* filename, 
				     double t0, 
				     double timespan,
				     int* status)
{
  AttitudeCatalog* ac=NULL;
  AttitudeFile* af=NULL;
  long entry=0; // Counter for the AttitudeEntry elements in the AttitudeCatalog.
  char msg[MAXMSG];

  do {  // beginning of ERROR handling loop

    // Read-in the attitude data from the FITS file and store 
    // them in the AttitudeCatalog.

    // Open the attitude file:
    headas_chat(5, "open attitude catalog file '%s' ...\n",
		filename);
    af = open_AttitudeFile(filename, READONLY, status);
    if (NULL==af) break;
    headas_chat(5, " attitude catalog contains %ld entries\n", 
		af->nrows);

    // Get memory for the AttitudeCatalog:
    ac = (AttitudeCatalog*)malloc(sizeof(AttitudeCatalog));
    if (NULL==ac) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory available to store the AttitudeCatalog!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }
    ac->nentries = 0;
    ac->current_entry = 0;
    ac->entry = (AttitudeEntry*)malloc(af->nrows*sizeof(AttitudeEntry));
    if (NULL==ac->entry) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory available to store the AttitudeCatalog!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }
    
    
    // Read all lines from attitude file subsequently.
    AttitudeFileEntry afe;
    for (af->row=0, entry=0; EXIT_SUCCESS==(*status); af->row++) {
      
      afe = read_AttitudeFileEntry(af, status);
      if (EXIT_SUCCESS!=*status) break;

      if (afe.time >= t0) {
	// Check if we are already at the end of the attitude file:
	if (af->row>=af->nrows) {
	  *status=EXIT_FAILURE;
	  sprintf(msg, "Error: Not enough attitude data available for the "
		  "specified period from %lf to %lf!", t0, timespan);
	  HD_ERROR_THROW(msg, *status);
	  break;
	}
	  
	// Calculate and store attitude data:
	ac->entry[entry].time = afe.time;
	// Rescale from degrees to radians:
	afe.rollang = afe.rollang*M_PI/180.;

	// Telescope pointing direction nz:
	ac->entry[entry].nz = 
	  unit_vector(afe.viewra*M_PI/180., afe.viewdec*M_PI/180.);

	entry++;
      }
	  
      // Check whether the end of the reading process 
      // (end of specified period) is reached.
      if (afe.time >= t0+timespan) {
	break;
      }
    }  // End of the attitude readout loop
    
    // Save the number of AttitudeEntry elements.
    ac->nentries = entry;


    // Loop over all AttitudeEntry elements in the AttitudeCatalog in
    // order to determine the telescope nx-direction (so far only the
    // nz direction is set).
    // Check the change of the telescope pointing direction between two subsequent
    // AttitudeEntry elements.
    Vector dnz = 
      vector_difference(ac->entry[1].nz, ac->entry[0].nz);
    if (sqrt(scalar_product(&dnz, &dnz))<1.e-7) { // 1.e-3
      // Change of the telescope axis is too small to be significant.
      Vector ny = {0., 1., 0.};
      ac->entry[0].nx=vector_product(ac->entry[0].nz, ny); // TODO
    } else {
      // nx = (nz_0 x nz_1) x nz_0
      ac->entry[0].nx=
	normalize_vector(vector_product(vector_product(ac->entry[0].nz,
						       ac->entry[1].nz),
					ac->entry[0].nz));
    }

    for (entry=1; entry<ac->nentries; entry++) {

      Vector dnz = 
	vector_difference(ac->entry[entry].nz, ac->entry[entry-1].nz);
      if (sqrt(scalar_product(&dnz, &dnz))<1.e-7) { // 1.e-3
	// Change of the telescope axis is too small to be significant.
	Vector ny = {0., 1., 0.};
	ac->entry[entry].nx=vector_product(ac->entry[entry].nz, ny); // TODO
      } else {
	printf("here");
	ac->entry[entry].nx=
	  normalize_vector(vector_difference(ac->entry[entry].nz,
					     ac->entry[entry-1].nz));
      }

    } // END of loop over all AttitudeEntry elements for the calculation of nx.

  } while (0); // End of error handling loop


  // --- clean up ---

  // Close FITS file
  if (NULL!=af) {
    if (af->fptr) fits_close_file(af->fptr, status);
    free(af);
  }

  if (EXIT_SUCCESS != *status) ac = NULL;
  return(ac);
}



AttitudeCatalog* getEntireAttitudeCatalog(const char* filename, int* status)
{
  AttitudeCatalog* ac=NULL;
  AttitudeFile* af=NULL;
  long entry=0; // Counter for the AttitudeEntry elements in the AttitudeCatalog.

  do {  // beginning of ERROR handling loop

    // Read-in the attitude data from the FITS file and store 
    // them in the AttitudeCatalog.

    // Open the attitude file:
    af = open_AttitudeFile(filename, READONLY, status);
    if (NULL==af) break;

    // Get memory for the AttitudeCatalog:
    ac = (AttitudeCatalog*)malloc(sizeof(AttitudeCatalog));
    if (NULL==ac) {
      *status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: not enough memory available to store the AttitudeCatalog!\n", 
		     *status);
      break;
    }
    ac->nentries = 0;
    ac->current_entry = 0;
    ac->entry = (AttitudeEntry*)malloc(af->nrows*sizeof(AttitudeEntry));
    if (NULL==ac->entry) {
      *status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: not enough memory available to store the AttitudeCatalog!\n", 
		     *status);
      break;
    }

    
    // Read all lines from attitude file subsequently.
    AttitudeFileEntry afe;
    for (af->row=0, entry=0; (EXIT_SUCCESS==*status)&&(af->row<af->nrows); af->row++) {
      
      afe = read_AttitudeFileEntry(af, status);
      if (EXIT_SUCCESS!=*status) break;

      // Calculate and store attitude data:
      ac->entry[entry].time = afe.time;
      // Rescale from degrees to radians:
      afe.rollang = afe.rollang*M_PI/180.;
      
      // Telescope pointing direction nz:
      ac->entry[entry].nz = unit_vector(afe.viewra*M_PI/180., afe.viewdec*M_PI/180.);

      entry++;
    }  // End of the attitude readout loop
    
    // Save the number of AttitudeEntry elements.
    ac->nentries = entry;


    // Loop over all AttitudeEntry elements in the AttitudeCatalog in
    // order to determine the telescope nx-direction (so far only the
    // nz direction is set).
    // Check the change of the telescope pointing direction between two subsequent
    // AttitudeEntry elements.
    Vector dnz = 
      vector_difference(ac->entry[1].nz, ac->entry[0].nz);
    if (sqrt(scalar_product(&dnz, &dnz))<1.e-7) { 
      // Change of the telescope axis is too small to be significant.
      Vector ny = {0., 1., 0.};
      ac->entry[0].nx=vector_product(ac->entry[0].nz, ny); // TODO
    } else {
      // nx = (nz_0 x nz_1) x nz_0
      ac->entry[0].nx=
	normalize_vector(vector_product(vector_product(ac->entry[0].nz, ac->entry[1].nz),
					ac->entry[0].nz));
    }

    for (entry=1; entry<ac->nentries; entry++) {

      Vector dnz = 
	vector_difference(ac->entry[entry].nz, ac->entry[entry-1].nz);
      if (sqrt(scalar_product(&dnz, &dnz))<1.e-7) { 
	// Change of the telescope axis is too small to be significant.
	Vector ny = {0., 1., 0.};
	ac->entry[entry].nx=vector_product(ac->entry[entry].nz, ny); // TODO
      } else {
	ac->entry[entry].nx=
	  normalize_vector(vector_difference(ac->entry[entry].nz,
					     ac->entry[entry-1].nz));
      }

    } // END of loop over all AttitudeEntry elements for the calculation of nx.

  } while (0); // End of error handling loop


  // --- clean up ---

  // Close FITS file
  if (NULL!=af) {
    if (af->fptr) fits_close_file(af->fptr, status);
    free(af);
  }

  if (EXIT_SUCCESS != *status) ac = NULL;
  return(ac);
}



void free_AttitudeCatalog(AttitudeCatalog* ac)
{
  if (NULL != ac) {
    if (NULL != ac->entry) {
      free(ac->entry);
    }
    free(ac);
  }
}



Vector getTelescopePointing(AttitudeCatalog* ac, double time, int* status)
{
  Vector nz = { .x = 0., .y = 0., .z = 0. };
  char msg[MAXMSG]; // Error message buffer.

  // Check if the requested time lies within the current time bin.
  while (time < ac->entry[ac->current_entry].time) {
    // Check if the beginning of the AttitudeCatalog is reached.
    if (ac->current_entry <= 0) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: no orbit entry available for time %lf!\n", time);
      HD_ERROR_THROW(msg, *status);
      return(nz);
    }
    // If not, go one step back.
    ac->current_entry--;
  }

  while (time > ac->entry[ac->current_entry+1].time) {
    // Check if the end of the AttitudeCatalog is reached.
    if (ac->current_entry >= ac->nentries-2) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: no orbit entry available for time %lf!\n", time);
      HD_ERROR_THROW(msg, *status);
      return(nz);
    }
    // If not, go one step back.
    ac->current_entry++;
  }
    
  // The requested time lies within the current time bin.
  // Interpolation:
  // TODO: replace this calculation by proper attitude interpolation.
  nz = interpolate_vec(ac->entry[ac->current_entry].nz, 
		       ac->entry[ac->current_entry].time, 
		       ac->entry[ac->current_entry].nz, 
		       ac->entry[ac->current_entry].time, 
		       time);
  normalize_vector_fast(&nz);

  return(nz);
}


