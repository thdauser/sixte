#include "attitudecatalog.h"


AttitudeCatalog* getAttitudeCatalog(int* const status)
{
  AttitudeCatalog* ac=(AttitudeCatalog*)malloc(sizeof(AttitudeCatalog));
  if (NULL==ac) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("memory allocation for AttitudeCatalog "
		   "failed!\n", *status);
    return(NULL);
  }
  
  // Initialize.
  ac->nentries = 0;
  ac->current_entry = 0;
  ac->entry = NULL;

  return(ac);
}


AttitudeCatalog* loadAttitudeCatalog(const char* filename, 
				     const double t0, 
				     const double timespan,
				     int* const status)
{
  AttitudeCatalog* ac=NULL;
  AttitudeFile* af=NULL;
  char msg[MAXMSG];

  do { // Beginning of ERROR handling loop

    // Read-in the attitude data from the FITS file and store 
    // them in the AttitudeCatalog.

    // Open the attitude file:
    headas_chat(5, "open attitude catalog file '%s' ...\n",
		filename);
    af = open_AttitudeFile(filename, READONLY, status);
    if (NULL==af) break;
    headas_chat(5, " attitude catalog contains %ld entries\n", 
		af->nrows);

    // Get a new AttitudeCatalog object.
    ac = getAttitudeCatalog(status);
    CHECK_STATUS_BREAK(*status);

    // Allocate memory for the entries in the catalog.
    ac->entry=(AttitudeEntry*)malloc(af->nrows*sizeof(AttitudeEntry));
    if (NULL==ac->entry) {
      *status = EXIT_FAILURE;
      HD_ERROR_THROW("not enough memory available to store "
		     "the AttitudeCatalog!\n", *status);
      break;
    }
    
    // Read all lines from attitude file subsequently.
    AttitudeFileEntry afe;
    for (af->row=0; af->row<af->nrows; af->row++) {
      
      afe = read_AttitudeFileEntry(af, status);
      if (EXIT_SUCCESS!=*status) break;

      // Calculate and store attitude data:
      ac->entry[af->row].time = afe.time;
      // Rescale from degrees to radians:
      afe.rollang = afe.rollang*M_PI/180.;

      // Telescope pointing direction nz:
      ac->entry[af->row].nz = 
	unit_vector(afe.viewra*M_PI/180., afe.viewdec*M_PI/180.);
      // Roll-Angle:
      ac->entry[af->row].roll_angle = afe.rollang;
    }
    if (EXIT_SUCCESS!=*status) break;
    // End of the attitude readout loop

    // Save the number of AttitudeEntry elements.
    ac->nentries = af->nrows;

    // Check if the required time interval for the simulation
    // is a subset of the time described by the attitude file.
    if (t0>0.) {
      if ((ac->entry[0].time > t0) || 
	  (ac->entry[ac->nentries-1].time < t0+timespan)) {
	*status=EXIT_FAILURE;
      sprintf(msg, "not enough attitude data available for the "
	      "specified period from %lf to %lf!", t0, t0+timespan);
      HD_ERROR_THROW(msg, *status);
      break;
      }
    }

    // Loop over all AttitudeEntry elements in the AttitudeCatalog in
    // order to determine the telescope nx-direction (so far only the
    // nz direction is set).
    // Check the change of the telescope pointing direction between two subsequent
    // AttitudeEntry elements.
    Vector dnz = vector_difference(ac->entry[1].nz, ac->entry[0].nz);
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

    long ii;
    for (ii=1; ii<ac->nentries; ii++) {

      Vector dnz = 
	vector_difference(ac->entry[ii].nz, ac->entry[ii-1].nz);
      if (sqrt(scalar_product(&dnz, &dnz))<1.e-7) { // 1.e-3
	// Change of the telescope axis is too small to be significant.
	Vector ny = {0., 1., 0.};
	ac->entry[ii].nx=vector_product(ac->entry[ii].nz, ny); // TODO
      } else {
	ac->entry[ii].nx=
	  normalize_vector(vector_difference(ac->entry[ii].nz,
					     ac->entry[ii-1].nz));
      }

    } // END of loop over all AttitudeEntry elements for the calculation of nx.

  } while (0); // End of error handling loop


  // --- Clean up ---

  // Close FITS file
  if (NULL!=af) {
    if (af->fptr) fits_close_file(af->fptr, status);
    free(af);
  }

  if (EXIT_SUCCESS != *status) ac = NULL;
  return(ac);
}


void freeAttitudeCatalog(AttitudeCatalog** const ac)
{
  if (NULL != (*ac)) {
    if (NULL != (*ac)->entry) {
      free((*ac)->entry);
    }
    free(*ac);
    *ac=NULL;
  }
}



Vector getTelescopePointing(AttitudeCatalog* const ac, 
			    const double time, 
			    int* const status)
{
  Vector nz = { .x = 0., .y = 0., .z = 0. };
  char msg[MAXMSG]; // Error message buffer.

  // Check if the requested time lies within the current time bin.
  while (time < ac->entry[ac->current_entry].time) {
    // Check if the beginning of the AttitudeCatalog is reached.
    if (ac->current_entry <= 0) {
      *status = EXIT_FAILURE;
      sprintf(msg, "no attitude entry available for time %lf!\n", time);
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
      sprintf(msg, "no attitude entry available for time %lf!\n", time);
      HD_ERROR_THROW(msg, *status);
      return(nz);
    }
    // If not, go one step further.
    ac->current_entry++;
  }
    
  // The requested time lies within the current time bin.
  // Interpolation:
  nz = interpolateCircleVector(ac->entry[ac->current_entry].nz, 
			       ac->entry[ac->current_entry+1].nz, 
			       (time-ac->entry[ac->current_entry].time)/
			       (ac->entry[ac->current_entry+1].time-
				ac->entry[ac->current_entry].time));

  return(nz);
}



double getRollAngle(AttitudeCatalog* ac, double time, int* status)
{
  char msg[MAXMSG]; // Error message buffer.

  // Check if the requested time lies within the current time bin.
  while (time < ac->entry[ac->current_entry].time) {
    // Check if the beginning of the AttitudeCatalog is reached.
    if (ac->current_entry <= 0) {
      *status = EXIT_FAILURE;
      sprintf(msg, "no attitude entry available for time %lf!\n", time);
      HD_ERROR_THROW(msg, *status);
      return(0.);
    }
    // If not, go one step back.
    ac->current_entry--;
  }

  while (time > ac->entry[ac->current_entry+1].time) {
    // Check if the end of the AttitudeCatalog is reached.
    if (ac->current_entry >= ac->nentries-2) {
      *status = EXIT_FAILURE;
      sprintf(msg, "no attitude entry available for time %lf!\n", time);
      HD_ERROR_THROW(msg, *status);
      return(0.);
    }
    // If not, go one step further.
    ac->current_entry++;
  }
    
  // The requested time lies within the current time bin.
  // Interpolation:
  double fraction = 
    (time-ac->entry[ac->current_entry].time)/
    (ac->entry[ac->current_entry+1].time-ac->entry[ac->current_entry].time);
  return(ac->entry[ac->current_entry  ].roll_angle*(1.-fraction) + 
	 ac->entry[ac->current_entry+1].roll_angle*    fraction );
}


AttitudeEntry defaultAttitudeEntry()
{
  AttitudeEntry ae;

  ae.time = 0.;

  ae.nz.x = 0.;
  ae.nz.y = 0.;
  ae.nz.z = 0.;

  ae.nx.x = 0.;
  ae.nx.y = 0.;
  ae.nx.z = 0.;

  ae.roll_angle = 0.;

  return(ae);
}


