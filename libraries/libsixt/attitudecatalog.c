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
  ac->entry    = NULL;

  return(ac);
}


AttitudeCatalog* loadAttitudeCatalog(const char* filename, 
				     int* const status)
{
  AttitudeCatalog* ac=NULL;
  AttitudeFile* af=NULL;

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
      SIXT_ERROR("not enough memory available to store "
		 "the AttitudeCatalog\n");
      break;
    }

    // Determine whether the rollangle alignment should refer
    // to the telescope's direction of motion or to the equatorial 
    // plane.
    int status2=EXIT_SUCCESS;
    char comment[MAXMSG], sbuffer[MAXMSG]={""}; // String buffers.
    fits_write_errmark();
    fits_read_key(af->fptr, TSTRING, "ALIGNMEN", &sbuffer, comment, &status2);
    fits_clear_errmark();
    // Check the value of the header keyword and set the alignment flag
    // appropriately.
    strtoupper(sbuffer);
    if ((0==strlen(sbuffer)) || (!strcmp(sbuffer, "EQUATOR"))) {
      ac->alignment=0;
    } else if (!strcmp(sbuffer, "MOTION")) {
      ac->alignment=1;
    } else {
      *status=EXIT_FAILURE;
      SIXT_ERROR("invalid value for telescope alignment flag");
      break;
    } 

    
    // Read all lines from attitude file subsequently.
    AttitudeFileEntry afe;
    for (af->row=0; af->row<af->nrows; af->row++) {
      
      afe = read_AttitudeFileEntry(af, status);
      if (EXIT_SUCCESS!=*status) break;

      // Calculate and store attitude data:
      ac->entry[af->row].time = afe.time;

      // Telescope pointing direction nz:
      ac->entry[af->row].nz = 
	unit_vector(afe.viewra*M_PI/180., afe.viewdec*M_PI/180.);

      // Roll-Angle:
      ac->entry[af->row].roll_angle = afe.rollang*M_PI/180.;
    }
    if (EXIT_SUCCESS!=*status) break;
    // End of the attitude readout loop

    // Save the number of AttitudeEntry elements.
    ac->nentries = af->nrows;


    // Determine the telescope nx-direction for all entries in the
    // AttitudeCatalog (so far only the nz direction is set).

    // TODO Check whether nx should explicitly be aligned parallel to 
    // the equatorial plane or to the direction of the telescope motion.

    // Special case for the first entry in the attitude catalog.
    // Determine the change of the telescope pointing direction 
    // between two subsequent AttitudeEntry elements.
    Vector dnz = vector_difference(ac->entry[1].nz, ac->entry[0].nz);
    if (sqrt(scalar_product(&dnz, &dnz))>1.e-7) {
      // The vector nx is pointing along the direction of the 
      // telescope motion.
      // nx = (nz_0 x nz_1) x nz_0
      ac->entry[0].nx=
	normalize_vector(vector_product(vector_product(ac->entry[0].nz,
						       ac->entry[1].nz),
					ac->entry[0].nz));
    } else {
      // Change of the telescope axis is too small to be significant.
      // Therefore the nx vector is selected to be aligned parallel 
      // to the equatorial plane.
      Vector ny = {0., 1., 0.};
      ac->entry[0].nx=vector_product(ac->entry[0].nz, ny); 
      // TODO If nz is pointing towards of one of the poles.
    }

    // Loop over all other (than the first) AttitudeEntry elements 
    // in the AttitudeCatalog.
    long ii;
    for (ii=1; ii<ac->nentries; ii++) {

      Vector dnz = 
	vector_difference(ac->entry[ii].nz, ac->entry[ii-1].nz);
      if (sqrt(scalar_product(&dnz, &dnz))>1.e-7) {
	// The vector nx is pointing along the direction of the 
	// telescope motion.
	ac->entry[ii].nx=
	  normalize_vector(vector_difference(ac->entry[ii].nz,
					     ac->entry[ii-1].nz));
      } else {
	// Change of the telescope axis is too small to be significant.
	// Therefore the nx vector is selected to be aligned parallel 
	// to the equatorial plane.
	Vector ny = {0., 1., 0.};
	ac->entry[ii].nx=vector_product(ac->entry[ii].nz, ny); // TODO
      }

    } 
    // END of loop over all AttitudeEntry elements for the calculation of nx.

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


static void setAttitudeCatalogCurrentEntry(AttitudeCatalog* const ac,
					   const double time,
					   int* const status)
{
  // Check if the requested time lies within the current time bin.
  while (time < ac->entry[ac->current_entry].time) {
    // Check if the beginning of the AttitudeCatalog is reached.
    if (ac->current_entry <= 0) {
      *status = EXIT_FAILURE;
      char msg[MAXMSG]; 
      sprintf(msg, "no attitude entry available for time %lf!\n", time);
      HD_ERROR_THROW(msg, *status);
      return;
    }
    // If not, go one step back.
    ac->current_entry--;
  }

  while (time > ac->entry[ac->current_entry+1].time) {
    // Check if the end of the AttitudeCatalog is reached.
    if (ac->current_entry >= ac->nentries-2) {
      *status = EXIT_FAILURE;
      char msg[MAXMSG]; 
      sprintf(msg, "no attitude entry available for time %lf!\n", time);
      HD_ERROR_THROW(msg, *status);
      return;
    }
    // If not, go one step further.
    ac->current_entry++;
  }
}


Vector getTelescopeNz(AttitudeCatalog* const ac, 
		      const double time, 
		      int* const status)
{
  Vector nz = {.x=0., .y=0., .z=0.};
 
  // Find the appropriate entry in the AttitudeCatalog for the 
  // requested time.
  setAttitudeCatalogCurrentEntry(ac, time, status);
  CHECK_STATUS_RET(*status,nz);
   
  // The requested time lies within the current time bin.
  // Interpolation:
  nz=interpolateCircleVector(ac->entry[ac->current_entry].nz, 
			     ac->entry[ac->current_entry+1].nz, 
			     (time-ac->entry[ac->current_entry].time)/
			     (ac->entry[ac->current_entry+1].time-
			      ac->entry[ac->current_entry].time));
  
  return(nz);
}


void getTelescopeAxes(AttitudeCatalog* const ac,
		      Vector* const nx,
		      Vector* const ny,
		      Vector* const nz,
		      const double time, 
		      int* const status)
{   
  // Determine the z vector (telescope pointing direction):
  *nz=getTelescopeNz(ac, time, status);
  CHECK_STATUS_VOID(*status);
  // After calling getTelescopeNz() the internal current entry pointer
  // within the attitude catalog is set to the appropriate time bin.

  // Check if this is a pointed observation.
  int pointed=0;
  Vector dnz;
  if (ac->nentries==1) {
    // There is only one entry in the AttitudeCatalog.
    pointed=1;
  } else {
    if (ac->current_entry>0) {
      dnz = vector_difference(ac->entry[ac->current_entry].nz, 
			      ac->entry[ac->current_entry-1].nz);
    } else {
      dnz = vector_difference(ac->entry[ac->current_entry+1].nz, 
			      ac->entry[ac->current_entry].nz);
    }
    if (scalar_product(&dnz, &dnz)<1.e-10) {
      pointed=1;
    }
  }

  // Check if the x1 vector should be aligned parallel to the equatorial 
  // plane or along the direction of motion of the telescope axis 
  // (neglecting the roll angle). For a pointed observation the alignment 
  // is done parallel to the equatorial plane.
  Vector x1;
  if ((1==pointed) || (0==ac->alignment)) {
    // Alignment parallel to the equatorial plane.
    Vector b = {0., 0., 1.};
    x1=vector_product(*nz, b);
    
    // CHECK if the telescope is pointing towards one of the poles.
    if (scalar_product(&x1, &x1)<1.e-10) {
      Vector c = {1., 0., 0.};
      x1 = vector_product(*nz, c);
    }

    // Normalize the vector in order to obtain a unit vector.
    x1=normalize_vector(x1);

  } else {
    // Alignment along the direction of motion of the telescope axis.
    x1=normalize_vector(dnz);
  }

  // Determine the y1 vector, which is perpendicular 
  // to the other 2 axes nz and x1:
  Vector y1=normalize_vector(vector_product(*nz, x1));

  // Take into account the roll angle.
  float roll_angle = getRollAngle(ac, time, status);
  CHECK_STATUS_VOID(*status);
  nx->x = x1.x * cos(roll_angle) + y1.x * sin(roll_angle);
  nx->y = x1.y * cos(roll_angle) + y1.y * sin(roll_angle);
  nx->z = x1.z * cos(roll_angle) + y1.z * sin(roll_angle);
  ny->x = - x1.x * sin(roll_angle) + y1.x * cos(roll_angle);
  ny->y = - x1.y * sin(roll_angle) + y1.y * cos(roll_angle);
  ny->z = - x1.z * sin(roll_angle) + y1.z * cos(roll_angle);
}


float getRollAngle(AttitudeCatalog* const ac, 
		   const double time, 
		   int* const status)
{
  // Find the appropriate entry in the AttitudeCatalog for the 
  // requested time.
  setAttitudeCatalogCurrentEntry(ac, time, status);
  CHECK_STATUS_RET(*status,0.);
    
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


