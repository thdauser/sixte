#include "attitude.h"


Attitude* getAttitude(int* const status)
{
  Attitude* ac=(Attitude*)malloc(sizeof(Attitude));
  if (NULL==ac) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for Attitude failed");
    return(NULL);
  }
  
  // Initialize.
  ac->entry    =NULL;
  ac->nentries =0;
  ac->currentry=0;
  ac->align    =ATTNX_NORTH;

  return(ac);
}


Attitude* loadAttitude(const char* filename, int* const status)
{
  Attitude* ac=NULL;
  AttitudeFile* af=NULL;

  do { // Beginning of ERROR handling loop

    // Read-in the attitude data from the FITS file and store 
    // them in the Attitude data structure.

    // Open the attitude file:
    headas_chat(5, "open attitude catalog file '%s' ...\n",
		filename);
    af=open_AttitudeFile(filename, READONLY, status);
    if (NULL==af) break;
    headas_chat(5, " attitude catalog contains %ld entries\n", 
		af->nrows);

    // Get a new Attitude object.
    ac=getAttitude(status);
    CHECK_STATUS_BREAK(*status);

    // Allocate memory for the entries in the catalog.
    ac->entry=(AttitudeEntry*)malloc(af->nrows*sizeof(AttitudeEntry));
    if (NULL==ac->entry) {
      *status=EXIT_FAILURE;
      SIXT_ERROR("not enough memory available to store the Attitude");
      break;
    }

    // Determine whether the roll angle alignment refers to the 
    // telescope's direction of motion or to the equatorial plane.
    int status2=EXIT_SUCCESS;
    char comment[MAXMSG], sbuffer[MAXMSG]; // String buffers.
    fits_write_errmark();
    strcpy(sbuffer, "");
    fits_read_key(af->fptr, TSTRING, "ALIGNMEN", sbuffer, comment, &status2);
    fits_clear_errmark();

    // Check the value of the header keyword and set the alignment flag
    // appropriately.
    strtoupper(sbuffer);
    if ((0==strlen(sbuffer)) || (!strcmp(sbuffer, "NORTH"))) {
      ac->align=ATTNX_NORTH;
    } else if (!strcmp(sbuffer, "MOTION")) {
      ac->align=ATTNX_MOTION;
    } else {
      *status=EXIT_FAILURE;
      SIXT_ERROR("invalid value for keyword ALIGNMEN in attitude file");
      break;
    } 
    
    // Read all lines from attitude file subsequently.
    AttitudeFileEntry afe;
    for (af->row=0; af->row<af->nrows; af->row++) {
      
      afe=read_AttitudeFileEntry(af, status);
      CHECK_STATUS_BREAK(*status);

      // Calculate and store attitude data:
      ac->entry[af->row].time=afe.time;

      // Telescope pointing direction nz:
      ac->entry[af->row].nz=
	unit_vector(afe.ra*M_PI/180., afe.dec*M_PI/180.);

      // Roll-Angle:
      ac->entry[af->row].roll_angle=afe.rollang*M_PI/180.;
    }
    CHECK_STATUS_BREAK(*status);
    // End of the attitude readout loop

    // Save the number of AttitudeEntry elements.
    ac->nentries=af->nrows;

  } while (0); // End of error handling loop


  // --- Clean up ---

  // Close FITS file
  if (NULL!=af) {
    if (af->fptr) fits_close_file(af->fptr, status);
    free(af);
  }

  return(ac);
}


void freeAttitude(Attitude** const ac)
{
  if (NULL!=(*ac)) {
    if (NULL!=(*ac)->entry) {
      free((*ac)->entry);
    }
    free(*ac);
    *ac=NULL;
  }
}


static void setAttitudeCurrEntry(Attitude* const ac,
				 const double time,
				 int* const status)
{
  // Check if this is a pointing attitude with only one data point.
  if (1==ac->nentries) {
    ac->currentry=0;
    return;
  }

  // Check if the requested time lies within the current time bin.
  while (time < ac->entry[ac->currentry].time) {
    // Check if the beginning of the Attitude is reached.
    if (ac->currentry <= 0) {
      *status=EXIT_FAILURE;
      char msg[MAXMSG]; 
      sprintf(msg, "no attitude entry available for time %lf", time);
      SIXT_ERROR(msg);
      return;
    }
    // If not, go one step back.
    ac->currentry--;
  }

  while (time > ac->entry[ac->currentry+1].time) {
    // Check if the end of the Attitude is reached.
    if (ac->currentry >= ac->nentries-2) {
      *status=EXIT_FAILURE;
      char msg[MAXMSG]; 
      sprintf(msg, "no attitude entry available for time %lf", time);
      SIXT_ERROR(msg);
      return;
    }
    // If not, go one step further.
    ac->currentry++;
  }
}


Vector getTelescopeNz(Attitude* const ac, 
		      const double time, 
		      int* const status)
{
  Vector nz={.x=0., .y=0., .z=0.};
 
  // Check if survey attitude.
  if (ac->nentries>1) {
    // Find the appropriate entry in the Attitude for the 
    // requested time.
    setAttitudeCurrEntry(ac, time, status);
    CHECK_STATUS_RET(*status,nz);
   
    // The requested time lies within the current time bin.
    // Interpolation:
    nz=interpolateCircleVector(ac->entry[ac->currentry].nz, 
			       ac->entry[ac->currentry+1].nz, 
			       (time-ac->entry[ac->currentry].time)/
			       (ac->entry[ac->currentry+1].time-
				ac->entry[ac->currentry].time));

  } else { // Pointing attitude.
    nz=ac->entry[0].nz;
  }
    
  return(nz);
}


void getTelescopeAxes(Attitude* const ac,
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
  if (1==ac->nentries) {
    // There is only one entry in the Attitude.
    pointed=1;
  } else {
    if (ac->currentry>0) {
      dnz=vector_difference(ac->entry[ac->currentry].nz, 
			    ac->entry[ac->currentry-1].nz);
    } else {
      dnz=vector_difference(ac->entry[ac->currentry+1].nz, 
			    ac->entry[ac->currentry].nz);
    }
    if (scalar_product(&dnz, &dnz)<1.e-10) {
      pointed=1;
    }
  }

  // Check if the x1 vector should be aligned perpendicular to the north 
  // direction or along the direction of motion of the telescope axis 
  // (neglecting the roll angle). For a pointed observation the alignment 
  // is done perpendicular to the north direction.
  Vector x1;
  if ((1==pointed) || (ATTNX_NORTH==ac->align)) {
    // Alignment perpendicular to the north direction.
    Vector b={0., 0., 1.};
    x1=vector_product(*nz, b);
    
    // CHECK if the telescope is pointing towards one of the poles.
    if (scalar_product(&x1, &x1)<1.e-20) {
      Vector c={1., 0., 0.};
      x1=vector_product(*nz, c);
    }

    // Normalize the vector in order to obtain a unit vector.
    x1=normalize_vector(x1);

  } else if (ATTNX_MOTION==ac->align) {
    // Alignment along the direction of motion of the telescope axis.
    Vector perpendicular=vector_product(*nz, dnz);
    x1=normalize_vector(vector_product(perpendicular, *nz));

  } else {
    *status=EXIT_FAILURE;
    SIXT_ERROR("invalid attitude alignment");
    return;
  }

  // Determine the y1 vector, which is perpendicular 
  // to the other two axes nz and x1:
  Vector y1=normalize_vector(vector_product(*nz, x1));

  // Take into account the roll angle.
  float roll_angle=getRollAngle(ac, time, status);
  CHECK_STATUS_VOID(*status);
  nx->x= x1.x * cos(roll_angle) + y1.x * sin(roll_angle);
  nx->y= x1.y * cos(roll_angle) + y1.y * sin(roll_angle);
  nx->z= x1.z * cos(roll_angle) + y1.z * sin(roll_angle);
  ny->x=-x1.x * sin(roll_angle) + y1.x * cos(roll_angle);
  ny->y=-x1.y * sin(roll_angle) + y1.y * cos(roll_angle);
  ny->z=-x1.z * sin(roll_angle) + y1.z * cos(roll_angle);
}


float getRollAngle(Attitude* const ac, 
		   const double time, 
		   int* const status)
{
  // Check if survey attitude.
  if (ac->nentries>1) {

    // Find the appropriate entry in the Attitude for the 
    // requested time.
    setAttitudeCurrEntry(ac, time, status);
    CHECK_STATUS_RET(*status,0.);
    
    // The requested time lies within the current time bin.
    // Interpolation:
    double fraction=
      (time-ac->entry[ac->currentry].time)/
      (ac->entry[ac->currentry+1].time-ac->entry[ac->currentry].time);
    return(ac->entry[ac->currentry  ].roll_angle*(1.-fraction) + 
	   ac->entry[ac->currentry+1].roll_angle*    fraction );

  } else { // Pointing attitude.
    return(ac->entry[0].roll_angle);
  }
}


AttitudeEntry defaultAttitudeEntry()
{
  AttitudeEntry ae;

  ae.time=0.0;

  ae.nz.x=0.0;
  ae.nz.y=0.0;
  ae.nz.z=0.0;

  ae.roll_angle=0.0;

  return(ae);
}

