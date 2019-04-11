/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

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
  ac->mjdref   =0.0;
  ac->tstart   =0.0;
  ac->tstop    =0.0;

  return(ac);
}

AttitudeEntry initializeAttitudeEntry ()
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

Attitude* loadAttitude(const char* const filename, int* const status)
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

    // Determine MJDREF, TSTART, and TSTOP.
    fits_read_key(af->fptr, TDOUBLE, "MJDREF", &ac->mjdref, comment, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'MJDREF' from attitude file '%s'",
	      filename);
      SIXT_ERROR(msg);
      break;
    }
    fits_read_key(af->fptr, TDOUBLE, "TSTART", &ac->tstart, comment, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'TSTART' from attitude file '%s'",
	      filename);
      SIXT_ERROR(msg);
      break;
    }
    fits_read_key(af->fptr, TDOUBLE, "TSTOP", &ac->tstop, comment, status);
    if (EXIT_SUCCESS!=*status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'TSTOP' from attitude file '%s'",
	      filename);
      SIXT_ERROR(msg);
      break;
    }

    // Check if TIMEZERO is present. If yes, it must be zero.
    double timezero;
    fits_write_errmark();
    status2=EXIT_SUCCESS;
    fits_read_key(af->fptr, TDOUBLE, "TIMEZERO", &timezero, comment, &status2);
    fits_clear_errmark();
    if (EXIT_SUCCESS==status2) {
      verifyTIMEZERO(timezero, status);
      CHECK_STATUS_BREAK(*status);
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
      dnz=vector_difference(*nz, ac->entry[ac->currentry-1].nz);
    } else {
      dnz=vector_difference(ac->entry[ac->currentry+1].nz, *nz);
    }
    if (scalar_product(&dnz, &dnz)<1.e-10) {
      pointed=1;
    }
  }

  // Check if the x1 vector should be aligned along the north direction
  // or along the direction of motion of the telescope axis
  // (neglecting the rotation according to the roll angle). For a pointed
  // observation the alignment is done along the north direction.
  Vector x1;
  if ((1==pointed) || (ATTNX_NORTH==ac->align)) {
    // Alignment along the north direction.

    // Check if the telescope is pointing towards one of the poles.
    Vector n1={1.0, 0.0, 0.0};
    Vector n2={0.0, 1.0, 0.0};
    if ((fabs(scalar_product(nz, &n1))<1.e-20) &&
	(fabs(scalar_product(nz, &n2))<1.e-20)) {
      x1.x=1.0;
      x1.y=0.0;
      x1.z=0.0;
    } else {
      // If not, align the x1 vector along the north direction.
      x1.x=0.0;
      x1.y=0.0;
      x1.z=1.0;
    }
  } else if (ATTNX_MOTION==ac->align) {
    // Alignment of nx along the direction of motion of the telescope axis.
    x1=normalize_vector(dnz);

  } else {
    *status=EXIT_FAILURE;
    SIXT_ERROR("invalid attitude alignment");
    return;
  }

  // Subtract the projection of the pointing direction
  // such that x1 and nz are perpendicular to each other.
  double scp=scalar_product(&x1, nz);
  x1.x-=scp*nz->x;
  x1.y-=scp*nz->y;
  x1.z-=scp*nz->z;

  // Normalize the vector in order to obtain a unit vector.
  x1=normalize_vector(x1);

  // Determine the y1 vector, which is perpendicular
  // to the other two axes nz and x1:
  Vector y1=normalize_vector(vector_product(*nz, x1));

  // Take into account the roll angle.
  float roll_angle=getRollAngle(ac, time, status);
  CHECK_STATUS_VOID(*status);
  double sinroll=sin(roll_angle);
  double cosroll=cos(roll_angle);
  nx->x= x1.x * cosroll + y1.x * sinroll;
  nx->y= x1.y * cosroll + y1.y * sinroll;
  nx->z= x1.z * cosroll + y1.z * sinroll;
  ny->x=-x1.x * sinroll + y1.x * cosroll;
  ny->y=-x1.y * sinroll + y1.y * cosroll;
  ny->z=-x1.z * sinroll + y1.z * cosroll;
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


AttitudeEntry* getAttitudeEntry(int* const status)
{
  AttitudeEntry* ae=(AttitudeEntry*)malloc(sizeof(AttitudeEntry));
  CHECK_NULL(ae, *status, "memory allocation for AttitudeEntry failed");

  ae->time=0.0;
  ae->nz.x=0.0;
  ae->nz.y=0.0;
  ae->nz.z=0.0;
  ae->roll_angle=0.0;

  return(ae);
}


Attitude* getPointingAttitude(const double mjdref,
			      const double tstart,
			      const double tstop,
			      const double ra,
			      const double dec,
			      const double roll_angle,
			      int* const status)
{
  // First allocate memory.
  Attitude* ac=getAttitude(status);
  CHECK_STATUS_RET(*status, ac);
  ac->entry=getAttitudeEntry(status);
  CHECK_STATUS_RET(*status, ac);
  ac->nentries=1;

  ac->mjdref=mjdref;
  ac->tstart=tstart;
  ac->tstop =tstop;
  ac->entry[0].time=tstart;
  ac->entry[0].nz=unit_vector(ra, dec);
  ac->entry[0].roll_angle=roll_angle;

  return(ac);
}


void checkAttitudeTimeCoverage(const Attitude* const ac,
			       const double mjdref,
			       const double tstart,
			       const double tstop,
			       int* const status)
{
  // Check that the values for MJDREF are equal.
   verifyMJDREF(mjdref, ac->mjdref, "in attitude file", status);
  CHECK_STATUS_VOID(*status);

  // Check the boundaries of the time interval.
  if ((ac->tstart > tstart) || (ac->tstop < tstop)) {
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "attitude does not cover the required period from %lf to %lf",
	    tstart, tstop);
    SIXT_ERROR(msg);
    return;
  }
}

void setWCScurrentPointing(const char* const filename, const Attitude* const ac,
			   Vector* const nz, struct wcsprm* wcs, int* const status)
{
  AttitudeFile* af=NULL;

  //open the attitude file:
  af=open_AttitudeFile(filename,READONLY,status);
  CHECK_STATUS_VOID(*status);
  AttitudeFileEntry afe;

  //get current AttitudeFileEntry
  af->row=ac->currentry;
  afe=read_AttitudeFileEntry(af,status);
  CHECK_STATUS_VOID(*status);
  double currRA=afe.ra;
  double currDEC=afe.dec;

  af->row=ac->currentry+1;
  afe=read_AttitudeFileEntry(af,status);
  CHECK_STATUS_VOID(*status);
  double nextRA=afe.ra;
  double nextDEC=afe.dec;

  if (NULL!=af) {
    if (af->fptr) fits_close_file(af->fptr,status);
    free(af);
  }

  //difference of current and succeeding attitude entries [deg]
  double diffRA=fabs(currRA-nextRA); //attitude-interval
  double diffDEC=fabs(currDEC-nextDEC);

  //get dec from actual current pointing (nz is approximated
   //and lies inbetween the attitude entries)
   //formula is from unit vector 'backwards'
  double DECnz=asin(nz->z)*180./M_PI; //could be the wrong value from asin
  if(fabs(currDEC-DECnz)>diffDEC){    //has to lie inbetween the attitude-interval
    DECnz=180.-DECnz;                 //if not, use the other solution
  }
  double RAnz=asin(nz->y/cos(DECnz*M_PI/180.))*180./M_PI;
  if(fabs(currRA-RAnz)>diffRA){
    RAnz=360.-RAnz;
  }

  //set wcs to actual current pointing
  wcs->crval[0]=RAnz;  //in deg
  wcs->crval[1]=DECnz;
}

void getCurrentVelocity(const char* const filename, const Attitude* const ac,
			double* vel_ra, double* vel_dec, double const att_start,
			double const att_stop, int* const status)
{
  AttitudeFile* af=NULL;

  //open the attitude file:
  af=open_AttitudeFile(filename,READONLY,status);
  CHECK_STATUS_VOID(*status);
  AttitudeFileEntry afe;

  //get current AttitudeFileEntry
  af->row=ac->currentry;
  afe=read_AttitudeFileEntry(af,status);
  CHECK_STATUS_VOID(*status);
  double currRA=afe.ra;
  double currDEC=afe.dec;

  af->row=ac->currentry+1;
  afe=read_AttitudeFileEntry(af,status);
  CHECK_STATUS_VOID(*status);
  double nextRA=afe.ra;
  double nextDEC=afe.dec;

  if (NULL!=af) {
    if (af->fptr) fits_close_file(af->fptr,status);
    free(af);
  }

  //difference of current and succeeding attitude entries [deg]
  double diffRA=fabs(currRA-nextRA); //attitude-interval
  double diffDEC=fabs(currDEC-nextDEC);

  //time interval between two current pointings
  double interval=att_stop-att_start;

  *vel_ra=diffRA/interval;
  *vel_dec=diffDEC/interval;
}


void convert_galLB2RAdec(double* world){
 /**  % following Carroll & Ostlie, Chapter 24.3, Coordinates obtained with HEASARC tools
  *      (taken from the Remeis ISISscripts, initially written by Moritz Böck)
  *
  *      expectiong l=world[0] and b=world[1]
  *      output     RA=world[1] and dec=world[1]
  */

	double l = world[0]/180.*M_PI;
	double b = world[1]/180.*M_PI;
	double cos_b = cos(b);
	double sin_b = sin(b);

	double cos_d_ngp = 0.8899880874849542;
	double sin_d_ngp = 0.4559837761750669;
	double l_ncp = 2.145566759798267518;

	world[0] = fmod((atan2(cos_b*sin(l_ncp - l) , cos_d_ngp*sin_b-sin_d_ngp*cos_b*cos(l_ncp - l))
			+3.3660332687500039)/M_PI*180. ,360.); // RA
	// need a positive modulo here
	if (world[0] < 0.0){
		world[0] += 360.;
	}
	world[1] =  asin( sin_d_ngp*sin_b + cos_d_ngp*cos_b*cos(l_ncp - l) )/M_PI*180.; // dec

}

// write an attitude to a FITS file
void write_attitude(Attitude* att, char* filename,
                     double mjdref,
                     int clobber, int* status) {

  CHECK_STATUS_VOID(*status);

  // open a FITS file
  fitsfile *fptr;
  fits_create_file_clobber(&fptr, filename, clobber, status);
  CHECK_STATUS_VOID(*status);


  // create an attitude extension
  char *ttype[]={"TIME","RA","DEC","ROLLANG"};
  char *tform[]={"1D","1D","1D","1D"};
  char *tunit[]={"s","deg","deg","deg"};
  long nrows = att->nentries;
  fits_create_tbl(fptr, BINARY_TBL, nrows, 4,
                  ttype, tform, tunit,
                  "ATT", status);

  // add header keywords
  fits_update_key(fptr, TDOUBLE, "MJDREF", &mjdref, NULL, status);
  fits_update_key(fptr, TDOUBLE, "TSTART", &(att->tstart), NULL, status);
  fits_update_key(fptr, TDOUBLE, "TSTOP", &(att->tstop), NULL, status);
  char *align="NORTH";
  fits_update_key(fptr, TSTRING, "ALIGNMEN", &align, NULL, status);
  char *orig="ECAP";
  fits_update_key(fptr, TDOUBLE, "ORIGIN", &orig, NULL, status);

  // write the columns
  // make empty arrays of time, ra, dec and angle
  double* time = (double*) malloc(nrows * sizeof(double));
  CHECK_MALLOC_VOID_STATUS(time, *status);
  double* ra = (double*) malloc(nrows * sizeof(double));
  CHECK_MALLOC_VOID_STATUS(ra, *status);
  double* dec = (double*) malloc(nrows * sizeof(double));
  CHECK_MALLOC_VOID_STATUS(dec, *status);
  double* roll = (double*) malloc(nrows * sizeof(double));
  CHECK_MALLOC_VOID_STATUS(roll, *status);

  // loop over attitude entries
  long ii;
  for (ii=0; ii<nrows; ii++) {
    AttitudeEntry entry = att->entry[ii];
    // time and rollis just copied
    time[ii] = entry.time;
    roll[ii] = entry.roll_angle;

    // ra and dec are calculated from nz and nx;
    calculate_ra_dec(entry.nz, &ra[ii], &dec[ii]);
    ra[ii] *= 180./M_PI;
    dec[ii] *= 180./M_PI;
  }

  // dump into the file
  fits_write_col(fptr, TDOUBLE, 1, 1, 1, nrows, time, status);
  fits_write_col(fptr, TDOUBLE, 2, 1, 1, nrows, ra, status);
  fits_write_col(fptr, TDOUBLE, 3, 1, 1, nrows, dec, status);
  fits_write_col(fptr, TDOUBLE, 4, 1, 1, nrows, roll, status);

  // free memory
  free(time);
  free(ra);
  free(dec);
  free(roll);

  fits_close_file(fptr, status);

  // check for FITS errors
  if (*status!=EXIT_SUCCESS) {
    char fits_err_msg[80];
    fits_get_errstatus(*status, fits_err_msg);
    SIXT_ERROR(fits_err_msg);
    fits_clear_errmark();
  }
}

// calculate a lissajous pattern from time 0 to 1, with Amplitude ampl
void lissajous_pattern(double ampl, double** t, double** x, double** y, int nt){

	  double afrac  = 3.0/4.;
	  double bfrac  = 4.0/4.;

	  double offset = M_PI/4.;

	  for (int ii=0; ii<nt; ii++){
		   (*t)[ii] = (1.0*ii) / (nt-1.0) ; //[0:1.0:#1000];
		   (*x)[ii] = sin(afrac*(*t)[ii]*2*M_PI*4+offset)*ampl;
		   (*y)[ii] = sin(bfrac*(*t)[ii]*2*M_PI*4)*ampl;
	  }

}

// ampl given in [arcsec]
Attitude* get_default_attitude_lissajous(double ampl, double ra0, double dec0,
		double tstart, double tstop, double mjdref, int* status){

	int nentries = 1000;  // for this setup 1000 points make a smooth curve

	// initialize the attitude file
	Attitude* ac = getAttitude(status);
    CHECK_MALLOC_RET_NULL(ac);

    ac->nentries = nentries;
	ac->mjdref = mjdref;
	ac->tstart = tstart;
	ac->tstop  = tstop;

    // Allocate memory for the entries in the catalog.
    ac->entry=(AttitudeEntry*)malloc(nentries*sizeof(AttitudeEntry));
    CHECK_MALLOC_RET_NULL(ac->entry);

    double* t = (double*) malloc(sizeof(double)*nentries);
    double* ra = (double*) malloc(sizeof(double)*nentries);
    double* dec = (double*) malloc(sizeof(double)*nentries);

    CHECK_MALLOC_RET_NULL_STATUS(t,*status);
    CHECK_MALLOC_RET_NULL_STATUS(ra,*status);
    CHECK_MALLOC_RET_NULL_STATUS(dec,*status);

    lissajous_pattern(ampl, &t, &ra, &dec, nentries);

    for (int ii=0; ii<nentries; ii++){
    	ac->entry[ii] = initializeAttitudeEntry ();

        ac->entry[ii].time = (tstop-tstart)*t[ii];
        ac->entry[ii].time += tstart;

        ac->entry[ii].nz = unit_vector(M_PI/180.*(ra[ii]+ra0), M_PI/180.*(dec[ii]+dec0));
    }

    free(t);
    free(ra);
    free(dec);

    return ac;
}
