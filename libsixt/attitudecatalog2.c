#include "attitudecatalog2.h"

AttCatalog* getAttCatalog(int* const status)
{
  AttCatalog* ac=(AttCatalog*)malloc(sizeof(AttCatalog));
  if (NULL==ac) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("memory allocation for AttitudeCatalog failed!\n", *status);
    return(NULL);
  }
  
  // Initialize.
  ac->nentries = 0;
  ac->current_entry = 0;
  ac->entry    = NULL;

  return(ac);
}

AttEntry initializeAttitudeEntry ()
{
   AttEntry ae;
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

AttCatalog* loadAttCatalog(const char* filename, int* const status)
{
  AttCatalog* ac=NULL;
  AttitudeFile* af=NULL;

  do{  //Beginning error handling loop

    //Open attitude file
    headas_chat(5, "open attitude file '%s'...\n", filename);
    af=open_AttitudeFile(filename, READONLY, status);
    if(NULL==af){
      *status = EXIT_FAILURE;
      SIXT_ERROR("Could not open Attitude File.");
      break;
    }

    //Get empty AttitudeCatalog-object.
    ac=getAttCatalog(status);
    CHECK_STATUS_BREAK(*status);

    //Memory-Allocation for AttitudeCatalog-entries:
    ac->entry=(AttEntry*)malloc(af->nrows*sizeof(AttEntry));
    if (NULL==ac->entry) {
      *status = EXIT_FAILURE;
      SIXT_ERROR("not enough memory available to store "
		 "the AttitudeCatalog\n");
      break;
    }

    //Determine rollangle alignment.
    //Read header keyword to string-buffer (sbuffer).
    int status2=EXIT_SUCCESS;
    char comment[MAXMSG], sbuffer[MAXMSG]={""};
    fits_write_errmark();
    fits_read_key(af->fptr, TSTRING, "ALIGNMEN", &sbuffer, comment, &status2);
    fits_clear_errmark();
    //Capital letters to compare keyword.
    strtoupper(sbuffer);
    //No header keyword or 'equator'-> alignment=0.
    //rollangle with respect to equatorial plane.
    if ((0==strlen(sbuffer)) || (!strcmp(sbuffer, "EQUATOR"))) {
      ac->alignment=0;
    //Header keyword 'motion'-> alignment=1.
    //rollangle with respect to telescope-motion-direction.
    } else if (!strcmp(sbuffer, "MOTION")) {
      ac->alignment=1;
    //Other header keyword-> invalid.
    } else {
      *status=EXIT_FAILURE;
      SIXT_ERROR("invalid value for telescope alignment flag");
      break;
    } 

    //Read data in attitude file subsequently.
    AttitudeFileEntry afe;

    for(af->row=0; af->row < af->nrows; af->row++ ){
      afe=read_AttitudeFileEntry(af, status);
      if (EXIT_SUCCESS!=*status) break;

      ac->entry[af->row].time=afe.time;

      //Telescope pointing direction:
      ac->entry[af->row].nz=unit_vector(afe.ra*M_PI/180., afe.dec*M_PI/180.);

      //Rollangle:
      ac->entry[af->row].roll_angle=afe.rollang*M_PI/180.;
    }
    if (EXIT_SUCCESS!=*status) break;

    //Number of AttitudeEntry-elements:
    ac->nentries=af->nrows;

    //Telescope-nx-direction (for all entries in the AttitudeCatalog)
    //If the change in telescope-pointing-direction nz is too small
    //then align nx parallel to the equatorial plane (perpendicular ny-nz-plane)

    //First entry:
    //Change of the telescope pointing direction:
    Vector dnz=vector_difference(ac->entry[1].nz, ac->entry[0].nz);
    //if change is not negligible
    if (sqrt(scalar_product(&dnz, &dnz))>1.e-10){
      //nx pointing along telescope motion (nx=normalize((nz_0xnz_1)xnz_0))
      ac->entry[0].nx=normalize_vector(vector_product
				      (vector_product(ac->entry[0].nz, ac->entry[1].nz),
				       ac->entry[0].nz ));
    }else{
    // change is negligible
      Vector ny={0., 1., 0.};
      //nx parallel to equatorial plane
      ac->entry[0].nx=vector_product(ac->entry[0].nz, ny);
    }
    
    //all other entries:
    long ii;
    for(ii=1; ii < ac->nentries; ii++){
      //Change of the telescope pointing direction:
      Vector dnz=vector_difference(ac->entry[ii].nz, ac->entry[ii-1].nz);
      //if change is not negligible
      if (sqrt(scalar_product(&dnz, &dnz))>1.e-10){
      //nx pointing along telescope motion (difference in nz equals direction of motion)
	ac->entry[ii].nx=normalize_vector(dnz);
      }else{
      // change is negligible
	Vector ny={0., 1., 0.};
	//nx parallel to equatorial plane
	ac->entry[ii].nx=vector_product(ac->entry[ii].nz, ny);
      }
    }

  }while(0); //End error handling loop

  //Cleaning-up
  
  //Close FITS-file
    if (NULL!=af) {
    if (af->fptr) fits_close_file(af->fptr, status);
    free(af);
  }

  if (EXIT_SUCCESS != *status) ac = NULL;
  return(ac);
}

static void setAttCatalogCurrentEntry(AttCatalog* const ac,
					   const double time,
					   int* const status)
{
  // Check if this is a pointing attitude with only one data point.
  if (1==ac->nentries) {
    ac->current_entry=0;
    return;
  }

  //Check if the requested time is < the current time.
  while (time < ac->entry[ac->current_entry].time){
    //Check if the beginning of the AttitudeCatalog is reached.
    if (ac->current_entry <= 0) {
      *status = EXIT_FAILURE;
      char msg[MAXMSG]; 
      sprintf(msg, "no attitude entry available for time %lf!\n", time);
      HD_ERROR_THROW(msg, *status);
      return;
    }
    //Go one step back until the requested time is >= current time.
    //The loop will stop at a entry whose time is < requested time.
    ac->current_entry--;
  }

  //Check if the requested time is > the current time.
  while (time > ac->entry[ac->current_entry+1].time) {
    //Check if the end of the AttitudeCatalog is reached.
    if (ac->current_entry >= ac->nentries-2) {
      *status = EXIT_FAILURE;
      char msg[MAXMSG]; 
      sprintf(msg, "no attitude entry available for time %lf!\n", time);
      HD_ERROR_THROW(msg, *status);
      return;
    }
    //Go one step further until requested time <= (current time +1).
    //The loop will stop at a entry whose time is < requested time.
    ac->current_entry++;
  }
}



Vector GetTelescopeNz (AttCatalog* ac, const double time, 
			int* const status)
{
  //Initialize to zero:
  Vector nz={.x=0., .y=0., .z=0.};

  //Check if attitude was loaded from file.
  if(ac->nentries > 1){
    //Find the entry (the closest one) for the requested time.
    setAttCatalogCurrentEntry(ac, time, status);
    CHECK_STATUS_RET(*status,nz);

    //Interpolate actual nz which lies inbetween the current_entry.nz and
    //current_entry+1.nz
    //the function determines the factor of the angle (between the two entries)
    //by taking the following ratio into account
    nz=interpolateCircleVector(ac->entry[ac->current_entry].nz,
			       ac->entry[ac->current_entry+1].nz,
			       (time - ac->entry[ac->current_entry].time)/
			       (ac->entry[ac->current_entry+1].time - 
				ac->entry[ac->current_entry].time));

  }else{ //Pointing attitude.
    nz=ac->entry[0].nz;
  }

  return(nz);
}

void freeAttCatalog(AttCatalog** const ac)
{
  if (NULL != (*ac)) {
    if (NULL != (*ac)->entry) {
      free((*ac)->entry);
    }
    free(*ac);
    *ac=NULL;
  }
}


void getTelAxes(AttCatalog* const ac,
		Vector* const nx, //no return value, BUT
		Vector* const ny, //pointers to structure of type Vector;
		Vector* const nz, //nx,ny,nz will be changed.
		const double time,
		int* const status)
{
  //Determine current telescope pointing direction
  *nz=GetTelescopeNz(ac, time, status);
  CHECK_STATUS_VOID(*status);
  //the current-entry-pointer points now to the entry closest
  //to requested time.

  //Check if this is a pointed observation.
  int pointobs=0;
  Vector dnz;
  if(ac->nentries==1){
    //only one entry-> pointed observation
    pointobs=1;
  }else{
    //more than one entry-> attitude-file
    //determine the difference in nz between the two
    //subsequent entries.
    if(ac->current_entry>0){
      //current entry is not the first
      dnz=vector_difference(ac->entry[ac->current_entry].nz,
			    ac->entry[ac->current_entry-1].nz);
    }else{
      //current entry is the first
      dnz=vector_difference(ac->entry[ac->current_entry+1].nz,
			    ac->entry[ac->current_entry].nz);
    }
    if(scalar_product(&dnz,&dnz)<1.e-10){
      //change in nz is negligible
      pointobs=1;
    }
  }

  //Check how the x1-vector should be aligned.
  //either parallel to the equatorial plane or in direction
  //of the telescope motion. pointed observation means
  //parallel to the equatorial-plane.
  Vector x1;
  if((pointobs==1) || (ac->alignment==0)){
    //parallel to the equatorial plane

    Vector z={0.,0.,1.};
    //Vector perpendicular to z-nz-plane,within equatorial plane
    x1=vector_product(*nz,z);
    //special case:check if the telescope is pointing towards one of the 
    //poles-> angle between z and nz is small-> x1 is small
    if(scalar_product(&x1,&x1)<1.e-10){
      Vector x={1.,0.,0.};
      x1=vector_product(*nz,x);
    }
    //make x1 a unit vector
    x1=normalize_vector(x1);

  }else{
    //along the telescope-motion-direction

    Vector perpendicular=vector_product(*nz,dnz);
    //vector-product of perpendicular-vector with nz-> resulting
    //vector in direction of motion
    x1=normalize_vector(vector_product(perpendicular,*nz));
  }

  //Determine the third axis perpendicular to nz and x1:
  Vector y1=normalize_vector(vector_product(*nz, x1));

  //Now take the roll-angle into account.
  //The x1-axis has already been set with respect to the
  //roll-angle definition.
  float roll_angle=GetRollAngle(ac, time, status);
  CHECK_STATUS_VOID(*status);

  //nx(ny) is the projection onto the x1- and y1-axis via the roll-angle.
  nx->x = x1.x * cos(roll_angle) + y1.x * sin(roll_angle);
  nx->y = x1.y * cos(roll_angle) + y1.y * sin(roll_angle);
  nx->z = x1.z * cos(roll_angle) + y1.z * sin(roll_angle);

  ny->x = - x1.x * sin(roll_angle) + y1.x * cos(roll_angle);
  ny->y = - x1.y * sin(roll_angle) + y1.y * cos(roll_angle);
  ny->z = - x1.z * sin(roll_angle) + y1.z * cos(roll_angle);
}


float GetRollAngle(AttCatalog* const ac,
		   const double time,
		   int* const status)
{
  //Check if pointing or survey attitude
  if(ac->nentries>1){
    //survey attitude(attitude-file)
    
    //Find the appropriate attitude-catalog-entry for the
    //requested time
    setAttCatalogCurrentEntry(ac,time,status);
    CHECK_STATUS_RET(*status,0.);

    //Calculate the fraction for the interpolation
    //(difference to requested time/time-difference in two subsequent
    //attitude-entries)
    double fraction=
      (time - ac->entry[ac->current_entry].time)/
      (ac->entry[ac->current_entry+1].time - ac->entry[ac->current_entry].time);

    //Interpolate current roll-angle
    //(current angle + fraction * difference in subsequent angles)
    return(ac->entry[ac->current_entry].roll_angle + fraction*
	  (ac->entry[ac->current_entry+1].roll_angle -
	   ac->entry[ac->current_entry].roll_angle));
    
  }else{
    //pointing attitude

    return(ac->entry[0].roll_angle);
  }

}
