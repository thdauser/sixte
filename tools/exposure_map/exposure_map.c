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
*/

#include "exposure_map.h"


void saveExpoMap(float** const map,
		 const char* const filename,
		 const long naxis1, const long naxis2,
		 struct wcsprm* const wcs,
		 const char clobber,
		 int* const status) 
{
  // 1d image buffer for storing in FITS image.
  float* map1d=NULL;

  // FITS file pointer for exposure map image.
  fitsfile* fptr=NULL;

  // String buffer for FITS header.
  char* headerstr=NULL;

  // Store the exposure map in a FITS file image.
  headas_chat(3, "\nstore exposure map in FITS image '%s' ...\n", 
	      filename);

  do { // Beginning of error handling loop.

    // Check if the file already exists.
    int exists;
    fits_file_exists(filename, &exists, status);
    CHECK_STATUS_VOID(*status);
    if (0!=exists) {
      if (0!=clobber) {
	// Delete the file.
	remove(filename);
      } else {
	// Throw an error.
	char msg[MAXMSG];
	sprintf(msg, "file '%s' already exists", filename);
	SIXT_ERROR(msg);
	*status=EXIT_FAILURE;
	return;
      }
    }

    // Create a new FITS-file (remove existing one before):
    fits_create_file(&fptr, filename, status);
    CHECK_STATUS_BREAK(*status);

    // Convert the exposure map to a 1d-array to store it in the FITS image.
    map1d=(float*)malloc(naxis1*naxis2*sizeof(float));
    CHECK_NULL_BREAK(map1d, *status, 
		    "memory allocation for 1d exposure map buffer failed");
    long x;
    for (x=0; x<naxis1; x++) {
      long y;
      for (y=0; y<naxis2; y++) {
	map1d[x + y*naxis1]=map[x][y];
      }
    }

    // Create an image in the FITS-file (primary HDU):
    long naxes[2]={ naxis1, naxis2 };
    fits_create_img(fptr, FLOAT_IMG, 2, naxes, status);
    //                               |-> naxis
    CHECK_STATUS_BREAK(*status);

    // Add header information about program parameters.
    // The second parameter "1" means that the headers are written
    // to the first extension.
    HDpar_stamp(fptr, 1, status);
    CHECK_STATUS_BREAK(*status);

    // Write WCS header keywords.
    int nkeyrec;
    if (0!=wcshdo(0, wcs, &nkeyrec, &headerstr)) {
      SIXT_ERROR("construction of WCS header failed");
      *status=EXIT_FAILURE;
      break;
    }
    char* strptr=headerstr;
    while (strlen(strptr)>0) {
      char strbuffer[81];
      strncpy(strbuffer, strptr, 80);
      strbuffer[80]='\0';
      fits_write_record(fptr, strbuffer, status);
      CHECK_STATUS_BREAK(*status);
      strptr+=80;
    }
    CHECK_STATUS_BREAK(*status);

    // Write the image to the file.
    long fpixel[2]={1, 1}; // Lower left corner.
    //              |--|--> FITS coordinates start at (1,1), NOT (0,0).
    // Upper right corner.
    long lpixel[2]={naxis1, naxis2}; 
    fits_write_subset(fptr, TFLOAT, fpixel, lpixel, map1d, status);
    CHECK_STATUS_BREAK(*status);

  } while(0); // End of error handling loop.

  // Close the exposure map FITS file.
  if(NULL!=fptr) fits_close_file(fptr, status);

  // Release memory.
  if (NULL!=map1d) {
    free(map1d);
  }
  if (NULL!=headerstr) {
    free(headerstr);
  }
}  

static char* parse_string(char *str){
	char *pch;
	pch = strtok (str,";");
	return pch;
}

static int parse_string2array(char *str, char ***strarray, int *status){

	const int MAXSTR = 9999;
	char buffer[MAXSTR];
	strcpy(buffer,str); // make sure we do not destroy the string

	char *pch = NULL;
	pch = parse_string(buffer);

	int ii=0;

	// allocate array
	*strarray = (char**) malloc(  1 * sizeof(char*) );
	CHECK_NULL_RET(*strarray,*status,"malloc failed",-1);
	*(strarray[ii]) = (char*) malloc (MAXSTR * sizeof(char));
	strcpy((*strarray)[ii],pch);

	while ( pch != NULL){
		ii++;
		pch = parse_string(NULL);
		if (pch != NULL){
			*strarray = (char**) realloc( *strarray,  (ii+1) * sizeof(char*) );
			CHECK_NULL_RET(*strarray,*status,"realloc failed",-1);
			(*strarray)[ii] = (char*) malloc (MAXSTR * sizeof(char));
			strcpy((*strarray)[ii],pch);
		}
	}

	return ii;
}


static void init_expoMap(struct Parameters par, float ***expoMap, int *status){

	// float ** expoMap = NULL;

	// initialize the exposure map
	(*expoMap)=(float**)malloc(par.ra_bins*sizeof(float*));
  	CHECK_NULL_VOID((*expoMap),*status,"Error when allocating memory for the Expsore Map");

	long ix;
	for (ix=0; ix<par.ra_bins; ix++) {
		(*expoMap)[ix]=(float*)malloc(par.dec_bins*sizeof(float));
	  	CHECK_NULL_VOID((*expoMap)[ix],*status,"Error when allocating memory for the Expsore Map");
	}

	// Clear the exposure map
	for (ix=0; ix<par.ra_bins; ix++) {
		// Clear the exposure map.
		long iy;
		for (iy=0; iy<par.dec_bins; iy++) {
			(*expoMap)[ix][iy]=0.0;
		}
	}

	// expoMapIn = &expoMap;

	CHECK_STATUS_VOID(status);

	return;
}

static void init_expo_wcs(struct Parameters par, struct wcsprm *wcs, int *status){

    // Set up the WCS data structure.
    if (0!=wcsini(1, 2, wcs)) {
      SIXT_ERROR("initalization of WCS data structure failed");
      *status=EXIT_FAILURE;
      return;
    }

    wcs->naxis=2;
    wcs->crpix[0]=par.ra_bins/2 +0.5;
    wcs->crpix[1]=par.dec_bins/2+0.5;
    wcs->crval[0]=0.5*(par.ra1 +par.ra2 )*180./M_PI;
    wcs->crval[1]=0.5*(par.dec1+par.dec2)*180./M_PI;
    wcs->cdelt[0]=(par.ra2 -par.ra1 )*180./M_PI/par.ra_bins;
    wcs->cdelt[1]=(par.dec2-par.dec1)*180./M_PI/par.dec_bins;
    strcpy(wcs->cunit[0], "deg");
    strcpy(wcs->cunit[1], "deg");
    if ((1==par.projection)||(3==par.projection)) {
      strcpy(wcs->ctype[0], "RA---AIT");
      strcpy(wcs->ctype[1], "DEC--AIT");
    } else if (2==par.projection) {
      strcpy(wcs->ctype[0], "RA---SIN");
      strcpy(wcs->ctype[1], "DEC--SIN");
    } else {
      SIXT_ERROR("projection type not supported");
      *status=EXIT_FAILURE;
      return;
    }

    return;
}

static FILE* get_progressout(struct Parameters par,int *status){
	char ucase_buffer[MAXFILENAME];
	strcpy(ucase_buffer, par.ProgressFile);
	strtoupper(ucase_buffer);

	FILE* progressfile=NULL;

	if (0!=strcmp(ucase_buffer, "STDOUT")) {
		progressfile=fopen(par.ProgressFile, "w+");
		char msg[MAXMSG];
		sprintf(msg, "could not open file '%s' for output of progress status",
				par.ProgressFile);
		CHECK_NULL_RET(progressfile, *status, msg,NULL);
	}
	return progressfile;
}

static Attitude* get_attitude(struct Parameters par, int *status){
	char ucase_buffer[MAXFILENAME];
	strcpy(ucase_buffer, par.Attitude);
	strtoupper(ucase_buffer);
	Attitude *ac = NULL;
	if ((strlen(par.Attitude)==0)||(0==strcmp(ucase_buffer, "NONE"))) {
		// Set up a simple pointing attitude.
		ac=getPointingAttitude(0., par.TSTART, par.TSTART+par.timespan,
				par.RA*M_PI/180., par.Dec*M_PI/180., status);
		CHECK_STATUS_RET(*status,NULL);

	} else {
		// Load the attitude from the given file.
		ac=loadAttitude(par.Attitude, status);
		CHECK_STATUS_RET(*status,NULL);

		// Check if the required time interval for the simulation
		// is a subset of the time described by the attitude file.
		// Note that MJDREF is assumed as the value from the attitude file.
		checkAttitudeTimeCoverage(ac, ac->mjdref, par.TSTART,
				par.TSTART+par.timespan, status);
		CHECK_STATUS_RET(*status,NULL);
	}
    CHECK_NULL_RET(ac,*status,"Failed to read Attitudefile",NULL);
	return ac;
	// END of setting up the attitude.
}

static Vignetting* get_vign(struct Parameters par, int *status){

	char ucase_buffer[MAXFILENAME];
	Vignetting* vignetting = NULL;
	if (strlen(par.Vignetting)>0) {
		strcpy(ucase_buffer, par.Vignetting);
		strtoupper(ucase_buffer);
		if (0!=strcmp(ucase_buffer, "NONE")) {
			vignetting=newVignetting(par.Vignetting, status);
			CHECK_STATUS_RET(*status,NULL);
		}
	}
	return vignetting;
}

static unsigned int init_progress(FILE *progressfile){
	if (NULL==progressfile) {
		headas_chat(2, "\r%.0lf %%", 0.);
		fflush(NULL);
	} else {
		rewind(progressfile);
		fprintf(progressfile, "%.2lf", 0.);
		fflush(progressfile);
	}
	return 0;
}

static double get_min_fov_align(GenInst* inst, struct Parameters par){
	const double fov_min_align=cos(inst->tel->fov_diameter*0.5); // / inst->tel->focal_length; //  cos(  fov/2.);
	double fov_diameter = acos(fov_min_align)*2;
    double field_min_align;

    if ((par.ra2-par.ra1 > M_PI/6.) ||
    		(par.dec2-par.dec1 > M_PI/6.)) {
    	// Actually -1 should be sufficient, but -2 is even safer.
    	field_min_align=-2.;
    } else {
    	field_min_align=
    			cos((sqrt(pow(par.ra2-par.ra1, 2.)+pow(par.dec2-par.dec1, 2.))+
    					fov_diameter)/2.);
    }
    return field_min_align;
}

static int file_exist (char *filename) {
  struct stat   buffer;
  return stat (filename, &buffer) ;
}

static void get_world_coords(int x, int y, struct wcsprm *wcs,
		double **world, double *phi, double *theta, int *status){

	double pixcrd[2]={ x+1., y+1. };
	double imgcrd[2];

	(*world) =(double*)malloc(2*sizeof(float));
	CHECK_NULL_VOID(*world,*status,"malloc failed");

	int status2=0;
	wcsp2s(wcs, 1, 2, pixcrd, imgcrd, phi, theta, *world, &status2);
	if (3==status2) {
		// Pixel does not correspond to valid world coordinates.
		SIXT_ERROR("Pixel does not correspond to valid world coordinates.");
		*status=EXIT_FAILURE;
		return;
	} else if (0!=status2) {
		SIXT_ERROR("projection failed");
		*status=EXIT_FAILURE;
		return;
	}
	return;
}

static int get_pixel_hit(GenInst **inst,int ninst, Vector pixel_skypos, double theta, double phi,
		struct Telescope telescope){


	// Calculate the off-axis angle ([rad]).
	double cos_theta=scalar_product(&telescope.nz, &pixel_skypos);
	// Avoid numerical problems with numbers slightly larger than 1.
	if ((cos_theta>1.0) && (cos_theta-1.0<1.e-10)) {
		cos_theta=1.0;
	}
	assert(cos_theta<=1.0);
	theta=acos(cos_theta);

	// Calculate the azimuthal angle ([rad]) of the source position.
	phi=atan2(scalar_product(&telescope.ny, &pixel_skypos),
			scalar_product(&telescope.nx, &pixel_skypos));


	double sinp, cosp;
#if defined( __APPLE__) && defined(__MACH__)
  __sincos(phi, &sinp, &cosp);
#else
  sincos(phi, &sinp, &cosp);
#endif

	double distance=inst[0]->tel->focal_length * tan(theta);
  	double posx = cosp*distance;
  	double posy = sinp*distance;

  	int rawx, rawy;
  	double rfx, rfy;

  	int ii;
  	for (ii=0;ii<ninst;ii++){
  		getGenDetAffectedPixel(inst[ii]->det->pixgrid, posx, posy,
  				&rawx,&rawy,&rfx,&rfy);
  		if (rawx!=-1 && rawy!=-1){
  			return 1;
  		}
  	}
  	// return 0 if no pixel is hit
	return 0;
}

static float get_single_expos_value(int x,int y, struct wcsprm *wcs,
		struct Telescope telescope, GenInst **inst, int ninst, int *status){

	// get world coordinates from the projection
	double *world, theta, phi;
	get_world_coords(x,y,wcs,&world,&theta,&phi,status);
	CHECK_STATUS_RET(*status,0.0);

	// Determine a unit vector for the calculated RA and Dec.
	Vector pixpos=unit_vector(world[0]*M_PI/180., world[1]*M_PI/180.);
	free(world);

	// instrument is the same
	const double fov_min_align=inst[0]->tel->fov_diameter*0.5;

	// double delta=acos(scalar_product(&telescope_nz, &pixpos));
//	double delta=acos(scalar_product(&telescope.nz, &pixpos));

	// Check if the current pixel lies within the FOV.
	//  (only a rough, conservative estimate to reduce computing power)
	if (check_fov(&pixpos, &telescope.nz, cos(fov_min_align))==0) {
		// Pixel lies inside the FOV!

		if (get_pixel_hit(inst,ninst,pixpos,theta,phi,telescope)==1){

			// Add the exposure time step weighted with the vignetting
			// factor for this particular off-axis angle at 1 keV.
			// Calculate the off-axis angle ([rad])
			return  acos(scalar_product(&telescope.nz, &pixpos));
		} else {
			return -1.0;
		}
	} else {
		return -1.0;
	}
}

static void get_geninst_all(GenInst **inst,xmlarray xmls,const unsigned int seed,int *status){
    int ii;
    for (ii=0;ii<xmls.n;ii++){
    // Load the instrument configuration
    	if (file_exist(xmls.xmlarray[ii])!=0){
    		printf("Cannot read XMLFile %s\n",xmls.xmlarray[ii]);
    		// sixt_error("exposure_map","Cannot read XMLFile");
    		*status = EXIT_FAILURE;
    		CHECK_STATUS_VOID(*status);
    	} else {
    		headas_chat(3, "Loading XMLFile %s ...\n",xmls.xmlarray[ii]);
    		inst[ii]=loadGenInst(xmls.xmlarray[ii], seed, status);
    		CHECK_STATUS_VOID(*status);

    	}
    }

}

static void release_expoMap(float *** expoMap, int n){
if (NULL!=*expoMap) {
	  long x;
	  for (x=0; x<n; x++) {
		  if (NULL!=(*expoMap)[x]) {
			  free((*expoMap)[x]);
		  }
	  }
	  free(*expoMap);
}
}

static int check_consist(GenInst **inst, int n){

	int ii;
	double buf;
	buf = inst[0]->tel->fov_diameter;
	for (ii=1;ii<n;ii++){
		if (inst[ii]->tel->fov_diameter != buf){
			printf("Error: XML Files not consistent!");
			return EXIT_FAILURE;
		}
	}
	return EXIT_SUCCESS;
}

int exposure_map_main()
{
  // Program parameters.
  struct Parameters par;

  // Error status.
  int status=EXIT_SUCCESS;

  // Register HEATOOL:
  set_toolname("exposure_map");
  set_toolversion("0.1");

  struct wcsprm wcs={ .flag=-1 };


  Vignetting* vignetting=NULL;
  FILE* progressfile=NULL;
  Attitude* ac=NULL;


  do { // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=exposure_map_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // init expo map
    float** expoMap=NULL;
    // Array for the calculation of the exposure map.
    init_expoMap(par,&expoMap,&status);
    CHECK_STATUS_BREAK(status);

    // see if we need the raw map (without vignetting)
    double rawMap = 0;
    float** rawExpoMap=NULL;
    if (strcmp(par.RawExposuremap,"NONE")!=0){
    	rawMap = 1;
        init_expoMap(par,&rawExpoMap,&status);
        CHECK_STATUS_BREAK(status);
    }


    // WCS data structure used for projection.
    init_expo_wcs(par,&wcs,&status);
    CHECK_STATUS_BREAK(status);

    // Initialize the random number generator. TODO: Do we need this?
    unsigned int seed=getSeed(par.seed);
    sixt_init_rng(seed, &status);
    CHECK_STATUS_BREAK(status);

    // Set the progress status output file.
    progressfile = get_progressout(par,&status);
    CHECK_STATUS_BREAK(status);

    // Set up the Attitude.
    ac = get_attitude(par,&status);
    CHECK_STATUS_BREAK(status);

    // Load the Vignetting data.
    vignetting = get_vign(par,&status);
    if (NULL==vignetting){
    	printf("Error: Vignetting file could not be loaded!!\n");
    }
    CHECK_STATUS_BREAK(status);

    // get xml array
    xmlarray xmls;
    xmls.n = -1;
    xmls.xmlarray = NULL;
    xmls.n = parse_string2array(par.XMLFile,&(xmls.xmlarray),&status);

    GenInst* inst[xmls.n];
    get_geninst_all(inst,xmls,seed,&status);
    CHECK_STATUS_BREAK(status);
    // check if xml files fir together
    CHECK_STATUS_BREAK(check_consist(inst,xmls.n));

    // Calculate the minimum cos-value for sources inside the FOV: (angle(x0,source) <= 1/2 * diameter)
    // (now we can choose the first one as we checked and all are the same)
    if (par.fov_diameter > 0){
    	inst[0]->tel->fov_diameter = par.fov_diameter; // we only use the [0] one in the following
    }
    double field_min_align = get_min_fov_align(inst[0], par);

    // ######## --- END of Initialization --- ######### //

    // --- Beginning of Exposure Map calculation
    headas_chat(3, "calculate the exposure map ...\n");

    // Simulation progress status (running from 0 to 100).
    unsigned int progress= init_progress(progressfile);
    
    // LOOP over the given time interval from TSTART to TSTART+timespan in steps of dt.
    // int intermaps=0;
    double time;
    for (time=par.TSTART; time<par.TSTART+par.timespan; time+=par.dt) {
      
   /**   // Check if an interim map should be saved now.
      if (par.intermaps>0) {
	if (time >= par.TSTART+intermaps*(par.timespan/par.intermaps)) {
	  // Construct the filename.
	  char filename[MAXFILENAME];
	  strncpy(filename, par.Exposuremap,
		  strlen(par.Exposuremap)-5);
	  filename[strlen(par.Exposuremap)-5]='\0';
	  char buffer[MAXFILENAME];
	  sprintf(buffer, "_%d.fits", intermaps);
	  strcat(filename, buffer);

	  // Save the interim map.
	  saveExpoMap(expoMap, filename, par.ra_bins, par.dec_bins,
		      &wcs, par.clobber, &status);
	  CHECK_STATUS_BREAK(status);

	  intermaps++;
	}
	}
	 **/
      // END of saving an interim map.


      // Determine the telescope pointing direction at the current time.
      struct Telescope telescope;
      getTelescopeAxes(ac, &telescope.nx, &telescope.ny, &telescope.nz,
    		     time, &status);
     CHECK_STATUS_BREAK(status);

      // Calculate the RA and DEC of the pointing direction.
      double telescope_ra, telescope_dec;
      calculate_ra_dec(telescope.nz, &telescope_ra, &telescope_dec);
      //printf("RA: %.4f Dec: %.4f \n",telescope_ra*180/M_PI, telescope_dec*180/M_PI);

      // Check if the specified field of the sky might be within the FOV.
      // Otherwise break this run and continue at the beginning of the loop 
      // with the next time step.
      Vector pixpos=unit_vector(0.5*(par.ra1+par.ra2), 0.5*(par.dec1+par.dec2));
      if (check_fov(&pixpos, &telescope.nz, field_min_align)!=0) {
    	  continue;
      }
     // printf("Current Pointing (t=%f):  %f %f %f\n",time,telescope.nz.x,telescope.nz.y,telescope.nz.z);

      // 2d Loop over the exposure map in order to determine all pixels that
      // are currently within the FOV.

      long x;
      float delta;
      for (x=0; x<par.ra_bins; x++) {
    	  long y;
    	  for (y=0; y<par.dec_bins; y++) {
    		  delta = get_single_expos_value(x,y,&wcs,
    				  telescope,inst,xmls.n,&status);
    		  if (delta>=0){
    			  expoMap[x][y]+=par.dt*get_Vignetting_Factor(vignetting, 1., delta, 0.);
    			  if (rawMap==1){
    				  rawExpoMap[x][y]+=par.dt;
    			  }
    		  }

    		  CHECK_STATUS_BREAK(status);
    	  }
     }

      // Program progress output.
      while((unsigned int)((time-par.TSTART)*100./par.timespan)>progress) {
    	  progress++;
    	  if (NULL==progressfile) {
    		  headas_chat(2, "\r%.0lf %%", progress*1.);
    		  fflush(NULL);
    	  } else {
    		  rewind(progressfile);
    		  fprintf(progressfile, "%.2lf", progress*1./100.);
    		  fflush(progressfile);
    	  }
      }
    } 


    CHECK_STATUS_BREAK(status);
    // END of LOOP over the specified time interval.


    // Progress output.
    if (NULL==progressfile) {
      headas_chat(2, "\r%.0lf %%\n", 100.);
      fflush(NULL);
    } else {
      rewind(progressfile);
      fprintf(progressfile, "%.2lf", 1.);
      fflush(progressfile);	
    }

    // END of generating the exposure map.



    // Store the exposure map in the output file.
    saveExpoMap(expoMap, par.Exposuremap, par.ra_bins, par.dec_bins,
		&wcs, par.clobber, &status);
    release_expoMap(&expoMap,par.ra_bins);

    // Write Raw Map?
    if (rawMap==1){
        saveExpoMap(rawExpoMap, par.RawExposuremap, par.ra_bins, par.dec_bins,
    		&wcs, par.clobber, &status);
        release_expoMap(&rawExpoMap,par.ra_bins);
        CHECK_STATUS_BREAK(status);
    }


    int ii;
    for (ii=0;ii<xmls.n;ii++){
    	destroyGenInst(&(inst[ii]),&status);
    }



  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");


  // Clean up the random number generator.
  sixt_destroy_rng();

  // Release memory.
  freeAttitude(&ac);
  destroyVignetting(&vignetting);
  wcsfree(&wcs);


  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int exposure_map_getpar(struct Parameters *par)
{
  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.
  
  query_simput_parameter_file_name("Attitude",&(par->Attitude), &status);
  query_simput_parameter_file_name("Vignetting",&(par->Vignetting), &status);
  query_simput_parameter_string("Exposuremap",&(par->Exposuremap), &status);
  query_simput_parameter_string("RawExposuremap",&(par->RawExposuremap), &status);


  query_simput_parameter_float("RA",&par->RA, &status);
  query_simput_parameter_float("Dec",&par->Dec, &status);



  // Read the diameter of the FOV (in arcmin).(Todo: do we really need this??)
  query_simput_parameter_float("fov_diameter",&(par->fov_diameter), &status);
  query_simput_parameter_file_name("XMLFile",&(par->XMLFile), &status);

  // Get the start time of the exposure map calculation.
  query_simput_parameter_double("TSTART",&par->TSTART, &status);
  // Get the timespan for the exposure map calculation.
  query_simput_parameter_double("timespan",&(par->timespan), &status);
  // Get the time step for the exposure map calculation.
  query_simput_parameter_double("dt",&(par->dt), &status);

  // Get the position of the desired section of the sky 
  // (right ascension and declination range).
  query_simput_parameter_double("ra1",&(par->ra1), &status);
  query_simput_parameter_double("ra2",&(par->ra2), &status);
  query_simput_parameter_double("dec1",&(par->dec1), &status);
  query_simput_parameter_double("dec2",&(par->dec2), &status);

  // Get the number of x- and y-bins for the exposure map.
  query_simput_parameter_int("ra_bins",&(par->ra_bins), &status);
  query_simput_parameter_int("dec_bins",&(par->dec_bins), &status);

  query_simput_parameter_int("projection",&(par->projection), &status);
  query_simput_parameter_int("intermaps",&(par->intermaps), &status);

  query_simput_parameter_int("seed",&(par->seed), &status);

  query_simput_parameter_string("ProgressFile",&(par->ProgressFile), &status);

  query_simput_parameter_bool("clobber",&par->clobber, &status);


  // Convert angles from [deg] to [rad].
  par->ra1 *=M_PI/180.;
  par->ra2 *=M_PI/180.;
  par->dec1*=M_PI/180.;
  par->dec2*=M_PI/180.;
  par->fov_diameter*=M_PI/180./60; // given in arcmin

  return(status);
}

