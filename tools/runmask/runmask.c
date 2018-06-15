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


   Copyright 2007 - 2018: Christian Schmid, Mirjam Oertel, FAU.
   Manuel Castro, National Institute for Space Research (INPE),
		 Brazil; under grant #2017/00968-6,
		 São Paulo Research Foundation (FAPESP).
*/

#include "runmask.h"

int runmask_main()
{

	const double timezero=0.0;

	// Program parameters
	struct Parameters par;

	// Instrument setup.
	GenInst* inst=NULL;

	//Coded-mask setup
	MaskSystem* mask_setup=NULL;	
	
	// Attitude.
	Attitude* ac=NULL;

	// Catalog of input X-ray sources.
	SourceCatalog* srccat=NULL;

	// Photon list file.
	PhotonFile* plf=NULL;

	// Error status.
	int status=EXIT_SUCCESS;  

	// Register HEATOOL.
	set_toolname("runmask");
	set_toolversion("0.05");
	
	do {

	// --- Initialization ---

	// Read the parameters using PIL
	status = runmask_getpar(&par);
	CHECK_STATUS_BREAK(status);

	// Initialize the random number generator.
	unsigned int seed = getSeed(par.Seed);
	sixt_init_rng(seed , &status);
	CHECK_STATUS_BREAK(status);

	// Determine the appropriate instrument XML definition file.
	char xml_filename[MAXFILENAME];
	sixt_get_XMLFile(xml_filename , par.XMLFile ,
			 par.Mission , par.Instrument , par.Mode ,
			 &status);
	CHECK_STATUS_BREAK(status);

	// Load the instrument configuration.
	inst = loadGenInst(xml_filename , seed , &status);
	CHECK_STATUS_BREAK(status);

	// Determine the appropriate advanced XML definition file.
	char xml_filename_adv[MAXFILENAME];
	sixt_get_XMLFile(xml_filename_adv , par.AdvXMLFile ,
			 par.Mission , par.Instrument , par.Mode ,
			 &status);
	CHECK_STATUS_BREAK(status);

	//Load coded-mask setup

	mask_setup = loadMaskSystem(xml_filename_adv,seed, &status);
	CHECK_STATUS_BREAK(status);

	} while (0);


//------------------------------------
//-------Execute each task 

	do{

	status = photogen(&par, inst);	
	CHECK_STATUS_BREAK(status);

	status = comaimg(&par, mask_setup);
	CHECK_STATUS_BREAK(status);
	
	status = comadet(&par, mask_setup);
	CHECK_STATUS_BREAK(status);

	status = comarecon(&par, mask_setup);
	CHECK_STATUS_BREAK(status);

	}while (0);

return status;

}


int runmask_getpar(struct Parameters* par)
{
	// String input buffer
	char* sbuffer = NULL;
	
	// Error status 
	int status = EXIT_SUCCESS;
	
	//Read all parameters via the aped_trad routines
	
	status = ape_trad_query_string("PhotonList" , &sbuffer);
	if (EXIT_SUCCESS != status) {
	   	SIXT_ERROR("failed reading the name of the photon list");
		return(status);
	}
	strcpy(par->PhotonList , sbuffer);
	free(sbuffer);

	status = ape_trad_query_string("ImpactList" , &sbuffer);
	if (EXIT_SUCCESS != status) {
	   	SIXT_ERROR("failed reading the name of the impact list");
		return(status);
	}
	strcpy(par->ImpactList , sbuffer);
	free(sbuffer);

	status = ape_trad_query_string("EventList" , &sbuffer);
	if (EXIT_SUCCESS != status) {
	   	SIXT_ERROR("failed reading the name of the event list");
		return(status);
	}
	strcpy(par->EventList , sbuffer);
	free(sbuffer);

	status = ape_trad_query_string("Image" , &sbuffer);
	if (EXIT_SUCCESS != status) {
	   	SIXT_ERROR("failed reading the name of the image file");
		return(status);
	}
	strcpy(par->Image , sbuffer);
	free(sbuffer);

	status = ape_trad_query_string("PositionList" , &sbuffer);
	if (EXIT_SUCCESS != status) {
	   	SIXT_ERROR("failed reading the name of the position list");
		return(status);
	}
	strcpy(par->PositionList , sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("Mission" , &sbuffer);
	if (EXIT_SUCCESS != status) {
		SIXT_ERROR("failed reading the name of the mission");
		return(status);
	}
	strcpy(par->Mission , sbuffer);
	free(sbuffer);

	status = ape_trad_query_string("Instrument" , &sbuffer);
	if (EXIT_SUCCESS != status) {
		SIXT_ERROR("failed reading the name of the instrument");	
		return(status);
	}
	strcpy(par->Instrument , sbuffer);
	free(sbuffer);

	status = ape_trad_query_string("Mode" , &sbuffer);
	if (EXIT_SUCCESS != status) {
		SIXT_ERROR("failed reading the name of the instrument mode");
		return(status);
	}
	strcpy(par->Mode , sbuffer);
	free(sbuffer);
	
	status=ape_trad_query_string("XMLFile" , &sbuffer);
  	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the XML file");
		return(status);
  	} 
	strcpy(par->XMLFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("AdvXMLFile" , &sbuffer);
  	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the advanced XML file");
		return(status);
  	} 
	strcpy(par->AdvXMLFile, sbuffer);
	free(sbuffer);
	
	status = ape_trad_query_string("Attitude" , &sbuffer);
	if (EXIT_SUCCESS != status) {
		 SIXT_ERROR("failed reading the name of the attitude");
		return(status);
	} 
	strcpy(par->Attitude , sbuffer);
	free(sbuffer);

	status = ape_trad_query_double("RA" , &par->RA);
	if (EXIT_SUCCESS != status) {
		SIXT_ERROR("failed reading the right ascension of the telescope pointing");
		return(status);
	} 

	status = ape_trad_query_double("DEC" , &par->DEC);
	if (EXIT_SUCCESS != status) {
		SIXT_ERROR("failed reading the declination of the telescope pointing");
		return(status);
	} 

	status = ape_trad_query_file_name("Simput" , &sbuffer);
	if (EXIT_SUCCESS != status) {
		SIXT_ERROR("failed reading the name of the SIMPUT file");
		return(status);
	} 
	strcpy(par->Simput , sbuffer);
	free(sbuffer);

	status = ape_trad_query_double("MJDREF" , &par->MJDREF);
	if (EXIT_SUCCESS != status) {
		SIXT_ERROR("failed reading MJDREF");
		return(status);
	} 

	status = ape_trad_query_double("dt" , &par->dt);
	if (EXIT_SUCCESS != status) {
		SIXT_ERROR("failed reading dt");
		return(status);
	 } 

	 status = ape_trad_query_double("TSTART" , &par->TSTART);
	 if (EXIT_SUCCESS != status) {
	 	SIXT_ERROR("failed reading TSTART");
		return(status);
	 } 

	status = ape_trad_query_double("Exposure", &par->Exposure);
	if (EXIT_SUCCESS != status) {
		SIXT_ERROR("failed reading the exposure time");
		return(status);
	 } 

	status = ape_trad_query_int("seed" , &par->Seed);
	if (EXIT_SUCCESS != status) {
		SIXT_ERROR("failed reading the seed for the random number generator");
		return(status);
	}

	status = ape_trad_query_bool("clobber" , &par->clobber);
	if (EXIT_SUCCESS != status) {
		SIXT_ERROR("failed reading the clobber parameter");
		return(status);
	}

  return(status);


}

int photogen(struct Parameters* par, GenInst* inst)
{


	const double timezero=0.0;

	// Attitude.
	Attitude* ac=NULL;

	// Catalog of input X-ray sources.
	SourceCatalog* srccat=NULL;

	// Photon list file.
	PhotonFile* plf=NULL;

	// Error status.
	int status=EXIT_SUCCESS;  



	//--------checking 
	printf("PhotonList= %s \n",par->PhotonList);
	printf("ImpactList= %s \n",par->ImpactList);
	printf("XMLFile   = %s \n",par->XMLFile);
	printf("EventList = %s \n",par->EventList);
	printf("Image	  = %s \n",par->Image);
	printf("PositionList= %s \n",par->PositionList);
	printf("advanced XML= %s \n",par->AdvXMLFile);


	do { // Beginning the ERROR HANDLING loop



	headas_chat(3 , "initialize photon generation ...\n");

	// Determine the photon list output file.
	char photonlist_filename[MAXFILENAME];
	strcpy(photonlist_filename , par->PhotonList);

	// Check

	printf("RA = %2.5f \n",par->RA);
	printf("DEC = %2.5f \n",par->DEC);
	printf("Exposure = %2.5f \n",par->Exposure);


	// Set up the Attitude.
	char ucase_buffer[MAXFILENAME];
	strcpy(ucase_buffer , par->Attitude);
	strtoupper(ucase_buffer);
	if (0 == strcmp(ucase_buffer , "NONE")) {
	// Set up a simple pointing attitude.
	ac = getPointingAttitude(par->MJDREF , par->TSTART , par->TSTART+par->Exposure ,
	                         par->RA*M_PI/180. , par->DEC*M_PI/180. , &status);
	CHECK_STATUS_BREAK(status);

	} else {
	// Load the attitude from the given file.
	ac = loadAttitude(par->Attitude , &status);
	CHECK_STATUS_BREAK(status);

	// Check if the required time interval for the simulation
	// is a subset of the period covered by the attitude file.
	checkAttitudeTimeCoverage(ac , par->MJDREF , par->TSTART ,
				 par->TSTART+par->Exposure, &status);
	CHECK_STATUS_BREAK(status);
	}
	// END of setting up the attitude.

	// Load the SIMPUT X-ray source catalog.
	srccat=loadSourceCatalog(par->Simput , inst->tel->arf , &status);
	CHECK_STATUS_BREAK(status);

// --- End of Initialization ---


// --- Photon Generation Process ---

    // Open the output photon list file.
	char telescop[MAXMSG]={""};
	char instrume[MAXMSG]={""};
	if (NULL!=inst->telescop) {
	strcpy(telescop, inst->telescop);
	}
	if (NULL!=inst->instrume) {
	strcpy(instrume, inst->instrume);
	}
	plf=openNewPhotonFile(photonlist_filename,
			  telescop, instrume, "Normal",
			  inst->tel->arf_filename, inst->det->rmf_filename,
			  par->MJDREF, timezero, par->TSTART, par->TSTART+par->Exposure,
			  par->clobber, &status);
	CHECK_STATUS_BREAK(status);

	// Set FITS header keywords.
	fits_update_key(plf->fptr, TSTRING, "ATTITUDE", par->Attitude,
		    "attitude file", &status);

	// Start the actual photon generation (after loading required data):
	headas_chat(3, "start photon generation process ...\n");

	// Loop over photon generation, till the time of the photon exceeds 
	// the requested exposure time.
	// Simulation progress status (running from 0 to 1000).
	int progress=0;
	do {

	// Photon generation.
	Photon ph;
	int isph=phgen(ac, &srccat, 1, 
		     par->TSTART, par->TSTART+par->Exposure, 
		     par->MJDREF, par->dt, 
		     inst->tel->fov_diameter, &ph, &status);
	CHECK_STATUS_BREAK(status);

	// If no photon has been generated, break the loop.
	if (0==isph) break;

	// Check if the photon still is within the requested exposre time.
	if (ph.time>par->TSTART+par->Exposure) break;

	// Write the photon to the output file.
	status=addPhoton2File(plf, &ph);
	CHECK_STATUS_BREAK(status);

	// Program progress output.
	while ((int)((ph.time-par->TSTART)*1000./par->Exposure)>progress) {
	progress++;
	headas_chat(2, "\r%.1lf %%", progress*1./10.);
	fflush(NULL);
	}

	} while(1); // END of photon generation loop.
	CHECK_STATUS_BREAK(status);

	// Progress output.
	headas_chat(2, "\r%.1lf %%\n", 100.);
	fflush(NULL);

// --- End of photon generation ---

	} while (0); // End of ERROR HANDLING loop


// --- Clean up ---

	headas_chat(3, "\ncleaning up ...\n");

	// Release memory.
	freePhotonFile(&plf, &status);
	freeSourceCatalog(&srccat, &status);
	freeAttitude(&ac);
	destroyGenInst(&inst, &status);

	// Clean up the random number generator.
	sixt_destroy_rng();

	if (EXIT_SUCCESS==status) {
		headas_chat(3, "photon generation finished successfully!\n\n");
		return(EXIT_SUCCESS);
	} else {
		return(EXIT_FAILURE);
	}

}

int comaimg(struct Parameters* par,  MaskSystem* mask_setup)
{

  struct Telescope telescope; //Telescope coordinate system
  PhotonFile* plf=NULL;
  ImpactFile* ilf=NULL;
  CodedMask* mask=NULL;
  Attitude* ac=NULL;

  double refxcrvl=0., refycrvl=0.;

  int status=EXIT_SUCCESS; // Error status.


  do { //Beginning of the ERROR handling loop (will at most be run once)

    headas_chat(3 , "initialize imaging ...\n");
     
    //Detector-setup:
    float x_det = mask_setup->x_det;
    float y_det = mask_setup->y_det;
    //Distance mask-detector:
    float distance = mask_setup->mask_distance;
    float det_pixelwidth= mask_setup->det_pixelwidth;
    double ra=par->RA;
    double dec=par->DEC;
    //Distance between the mask plane and the collimator top plane
    float coll_distance = mask_setup->mask_distance - mask_setup->collimator_height;
    //Width of one detector-element that is surrounded by the collimator
    float det_width = mask_setup->det_width; 
    // Initialize the random number generator.
    sixt_init_rng(getSeed(-1),&status);
    CHECK_STATUS_BREAK(status);

    //Open the FITS file with the input photon list.
    plf=openPhotonFile(par->PhotonList,READONLY,&status);
    if (EXIT_SUCCESS!=status) break;

    //Load the coded mask from the file.
    mask=getCodedMaskFromFile(mask_setup->mask_pattern,&status);
    if (EXIT_SUCCESS!=status) break;


    //Calculate min cos-value for sources inside FOV.
    float phi_max=det_phi_max(distance,mask_setup->x_mask,mask_setup->y_mask,x_det,y_det);
    float fov_min_align = cos(phi_max);
    //Set up the telescope attitude:

    //Buffer for attitude-filename.
    char att_buffer[MAXFILENAME];
    //Copy attitude-filename from par-file to buffer.
    strcpy(att_buffer, par->Attitude);
    //Make all letters capital to compare with any spelling of 'none'.
    strtoupper(att_buffer);

    if (0==strcmp(att_buffer, "NONE")){
      //Set up attitude(no attitude-file to use).

      //Memory-allocation:

      //Allocates memory for struct Attitude.
      //Initializes nentries, current_entry to 0;entry to NULL.
      ac=getAttitude(&status);
      CHECK_STATUS_BREAK(status);

      //Allocates memory for struct AttitudeEntry.
      ac->entry=(AttitudeEntry*)malloc(sizeof(AttitudeEntry));
      if (NULL==ac->entry) {
	status=EXIT_FAILURE;
	SIXT_ERROR("memory allocation for AttitudeEntry failed");
	break;
      }

      //Set the values of the AttitudeEntry.

      ac->nentries=1;
      ac->entry[0]=initializeAttitudeEntry();  
      //ac->entry[0].time=0;

      //Telescope pointing direction:
      ac->entry[0].nz=normalize_vector(unit_vector(par->RA*M_PI/180.0,par->DEC* M_PI/180.0));
      //Unit-vector in z-direction:
      Vector vz = {0.,0.,1.};
      //Vector perpendicular to nz-vz-plane:
      ac->entry[0].nx=vector_product(vz,ac->entry[0].nz);
      //initialize telescope coordinate sytem
      telescope.nz=ac->entry[0].nz;
      telescope.nx=ac->entry[0].nx;

    } else {
      //Load the attitude from file.
      ac=loadAttitude(par->Attitude, &status);
      CHECK_STATUS_BREAK(status);

      Vector initializey={0.,0.,0.};
      telescope.ny =initializey;
      telescope.nx=ac->entry[0].nx;

      //Check if the required time interval for the simulation
      //lies within the time interval in the attitude-file
      /* if((ac->entry[0].time > 0) || ((ac->entry[ac->nentries-1].time) < (par.Exposure + par.Timezero))){
	status=EXIT_FAILURE;
	SIXT_ERROR("attitude data does not cover the specified period");
	break;
	}*/
      
    }//END of setting up the attitude.


    // Read header keywords.
    char telescop[MAXMSG], instrume[MAXMSG], filter[MAXMSG];
    char ancrfile[MAXMSG], respfile[MAXMSG];
    char comment[MAXMSG];
    fits_read_key(plf->fptr, TSTRING, "TELESCOP", &telescop, comment, &status);
    fits_read_key(plf->fptr, TSTRING, "INSTRUME", &instrume, comment, &status);
    fits_read_key(plf->fptr, TSTRING, "FILTER", &filter, comment, &status);
    fits_read_key(plf->fptr, TSTRING, "ANCRFILE", &ancrfile, comment, &status);
    fits_read_key(plf->fptr, TSTRING, "RESPFILE", &respfile, comment, &status);
    CHECK_STATUS_BREAK(status);

    double mjdref, timezero, tstart, tstop;
    fits_read_key(plf->fptr, TDOUBLE, "MJDREF", &mjdref, comment, &status);
    fits_read_key(plf->fptr, TDOUBLE, "TIMEZERO", &timezero, comment, &status);
    fits_read_key(plf->fptr, TDOUBLE, "TSTART", &tstart, comment, &status);
    fits_read_key(plf->fptr, TDOUBLE, "TSTOP", &tstop, comment, &status);
    CHECK_STATUS_BREAK(status);

    // Create a new FITS file for the output of the impact list.
    ilf=openNewImpactFile(par->ImpactList, 
			  telescop, instrume, filter,
			  ancrfile, respfile,
			  mjdref, timezero, tstart, tstop,
			  par->clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Write WCS header keywords.
    fits_update_key(ilf->fptr, TDOUBLE, "REFXCRVL", &refxcrvl, "", &status);
    fits_update_key(ilf->fptr, TDOUBLE, "REFYCRVL", &refycrvl, "", &status);
    CHECK_STATUS_BREAK(status);

    //initialization of wcs parameter structure for determining the impact position
    struct wcsprm wcs = {
      .flag=-1
    }; //flag has to be set only at 1st init
    if (0!=wcsini(1, 2, &wcs)) {
      SIXT_ERROR("initalization of WCS data structure failed");
      status=EXIT_FAILURE;
      break;
    }
    wcs.naxis=2;
    wcs.crpix[0]=(x_det/det_pixelwidth)/2;
    wcs.crpix[1]=(y_det/det_pixelwidth)/2;
    wcs.crval[0]=ra;  //in deg
    wcs.crval[1]=dec;
    wcs.cdelt[0]=atan(det_pixelwidth/distance)*180./M_PI;
    wcs.cdelt[1]=atan(det_pixelwidth/distance)*180./M_PI;
    strcpy(wcs.cunit[0], "deg");
    strcpy(wcs.cunit[1], "deg");
    strcpy(wcs.ctype[0], "RA---TAN");
    strcpy(wcs.ctype[1], "DEC--TAN");

    //initialization of wcs parameter structure for determining the impact position at collimator height
    //for protoMirax

    struct wcsprm wcs2 = {
      .flag=-1
    }; //flag has to be set only at 1st init
    if (0!=wcsini(1, 2, &wcs2)) {
      SIXT_ERROR("initalization of WCS2 data structure failed");
      status=EXIT_FAILURE;
      break;
    }
    wcs2.naxis=2;
    wcs2.crpix[0]=(x_det/det_pixelwidth)/2;  //in pixels
    wcs2.crpix[1]=(y_det/det_pixelwidth)/2;
    wcs2.crval[0]=ra;  //in deg
    wcs2.crval[1]=dec;
    wcs2.cdelt[0]=atan(det_pixelwidth/coll_distance)*180./M_PI;
    wcs2.cdelt[1]=atan(det_pixelwidth/coll_distance)*180./M_PI;
    strcpy(wcs2.cunit[0], "deg");
    strcpy(wcs2.cunit[1], "deg");
    strcpy(wcs2.ctype[0], "RA---TAN");
    strcpy(wcs2.ctype[1], "DEC--TAN");
    // --- END of Initialization ---

    // --- Beginning of Imaging Process ---

    //Beginning of actual simulation (after loading required data).
    headas_chat(3, "start imaging process ...\n");

    //SCAN PHOTON LIST (starts with first entry, plf initialized to 0 at beginning)
    while (plf->row < plf->nrows) {
         
      //Read an entry from the photon list:
      Photon photon={.time= 0.0};
      status=PhotonFile_getNextRow(plf,&photon);
      if (EXIT_SUCCESS!=status) break;

      //Check whether photon-list-entry is within requested time interval.
      if (photon.time > (par->Exposure + par->Timezero))
	break;

      //Determine unit vector in photon-direction
      Vector phodir=normalize_vector(unit_vector(photon.ra,photon.dec));

      //Determine current telescope pointing direction.
      telescope.nz=getTelescopeNz(ac,photon.time,&status);
      CHECK_STATUS_BREAK(status);
      
      //Check whether photon is inside FOV:
      //Compare photon direction to direction of telescope axis
      if (check_fov(&phodir,&telescope.nz,fov_min_align)==0){
	//Photon is inside fov

	getTelescopeAxes(ac,&telescope.nx,&telescope.ny,&telescope.nz,photon.time,&status);

	if (0!=strcmp(att_buffer, "NONE")){ //attitude-file is used
	  //set wcs according to new pointing (attitude could have changed)
	  setWCScurrentPointing(par->Attitude,ac,&telescope.nz,&wcs,&status);
	}

	//Determine photon impact position on detector in [m].
	struct Point2d position;
	
	//first:impact position in mask-plane (transparent pixels only);
	//if photon then hits the detector (and not the walls), return value is 1, 0 else

	// If flag == 0 ---> mirax, flag == 1 ---> protoMirax

	int reval;	
	if (mask_setup->flag ==1 ){

	reval=getImpactPos_protoMirax(&wcs,&wcs2,&position,mask,photon.ra*180./M_PI,
					  photon.dec*180./M_PI,det_pixelwidth,det_width,
					  x_det,y_det,mask_setup->wall_thickness,&status);

	} else {
	reval=getImpactPos_wcs(&wcs,&position,mask,photon.ra*180./M_PI,
				photon.dec*180./M_PI,det_pixelwidth,x_det,y_det,&status);
	CHECK_STATUS_BREAK(status);
	}


	if (reval == 1){
	 
	  //Create new impact.
	  Impact impact;
	  impact.time       = photon.time;
	  impact.energy     = photon.energy;
	  impact.position.x = position.x;
	  impact.position.y = position.y;
	  impact.ph_id      = photon.ph_id;
	  impact.src_id     = photon.src_id;

	  //Write to output-file.
	  addImpact2File(ilf,&impact,&status);
	  CHECK_STATUS_BREAK(status);
	} // END of photon hits the detector.
      } // END of photon  inside fov.
    } // END of scanning LOOP over the photon list.
    if (EXIT_SUCCESS!=status) break;
  
    // Release memory.
    destroyCodedMask(&mask);
    wcsfree(&wcs);
  } while(0);  // END of the error handling loop.

  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the FITS files.
  freeImpactFile(&ilf,&status);
  freePhotonFile(&plf,&status);
  freeAttitude(&ac);

  // Clean up the random number generator.
  sixt_destroy_rng();

  if (EXIT_SUCCESS==status) headas_chat(3, "comaimg finished successfully!\n\n");
  

return(status);
}

int comadet(struct Parameters* par,  MaskSystem* mask_setup)
{


  ImpactFile* ilf=NULL;
  CoMaDetector* detector=NULL;

  //Error status.
  int status=EXIT_SUCCESS;


  do {  //Beginning of the ERROR handling loop (will at most be run once)

    headas_chat(3 , "initialize comadet ...\n");
   
  //Open the impact list FITS file.
    ilf=openImpactFile(par->ImpactList, READONLY, &status);
    CHECK_STATUS_RET(status, status);
    
    //Set the event list template file:
    strcpy(par->EventListTemplate, SIXT_DATA_PATH);
    strcat(par->EventListTemplate, "/templates/coma.eventlist.tpl");

    //DETECTOR setup.
    //initializes from par-file
    struct CoMaDetectorParameters cdp = {
      .pixels = 
      { .xwidth = mask_setup->width,          //detector-width [pixel]
	.ywidth = mask_setup->width,
	.DCU_length = mask_setup->DCU_length, //length of DCU [m]
	.DCU_gap = mask_setup->DCU_gap,       //gap between 2 DCU's [m]
	.DCA_gap = mask_setup->DCA_gap,       //gap between 2 DCA's [m]
	.xpixelwidth = mask_setup->pixelwidth,//width of one pixel [m]
	.ypixelwidth = mask_setup->pixelwidth 
      },
      .eventfile_filename=par->EventList,
      .eventfile_template=par->EventListTemplate
    };

    //create new CoMaDetector-object:
    detector=getCoMaDetector(&cdp, &status);
    CHECK_STATUS_RET(status, status);

    //END of DETECTOR CONFIGURATION SETUP    

    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---


    //Beginning of actual detector simulation (after loading required data):
    headas_chat(5, "start detection process ...\n");
    Impact impact;

    //Loop over all impacts in the FITS file.
    while (ilf->row<ilf->nrows) {

      getNextImpactFromFile(ilf, &impact, &status);
      CHECK_STATUS_RET(status, status);

      //detection routine:determines affected pixel and adds new event
      //to the event-file
	
     //flag if protomirax=1, mirax =0
	
      if(mask_setup->flag==1){
	status=addImpact2CoMaDetector_protoMirax(detector, &impact);
	CHECK_STATUS_RET(status, status);
      }else{
	status=addImpact2CoMaDetector(detector, &impact);
	CHECK_STATUS_RET(status, status);
      }
      
    } //END of scanning the impact list.
    CHECK_STATUS_RET(status, status);


  } while(0);  // END of the error handling loop.


  // --- cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  //Free the Coded Mask Detector.
  freeCoMaDetector(detector);

  //Close the FITS files.
  freeImpactFile(&ilf, &status);

  if (EXIT_SUCCESS==status) headas_chat(3, "comadet finished successfully!\n\n");
  return(status);


}

int comarecon(struct Parameters* par, MaskSystem* mask_setup)
{


  CoMaEventFile* eventfile=NULL;
  SquarePixels* detector_pixels=NULL;
  CodedMask* mask=NULL; 
  SourceImage* sky_pixels=NULL;
  ReconArray* recon=NULL;
  MaskShadow* mask_shadow=NULL;
  PixPositionList* position_list=NULL;
  double* median_list=NULL; //temp array of all background pix for determination of median 
  ReadEvent* ea=NULL;
  ReadEvent* ear=NULL;
  double* ReconImage1d=NULL;
  double* EventImage1d=NULL;
  fftw_complex* fftReconArray=NULL;
  fftw_complex* fftEventArray=NULL;
  fftw_complex* Multiply = NULL;
  fftw_complex* fftInvMultiply=NULL;

  int status=EXIT_SUCCESS; // Error status.

  int ii, jj; //counts
  int xdiff, ydiff; //difference in size of mask and det plane
  int PixAmount;



  do {  // Beginning of the ERROR handling loop (will at most be run once)


    // Open the event file.
    eventfile=openCoMaEventFile(par->EventList, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Load the coded mask from the file.
    mask=getCodedMaskFromFile(mask_setup->mask_pattern, &status);
    CHECK_STATUS_BREAK(status);
    
    // DETECTOR setup.
    double ra=par->RA;
    double dec=par->DEC;
    double distance=mask_setup->mask_distance;

    double xdetsize_beforeRepix=mask_setup->width; //in the case of Repix the detector_pixels are overwritten, 
    double ydetsize_beforeRepix=mask_setup->width; //but some functions need the old smaller size (of the det in pixels)
    double xpixelsize_beforeRepix=mask_setup->pixelwidth; //original detector-pixelsize
    double ypixelsize_beforeRepix=mask_setup->pixelwidth;

    struct SquarePixelsParameters spp = {
      .xwidth = mask_setup->width,
      .ywidth = mask_setup->width,
      .xpixelwidth = mask_setup->pixelwidth,
      .ypixelwidth = mask_setup->pixelwidth,
      .DCU_length=mask_setup->DCU_length,
      .DCU_gap=mask_setup->DCU_gap,
      .DCA_gap=mask_setup->DCA_gap
    };
    detector_pixels=newSquarePixels(&spp, &status);
    CHECK_STATUS_BREAK(status);
    // END of DETECTOR CONFIGURATION SETUP    


    // SKY IMAGE setup.
       if(mask_setup->repixsize==0.){
	 PixAmount=4;

	 float delta = atan(mask_setup->pixelwidth/distance);
	 struct SourceImageParameters sip = {
	   .naxis1 = 2*(mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth),
	   .naxis2 = 2*(mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth),
	   .crpix1 = (float)(sip.naxis1/2.+1), //even axes:ends to .5;odd axes (not possible,*2):ends to .0 (therefore: axis/2+0.5)
	   .crpix2 = (float)(sip.naxis2/2.+1), //further +0.5,since left boarder of 1st pix in sky image is 0.5(FITS)
	   .cdelt1 = delta,
	   .cdelt2 = delta,
	   .crval1 = ra*M_PI/180.,
	   .crval2 = dec*M_PI/180.
	};
       sky_pixels=getEmptySourceImage(&sip, &status);
       CHECK_STATUS_BREAK(status);
       }else{//only works, if new smaller pixel fit without remainder in former big ones!
	 //EventArray has been built with original pix-size, in order to distribute the events correctly
	 //from now on, the new smaller size is used -> also for the later called ReconArray
	 
	 detector_pixels->xwidth=(detector_pixels->xwidth*detector_pixels->xpixelwidth)/mask_setup->repixsize;
	 detector_pixels->ywidth=(detector_pixels->ywidth*detector_pixels->ypixelwidth)/mask_setup->repixsize;
	 detector_pixels->xpixelwidth=mask_setup->repixsize;
	 detector_pixels->ypixelwidth=mask_setup->repixsize;
	 
	 PixAmount=24;

	 //get source image with correct axes, since xpixelwidth has changed
	 float delta = atan(mask_setup->repixsize/distance);
	 struct SourceImageParameters sip = {
	   .naxis1 = 2*(mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth),
	   .naxis2 = 2*(mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth),
	   //due to repix the axes get even numbers -> have to be shifted half a pixel
	   //since one former pixel (for MIRAX) corresponds to 6 now -> shift: 4.0
	   .crpix1 = (float)(sip.naxis1/2.+1),
	   .crpix2 = (float)(sip.naxis2/2.+1),
	   .cdelt1 = delta,
	   .cdelt2 = delta,
	   .crval1 = ra*M_PI/180.,
	   .crval2 = dec*M_PI/180.
	 };
	sky_pixels=getEmptySourceImage(&sip, &status);
        CHECK_STATUS_BREAK(status);

	}
    // END of SKY IMAGE CONFIGURATION SETUP 

    //SOURCE-POSITION DETERMINATION:
    double pixval=0.;
    int threshold=0; //for first source threshold is '0' after that: '1' if still bright enough, '2' else
    //get empty PixPositionList structure (contains pointer to current PixPosition-element
    //and count for found sources)
    position_list=getPixPositionList(sky_pixels);
    //memory-allocation for median_list
    median_list=getMedian_list(sky_pixels, &status);

    //initialization of wcs parameter structure for deteremining ra/dec of source
    struct wcsprm wcs = {
      .flag=-1
    }; //flag has to be set only at 1st init
    if (0!=wcsini(1, 2, &wcs)) {
      SIXT_ERROR("initalization of WCS data structure failed");
      status=EXIT_FAILURE;
      break;
    }
    wcs.naxis=2;
    wcs.crpix[0]=sky_pixels->crpix1;
    wcs.crpix[1]=sky_pixels->crpix2; 
    wcs.crval[0]=sky_pixels->crval1*180./M_PI;
    wcs.crval[1]=sky_pixels->crval2*180./M_PI;
    wcs.cdelt[0]=sky_pixels->cdelt1*180./M_PI; //in deg
    wcs.cdelt[1]=sky_pixels->cdelt2*180./M_PI;

    //initialization of wcs parameter structure for getting mask shadow
    struct wcsprm wcs2 = {
      .flag=-1
    }; //flag has to be set only at 1st init
    if (0!=wcsini(1, 2, &wcs2)) {
      SIXT_ERROR("initalization of WCS data structure failed");
      status=EXIT_FAILURE;
      break;
    }
    wcs2.naxis=2;
    wcs2.crpix[0]=(detector_pixels->xwidth)/2;
    wcs2.crpix[1]=(detector_pixels->ywidth)/2;
    wcs2.crval[0]=ra;  //in deg
    wcs2.crval[1]=dec;
    wcs2.cdelt[0]=atan(detector_pixels->xpixelwidth/distance)*180./M_PI; //in deg
    wcs2.cdelt[1]=atan(detector_pixels->ypixelwidth/distance)*180./M_PI;


    //telescope coordinate system     //TODO: USE ATTITUDE
       //telescope pointing direction
       //Vector nz=normalize_vector(unit_vector(ra*M_PI/180.0,dec* M_PI/180.0));
       //unit-vector in z-direction:
       //Vector vz = {0.,0.,1.};

       // Vector nx= normalize_vector(vector_product(nz,vz));   
       // Vector ny= normalize_vector(vector_product(nz,nx)); 
    
    // --- END of Initialization ---


    // --- Beginning of Reconstruction Process ---

    // Beginning of actual detector simulation (after loading required data):
    headas_chat(3, "start image reconstruction process ...\n");

    //Get empty event array object (type: ReadEvent)
       //sizes are equal to those of the ReconArray -> mask-width but detector-pixelsize
       //mask has to be >= detector
    int ea_size1=2*(mask->naxis1*mask->cdelt1/xpixelsize_beforeRepix);
    int ea_size2=2*(mask->naxis2*mask->cdelt2/ypixelsize_beforeRepix);
    ea=getEventArray(ea_size1,ea_size2,&status);
     
       //detector size <= mask size
       //if mask size = det size -> shift is zero
       xdiff=((ea->naxis1)/2-xdetsize_beforeRepix)/2;
       ydiff=((ea->naxis2)/2-ydetsize_beforeRepix)/2;

       // Loop over all events in the FITS file.
       while (0==EventListEOF(&eventfile->generic)) {

	 status=readEventList_nextRow(eventfile, ea);
	 CHECK_STATUS_BREAK(status);

	 //Get the 2d-EventArray, padded to the upper right corner
	 ea->EventArray[ea->rawx+ea->naxis1/2+xdiff][ea->rawy+ea->naxis2/2+ydiff]+=ea->charge;

       } // END of scanning the event list.
       CHECK_STATUS_BREAK(status);
       
       //createTestImg(&ea->EventArray,2,detector_pixels->xwidth,detector_pixels->ywidth,
       // ea->naxis1/2+xdiff,ea->naxis2/2+ydiff,"eventArray.fits",&status);
       

       if(mask_setup->repixsize!=0.){
	 double pixelwidth_big=xpixelsize_beforeRepix;//width of the former EventArray pixels (the real det-pix-size)

	 //new pointer for EventArray with more entries, since smaller pix-size (EventArrayRepix)
	 int ear_size1=2*(mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth);
	 int ear_size2=2*(mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth);

	 ear=getEventArray(ear_size1, ear_size2, &status);

	 repixNoReminder(ea,ear,2,ea->naxis1,ea->naxis2,pixelwidth_big,detector_pixels->xpixelwidth);

	 ReadEvent* ea_temp=NULL;
	 ea_temp=ea;
	 FreeEventArray(ea_temp);
	 ea=ear; //set pointer to EventArray-data to just re-pixeled array

	 xdiff=((ea->naxis1)/2-detector_pixels->xwidth)/2;
	 ydiff=((ea->naxis2)/2-detector_pixels->ywidth)/2;
       }//end re-pixel EventArray to smaller size given by RePixSize

       //Get the reconstruction array:
       //type: 1: balanced cross correlation (rnd pattern), 2: MURA
       /*int type;
       if(par.protoMirax == 1){
	 type=2;
       }else{
	 type=1;
	 }*/ //TODO
       
       recon=getReconArray(mask,2,detector_pixels,&status);
       int Size1 = recon->naxis1;
       int Size2 = recon->naxis2;

        //Get the 1d image of the reconstruction array -> needed by FFTW
       ReconImage1d=SaveReconArray1d(recon, &status);  
      
       //perform a fft with the ReconArray       
       fftReconArray=FFTOfArray_1d(ReconImage1d, Size1, Size2, -1);

       //get repixeled mask from ReconArray, which is needed later for building the mask shadow during IROS
       //basic constructor for both,the whole re-pixeled mask&/shadow element
       mask_shadow=getMaskShadowElement(Size1/2, Size2/2, Size1, Size2, &status); 
       //gets re-pixeled mask as big as EventArray with values betw. 0...1
       getMaskRepix(recon, mask_shadow, 1);
  
       do{ //search for sources as long as pixval is above certain value
	 //run as long as threshold==1

	 //Get the 1d image of the event array -> needed by FFTW
	 EventImage1d=SaveEventArray1d(ea, &status);
	 //Check whether the ReconArray and the EventArray have the same size
	 if ((recon->naxis1 != ea->naxis1) || (recon->naxis2 != ea->naxis2)){
	   printf ("Error: ReconArrray and EventArray must have the same size!\n");
	   break;
	 }
   
	 //perform a fft with the EventArray       
	 fftEventArray=FFTOfArray_1d(EventImage1d, Size1, Size2, -1);
    
       //multiply fftEventArray with komplex conjugate of fftReconArray
       //Re-part: E(Re)*R(Re)+E(Im)*R(Im); Im-part: E(Re)*R(Im)-E(Im)*R(Re)       
       Multiply=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Size1*Size2));
       for(ii=0; ii<Size1; ii++){
	 for(jj=0; jj<Size2; jj++){
	   Multiply[ii+Size1*jj][0]=fftEventArray[ii+Size1*jj][0]*fftReconArray[ii+Size1*jj][0]
	     +fftEventArray[ii+Size1*jj][1]*(fftReconArray[ii+Size1*jj][1]);
	   Multiply[ii+Size1*jj][1]=-fftEventArray[ii+Size1*jj][0]*(fftReconArray[ii+Size1*jj][1])
	     +fftEventArray[ii+Size1*jj][1]*fftReconArray[ii+Size1*jj][0];
	 }
       }

       //Inverse FFT of Multilpy which already is of type fftw_complex       
       fftInvMultiply=FFTOfArray(Multiply, Size1, Size2, +1);

       //save real part of inverse fft in sky image
       for(ii=0; ii<Size1; ii++){
	 for(jj=0; jj<Size2; jj++){
	   sky_pixels->pixel[ii][jj]=fftInvMultiply[ii+Size1*jj][0]/(Size1*Size2);
	 }
       }

       //for testing:
       char name_image[MAXFILENAME];
       sprintf(name_image,"image_%lu", position_list->entryCount);

       // Write the reconstructed source function to the output FITS file.
       if(position_list->entryCount <=2){
       saveSourceImage(sky_pixels, name_image, &status);
       CHECK_STATUS_BREAK(status);
       }
   
       //finds current brightest pixel coordinates and saves PixPosition; returns current brightest pixval
       pixval=findBrightestPix(threshold, PixAmount, sky_pixels, pixval, position_list, &wcs, &status);
       threshold=getThresholdForSources(pixval, position_list, sky_pixels, median_list, par->Sigma);

       //get mask shadow for current source
       getMaskShadow2(mask_shadow,&wcs2,position_list,sky_pixels->crpix1,sky_pixels->crpix2,detector_pixels,Size1/2,Size2/2,1,&status);
       double norm=getNormalization2(mask_shadow, ea, detector_pixels, xdiff, ydiff);
       
       //new event array: method two
       if(norm>1.){

	 for(ii=0; ii<detector_pixels->xwidth; ii++){
	   for(jj=0; jj<detector_pixels->ywidth; jj++){
	     if(ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff]!=0.){
	       ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff]-=norm*mask_shadow->shadow[ii][jj];
	       if(ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff]<0.){
		 ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff]=0.;
	       }
	     }
	     
	   }
	 }

       }else{
	 threshold=2; 
       }
       FreeEventArray1d(EventImage1d);
       fftw_free(fftEventArray);
       fftw_free(fftInvMultiply);
        }while(threshold==1);

    //create FITS-file with all pix-coordinates
     savePositionList(position_list, par->PositionList, &status);

  // --- END of Reconstruction Process ---

  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Free the detector and sky image pixels.
  
    //set detector_pixels to original size again to be able to call destroy-fct from 'squarepixels.c'
  detector_pixels->xwidth=xdetsize_beforeRepix;
  destroySquarePixels(&detector_pixels);

  destroyCodedMask(&mask);
  FreeReconArray(&recon);
  FreeReconArray1d(ReconImage1d);
  FreeEventArray(ea);
  FreePixPositionList(position_list);
  FreeMaskShadow(mask_shadow,Size1);
  fftw_free(fftReconArray); 
  wcsfree(&wcs);
  wcsfree(&wcs2);
  free_SourceImage(sky_pixels);



  } while(0);  // END of the error handling loop.


  // Close the FITS files.
  status=closeCoMaEventFile(eventfile);
  free(eventfile);

  if (EXIT_SUCCESS==status) headas_chat(3, "comarecon finished successfully!\n\n");
  return(status);



}
