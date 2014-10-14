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


   Copyright 2007-2014 Christian Schmid, Mirjam Oertel, FAU
*/

#include "comaimg.h"

/////////////////////////////////////////////////////////////
//IMAGING: incoming photons are passing transparent mask-  //
//pixels and are hitting the detector (if not the walls)!  //
//Input:photon-list(t,E,RA,DEC,PH_ID,SRC_ID),mask-file,    //
//      pointing-attitude(DEC, RA)                         //
//      or attitude-file(DEC, RA for different times),     //
//      mask-width,det-width(in m) for FOV,distance        //
//Output:impact-list(x,y [m],t,E)                          //
/////////////////////////////////////////////////////////////

/** Main procedure. */
int comaimg_main() {
  struct Parameters par;

  struct Telescope telescope; //Telescope coordinate system
  PhotonFile* plf=NULL;
  ImpactFile* ilf=NULL;
  CodedMask* mask=NULL;
  Attitude* ac=NULL;

  double refxcrvl=0., refycrvl=0.;

  int status=EXIT_SUCCESS; // Error status.

  // Register HEATOOL:
  set_toolname("comaimg");
  set_toolversion("0.01");


  do { //Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    //Read parameters using PIL library.
    status=comaimg_getpar(&par);
    if (EXIT_SUCCESS!=status) break;
    
    //Detector-setup:
    float x_det = par.x_det;
    float y_det = par.y_det;
    //Distance mask-detector:
    float distance = par.MaskDistance;
    float det_pixelwidth= par.det_pixelwidth;
    double ra=par.RA;
    double dec=par.DEC;

    // Initialize the random number generator.
    sixt_init_rng(getSeed(-1),&status);
    CHECK_STATUS_BREAK(status);

    //Open the FITS file with the input photon list.
    plf=openPhotonFile(par.PhotonList,READONLY,&status);
    if (EXIT_SUCCESS!=status) break;

    //Load the coded mask from the file.
    mask=getCodedMaskFromFile(par.Mask,&status);
    if (EXIT_SUCCESS!=status) break;

    //Calculate min cos-value for sources inside FOV.
    const float fov_min_align = cos(det_phi_max(par.MaskDistance,
						par.x_mask, par.y_mask,
						x_det, y_det));
    //Set up the telescope attitude:

    //Buffer for attitude-filename.
    char att_buffer[MAXFILENAME];
    //Copy attitude-filename from par-file to buffer.
    strcpy(att_buffer, par.Attitude);
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
      ac->entry[0].nz=normalize_vector(unit_vector(par.RA*M_PI/180.0,par.DEC* M_PI/180.0));
      //Unit-vector in z-direction:
      Vector vz = {0.,0.,1.};
      //Vector perpendicular to nz-vz-plane:
      ac->entry[0].nx=vector_product(vz,ac->entry[0].nz);
      //initialize telescope coordinate sytem
      telescope.nz=ac->entry[0].nz;
      telescope.nx=ac->entry[0].nx;

    } else {
      //Load the attitude from file.
      ac=loadAttitude(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);

      Vector initializey={0.,0.,0.};
      telescope.ny =initializey;
      telescope.nx=ac->entry[0].nx;

      //Check if the required time interval for the simulation
      //lies within the time interval in the attitude-file
      if((ac->entry[0].time > 0) || ((ac->entry[ac->nentries-1].time) < (par.Exposure + par.Timezero))){
	status=EXIT_FAILURE;
	SIXT_ERROR("attitude data does not cover the specified period");
	break;
      }
      
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
    ilf=openNewImpactFile(par.ImpactList, 
			  telescop, instrume, filter,
			  ancrfile, respfile,
			  mjdref, timezero, tstart, tstop,
			  0, &status);
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
      if (photon.time > (par.Exposure + par.Timezero))
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
	  setWCScurrentPointing(par.Attitude,ac,&telescope.nz,&wcs,&status);
	}

	//Determine photon impact position on detector in [m].
	struct Point2d position;
	
	//first:impact position in mask-plane (transparent pixels only);
	//if photon then hits the detector (and not the walls), return value is 1, 0 else
	int reval=getImpactPos2(&wcs,&position,mask,photon.ra*180./M_PI,
				photon.dec*180./M_PI,det_pixelwidth,x_det,y_det,&status);

	CHECK_STATUS_BREAK(status);

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

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int comaimg_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  // Get the filename of the input photon list (FITS file).
  status=ape_trad_query_string("PhotonList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the filename of the photon list");
    return(status);
  }
  strcpy(par->PhotonList, sbuffer);
  free(sbuffer);
  
  // Get the filename of the coded mask file (FITS image file).
  status=ape_trad_query_string("Mask", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the filename of the coded mask");
    return(status);
  }
  strcpy(par->Mask, sbuffer);
  free(sbuffer);

  // Get the filename of the output impact list file (FITS file).
  status=ape_trad_query_string("ImpactList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the filename of the impact list");
    return(status);
  }
  strcpy(par->ImpactList, sbuffer);
  free(sbuffer);

  //Read distance between the detector and the mask plane [m].
  status=ape_trad_query_float("MaskDistance", &par->MaskDistance);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the distance between the mask and detection plane");
    return(status);
  }

  // Get the filename of the telescope-attitude input file (FITS file).
  status=ape_trad_query_string("Attitude", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the filename of the telescope attitude");
    return(status);
  }
  strcpy(par->Attitude, sbuffer);
  free(sbuffer);

  //Read width of the mask [m].
  status=ape_trad_query_float("x_mask", &par->x_mask);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading width of the mask");
    return(status);
  }

  //Read depth of the mask [m].
  status=ape_trad_query_float("y_mask", &par->y_mask);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the depth of the mask");
    return(status);
  }
 
  //Read width of the detector [m].
  status=ape_trad_query_float("x_det", &par->x_det);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the width of the detector");
    return(status);
  }

  //Read depth of the detector [m].
  status=ape_trad_query_float("y_det", &par->y_det);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the depth of the detector");
    return(status);
  }

  //Read width of one detector pixel[m].
  status=ape_trad_query_float("det_pixelwidth", &par->det_pixelwidth);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the width of one detector pixel");
    return(status);
  }
 
  status=ape_trad_query_double("RA", &par->RA);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the right ascension of the telescope");
    return(status);
  } 

 status=ape_trad_query_double("DEC", &par->DEC);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the declination of the telescope");
    return(status);
  } 

  //Read time-offset for simulated intervall [s].
  status=ape_trad_query_double("Timezero", &par->Timezero);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the time-offset for simulated intervall");
    return(status);
  }

  status=ape_trad_query_double("Exposure", &par->Exposure);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the exposure time");
    return(status);
  } 

  return(status);
}
