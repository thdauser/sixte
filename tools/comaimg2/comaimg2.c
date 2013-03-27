#include "comaimg2.h"

/////////////////////////////////////////////////////////////
//simulates detection process. All incoming photons are    //
//passing the mask and hitting the detector!               //
//Input:photon-list(t,E,RA,DEC,PH_ID,SRC_ID),mask-file,    //
//      pointing-attitude(DEC, RA)                         //
//      or attitude-file(DEC, RA for different times),     //
//      mask-width,det-width(in m) for FOV,distance        //
//Output:impact-list(x-and y-value on detection-plane)     //
/////////////////////////////////////////////////////////////

/** Main procedure. */
int comaimg2_main() {
  struct Parameters par;

  struct Telescope telescope; //Telescope coordinate system
  PhotonListFile* plf=NULL;
  ImpactListFile* ilf=NULL;
  CodedMask* mask=NULL;
  AttCatalog* ac=NULL;

  double refxcrvl=0., refycrvl=0.;

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("comaimg2");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read parameters using PIL library.
    status=comaimg2_getpar(&par);
    if (EXIT_SUCCESS!=status) break;

    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Open the FITS file with the input photon list.
    plf=openPhotonListFile(par.PhotonList, READONLY, &status);
    if (EXIT_SUCCESS!=status) break;

    // Load the coded mask from the file.
    mask=getCodedMaskFromFile(par.Mask, &status);
    if (EXIT_SUCCESS!=status) break;

    //Calculate min cos-value for sources inside FOV.
    const float fov_min_align = cos(det_phi_max(par.MaskDistance,
						par.x_mask, par.y_mask,
						par.x_det, par.y_det));
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

      //Allocates memory for struct AttitudeCatalog.
      //Initializes nentries, current_entry to 0;entry to NULL.
      ac=getAttCatalog(&status);
      CHECK_STATUS_BREAK(status);

      //Allocates memory for struct AttitudeEntry.
      ac->entry=(AttEntry*)malloc(sizeof(AttEntry));
      if (NULL==ac->entry) {
	status=EXIT_FAILURE;
	SIXT_ERROR("memory allocation for AttitudeEntry failed");
	break;
      }

      //Set the values of the AttitudeCatalog-Entry.

      ac->nentries=1;
      ac->entry[0]=initializeAttitudeEntry();  
      //ac->entry[0].time=0;

      //Telescope pointing direction:
      ac->entry[0].nz=unit_vector(par.RA*M_PI/180.0,par.DEC* M_PI/180.0);
      //Unit-vector in z-direction:
      Vector vz = {0.,0.,1.};
      //Vector perpendicular to nz-vz-plane:
      ac->entry[0].nx=vector_product(vz, ac->entry[0].nz);
      //initialize telescope coordinate sytem
      telescope.nz=ac->entry[0].nz;
      telescope.nx=ac->entry[0].nx;

    } else {
      //Load the attitude from file.
      ac=loadAttCatalog(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);

      Vector initializey={0.,0.,0.};
      telescope.ny =initializey;

      //Check if the required time interval for the simulation
      //lies within the time interval in the attitude-file
      if((ac->entry[0].time > 0) || ((ac->entry[ac->nentries-1].time) < (par.Exposure +par.Timezero))){
	status=EXIT_FAILURE;
	SIXT_ERROR("attitude data does not cover the specified period");
	break;
      }
      
    }//END of setting up the attitude.
    

    //Distance mask-detector:
    float distance = par.MaskDistance;
    //Detector- dimensions:
    float x_det = par.x_det;
    float y_det = par.y_det;

    //Create a new FITS-file for the impact-list (output).
    ilf = openNewImpactListFile(par.ImpactList, 1, &status);
    CHECK_STATUS_BREAK(status);

    //Write WCS header keywords.
    fits_update_key(ilf->fptr, TDOUBLE, "REFXCRVL", &refxcrvl, "", &status);
    fits_update_key(ilf->fptr, TDOUBLE, "REFYCRVL", &refycrvl, "", &status);
    CHECK_STATUS_BREAK(status);
  
    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    //Beginning of actual simulation (after loading required data).
    headas_chat(3, "start imaging process ...\n");

    //SCAN PHOTON LIST (starts with first entry, plf initialized to 0 at beginning)
    while (plf->row < plf->nrows) {
         
      //Read an entry from the photon list:
      Photon photon={.time= 0.0};
      status=PhotonListFile_getNextRow(plf, &photon);
      if (EXIT_SUCCESS!=status) break;

      //Check whether photon-list-entry is within requested time interval.
      if (photon.time > (par.Exposure + par.Timezero))
	break;

      //Determine unit vector in photon-direction
      Vector phodir=unit_vector(photon.ra, photon.dec);

      //Determine current telescope pointing direction.
      telescope.nz=GetTelescopeNz(ac, photon.time, &status);
      CHECK_STATUS_BREAK(status);
      
    
      //Check whether photon is inside FOV:
      //Compare photon direction to direction of telescope axis
      if (check_fov(&phodir, &telescope.nz, fov_min_align)==0){
	//Photon is inside fov

	telescope.nx=ac->entry[0].nx;
	getTelAxes(ac,&telescope.nx,&telescope.ny,&telescope.nz,photon.time,&status);

	//Determine photon impact position on detector in [m].
	struct Point2d position;
	
	//function determines impact position on mask first;
	//checks wheather pixel is transparent;
	//if photon then hits the detector, return value is 1, 0 else
       	int reval = getImpactPos(&position, &phodir,
				 mask, &telescope,
				 distance, x_det, y_det, &status);
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
	  addImpact2File(ilf, &impact, &status);
	  CHECK_STATUS_BREAK(status);
	} // END of photon hits the detector.
      } // END of photon  inside fov.
    } // END of scanning LOOP over the photon list.
    if (EXIT_SUCCESS!=status) break;
      
  } while(0);  // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the FITS files.
  freeImpactListFile(&ilf, &status);
  freePhotonListFile(&plf, &status);
  freeAttCatalog(&ac);

  // Release memory.
  destroyCodedMask(&mask);

  // release HEADAS random number generator
  HDmtFree();

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int comaimg2_getpar(struct Parameters* par)
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
