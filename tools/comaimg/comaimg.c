#include "comaimg.h"


////////////////////////////////////
/** Main procedure. */
int comaimg_main() {
  struct Parameters par;

  AttitudeCatalog* ac=NULL;
  struct Telescope telescope; // Telescope data.
  PhotonListFile* plf=NULL;
  ImpactListFile* ilf=NULL;
  double refxcrvl=0., refycrvl=0.;
  CodedMask* mask=NULL;

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("comaimg");
  set_toolversion("0.02");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read parameters using PIL library.
    if ((status=comaimg_getpar(&par))) break;

    float focal_length=par.MaskDistance;

    // Calculate the minimum cos-value for sources inside the FOV: 
    // (angle(x0,source) <= 1/2 * diameter)
    const double fov_min_align=cos(M_PI/3.);
    
    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Open the FITS file with the input photon list:
    plf=openPhotonListFile(par.PhotonList, READONLY, &status);
    CHECK_STATUS_BREAK(status);


    // Set up the Attitude.
    char ucase_buffer[MAXFILENAME];
    strcpy(ucase_buffer, par.Attitude);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer, "NONE")) {
      // Set up a simple pointing attitude.

      // First allocate memory.
      ac=getAttitudeCatalog(&status);
      CHECK_STATUS_BREAK(status);

      ac->entry=(AttitudeEntry*)malloc(sizeof(AttitudeEntry));
      if (NULL==ac->entry) {
	status=EXIT_FAILURE;
	SIXT_ERROR("memory allocation for AttitudeCatalog failed");
	break;
      }

      // Set the values of the entries.
      ac->nentries=1;
      ac->entry[0]=defaultAttitudeEntry();
      ac->entry[0].time=0.;
      ac->entry[0].nz=unit_vector(par.RA*M_PI/180., par.Dec*M_PI/180.);

      Vector vz={0., 0., 1.};
      ac->entry[0].nx=vector_product(vz, ac->entry[0].nz);

    } else {
      // Load the attitude from the given file.
      ac=loadAttitudeCatalog(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);

      // Check if the required time interval for the simulation
      // is a subset of the time described by the attitude file.
      if ((ac->entry[0].time > 0.) || 
	  (ac->entry[ac->nentries-1].time < par.Exposure)) {
	status=EXIT_FAILURE;
	SIXT_ERROR("attitude data does not cover the specified period");
	break;
      }
    }
    // END of setting up the attitude.


    // Load the coded mask from the file.
    mask=getCodedMaskFromFile(par.Mask, &status);
    CHECK_STATUS_BREAK(status);

    // Create a new FITS file for the output of the impact list.
    ilf=openNewImpactListFile(par.ImpactList, 0, &status);
    CHECK_STATUS_BREAK(status);

    // Write WCS header keywords.
    fits_update_key(ilf->fptr, TDOUBLE, "REFXCRVL", &refxcrvl, "", &status);
    fits_update_key(ilf->fptr, TDOUBLE, "REFYCRVL", &refycrvl, "", &status);
    CHECK_STATUS_BREAK(status);
    
    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    // Beginning of actual simulation (after loading required data):
    headas_chat(3, "start imaging process ...\n");

    // LOOP over all timesteps given the specified timespan from t0 to t0+timespan
    long attitude_counter=0;  // counter for AttitudeCatalog

    // SCAN PHOTON LIST 
    while (plf->row < plf->nrows) {
   
      Photon photon={.time=0.};
      
      // Read an entry from the photon list:
      status=PhotonListFile_getNextRow(plf, &photon);
      CHECK_STATUS_BREAK(status);

      // Check whether we are within the requested time interval.
      if (photon.time > par.Exposure) break;

      // Determine the unit vector pointing in the direction of the photon.
      Vector phodir=unit_vector(photon.ra, photon.dec);
   
      // Determine telescope pointing direction at the current time.
      telescope.nz=getTelescopeNz(ac, photon.time, &status);
      CHECK_STATUS_BREAK(status);

      // Check whether the photon is inside the FOV:
      // Compare the photon direction to the unit vector specifiing the 
      // direction of the telescope axis:
      if (check_fov(&phodir, &telescope.nz, fov_min_align)==0) {
	// Photon is inside the FOV!
	
	// Determine telescope data like direction etc. (attitude).
	// The telescope coordinate system consists of a nx, ny, and nz axis.
	// The nz axis is perpendicular to the detector plane and pointing along
	// the telescope direction. The nx axis is align along the detector 
	// x-direction, which is identical to the detector COLUMN.
	// The ny axis ix pointing along the y-direction of the detector,
	// which is also referred to as ROW.

	// Determine the current nx: perpendicular to telescope axis nz
	// and in the direction of the satellite motion.
	telescope.nx = 
	  normalize_vector(interpolate_vec(ac->entry[attitude_counter].nx, 
					   ac->entry[attitude_counter].time, 
					   ac->entry[attitude_counter+1].nx, 
					   ac->entry[attitude_counter+1].time, 
					   photon.time));
	
	// Remove the component along the vertical direction nz 
	// (nx must be perpendicular to nz!):
	double scp = scalar_product(&telescope.nz, &telescope.nx);
	telescope.nx.x -= scp*telescope.nz.x;
	telescope.nx.y -= scp*telescope.nz.y;
	telescope.nx.z -= scp*telescope.nz.z;
	telescope.nx = normalize_vector(telescope.nx);

	// The third axis of the coordinate system ny is perpendicular 
	// to telescope axis nz and nx:
	telescope.ny=normalize_vector(vector_product(telescope.nz, telescope.nx));
	
	// Determine the photon impact position on the detector (in [m]):
	struct Point2d position;  

	// Convolution with PSF:
	// Function returns 0, if the photon does not fall on the detector. 
	// If it hits the detector, the return value is 1.
	int retval=
	  getCodedMaskImpactPos(&position, &photon, mask, 
				&telescope, focal_length, &status);
	CHECK_STATUS_BREAK(status);
	if (1==retval) {
	  // Check whether the photon hits the detector within the FOV. 
	  // (Due to the effects of the mirrors it might have been scattered over 
	  // the edge of the FOV, although the source is inside the FOV.)
	  //if (sqrt(pow(position.x,2.)+pow(position.y,2.)) < 
	  //    tan(telescope.fov_diameter)*telescope.focal_length) {
	  
	  // New impact.
	  Impact impact;
	  impact.time  =photon.time;
	  impact.energy=photon.energy;
	  impact.position.x=position.x;
	  impact.position.y=position.y;
	  impact.ph_id     =photon.ph_id;
	  impact.src_id    =photon.src_id;

	  // Write the impact to the output file.
	  addImpact2File(ilf, &impact, &status);
	  CHECK_STATUS_BREAK(status);

	  //}
	} // END getCodedMaskImpactPos(...)
      } // End of FOV check.
    } // END of scanning LOOP over the photon list.
    CHECK_STATUS_BREAK(status);
      
  } while(0);  // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // release HEADAS random number generator
  HDmtFree();

  // Close the FITS files.
  freeImpactListFile(&ilf, &status);
  freePhotonListFile(&plf, &status);

  freeAttitudeCatalog(&ac);
  destroyCodedMask(&mask);

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
  
  // Get the filename of the Coded Mask file (FITS image file).
  status=ape_trad_query_string("Mask", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the filename of the coded mask");
    return(status);
  }
  strcpy(par->Mask, sbuffer);
  free(sbuffer);

  // Get the filename of the impact list file (FITS output file).
  status=ape_trad_query_string("ImpactList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the filename of the impact list");
    return(status);
  }
  strcpy(par->ImpactList, sbuffer);
  free(sbuffer);

  // Read the distance between the coded mask and the detector plane [m].
  status=ape_trad_query_double("MaskDistance", &par->MaskDistance);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the distance between the mask and the detector");
    return(status);

  }

  status=ape_trad_query_string("Attitude", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the attitude file");
    return(status);
  } 
  strcpy(par->Attitude, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("RA", &par->RA);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the right ascension of the telescope pointing");
    return(status);
  } 

  status=ape_trad_query_float("Dec", &par->Dec);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the declination of the telescope pointing");
    return(status);
  } 

  status=ape_trad_query_double("Exposure", &par->Exposure);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the exposure time");
    return(status);
  } 

  return(status);
}



