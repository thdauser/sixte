#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "comaimg.h"


////////////////////////////////////
/** Main procedure. */
int comaimg_main() {
  struct Parameters parameters;

  AttitudeCatalog* attitudecatalog=NULL;
  struct Telescope telescope; // Telescope data.
  PhotonListFile photonlistfile;
  ImpactListFile impactlistfile;
  double refxcrvl=0., refycrvl=0.;
  CodedMask* mask=NULL;

  int status=EXIT_SUCCESS; // Error status.
  char msg[MAXMSG];


  // Register HEATOOL:
  set_toolname("comaimg");
  set_toolversion("0.01");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read parameters using PIL library.
    if ((status=comaimg_getpar(&parameters))) break;


    // Calculate the minimum cos-value for sources inside the FOV: 
    // (angle(x0,source) <= 1/2 * diameter)
    const double fov_min_align = cos(M_PI/3.);
    
    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Open the FITS file with the input photon list:
    status = openPhotonListFile(&photonlistfile, parameters.photonlist_filename, 
				READONLY);
    if (EXIT_SUCCESS!=status) break;

    // Open the attitude file specified in the header keywords of the photon list.
    char comment[MAXMSG];
    if (fits_read_key(photonlistfile.fptr, TSTRING, "ATTITUDE", 
		      &parameters.attitude_filename, comment, &status)) break;
    if (NULL==(attitudecatalog=getEntireAttitudeCatalog(parameters.attitude_filename,
							&status))) break;

    // Load the coded mask from the file.
    mask = getCodedMaskFromFile(parameters.mask_filename, &status);
    if(EXIT_SUCCESS!=status) break;

    // Create a new FITS file for the output of the impact list.
    status = openNewImpactListFile(&impactlistfile, parameters.impactlist_filename, 
				   parameters.impactlist_template);
    if (EXIT_SUCCESS!=status) break;
    // Write WCS header keywords.
    if (fits_update_key(impactlistfile.fptr, TDOUBLE, "REFXCRVL", &refxcrvl, "", &status)) break;
    if (fits_update_key(impactlistfile.fptr, TDOUBLE, "REFYCRVL", &refycrvl, "", &status)) break;
    // Add attitude filename.
    if (fits_update_key(impactlistfile.fptr, TSTRING, "ATTITUDE", parameters.attitude_filename,
		       "name of the attitude FITS file", &status)) break;
    
    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    // Beginning of actual simulation (after loading required data):
    headas_chat(5, "start imaging process ...\n");

    // LOOP over all timesteps given the specified timespan from t0 to t0+timespan
    long attitude_counter=0;  // counter for AttitudeCatalog

    // SCAN PHOTON LIST    
    for(photonlistfile.row=0; (photonlistfile.row<photonlistfile.nrows)&&(EXIT_SUCCESS==status); 
	photonlistfile.row++) {
      
      // Read an entry from the photon list:
      int anynul = 0;
      Photon photon = { .time=0., .energy=0., .ra=0., .dec=0. };
      fits_read_col(photonlistfile.fptr, TDOUBLE, photonlistfile.ctime, 
		    photonlistfile.row+1, 1, 1, &photon.time, &photon.time, &anynul, &status);
      fits_read_col(photonlistfile.fptr, TFLOAT, photonlistfile.cenergy, 
		    photonlistfile.row+1, 1, 1, &photon.energy, &photon.energy, &anynul, &status);
      fits_read_col(photonlistfile.fptr, TDOUBLE, photonlistfile.cra, 
		    photonlistfile.row+1, 1, 1, &photon.ra, &photon.ra, &anynul, &status);
      fits_read_col(photonlistfile.fptr, TDOUBLE, photonlistfile.cdec, 
		    photonlistfile.row+1, 1, 1, &photon.dec, &photon.dec, &anynul, &status);
      if (status!=EXIT_SUCCESS) break;

      // Rescale from [deg] -> [rad]
      photon.ra  = photon.ra *M_PI/180.;
      photon.dec = photon.dec*M_PI/180.;
      // Determine a unit vector pointing in the direction of the photon.
      photon.direction = unit_vector(photon.ra, photon.dec);


      // Get the last attitude entry before 'photon.time'
      // (in order to interpolate the attitude at this time between 
      // the neighboring calculated values):
      for( ; attitude_counter<attitudecatalog->nentries-1; attitude_counter++) {
	if(attitudecatalog->entry[attitude_counter+1].time>photon.time) {
	  break;
	}
      }
      if(fabs(attitudecatalog->entry[attitude_counter].time-photon.time)>600.) { 
	// no entry within 10 minutes !!
	status = EXIT_FAILURE;
	sprintf(msg, "Error: no adequate orbit entry for time %lf!\n", photon.time);
	HD_ERROR_THROW(msg,status);
	break;
      }

   
      // Check whether the photon is inside the FOV:
      // First determine telescope pointing direction at the current time.
      // TODO: replace this calculation by proper attitude interpolation.
      telescope.nz = 
	normalize_vector(interpolate_vec(attitudecatalog->entry[attitude_counter].nz, 
					 attitudecatalog->entry[attitude_counter].time, 
					 attitudecatalog->entry[attitude_counter+1].nz, 
					 attitudecatalog->entry[attitude_counter+1].time, 
					 photon.time));

      // Compare the photon direction to the unit vector specifiing the 
      // direction of the telescope axis:
      if (check_fov(&photon.direction, &telescope.nz, fov_min_align)==0) {
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
	  normalize_vector(interpolate_vec(attitudecatalog->entry[attitude_counter].nx, 
					   attitudecatalog->entry[attitude_counter].time, 
					   attitudecatalog->entry[attitude_counter+1].nx, 
					   attitudecatalog->entry[attitude_counter+1].time, 
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
	if (1==getCodedMaskImpactPos(&position, &photon, mask, &telescope)) {
	  // Check whether the photon hits the detector within the FOV. 
	  // (Due to the effects of the mirrors it might have been scattered over 
	  // the edge of the FOV, although the source is inside the FOV.)
	  //if (sqrt(pow(position.x,2.)+pow(position.y,2.)) < 
	  //    tan(telescope.fov_diameter)*telescope.focal_length) {
	    
	  // Insert the impact position with the photon data into the impact list:
	  fits_insert_rows(impactlistfile.fptr, impactlistfile.row++, 1, &status);
	  fits_write_col(impactlistfile.fptr, TDOUBLE, impactlistfile.ctime, 
			 impactlistfile.row, 1, 1, &photon.time, &status);
	  fits_write_col(impactlistfile.fptr, TFLOAT, impactlistfile.cenergy, 
			 impactlistfile.row, 1, 1, &photon.energy, &status);
	  fits_write_col(impactlistfile.fptr, TDOUBLE, impactlistfile.cx, 
			 impactlistfile.row, 1, 1, &position.x, &status);
	  fits_write_col(impactlistfile.fptr, TDOUBLE, impactlistfile.cy, 
			 impactlistfile.row, 1, 1, &position.y, &status);
	  impactlistfile.nrows++;
	  //}
	} // END get_psf_pos(...)
      } // End of FOV check.
    } // END of scanning LOOP over the photon list.

  } while(0);  // END of the error handling loop.


  // --- cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // release HEADAS random number generator
  HDmtFree();

  // Close the FITS files.
  status += closeImpactListFile(&impactlistfile);
  status += closePhotonListFile(&photonlistfile);

  free_AttitudeCatalog(attitudecatalog);
  freeCodedMask(mask);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully!\n\n");
  return(status);
}



int comaimg_getpar(struct Parameters* parameters)
{
  int status=EXIT_SUCCESS; // Error status.

  // Get the filename of the input photon list (FITS file).
  if ((status = PILGetFname("photonlist_filename", parameters->photonlist_filename))) {
    HD_ERROR_THROW("Error reading the filename of the photon list!\n", status);
  }
  
  // Get the filename of the Coded Mask file (FITS image file).
  else if ((status = PILGetFname("mask_filename", parameters->mask_filename))) {
    HD_ERROR_THROW("Error reading the filename of the Coded Mask!\n", status);
  }

  // Get the filename of the impact list file (FITS output file).
  else if ((status = PILGetFname("impactlist_filename", parameters->impactlist_filename))) {
    HD_ERROR_THROW("Error reading the filename of the impact list output file!\n", status);
  }

  // Read the distance between the coded mask and the detector plane [m].
  else if ((status = PILGetReal("mask_distance", &parameters->mask_distance))) {
    HD_ERROR_THROW("Error reading the distance between the mask and the detector plane!\n", 
		   status);
  }
  if (EXIT_SUCCESS!=status) return(status);

  // Get the name of the FITS template directory.
  // First try to read it from the environment variable.
  // If the variable does not exist, read it from the PIL.
  char* buffer;
  if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
    strcpy(parameters->impactlist_template, buffer);
  } else {
    if ((status = PILGetFname("fits_templates", parameters->impactlist_template))) {
      HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);
      
    }
  }
  // Set the impact list template file:
  strcat(parameters->impactlist_template, "/impactlist.tpl");

  return(status);
}



