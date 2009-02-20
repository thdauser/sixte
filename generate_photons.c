#include "generate_photons.h"



////////////////////////////
int generate_photons_getpar(
			    char orbit_filename[],
			    char attitude_filename[],
			    char photonlist_filename[]
			    )
{
  char msg[MAXMSG];           // error message buffer
  int status = EXIT_SUCCESS;  // error status flag

  // Get the filename of the Orbit file (FITS file):
  if ((status = PILGetFname("orbit_filename", orbit_filename))) {
    sprintf(msg, "Error reading the filename of the orbit file!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the filename of the Attitude file (FITS file):
  else if ((status = PILGetFname("attitude_filename", attitude_filename))) {
    sprintf(msg, "Error reading the filename of the attitude file!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the filename of the Photon-List file (FITS file):
  else if ((status = PILGetFname("photonlist_filename", photonlist_filename))) {
    sprintf(msg, "Error reading the filename of the output file for "
	    "the photon list!\n");
    HD_ERROR_THROW(msg,status);
  }

  return(status);
}




//////////////////////////
int generate_photons_main() 
{
  // Names of several input and output files:
  char orbit_filename[FILENAME_LENGTH];      // input: orbit
  char attitude_filename[FILENAME_LENGTH];   // input: attitude
  char photonlist_filename[FILENAME_LENGTH]; // output: photon list

  // Number of entries in the orbit list ( <= orbit_nrows)
  long sat_nentries;                      

  char msg[MAXMSG];           // error message buffer
  int status = EXIT_SUCCESS;  // error status flag


  do {  // Beginning of ERROR HANDLING Loop.

    if ((status = generate_photons_getpar(orbit_filename, attitude_filename,
					  photonlist_filename))) 
      break;


  } while(0); // END of ERROR HANDLING Loop.



  // --- clean up ---

  return(status);
}





