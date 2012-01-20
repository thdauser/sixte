#include "phovign.h"


////////////////////////////////////
/** Main procedure. */
int phovign_main() {
  struct Parameters par;

  AttitudeCatalog* ac=NULL;
  PhotonListFile* plf=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL:
  set_toolname("phovign");
  set_toolversion("0.01");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read parameters using PIL library.
    if ((status=phovign_getpar(&par))) break;

    headas_chat(3, "initialize ...\n");

    // Determine the photon list file name.
    char photonlist_filename[MAXFILENAME];
    strcpy(photonlist_filename, par.PhotonList);

    // Determine the random number seed.
    int seed;
    if (-1!=par.Seed) {
      seed = par.Seed;
    } else {
      // Determine the seed from the system clock.
      seed = (int)time(NULL);
    }

    // Initialize HEADAS random number generator.
    HDmtInit(seed);

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
	status = EXIT_FAILURE;
	SIXT_ERROR("memory allocation for AttitudeCatalog failed");
	break;
      }

      // Set the values of the entries.
      ac->nentries=1;
      ac->entry[0] = defaultAttitudeEntry();
      ac->entry[0].time = par.TIMEZERO;
      ac->entry[0].nz = unit_vector(par.RA*M_PI/180., par.Dec*M_PI/180.);

      Vector vz = {0., 0., 1.};
      ac->entry[0].nx = vector_product(vz, ac->entry[0].nz);

    } else {
      // Load the attitude from the given file.
      ac=loadAttitudeCatalog(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);

      // Check if the required time interval for the simulation
      // is a subset of the time described by the attitude file.
      if ((ac->entry[0].time > par.TIMEZERO) || 
	  (ac->entry[ac->nentries-1].time < par.TIMEZERO+par.Exposure)) {
	status=EXIT_FAILURE;
	char msg[MAXMSG];
	sprintf(msg, "attitude data does not cover the "
		"specified period from %lf to %lf!", 
		par.TIMEZERO, par.TIMEZERO+par.Exposure);
	SIXT_ERROR(msg);
	break;
      }
    }
    // END of setting up the attitude.
    
    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    headas_chat(3, "apply vignetting ...\n");

    // Open the photon list file.
    plf=openPhotonListFile(photonlist_filename, READWRITE, &status);
    CHECK_STATUS_BREAK(status);

    // Scan the entire photon list.
    int progress=0;  
    while (plf->row < plf->nrows) {

      Photon photon={.time=0.};
      
      // Read an entry from the photon list:
      status=PhotonListFile_getNextRow(plf, &photon);
      CHECK_STATUS_BREAK(status);

      // Check whether we are still within the requested time interval.
      if (photon.time < par.TIMEZERO) continue;
      if (photon.time > par.TIMEZERO+par.Exposure) break;

      // Apply the vignetting.
      // Compare the photon direction to the direction of the telescope axis.
      Vector nz=getTelescopeNz(ac, photon.time, &status);
      CHECK_STATUS_BREAK(status);
      Vector phodir=unit_vector(photon.ra, photon.dec);

      double cosine=nz.x*phodir.x + nz.y*phodir.y + nz.z*phodir.z;
      double rnd=sixt_get_random_number();

      // Delete the photon from the FITS file.
      if (cosine<rnd) {
	fits_delete_rows(plf->fptr, plf->row, 1, &status);
	CHECK_STATUS_BREAK(status);

	plf->nrows--;
	plf->row--;
      }

      // Program progress output.
      while ((int)(plf->row*1000./plf->nrows)>progress) {
	progress++;
	headas_chat(2, "\r%.1lf %%", progress*1./10.);
	fflush(NULL);
      }

    }; // END of photon processing loop.
    CHECK_STATUS_BREAK(status);

    // Progress output.
    headas_chat(2, "\r%.1lf %%\n", 100.);
    fflush(NULL);

    // --- END of imaging process ---

  } while(0); // END of the error handling loop.


  // --- cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Close the FITS files.
  freePhotonListFile(&plf, &status);
  freeAttitudeCatalog(&ac);

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");

  return(status);
}


int phovign_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("PhotonList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the photon list!\n", status);
    return(status);
  } 
  strcpy(par->PhotonList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Attitude", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the attitude!\n", status);
    return(status);
  } 
  strcpy(par->Attitude, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("RA", &par->RA);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the right ascension of the "
		   "telescope pointing!\n", status);
    return(status);
  } 

  status=ape_trad_query_float("Dec", &par->Dec);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the declination of the "
		   "telescope pointing!\n", status);
    return(status);
  } 

  status=ape_trad_query_double("MJDREF", &par->MJDREF);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading MJDREF!\n", status);
    return(status);
  } 

  status=ape_trad_query_double("TIMEZERO", &par->TIMEZERO);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading TIMEZERO!\n", status);
    return(status);
  } 

  status=ape_trad_query_double("Exposure", &par->Exposure);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the exposure time!\n", status);
    return(status);
  } 

  status=ape_trad_query_int("seed", &par->Seed);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the seed for the random "
		   "number generator!\n", status);
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }

  return(status);
}

