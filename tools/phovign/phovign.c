#include "phovign.h"


////////////////////////////////////
/** Main procedure. */
int phovign_main() {
  struct Parameters par;

  AttitudeCatalog* ac=NULL;
  PhotonListFile* plif=NULL;
  PhotonListFile* plof=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL:
  set_toolname("phovign");
  set_toolversion("0.03");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read parameters using PIL library.
    if ((status=phovign_getpar(&par))) break;

    headas_chat(3, "initialize ...\n");

    // Determine the input photon list file name.
    char inputlist_filename[MAXFILENAME];
    strcpy(inputlist_filename, par.InputList);

    // Determine the output photon list file name.
    char outputlist_filename[MAXFILENAME];
    strcpy(outputlist_filename, par.OutputList);

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
      ac->entry[0].time = 0.;
      ac->entry[0].nz = unit_vector(par.RA*M_PI/180., par.Dec*M_PI/180.);

      Vector vz = {0., 0., 1.};
      ac->entry[0].nx = vector_product(vz, ac->entry[0].nz);

    } else {
      // Load the attitude from the given file.
      ac=loadAttitudeCatalog(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);
    }
    // END of setting up the attitude.
    
    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    headas_chat(3, "apply vignetting ...\n");

    // Open the photon list files.
    plif=openPhotonListFile(inputlist_filename, READONLY, &status);
    CHECK_STATUS_BREAK(status);
    plof=openNewPhotonListFile(outputlist_filename, par.clobber,
			       &status);
    CHECK_STATUS_BREAK(status);

    // Scan the entire photon list.
    int progress=0;  
    while (plif->row < plif->nrows) {

      Photon photon={.time=0.};
      
      // Read an entry from the photon list:
      status=PhotonListFile_getNextRow(plif, &photon);
      CHECK_STATUS_BREAK(status);

      // Apply the vignetting.
      // Compare the photon direction to the direction of the telescope axis.
      Vector nz=getTelescopeNz(ac, photon.time, &status);
      CHECK_STATUS_BREAK(status);
      Vector phodir=unit_vector(photon.ra, photon.dec);

      double cosine=nz.x*phodir.x + nz.y*phodir.y + nz.z*phodir.z;
      double rnd=sixt_get_random_number();

      // Delete the photon from the FITS file.
      if (cosine>=rnd) {
	// Write the photon to the output file.
	status=addPhoton2File(plof, &photon);
	CHECK_STATUS_BREAK(status);
      }

      // Program progress output.
      while ((int)(plif->row*1000./plif->nrows)>progress) {
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
  freePhotonListFile(&plof, &status);
  freePhotonListFile(&plif, &status);
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

  status=ape_trad_query_file_name("InputList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("could not read the name of the input photon list");
    return(status);
  } 
  strcpy(par->InputList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("OutputList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("could not read the name of the output photon list");
    return(status);
  } 
  strcpy(par->OutputList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Attitude", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("could not read the name of the attitude");
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

