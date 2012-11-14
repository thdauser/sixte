#include "wfmphovign.h"


////////////////////////////////////
/** Main procedure. */
int wfmphovign_main() {
  struct Parameters par;

  AttitudeCatalog* ac=NULL;
  PhotonListFile* plif=NULL;

  fitsfile* ofptr=NULL;
  long onrows=0;
  int cenergy, cra, cdec;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL:
  set_toolname("wfmphovign");
  set_toolversion("0.05");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    // Read parameters using PIL library.
    if ((status=wfmphovign_getpar(&par))) break;

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

    // Open the input photon list file.
    plif=openPhotonListFile(inputlist_filename, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Open the output photon list file.
    // Check if the file already exists.
    int exists;
    fits_file_exists(outputlist_filename, &exists, &status);
    CHECK_STATUS_BREAK(status);
    if (0!=exists) {
      if (0!=par.clobber) {
	// Delete the file.
	remove(outputlist_filename);
      } else {
	// Throw an error.
	char msg[MAXMSG];
	sprintf(msg, "file '%s' already exists", outputlist_filename);
	SIXT_ERROR(msg);
	status=EXIT_FAILURE;
	break;
      }
    }

    // Create a new empty photon list FITS file.
    fits_create_file(&ofptr, outputlist_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Create the binary table for the photon list.
    char *ttype[] = { "ENERGY", "RA", "DEC" };
    char *tform[] = { "E", "E", "E" };
    char *tunit[] = { "keV", "deg", "deg" };
    fits_create_tbl(ofptr, BINARY_TBL, 0, 3, ttype, tform, tunit, 
		    "PHOTONS", &status);
    CHECK_STATUS_BREAK(status);

    // Determine the column numbers in the output file.
    fits_get_colnum(ofptr, CASEINSEN, "ENERGY", &cenergy, &status); 
    fits_get_colnum(ofptr, CASEINSEN, "RA", &cra, &status);
    fits_get_colnum(ofptr, CASEINSEN, "DEC", &cdec, &status);
    CHECK_STATUS_BREAK(status);

    // Add header information about program parameters.
    // The second parameter "1" means that the headers are written
    // to the first extension.
    HDpar_stamp(ofptr, 1, &status);
    CHECK_STATUS_BREAK(status);

    // Move back to the second HDU.
    int hdu_type;
    fits_movabs_hdu(ofptr, 2, &hdu_type, &status);
    CHECK_STATUS_BREAK(status);
    
    // Copy FITS header keywords.
    char attitude[MAXMSG], comment[MAXMSG];
    fits_read_key(plif->fptr, TSTRING, "ATTITUDE", 
		  attitude, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_update_key(ofptr, TSTRING, "ATTITUDE",
		    attitude, comment, &status);
    CHECK_STATUS_BREAK(status);

    double dbuffer=0.;
    fits_read_key(plif->fptr, TDOUBLE, "MJDREF", 
		  &dbuffer, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_update_key(ofptr, TDOUBLE, "MJDREF",
		    &dbuffer, comment, &status);
    CHECK_STATUS_BREAK(status);

    fits_read_key(plif->fptr, TDOUBLE, "TIMEZERO", 
		  &dbuffer, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_update_key(ofptr, TDOUBLE, "TIMEZERO", 
		    &dbuffer, comment, &status);
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
      double rnd=sixt_get_random_number(&status);
      CHECK_STATUS_BREAK(status);

      // Delete the photon from the FITS file.
      if (cosine>=rnd) {
	// Append the photon to the output file.

	// Insert a new, empty row to the table:
	fits_insert_rows(ofptr, onrows, 1, &status);
	CHECK_STATUS_BREAK(status);
	onrows++;
	
	// Store the data in the FITS file.
	fits_write_col(ofptr, TFLOAT, cenergy, 
		       onrows, 1, 1, &photon.energy, &status);
	float fbuffer=photon.ra * 180./M_PI;
	fits_write_col(ofptr, TFLOAT, cra, 
		       onrows, 1, 1, &fbuffer, &status);
	fbuffer=photon.dec * 180./M_PI;
	fits_write_col(ofptr, TFLOAT, cdec, 
		       onrows, 1, 1, &fbuffer, &status);
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
  if (NULL!=ofptr) {
    fits_close_file(ofptr, &status);
    ofptr=NULL;
  }
  freePhotonListFile(&plif, &status);
  freeAttitudeCatalog(&ac);

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");

  return(status);
}


int wfmphovign_getpar(struct Parameters* par)
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

