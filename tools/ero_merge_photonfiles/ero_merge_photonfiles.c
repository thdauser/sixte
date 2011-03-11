#include "ero_merge_photonfiles.h"


int ero_merge_photonfiles_getpar(struct Parameters* parameters)
{
  int status = EXIT_SUCCESS;

  if ((status = PILGetInt("n_inputfiles", &parameters->n_inputfiles))) {
    HD_ERROR_THROW("Error reading the number of input files!\n", status);
  }

  else if ((status = PILGetFname("input_prefix", parameters->input_prefix))) {
    HD_ERROR_THROW("Error reading the prefix for the input files!\n", status);
  }

  else if ((status = PILGetFname("output_filename", parameters->output_filename))) {
    HD_ERROR_THROW("Error reading the name of the output file!\n", status);
  }

  // Get the name of the FITS template directory.
  // First try to read it from the environment variable.
  // If the variable does not exist, read it from the PIL.
  else { 
    char* buffer;
    if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
      strcpy(parameters->photonlist_template, buffer);
    } else {
      if ((status = PILGetFname("fits_templates", parameters->photonlist_template))) {
	HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);
      }
    }
  }
  if (EXIT_SUCCESS!=status) return(status);
  // Set the event list template file for eROSITA:
  strcat(parameters->photonlist_template, "/photonlist.tpl");

  return(status);
}



int ero_merge_photonfiles_main() {
  struct Parameters parameters;
  // Array of input photonlist files.
  PhotonListFile* inputfiles[MAX_N_INPUTFILES];
  int fileidx;
  // Output (merged) photonlist file.
  PhotonListFile* outputfile=NULL;
  // Buffer for the photons.
  Photon photons[MAX_N_INPUTFILES];

  int status = EXIT_SUCCESS;

  // Register HEATOOL
  set_toolname("ero_merge_photonfiles");
  set_toolversion("0.01");

  do { // ERROR handling loop

    // Read parameters by PIL:
    status = ero_merge_photonfiles_getpar(&parameters);
    if (EXIT_SUCCESS!=status) break;

    // Open the INPUT event files:
    char filename[MAXMSG];
    for (fileidx=0; fileidx<parameters.n_inputfiles; fileidx++) {
      sprintf(filename, "%s%d.fits", 
	      parameters.input_prefix, fileidx);
      inputfiles[fileidx]=openPhotonListFile(filename, READONLY, &status);
      if (EXIT_SUCCESS!=status) break;
    }
    if (EXIT_SUCCESS!=status) break;

    // Create and open a new output (merged) event file:
    outputfile=openNewPhotonListFile(parameters.output_filename, 
				     parameters.photonlist_template,
				     &status);
    if (EXIT_SUCCESS!=status) break;


    // Copy header keywords.
    // Read the keywords from the first PhotonListFile and write
    // them to the output file.
    struct HKeys {
      double refxcrvl;
      double refycrvl;
      char attitude[MAXMSG];
    } hkeys;

    // Read from the first input file.
    char comment[MAXMSG]; // String buffer.
    if (fits_read_key(inputfiles[0]->fptr, TDOUBLE, "REFXCRVL", 
		      &hkeys.refxcrvl, comment, &status)) break;    
    if (fits_read_key(inputfiles[0]->fptr, TDOUBLE, "REFYCRVL", 
		      &hkeys.refycrvl, comment, &status)) break;    
    if (fits_read_key(inputfiles[0]->fptr, TSTRING, "ATTITUDE", 
		      hkeys.attitude, comment, &status)) break;

    // Write to output file.
    if (fits_update_key(outputfile->fptr, TDOUBLE, "REFXCRVL", 
			&hkeys.refxcrvl, "", &status)) break;
    if (fits_update_key(outputfile->fptr, TDOUBLE, "REFYCRVL", 
			&hkeys.refycrvl, "", &status)) break;
    if (fits_update_key(outputfile->fptr, TSTRING, "ATTITUDE", 
			hkeys.attitude, "name of the attitude FITS file", 
			&status)) break;
    // END of copying header keywords.


    // Transfer all photons from the input files to the output file.
    int eof[MAX_N_INPUTFILES];
    int sum_eof= 0;
    // Read the first entry from each photon list file:
    for (fileidx=0; fileidx<parameters.n_inputfiles; fileidx++) {
      if (inputfiles[fileidx]->nrows>0) {
	status = PhotonListFile_getNextRow(inputfiles[fileidx], 
					   &(photons[fileidx]));
	if (status!=EXIT_SUCCESS) break;
	eof[fileidx]=0;
      } else {
	eof[fileidx]=1;
	sum_eof++;
      }
    }
    // Repeat this loop as long as at least one the input event 
    // files contains furhter un-read lines.
    int minidx;
    while(sum_eof<parameters.n_inputfiles) {

      // Find the photon with the lowest time value.
      for (minidx=0; minidx<parameters.n_inputfiles; minidx++) {
	if (0==eof[minidx]) break;
      }
      for (fileidx=minidx+1; fileidx<parameters.n_inputfiles; fileidx++) {
	if (photons[fileidx].time < photons[minidx].time) {
	  minidx = fileidx;
	}
      }

      // Add the photon to the output file.
      status = addPhoton2File(outputfile, &photons[minidx]);
      if (EXIT_SUCCESS!=status) break;

      // Read new photon from the respective input file.
      if (inputfiles[minidx]->row<inputfiles[minidx]->nrows) {
	status = PhotonListFile_getNextRow(inputfiles[minidx], 
					   &(photons[minidx]));
	if (status!=EXIT_SUCCESS) break;
      } else {
	eof[minidx]=1;
	sum_eof++;
      }

    }
    if (status!=EXIT_SUCCESS) break;
    // END of event transfer loop.
        
  } while(0); // End of error handling loop


  // --- Clean Up ---
  
  // Close the event files:
  for (fileidx=0; fileidx<parameters.n_inputfiles; fileidx++) {
    freePhotonListFile(&(inputfiles[fileidx]), &status);
  }
  freePhotonListFile(&outputfile, &status);

  return(status);
}


