#include "photonmerge.h"


int ero_merge_photonfiles_getpar(struct Parameters* par)
{
  char* sbuffer=NULL;

  int status = EXIT_SUCCESS;

  status=ape_trad_query_int("NInputFiles", &par->NInputFiles);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the number of input files!\n", status);
    return(status);
  }

  status=ape_trad_query_string("InputPrefix", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the prefix of the input files!\n", status);
    return(status);
  } 
  strcpy(par->InputPrefix, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("OutputFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the outpuf file!\n", status);
    return(status);
  } 
  strcpy(par->OutputFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }

  return(status);
}



int photonmerge_main() {
  struct Parameters par;
  // Array of input photonlist files.
  PhotonListFile* inputfiles[MAX_N_INPUTFILES];
  int fileidx;
  // Output (merged) photonlist file.
  PhotonListFile* outputfile=NULL;
  // Buffer for the photons.
  Photon photons[MAX_N_INPUTFILES];

  int status = EXIT_SUCCESS;

  // Register HEATOOL
  set_toolname("photonmerge");
  set_toolversion("0.01");

  do { // ERROR handling loop

    // Initialize with NULL.
    for (fileidx=0; fileidx<MAX_N_INPUTFILES; fileidx++) {
      inputfiles[fileidx]=NULL;
    }

    // Read parameters by PIL:
    status = ero_merge_photonfiles_getpar(&par);
    if (EXIT_SUCCESS!=status) break;

    // Open the INPUT event files:
    char filename[MAXMSG];
    for (fileidx=0; fileidx<par.NInputFiles; fileidx++) {
      sprintf(filename, "%s%d.fits", 
	      par.InputPrefix, fileidx);
      inputfiles[fileidx]=openPhotonListFile(filename, READONLY, &status);
      if (EXIT_SUCCESS!=status) break;
    }
    if (EXIT_SUCCESS!=status) break;

    // Create and open a new output photon file:
    outputfile=openNewPhotonListFile(par.OutputFile, &status);
    if (EXIT_SUCCESS!=status) break;


    // Copy header keywords.
    // Read the keywords from the first PhotonListFile and write
    // them to the output file.
    struct HKeys {
      double mjdref;
      double timezero;
    } hkeys;

    // Read from the first input file.
    char comment[MAXMSG]; // String buffer.
    if (fits_read_key(inputfiles[0]->fptr, TDOUBLE, "MJDREF", 
		      &hkeys.mjdref, comment, &status)) break;    
    if (fits_read_key(inputfiles[0]->fptr, TDOUBLE, "TIMEZERO", 
		      &hkeys.timezero, comment, &status)) break;    

    // Write to output file.
    if (fits_update_key(outputfile->fptr, TDOUBLE, "MJDREF", 
			&hkeys.mjdref, "", &status)) break;
    if (fits_update_key(outputfile->fptr, TDOUBLE, "TIMEZERO", 
			&hkeys.timezero, "", &status)) break;
    // END of copying header keywords.


    // Transfer all photons from the input files to the output file.
    int eof[MAX_N_INPUTFILES];
    int sum_eof= 0;
    // Read the first entry from each photon list file:
    for (fileidx=0; fileidx<par.NInputFiles; fileidx++) {
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
    // files contains further un-read lines.
    int minidx;
    while(sum_eof<par.NInputFiles) {

      // Find the photon with the lowest time value.
      for (minidx=0; minidx<par.NInputFiles; minidx++) {
	if (0==eof[minidx]) break;
      }
      for (fileidx=minidx+1; fileidx<par.NInputFiles; fileidx++) {
	if (0==eof[fileidx]) {
	  if (photons[fileidx].time < photons[minidx].time) {
	    minidx = fileidx;
	  }
	}
      }

      // Reset the photon ID.
      photons[minidx].ph_id = 0;

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
  for (fileidx=0; fileidx<par.NInputFiles; fileidx++) {
    freePhotonListFile(&(inputfiles[fileidx]), &status);
  }
  freePhotonListFile(&outputfile, &status);

  return(status);
}


