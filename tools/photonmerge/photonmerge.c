#include "photonmerge.h"


int ero_merge_photonfiles_getpar(struct Parameters* par)
{
  char* sbuffer=NULL;

  int status=EXIT_SUCCESS;

  status=ape_trad_query_int("NInputFiles", &par->NInputFiles);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the number of input files");
    return(status);
  }

  status=ape_trad_query_string("InputPrefix", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the input file prefix");
    return(status);
  } 
  strcpy(par->InputPrefix, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("OutputFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the outpuf file");
    return(status);
  } 
  strcpy(par->OutputFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
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

  int status=EXIT_SUCCESS;

  // Register HEATOOL
  set_toolname("photonmerge");
  set_toolversion("0.01");

  do { // ERROR handling loop

    // Initialize with NULL.
    for (fileidx=0; fileidx<MAX_N_INPUTFILES; fileidx++) {
      inputfiles[fileidx]=NULL;
    }

    // Read parameters by PIL:
    status=ero_merge_photonfiles_getpar(&par);
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
    char telescop[MAXMSG], instrume[MAXMSG], filter[MAXMSG];
    char ancrfile[MAXMSG], respfile[MAXMSG];
    char comment[MAXMSG];
    fits_read_key(inputfiles[0]->fptr, TSTRING, "TELESCOP",
		  &telescop, comment, &status);
    fits_read_key(inputfiles[0]->fptr, TSTRING, "INSTRUME",
		  &instrume, comment, &status);
    fits_read_key(inputfiles[0]->fptr, TSTRING, "FILTER",
		  &filter, comment, &status);
    fits_read_key(inputfiles[0]->fptr, TSTRING, "ANCRFILE", 
		  &ancrfile, comment, &status);
    fits_read_key(inputfiles[0]->fptr, TSTRING, "RESPFILE", 
		  &respfile, comment, &status);
    CHECK_STATUS_BREAK(status);

    double mjdref, timezero, tstart, tstop;
    fits_read_key(inputfiles[0]->fptr, TDOUBLE, "MJDREF", 
		  &mjdref, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(inputfiles[0]->fptr, TDOUBLE, "TIMEZERO", 
		  &timezero, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(inputfiles[0]->fptr, TDOUBLE, "TSTART", 
		  &tstart, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(inputfiles[0]->fptr, TDOUBLE, "TSTOP", 
		  &tstop, comment, &status);
    CHECK_STATUS_BREAK(status);
    outputfile=openNewPhotonListFile(par.OutputFile, 
				     telescop, instrume, filter,
				     ancrfile, respfile,
				     mjdref, timezero, tstart, tstop,
				     par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Copy header keywords.
    char attitude[MAXMSG];
    fits_read_key(inputfiles[0]->fptr, TSTRING, "ATTITUDE", 
		  attitude, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_update_key(outputfile->fptr, TSTRING, "ATTITUDE",
		    attitude, comment, &status);
    CHECK_STATUS_BREAK(status);


    // Transfer all photons from the input files to the output file.
    int eof[MAX_N_INPUTFILES];
    int sum_eof= 0;
    // Read the first entry from each photon list file:
    for (fileidx=0; fileidx<par.NInputFiles; fileidx++) {
      if (inputfiles[fileidx]->nrows>0) {
	status=PhotonListFile_getNextRow(inputfiles[fileidx], 
					 &(photons[fileidx]));
	CHECK_STATUS_BREAK(status);
	eof[fileidx]=0;
      } else {
	eof[fileidx]=1;
	sum_eof++;
      }
    }
    CHECK_STATUS_BREAK(status);
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
      status=addPhoton2File(outputfile, &photons[minidx]);
      CHECK_STATUS_BREAK(status);

      // Read new photon from the respective input file.
      if (inputfiles[minidx]->row<inputfiles[minidx]->nrows) {
	status=PhotonListFile_getNextRow(inputfiles[minidx], 
					 &(photons[minidx]));
	CHECK_STATUS_BREAK(status);
      } else {
	eof[minidx]=1;
	sum_eof++;
      }

    }
    CHECK_STATUS_BREAK(status);
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


