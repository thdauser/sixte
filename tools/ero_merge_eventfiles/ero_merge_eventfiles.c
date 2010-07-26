#include "ero_merge_eventfiles.h"


int ero_merge_eventfiles_getpar(struct Parameters* parameters)
{
  int status = EXIT_SUCCESS;

  if ((status = PILGetFname("input_prefix", parameters->input_prefix))) {
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
      strcpy(parameters->eventlist_template, buffer);
    } else {
      if ((status = PILGetFname("fits_templates", parameters->eventlist_template))) {
	HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);
      }
    }
  }
  if (EXIT_SUCCESS!=status) return(status);
  // Set the event list template file for eROSITA:
  strcat(parameters->eventlist_template, "/erosita.eventlist.tpl");

  return(status);
}



int ero_merge_eventfiles_main() {
  struct Parameters parameters;
  // Array of input event files.
  eROSITAEventFile inputfiles[7];
  int filecounter;
  // Output (merged) event file.
  eROSITAEventFile outputfile;

  int status = EXIT_SUCCESS;

  // Register HEATOOL
  set_toolname("ero_merge_eventfiles");
  set_toolversion("0.01");

  do { // ERROR handling loop

    // Read parameters by PIL:
    status = ero_merge_eventfiles_getpar(&parameters);
    if (EXIT_SUCCESS!=status) break;

    // Open the INPUT event files:
    char filename[MAXMSG];
    for (filecounter=0; filecounter<7; filecounter++) {
      sprintf(filename, "%s%d.fits", 
	      parameters.input_prefix, filecounter);
      status = openeROSITAEventFile(&inputfiles[filecounter], 
				    filename, 
				    READONLY);
      if (EXIT_SUCCESS!=status) break;
    }
    if (EXIT_SUCCESS!=status) break;

    // Create and open a new output (merged) event file:
    status = openNeweROSITAEventFile(&outputfile, 
				     parameters.output_filename, 
				     parameters.eventlist_template);
    if (EXIT_SUCCESS!=status) break;


    // Copy header keywords.
    // Read the keywords from the first event file (for telescope 0)
    // and write them to the common (merged) event file.
    headas_chat(5, "Copy header keywords ...\n");
    struct HKeys {
      char attitude[MAXMSG];
      
      // Sky WCS keywords.
      double tcrvlx, tcdltx, tcrpxx;
      double tcrvly, tcdlty, tcrpxy;
      // Detector WCS keywords.
      double tcdltrawx, tcdltrawy;

    } hkeys;
    // Read from the first input event file and write to the output file.
    char comment[MAXMSG];
    char keyword[MAXMSG];
    // Attitude.
    if (fits_read_key(inputfiles[0].generic.fptr, TSTRING, "ATTITUDE", 
		      hkeys.attitude, comment, &status)) break;
    if (fits_update_key(outputfile.generic.fptr, TSTRING, "ATTITUDE", 
			hkeys.attitude, comment, &status)) break;
    // The sky "X" column.
    sprintf(keyword, "TCRVL%d", inputfiles[0].cskyx);
    if (fits_read_key(inputfiles[0].generic.fptr, TDOUBLE, keyword, &hkeys.tcrvlx, 
		      comment, &status)) break; 
    if (fits_update_key(outputfile.generic.fptr, TDOUBLE, keyword, &hkeys.tcrvlx, 
			comment, &status)) break;    
    sprintf(keyword, "TCDLT%d", inputfiles[0].cskyx);
    if (fits_read_key(inputfiles[0].generic.fptr, TDOUBLE, keyword, &hkeys.tcdltx, 
		      comment, &status)) break;
    if (fits_update_key(outputfile.generic.fptr, TDOUBLE, keyword, &hkeys.tcdltx, 
			comment, &status)) break;
    sprintf(keyword, "TCRPX%d", inputfiles[0].cskyx);
    if (fits_read_key(inputfiles[0].generic.fptr, TDOUBLE, keyword, &hkeys.tcrpxx, 
		      comment, &status)) break;
    if (fits_update_key(outputfile.generic.fptr, TDOUBLE, keyword, &hkeys.tcrpxx, 
			comment, &status)) break;
    // The sky "Y" column.
    sprintf(keyword, "TCRVL%d", inputfiles[0].cskyy);
    if (fits_read_key(inputfiles[0].generic.fptr, TDOUBLE, keyword, &hkeys.tcrvly, 
		      comment, &status)) break;
    if (fits_update_key(outputfile.generic.fptr, TDOUBLE, keyword, &hkeys.tcrvly, 
			comment, &status)) break;
    sprintf(keyword, "TCDLT%d", inputfiles[0].cskyy);
    if (fits_read_key(inputfiles[0].generic.fptr, TDOUBLE, keyword, &hkeys.tcdlty,
		      comment, &status)) break;
    if (fits_update_key(outputfile.generic.fptr, TDOUBLE, keyword, &hkeys.tcdlty,
			comment, &status)) break;
    sprintf(keyword, "TCRPX%d", inputfiles[0].cskyy);
    if (fits_read_key(inputfiles[0].generic.fptr, TDOUBLE, keyword, &hkeys.tcrpxy,
		      comment, &status)) break;
    if (fits_update_key(outputfile.generic.fptr, TDOUBLE, keyword, &hkeys.tcrpxy,
			comment, &status)) break;

    // Detector coordinates.
    sprintf(keyword, "TCDLT%d", inputfiles[0].crawx);
    if (fits_read_key(inputfiles[0].generic.fptr, TDOUBLE, keyword, &hkeys.tcdltrawx,
		      comment, &status)) break;
    if (fits_update_key(outputfile.generic.fptr, TDOUBLE, keyword, &hkeys.tcdltrawx,
			comment, &status)) break;

    sprintf(keyword, "TCDLT%d", inputfiles[0].crawy);
    if (fits_read_key(inputfiles[0].generic.fptr, TDOUBLE, keyword, &hkeys.tcdltrawy,
		      comment, &status)) break;
    if (fits_update_key(outputfile.generic.fptr, TDOUBLE, keyword, &hkeys.tcdltrawy,
			comment, &status)) break;
    // END of copying header keywords.


    // Transfer all events from the input files to the output file.
    headas_chat(5, "Merge events to 1 file ...\n");
    long ntotal_lines=0;
    for (filecounter=0; filecounter<7; filecounter++) {
      inputfiles[filecounter].generic.row = 1;
      ntotal_lines += inputfiles[filecounter].generic.nrows;
    }
    headas_chat(5, "Total number of events %ld\n", ntotal_lines);

    // Buffer for the events.
    eROSITAEvent event;
    // Current frame number.
    long frame=0, min_next_frame=0;
    filecounter=0;
    int eof[7] = { 0, 0, 0, 0, 0, 0, 0 };
    int sum_eof= 0;
    // Repeat this loop as long as at least one the input event files contains
    // furhter un-read lines.
    while(sum_eof<7) {

      // Check if the end of the current file is reached.
      if (inputfiles[filecounter].generic.row>inputfiles[filecounter].generic.nrows) {
	if (0==eof[filecounter]) {
	  eof[filecounter]=1;
	  sum_eof++;
	  if (7==sum_eof) break;
	}
	// Increase the file counter.
	filecounter++;
	if (filecounter>=7) {
	  filecounter=0;
	  assert(min_next_frame>frame);
	  frame=min_next_frame;
	}
      } else {
	// Read the current event from this file.
	status=eROSITAEventFile_getRow(&inputfiles[filecounter], 
				       &event, 
				       inputfiles[filecounter].generic.row);
	if (status!=EXIT_SUCCESS) break;
	
	if (event.frame > frame) {
	  if ((min_next_frame==frame) || (event.frame < min_next_frame)) {
	    min_next_frame = event.frame;
	  }
	  // Move to next file.
	  filecounter++;
	  if (filecounter>=7) {
	    filecounter=0;
	    assert(min_next_frame>frame);
	    frame=min_next_frame;
	  }
	} else {
	  // Assign the right CCD number to the event.
	  event.ccdnr = filecounter+1;

	  // Add the event to the output file
	  status=addeROSITAEvent2File(&outputfile, &event);
	  if (status!=EXIT_SUCCESS) break;
	  inputfiles[filecounter].generic.row++;

	  // Status output.
	  if (0==outputfile.generic.nrows%1000) {
	    headas_printf("\revent %d/%d (%.1lf%%) ", 
			  outputfile.generic.nrows, ntotal_lines,
			  outputfile.generic.nrows*100./ntotal_lines);
	    fflush(NULL);
	  }
	}
      }
    }
    if (status!=EXIT_SUCCESS) break;
    // END of event transfer loop.
        
  } while(0); // End of error handling loop


  // --- Clean Up ---
  
  // Close the event files:
  for (filecounter=0; filecounter<7; filecounter++) {
    closeeROSITAEventFile(&inputfiles[filecounter]);
  }
  closeeROSITAEventFile(&outputfile);

  return(status);
}


