//////////////////////////////////////////////////////////////////////////////////////
// This application evaluates the event list and creates spectrum data from the     //
// individual photon events.                                                        //
//////////////////////////////////////////////////////////////////////////////////////
//
// @author      Christian Schmid
// @date        2008/07
// @param       eventlist - filename of the FITS file containing the event list
//
//////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "event_list.h"
#include "global_constants.h"
#include "detector.h"

#define TOOLSUB plot_spectrum_main
#include "headas_main.c"



// reads the program parameters using PIL
int plot_spectrum_getpar(char inputfile[], char rmffile[]);

// does the actual work: open FITS file, read eventlist and create detector frame plots
int plot_spectrum_work(const char inputfile[], const char rmffile[]);



// main program
int plot_spectrum_main() {
  char eventlist_filename[FILENAME_LENGTH]; // filename of the eventlist FITS file
  char rmffile[FILENAME_LENGTH];            

  int status=0;                             // error status

  
  // register HEATOOL
  set_toolname("plot_spectrum");
  set_toolversion("0.01");

  
  // read parameters using PIL
  status = plot_spectrum_getpar(eventlist_filename, rmffile);

  if (!status) {
    // call the routine which performs the actual work: load eventlist from FITS file,
    // create detector frames and plot them
    status = plot_spectrum_work(eventlist_filename, rmffile);

    headas_chat(5, "finished\n");
  }  
  
  return(status);
}




////////////////////////////////////////////////////////////////////////////////////
// reads the program parameters using PIL
int plot_spectrum_getpar(
			 char eventlist_filename[],   // filename of FITS file containing the eventlist
			 char rmffile[]               // filename of the RMF file 
			                              // (containing the energy bins in the EBOUNDS table)
			 )
{
  int status=0;           // error status
  char msg[MAXMSG];       // error message

  if ((status = PILGetFname("eventlist", eventlist_filename))) {
    sprintf(msg, "Error reading the 'eventlist' parameter!");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetFname("rmffile", rmffile))) {
    sprintf(msg, "Error reading the name of the RMF file");
    HD_ERROR_THROW(msg,status);
  }

  return(status);
}




////////////////////////////////////////////////////////////////////////////////////
// does the actual work: open FITS file, read eventlist and create detector frame plots
int plot_spectrum_work(
		       const char eventlist_filename[], // filename of the FITS file containing the eventlist
		       const char rmffile[]
		       )
{
  fitsfile *eventlist_fptr=NULL;         // fitsfile pointer to input FITS file containing the eventlist
  long eventlist_nrows, eventlist_row;   // number of rows and actual row in the eventlist FITS file
  struct Ebounds ebounds;                // energy bins of PHA channels
  int N_channels;                        // number of detector PHA channels
  float *spectrum=NULL;

  int status=0;                          // error status
  char msg[MAXMSG];                      // buffer for error messages


  do {     // error handling loop (only run once)

    // open eventlist FITS file
    headas_chat(5,"open eventlist '%s' ...\n", eventlist_filename);
    int eventlist_hdunum, eventlist_hdutype;
    if (fits_open_file(&eventlist_fptr, eventlist_filename, READONLY, &status)) break;

    // get the HDU number
    if (fits_get_hdu_num(eventlist_fptr, &eventlist_hdunum)==1) {
      // this is the primary array
      // try to move to the first extension and see if it is a table
      if (fits_movabs_hdu(eventlist_fptr, 2, &eventlist_hdutype, &status)) break;
    } else {
      // get the HDU type
      if (fits_get_hdu_type(eventlist_fptr, &eventlist_hdutype, &status)) break;
    }

    // image HDU results in an error message
    if (eventlist_hdutype==IMAGE_HDU) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: FITS extension in eventlist file '%s' is not a table but an image (HDU number: %d)\n", eventlist_filename, eventlist_hdunum);
      HD_ERROR_THROW(msg,status);
      break;
    }
    
    // determine the width of the detector (number of pixels) from the header information
    char comment[MAXMSG];   // input buffer for header comment
    if (fits_read_key(eventlist_fptr, TLONG, "DETCHANS", &N_channels, comment, &status)) break;

    // get memory for the detector array
    spectrum = (float *)malloc(N_channels * sizeof(float));
    if (!spectrum) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory to store spectrum array!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }

    // clear spectrum
    long count;
    for (count=0; count<N_channels; count++) {
      spectrum[count] = 0.;
    }


    // determine number of rows in the eventlist (i.e. number of listed sources)
    if (fits_get_num_rows(eventlist_fptr, &eventlist_nrows, &status)) break;
	

    double time;     // time of the event 
    int xi, yi;
    long pha;        // and corresponding PHA value
    int grade;       // the grade is determined according to the event pattern
    long frame, patnum, patid;

    // loop over all events in the list
    headas_chat(5, "processing events ...\n");
    for (eventlist_row=0; eventlist_row<eventlist_nrows; eventlist_row++) {
      if(get_eventtbl_row(eventlist_fptr, eventlist_row, &time, &pha, &grade, 
			  &xi, &yi, &frame, &patnum, &patid, &status)) break;

      // add photon to spectrum
      if ((pha >= 0) && (pha < N_channels)) {
	spectrum[pha] += 1./eventlist_nrows;
      }

    }  // end of loop over all events


    // rescale spectrum according to bin width:
    // therefore first get the energy bins of the PHA channels
    if ((status=get_ebounds(&ebounds, &N_channels, rmffile))!=EXIT_SUCCESS) break;
    // now the rescaling can be performed
    for (pha=0; pha<N_channels; pha++) {
      spectrum[pha] = spectrum[pha]/(ebounds.row[pha].E_max -ebounds.row[pha].E_min);
    }


    // plot spectrum to stdout
    for (count=0; count<N_channels; count++) {
      printf("%ld\t%lf\n", count, spectrum[count]);
    }

  } while (0);   // end of error handling loop
 
  // clean up:
  headas_chat(5, "cleaning up ...\n");
  
  // release memory
  free_ebounds(ebounds);
  if (spectrum) {
    free(spectrum);
  }

  // close FITS file
  if(eventlist_fptr) fits_close_file(eventlist_fptr, &status);

  return(status);
}


