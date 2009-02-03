//////////////////////////////////////////////////////////////////////////////////////
// This application makes detector images from an eventfile.                        //
//////////////////////////////////////////////////////////////////////////////////////
//
// @author       Christian Schmid
// @date         2008/04
// @param        eventlist - filename of the FITS file containing the eventlist
//                           (width of the detector is given in the FITS header)
// @param        plotpixelwidth - width of a detector pixel in the output .png-file
//
//////////////////////////////////////////////////////////////////////////////////////

#define DET_IMAGES 1  // plot detector images
//#define LIGHTCURVE 1  // plot lightcurve

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <png.h>
#include <math.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "event_list.h"
#include "imglib.h"

#define TOOLSUB plot_eventlist_main
#include "headas_main.c"

#define FILENAME_LENGTH 128 // maximum length of filenames
#define MAXMSG  256         // maximum length of an output/error message


// reads the program parameters using PIL
int plot_eventlist_getpar(char eventlist_filename[], int *plotpixelwidth, 
			  double *integration_time, double *timespan);

// does the actual work: open FITS file, read eventlist and create detector frame plots
int plot_eventlist_work(const char eventlist_filename[], const int plotpixelwidth, 
			double integration_time, double timespan);

// plots the content of a square array to a png file
int plot_array(double **array, int width, int plotpixelwidth, char filename[]);

// clears the given array
void clear_array(double **array, int width);


// main program
int plot_eventlist_main() {
  char eventlist_filename[FILENAME_LENGTH]; // filename of the eventlist FITS file
  int plotpixelwidth;                       // width of a pixel in the image file
  // "integration time" for this output program (can be bigger than the  
  // detector integration time for the measurement)
  double integration_time;                  
  double timespan;                          // period for the event list evaluation

  int status=0;                             // error status

  
  // register HEATOOL
  set_toolname("plot_eventlist");
  set_toolversion("0.01");

  
  // read parameters using PIL
  status = plot_eventlist_getpar(eventlist_filename, &plotpixelwidth, 
				 &integration_time, &timespan);

  if (!status) {
    // call the routine which performs the actual work: load eventlist from FITS file,
    // create detector frames and plot them
    status = plot_eventlist_work(eventlist_filename, plotpixelwidth, 
				 integration_time, timespan);

    headas_chat(5, "finished\n");
  }  
  
  return(status);
}




////////////////////////////////////////////////////////////////////////////////////
// reads the program parameters using PIL
int plot_eventlist_getpar(
			  // filename of FITS file containing the eventlist
			  char eventlist_filename[],   
			  // width of a detector pixel in the output .png-file
			  int *plotpixelwidth,         
			  double *integration_time,
			  double *timespan
			  )
{
  int status=0;           // error status
  char msg[MAXMSG];       // error message

  if ((status = PILGetFname("eventlist", eventlist_filename))) {
    sprintf(msg, "Error reading the 'eventlist' parameter!");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetInt("plotpixelwidth", plotpixelwidth))) {
    sprintf(msg, "Error reading the 'plotpixelwidth' parameter!");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetReal("integration_time", integration_time))) {
    sprintf(msg, "Error reading the integration time!");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetReal("timespan", timespan))) {
    sprintf(msg, "Error reading the timespan!");
    HD_ERROR_THROW(msg,status);
  }

  return(status);
}




////////////////////////////////////////////////////////////////////////////////////
// does the actual work: open FITS file, read eventlist and create detector frame plots
int plot_eventlist_work(
			// filename of the FITS file containing the eventlist
			const char input_filename[], 
			// width of a detector pixel in the output .png-file
			const int plotpixelwidth,        
			const double integration_time,
			const double timespan
			)
{ 
  struct Event_List_File event_list_file;
  char image_filename[FILENAME_LENGTH];    // filename of the output image file
  // filename format of the output image file (.png)
  const char image_filename_f[]="det%lf.png";  

  int x;                 // detector pixel coordinate counters
  int det_width;         // width of the detector (number of pixels)
  double **det=NULL;     // detector array
  long lightcurve_bin;   // current bin of the lightcurve

  int status=0;           // error status
  char msg[MAXMSG];       // buffer for error messages


  do {     // error handling loop (only run once)

    // open event list FITS file
    strcpy(event_list_file.filename, input_filename);
    headas_chat(5,"open event list '%s' ...\n", event_list_file.filename);

    if (fits_open_table(&event_list_file.fptr, event_list_file.filename, READONLY, &status)) break;

    // get the HDU type
    int hdutype;
    if (fits_get_hdu_type(event_list_file.fptr, &hdutype, &status)) break;

    // image HDU results in an error message
    if (hdutype==IMAGE_HDU) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: no table extension available in event list "
	      "FITS file '%s'!\n", event_list_file.filename);
      HD_ERROR_THROW(msg,status);
      break;
    }
   
    
    // Determine the width of the detector (number of pixels) 
    // from the header information
    char comment[MAXMSG];   // input buffer for header comment
    if (fits_read_key(event_list_file.fptr, TINT, "DETWIDTH", &det_width,
		      comment, &status)) 
      break;

#ifdef DET_IMAGES
    // get memory for the detector array
    det = (double **)malloc(det_width*sizeof(double*));
    if (!det) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory to store detector array!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
    for (x=0;x<det_width;x++) {
      det[x] = (double *)malloc(det_width*sizeof(double));
      if (!det[x]) {
	status=EXIT_FAILURE;
	sprintf(msg, "Error: not enough memory to store detector array!\n");
	HD_ERROR_THROW(msg,status);
	break;
      }
    }
#endif

    // determine number of rows in the event list
    if (fits_get_num_rows(event_list_file.fptr, &event_list_file.nrows, &status)) 
      break;
	

    struct Event event; // input buffer
    double last_time;   // time of the former event in the list
    double t0;        // time of the first photon event
    int first=1;      // flag to determine, if the handled frame is the very first

    event.time = 0.; t0 = 0.;

    // loop over all events in the list
    headas_chat(5, "processing events ...\n");
    for (event_list_file.row=0; 
	 (event_list_file.row<event_list_file.nrows)&&(event.time-t0<timespan); 
	 event_list_file.row++) {
      if(get_eventtbl_row(event_list_file, &event, &status)) break;

      while ((event.time-last_time >= integration_time)||(first)) {
	if(!first) {

#ifdef DET_IMAGES
	  sprintf(image_filename, image_filename_f, last_time);
	  headas_chat(1, "%s,", image_filename);

	  // plot the current detector image
	  sprintf(image_filename, image_filename_f, last_time);
	  plot_array(det, det_width, plotpixelwidth, image_filename);
#endif

#ifdef LIGHTCURVE
	  // plot the current lightcurve bin
	  headas_printf("%lf %ld\n", last_time, lightcurve_bin);
#endif
	  // update variables
	  last_time += integration_time; //= time;    // save time of the actual frame
	} else {
	  t0 = event.time;
	  first=0;             // now this is not the first frame any more
	  last_time = event.time;
	}
		 
#ifdef DET_IMAGES
	// clear the detector image
	clear_array(det, det_width);
#endif

#ifdef LIGHTCURVE
	lightcurve_bin = 0;
#endif
      }

#ifdef DET_IMAGES
      // add event to detector array
      det[event.xi][event.yi] += 1.; // +=;
#endif

#ifdef LIGHTCURVE
      lightcurve_bin++;
#endif
    } // end of loop over all events



#ifdef DET_IMAGES
    // create a final plot
    sprintf(image_filename, image_filename_f, event.time);
    headas_chat(1, "%s,", image_filename);
    plot_array(det, det_width, plotpixelwidth, image_filename);

    // plot detector array to stdout
    //    for (xi=178;xi<207 /*det_width*/;xi++) {
    //      for (yi=178;yi<207 /*det_width*/;yi++) {
    //	headas_chat(1, "%d %d %lf\n", xi, yi, det[xi][yi]);
    //      }
    //      headas_chat(1, "\n");
    //    }
#endif


#ifdef LIGHTCURVE
    // final output
    headas_printf("%lf %ld\n", last_time, lightcurve_bin);
#endif

  } while (0);   // end of error handling loop
 
  // clean up:
  headas_chat(5, "cleaning up ...\n");
  
#ifdef DET_IMAGES
  // release memory of detector array
  if (det) {
    for (x=0;x<det_width;x++) {
      if (det[x]) free(det[x]);
    }
    free(det);
  }
#endif

  // close FITS file
  if(event_list_file.fptr) fits_close_file(event_list_file.fptr, &status);

  return(status);
}




////////////////////////////////////////////////////////////////////////////////////
// Function plots the content of a square array to a png file.
int plot_array(double **array, int width, int plotpixelwidth, char filename[]) {
  int x, y;            // counters
  double maximum = 0.; // value of the brightest point in the array

  // variables for png-access:
  bitmap_t bitmap;
  pixel_t pixel;

  // initialize image:
  bitmap.width = width*plotpixelwidth;
  bitmap.height = width*plotpixelwidth;
  bitmap.name = filename;

  if(bmp_init(&bitmap) !=IMGLIB_OK){
    printf("Bitmap initialization failed\n");
    return(IMGLIB_ERROR);
  }

  
  // find the brightest pixel in the array (to get the right scale of the picture):
  for(x=0; x < width; x++) {
    for(y=0; y < width; y++) {
      if(array[x][y] > maximum) {
        maximum = array[x][y];
      }
    }
  }  


  // create the png-plot:
  for(x=0; x < bitmap.height; x++) {
    for(y=0; y < bitmap.width; y++) {

      if (array[x/plotpixelwidth][y/plotpixelwidth] > 0.) {
	// photon event
	pixel.r = //(int) (255. *(array[x/plotpixelwidth][y/plotpixelwidth]/maximum));
	          255 - (int) (255. *(array[x/plotpixelwidth][y/plotpixelwidth]/maximum));
	pixel.g = 255 - (int) (255. *(array[x/plotpixelwidth][y/plotpixelwidth]/maximum));
	pixel.b = 255 - (int) (255. *(array[x/plotpixelwidth][y/plotpixelwidth]/maximum));
      } else if ((int)(sqrt(pow((double)(x-width/2*plotpixelwidth),2.) + pow((double)(y-width/2*plotpixelwidth),2.))) 
		 >= width/2*plotpixelwidth) {
	pixel.r = 0;    // border of the FOV:
	pixel.g = 0;    // black
	pixel.b = 0;
      } else {
	pixel.r = 255;  // background pixel:
	pixel.g = 255;  // white
	pixel.b = 255;
      }

      if(bmp_set_pixel_rgb(&bitmap, &pixel, x, y)){
        printf("Failed to set RGB of pixel (%d,%d)\n",x,y);
        return(IMGLIB_ERROR);
      }
    }
  }

  if(png_write(&bitmap)!=IMGLIB_OK){
    printf("PNG output failed\n");
    return(IMGLIB_ERROR);
  }

  return(IMGLIB_OK);
}



// Function clears the given array:
void clear_array(double **array, int width) {
  int x,y;

  for (x=0; x<width; x++) {
    for (y=0; y<width; y++) {
      array[x][y] = 0.0;
    }
  }
}

