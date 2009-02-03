/////////////////////////////////////////////////////////////////////////////////////////
// This application creates psf images from the psf event list.                        //
/////////////////////////////////////////////////////////////////////////////////////////
//
// @author            Christian Schmid
// @date              2008/04
// @param             psf - filename of the ASCII file containing the psf photon event list
// @param             detwidth - with of the detector in pixels                        
// @param             plotpixelwidth - width of a detector pixel in the output .png-file
//
/////////////////////////////////////////////////////////////////////////////////////////

#include <png.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "imglib.h"
#include "strftcpy.h"

#define TOOLSUB plot_psf_main
#include "headas_main.c"

#define FILENAME_LENGTH 128 // maximum length of filenames
#define LINELENGTH 1024     // maximum linelength in the ASCII file
#define MAXMSG  256         // maximum length of an output/error message


// reads the program parameters using PIL
int plot_psf_getpar(char psf_filename[], int *detwidth, int *plotpixelwidth);

// does the actual work: open ASCII file, read photon eventlist and create psf detector image
int plot_psf_work(const char psf_filename[], const int detwidth, const int plotpixelwidth);

// plots the content of a square array to a png file
int plot_array(int **array, int width, int plotpixelwidth, char filename[]);

// clears the given array
void clear_array(int **array, int width);


// main program
int plot_psf_main() {
  int detwidth;                             // width of the detector in pixels
  int plotpixelwidth;                       // width of a pixel in the image file
  char psf_filename[FILENAME_LENGTH];       // filename of the psf ASCII file

  int status=0;                             // error status

  
  // register HEATOOL
  set_toolname("plot_psf");
  set_toolversion("0.01");

  
  // read parameters using PIL
  status = plot_psf_getpar(psf_filename, &detwidth, &plotpixelwidth);

  if (!status) {
    // call the routine which performs the actual work: load psf data from ASCII file,
    // and create a picture
    status = plot_psf_work(psf_filename, detwidth, plotpixelwidth);

    headas_chat(5, "finished\n");
  }  
  
  return(status);
}




////////////////////////////////////////////////////////////////////////////////////
// reads the program parameters using PIL
int plot_psf_getpar(
		    char psf_filename[],         // filename of psf ASCII file containing the eventlist
		    int *detwidth,               // width of the detector in pixels
		    int *plotpixelwidth          // width of a detector pixel in the output .png-file
		    )
{
  int status=0;           // error status
  char msg[MAXMSG];       // error message

  if ((status = PILGetFname("psf", psf_filename))) {
    sprintf(msg, "Error reading the filename of the psf data file");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetInt("detwidth", detwidth))) {
    sprintf(msg, "Error reading the detector width");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetInt("plotpixelwidth", plotpixelwidth))) {
    sprintf(msg, "Error reading the 'plotpixelwidth' parameter");
    HD_ERROR_THROW(msg,status);
  }

  return(status);
}




////////////////////////////////////////////////////////////////////////////////////
// does the actual work: open FITS file, read eventlist and create detector frame plots
int plot_psf_work(
		  const char psf_filename[],       // filename of the psf ASCII file
		  const int detwidth,              // width of the detector in pixels
		  const int plotpixelwidth         // width of a detector pixel in the output .png-file
		  )
{
  FILE *psf_fptr=NULL;                     // ASCII file pointer to input psf photon event list
  char image_filename[FILENAME_LENGTH];    // filename of the output image file
  const char image_filename_f[]="psf/psf.png";   // filename format of the output image file (.png)

  char line[LINELENGTH];                   // input buffer for ASCII file
  char *cbuffer[1]={NULL};                 // string buffer for ASCII input
  int xi, yi;                              // detector pixel coordinate counters
  int **det=NULL;                          // detector array

  int status=0;                            // error status
  char msg[MAXMSG];                        // buffer for error messages


  do {     // error handling loop (only run once)

    // get memory for the detector array
    det = (int **)malloc(detwidth*sizeof(int*));
    if (!det) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory to store detector array!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
    for (xi=0;xi<detwidth;xi++) {
      det[xi] = (int *)malloc(detwidth*sizeof(int));
      if (!det[xi]) {
	status=EXIT_FAILURE;
	sprintf(msg, "Error: not enough memory to store detector array!\n");
	HD_ERROR_THROW(msg,status);
	break;
      }
    }

    // clear the detector image
    clear_array(det, detwidth);



    // open psf data file (ASCII)
    headas_chat(5,"open psf file '%s' ...\n", psf_filename);
    if (!(psf_fptr=fopen(psf_filename, "r+"))) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error opening the PSF data '%s' (ASCII)\n", psf_filename);
      HD_ERROR_THROW(msg,status);
      break;
    }


    // get memory for char input buffer
    cbuffer[0]=malloc(21*sizeof(char));
    if(!cbuffer[0]) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error allocating memory");
      HD_ERROR_THROW(msg,status);
      break;
    }

    // loop over all events in the list
    headas_chat(5, "processing photon events ...\n");
    while (fgets(line, LINELENGTH, psf_fptr)) {
      // parse ASCII line
      strftcpy(cbuffer[0],line,24,3);
      xi = atoi(cbuffer[0]);    // convert to integer number

      strftcpy(cbuffer[0],line,31,3);
      yi = atoi(cbuffer[0]);


      // add event to detector array
      det[xi][yi]++;
    } // end of loop over all events


    // plot the psf detector image
    sprintf(image_filename, image_filename_f);
    plot_array(det, detwidth, plotpixelwidth, image_filename);


    for(xi=309; xi<384; xi++) {
      for(yi=155; yi<230; yi++) {
	headas_chat(1, "%d %d %d\n", xi, yi, det[xi][yi]);
      }
      headas_chat(1, "\n");
    }


  } while (0);   // end of error handling loop
 


  // clean up:
  headas_chat(5, "cleaning up ...\n");

  // free memory of input buffer
  if(cbuffer[0]) free(cbuffer[0]);
  
  // release memory of detector array
  if (det) {
    for (xi=0;xi<detwidth;xi++) {
      if (det[xi]) free(det[xi]);
    }
    free(det);
  }

  // close psf file
  if(psf_fptr) fclose(psf_fptr);

  return(status);
}




////////////////////////////////////////////////////////////////////////////////////
// Function plots the content of a square array to a png file.
int plot_array(
	       int **array, 
	       int width, 
	       int plotpixelwidth, 
	       char filename[]
	       ) 
{
  int x, y;            // counters
  int maximum=0;       // value of the brightest point on the detector

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


  // find the brightest pixel in the array (to scale the brightness of the picture):
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
      pixel.r = (unsigned char)(255.*(double)array[x/plotpixelwidth][y/plotpixelwidth]/maximum);
      pixel.g = 0;
      pixel.b = 0;

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



////////////////////////////////////////////////////////
// Function clears the given array:
void clear_array(
		 int **array, 
		 int width
		 ) 
{
  int x,y;

  for (x=0; x<width; x++) {
    for (y=0; y<width; y++) {
      array[x][y] = 0;
    }
  }
}

