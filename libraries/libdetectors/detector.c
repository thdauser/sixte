#include "detector.h"



////////////////////////////////////////////////////////////////////////
Detector* get_Detector(int* status)
{
  Detector* detector=NULL;

  do { // Outer ERROR handling loop.

    headas_chat(5, "allocate memory for detector data structure ...\n");

    // Allocate memory for the detector:
    detector = (Detector*)malloc(sizeof(Detector));
    if (NULL==detector) {
      *status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: not enough memory available to store "
		     "the detector array!\n", *status);
      break;
    } else { // Memory has been allocated successfully.
      detector->pixel=NULL;
      detector->rmf=NULL;
      detector->specific=NULL;
    }

  } while (0); // END of Error handling loop.

  // Return a pointer to the newly created detector data structure.
  return(detector);
}




////////////////////////////////////////////////////////////////////////
int get_DetectorPixels(Detector* detector, int* status)
{

  do { // Outer ERROR handling loop

    headas_chat(5, "allocate memory for detector pixel array ...\n");

    // Allocate memory for the detector pixel array:
    int count;
    detector->pixel = (struct Pixel **) malloc(detector->width*sizeof(struct Pixel *));
    if (detector->pixel) {
      for (count=0; (count<detector->width)&&(*status==EXIT_SUCCESS); count++) {
	detector->pixel[count] = (struct Pixel *)
	  malloc(detector->width * sizeof(struct Pixel));
	
	if (!(detector->pixel[count])) {
	  *status = EXIT_FAILURE;
	  break;
	}
      }
    } else { *status = EXIT_FAILURE; }

    // Check if an error has occurred during memory allocation:
    if (*status==EXIT_FAILURE) {
      HD_ERROR_THROW("Error: not enough memory available to store "
		     "the detector array!\n", *status);
      break;
    }

    // Clear the detector array (at the beginning there are no charges).
    clear_detector(detector);

  } while (0); // END of Error handling loop

  return(*status);
}





////////////////////////////////////////////////////////////////////////////////
// Searches for the minimum distance in an array with 4 entries and returns the 
// corresponding index.
inline int min_dist(double array[], int directions) 
{
  int count, index=0;
  double minimum=array[0];

  for (count=1; count<directions; count++) {
    if ( (minimum < 0.) ||
	 ((array[count]<=minimum)&&(array[count]>=0.)) ) {
      minimum = array[count];
      index = count;
    }
  }

  return(index);
}





///////////////////////////////////////////////////////////////
int detector_assign_rsp(Detector *detector, char *filename) {
  fitsfile* fptr=NULL;

  int status=EXIT_SUCCESS;

  // Allocate memory:
  detector->rmf = (struct RMF*)malloc(sizeof(struct RMF));
  if (NULL==detector->rmf) {
    status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: could not allocate memory for RMF!\n", status);
    return(status);
  }

  // Load the RMF from the FITS file using the HEAdas RMF access routines
  // (part of libhdsp).
  fits_open_file(&fptr, filename, READONLY, &status);
  if (status != EXIT_SUCCESS) return(status);
  
  // Read the SPECRESP MATRIX or MATRIX extension:
  if ((status=ReadRMFMatrix(fptr, 0, detector->rmf)) != EXIT_SUCCESS) return(status);

  // Print some information:
  headas_chat(5, "RMF loaded with %ld energy bins and %ld channels\n",
	      detector->rmf->NumberEnergyBins, detector->rmf->NumberChannels);

#ifdef NORMALIZE_RMF
  // Normalize the RMF:
  headas_printf("### Warning: RMF is explicitly renormalized! ###\n");
  NormalizeRMF(detector->rmf);
#else
  // Check if the RSP file contains Matrix Rows with a sum of more than 1.
  // In that case the RSP probably also contains the mirror ARF, what should 
  // normally not be the case for this simulation.
  long chancount, bincount;
  double maxsum = 0.;
  for (bincount=0; bincount<detector->rmf->NumberEnergyBins ; bincount++) {
    double sum = 0.;
    for (chancount=0; chancount<detector->rmf->NumberChannels; chancount++) {
      sum += ReturnRMFElement(detector->rmf, chancount, bincount);
    }
    if (sum > maxsum) {
      maxsum = sum;
    }
  }
  if (maxsum > 1.) {
    headas_printf("### Warning: RSP probably contains mirror ARF (row-sum = %lf)! ###\n", maxsum);
  }
#endif

  // Read the EBOUNDS extension:
  if ((status=ReadRMFEbounds(fptr, 0, detector->rmf)) !=EXIT_SUCCESS) return(status);

  fits_close_file(fptr, &status);
  
  return(status);
}




////////////////////////////////////////////////////////////////
inline void clear_detector(Detector* detector) {
  int y;

  for (y=0; y<detector->width; y++) {
    clear_detector_line(detector, y);
  }
}



////////////////////////////////////////////////////////////////
inline void clear_detector_line(Detector* detector, int line) {
  int x;

  for (x=0; x<detector->width; x++) {
    detector->pixel[x][line].charge  = 0.;
    detector->pixel[x][line].arrival = -1.e-16; // TODO detector->dead_time;
  }
}




////////////////////////////////////////////////////////////////
int detector_active(
		    int x, int y,
		    Detector* detector,
		    double time
		    )
{
  if (detector->type == WFI) {
    if ((y==detector->readout_line)||(y==detector->width -detector->readout_line -1)) {
      // If we are at the beginning of the readout interval of the regarded line
      // within the clear time, the photon is not measured.
      if (time - detector->readout_time < detector->clear_time) {
	// The specified line is cleared at the moment, so no photon can 
	// be detected.
	return(0);
      }
    }

    // The specified detector pixel is active at the moment.
    return(1);

  } else {  // Default: detector unknown.
    // The specified detector pixel is active at the moment.
    return(1);
  }
}




////////////////////////////////////////////////////////////////
int htrs_detector_active(
			 int x, int y,
			 Detector* detector,
			 double time
			 )
{
  if (time - detector->pixel[x][y].arrival < detector->dead_time) {
    return(0);  // Pixel is INactive!
  } else {
    return(1);  // Pixel is active!
  }
}






////////////////////////////////////////////////////
long get_channel(
		 float energy, 
		 Detector* detector
		 )
{
  // Check if the charge is outside the range of the energy bins defined
  // in the EBOUNDS table. In that case the return value of this function is '-1'.
  if (detector->rmf->ChannelLowEnergy[0] > energy) {
    return(0); // TODO
  } else if (detector->rmf->ChannelHighEnergy[detector->rmf->NumberChannels-1] < energy) {
    return(detector->rmf->NumberChannels - 1 + detector->rmf->FirstChannel);
  }
  

  // Perform a binary search to obtain the detector PHA channel 
  // that corresponds to the given detector charge.
  long min, max, row;
  min = 0;
  max = detector->rmf->NumberChannels-1;
  while (max-min > 1) {
    row = (long)(0.5*(min+max));
    if (detector->rmf->ChannelHighEnergy[row] < energy) {
      min = row;
    } else {
      max = row;
    }
  }
  // Take the final decision wheter max or min is right:
  if (detector->rmf->ChannelLowEnergy[max] < energy) {
    row = max;
  } else {
    row = min;
  }
  
  // Return the PHA channel.
  return(row + detector->rmf->FirstChannel);
}





////////////////////////////////////////////////////
float get_energy(
		 long channel, 
		 Detector* detector
		 )
{
  channel -= detector->rmf->FirstChannel;
  if ((channel < 0) || (channel >= detector->rmf->NumberChannels)) {
    return(-1.);
  }

  // Return the mean of the energy that corresponds to the specified PHA channel
  // according to the EBOUNDS table.
  return(detector->rmf->ChannelLowEnergy[channel] +
	 get_random_number()*(detector->rmf->ChannelHighEnergy[channel]-
			      detector->rmf->ChannelLowEnergy[channel]));
}




///////////////////////////////////////////
inline double linear_function(double x, double m, double t)
{
  return(m*x + t);
}





///////////////////////////////////////////
inline int htrs_get_line(
				struct Point2d point, 
				double m, double dt, 
				Detector* detector
				)
{
  double dy = point.y - m * point.x;
  int line = (int)(dy/dt + detector->width + 1.) -1;

  // Check whether the point lies below the lowest or above the highest allowed line:
  if ((line < 0) || (line>= 2*detector->width)) line = INVALID_PIXEL;

  return(line);
}



///////////////////////////////////////////
inline double htrs_distance2line(
					struct Point2d point,
					double m, double t
					)
{
  return(sqrt( pow(t-point.y+m*point.x, 2.) / (1+pow(m,2.)) ));
}



///////////////////////////////////////////
inline void htrs_get_lines(
				  struct Point2d point, 
				  Detector* detector, 
				  int* l
				  )
{
  l[0] = htrs_get_line(point,        0.,    detector->h, detector);
  l[1] = htrs_get_line(point, -sqrt(3.), 2* detector->h, detector);
  l[2] = htrs_get_line(point,  sqrt(3.), 2* detector->h, detector);
}




///////////////////////////////////////////
int htrs_get_pixel(
		   Detector* detector, 
		   struct Point2d point,
		   int* x, int* y,
		   double* fraction
		   )
{
  int npixels = 0;
  int l[3];

  htrs_get_lines(point, detector, l);

  // Store the pixel coordinates:
  int pixel = htrs_get_lines2pixel(l, detector);
  struct Point2i pixel_coordinates = htrs_get_pixel2icoordinates(pixel, detector);
  x[0] = pixel_coordinates.x;
  y[0] = pixel_coordinates.y;

  
  // Check for split events:
  int dl[3][6] =  {
    {1, 0, 0, -1, 0, 0},
    {0, 1, 0, 0, -1, 0},
    {0, 0, -1, 0, 0, 1}
  };

  // Distances to neighbouring pixel segments (equilateral triangles, CAN possibly
  // belong to the same pixel).
  double distances[6] = {
    // upper
    htrs_distance2line(point,        0., (l[0]+1 -detector->width)*   detector->h),
    // upper right
    htrs_distance2line(point, -sqrt(3.), (l[1]+1-detector->width)*2.*detector->h),
    // lower right
    htrs_distance2line(point,  sqrt(3.), (l[2]  -detector->width)*2.*detector->h),
    // lower
    htrs_distance2line(point,        0., (l[0]  -detector->width)*   detector->h),
    // lower left
    htrs_distance2line(point, -sqrt(3.), (l[1]  -detector->width)*2.*detector->h),
    // upper left
    htrs_distance2line(point,  sqrt(3.), (l[2]+1-detector->width)*2.*detector->h),
  };


  // Find the closest distance to the nearest neighbouring pixel.
  int count, mindist, secpixel;
  double minimum;
  for(count=0; count<3; count++) {
    mindist = min_dist(distances, 6);
    minimum = distances[mindist];
    distances[mindist] = -1.;

    int k[3] = {l[0]+dl[0][mindist], l[1]+dl[1][mindist], l[2]+dl[2][mindist]};
    secpixel = htrs_get_lines2pixel(k, detector);
    
    if (secpixel != pixel) break;
  }

  if ((minimum > detector->ccsize) || (secpixel == pixel)) {
    // Single event!
    npixels=1;
    fraction[0] = 1.;

  } else {
    // Double event!
    npixels=2;

    pixel_coordinates = htrs_get_pixel2icoordinates(secpixel, detector);
    x[1] = pixel_coordinates.x;
    y[1] = pixel_coordinates.y;
      
    double mindistgauss = gaussint(distances[mindist]/detector->ccsigma);

    fraction[0] = 1. - mindistgauss;
    fraction[1] =      mindistgauss;
  }

  return(npixels);
}





///////////////////////////////////////////
inline int htrs_get_lines2pixel(int* l, Detector* detector)
{
  if ((l[0]<0)||(l[0]>=2*detector->width)||
      (l[1]<0)||(l[1]>=2*detector->width)||
      (l[2]<0)||(l[2]>=2*detector->width)) {
    return(INVALID_PIXEL);
  } else {
    return(detector->htrs_lines2pixel[l[0]][l[1]][l[2]]);
  }
}





///////////////////////////////////////////
inline struct Point2i htrs_get_pixel2icoordinates(int pixel, 
							 Detector* detector)
{
  if (pixel != INVALID_PIXEL) {
    return(detector->htrs_pixel2icoordinates[pixel]);
  } else {
    struct Point2i point2i = {-1, -1};
    return(point2i);
  }
}








///////////////////////////////////////////
int htrs_init_Detector(Detector* detector)
{
  //  Detector* detector=NULL;
  struct Point2d* centers = NULL;

  int status=EXIT_SUCCESS;
  char msg[MAXMSG];        // buffer for error output messages

  // Determine hexagonal pixel dimensions:
  detector->h = detector->pixelwidth/2.;      
  detector->a = 2. * detector->h / sqrt(3.);

  do { // Error handling loop 
    
    // Allocate memory and set the relation between the two different 
    // numbering arrays of the pixels in the hexagonal structure.
    // (linear numbering starting at the left bottom and 2D numbering)
    detector->htrs_pixel2icoordinates = 
      (struct Point2i*)malloc(HTRS_N_PIXELS * sizeof(struct Point2i));
    if (detector->htrs_pixel2icoordinates == NULL) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error: Not enough memory available for HTRS initialization!\n");
      HD_ERROR_THROW(msg, status);
      break;
    }

    detector->htrs_icoordinates2pixel = 
      (int**)malloc(detector->width * sizeof(int*));
    if (detector->htrs_icoordinates2pixel == NULL) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error: Not enough memory available for HTRS initialization!\n");
      HD_ERROR_THROW(msg, status);
      break;
    }

    int xi, yi, pixel=0;
    for(xi=0; xi<detector->width; xi++) {
      detector->htrs_icoordinates2pixel[xi] = 
	(int*)malloc(detector->width * sizeof(int));
      if (detector->htrs_icoordinates2pixel[xi] == NULL) {
	status = EXIT_FAILURE;
	sprintf(msg, "Error: Not enough memory available for HTRS "
		"initialization!\n");
	HD_ERROR_THROW(msg, status);
	break;
      }

      for(yi=0 ; yi<detector->width; yi++) {
	detector->htrs_icoordinates2pixel[xi][yi] = INVALID_PIXEL;
      }
    }

    
    // Fill the 2 different arrays for the conversion
    //     pixel index <-> coordinates.
    detector->htrs_icoordinates2pixel[0][1] = 22;
    detector->htrs_icoordinates2pixel[0][2] = 21;
    detector->htrs_icoordinates2pixel[0][3] = 20;
    detector->htrs_icoordinates2pixel[0][4] = 37;
    
    detector->htrs_icoordinates2pixel[1][1] = 23;
    detector->htrs_icoordinates2pixel[1][2] = 9;
    detector->htrs_icoordinates2pixel[1][3] = 8;
    detector->htrs_icoordinates2pixel[1][4] = 19;
    detector->htrs_icoordinates2pixel[1][5] = 36;

    detector->htrs_icoordinates2pixel[2][0] = 24;
    detector->htrs_icoordinates2pixel[2][1] = 10;
    detector->htrs_icoordinates2pixel[2][2] = 2;
    detector->htrs_icoordinates2pixel[2][3] = 7;
    detector->htrs_icoordinates2pixel[2][4] = 18;
    detector->htrs_icoordinates2pixel[2][5] = 35;

    detector->htrs_icoordinates2pixel[3][0] = 25;
    detector->htrs_icoordinates2pixel[3][1] = 11;
    detector->htrs_icoordinates2pixel[3][2] = 3;
    detector->htrs_icoordinates2pixel[3][3] = 1;
    detector->htrs_icoordinates2pixel[3][4] = 6;
    detector->htrs_icoordinates2pixel[3][5] = 17;
    detector->htrs_icoordinates2pixel[3][6] = 34;

    detector->htrs_icoordinates2pixel[4][0] = 26;
    detector->htrs_icoordinates2pixel[4][1] = 12;
    detector->htrs_icoordinates2pixel[4][2] = 4;
    detector->htrs_icoordinates2pixel[4][3] = 5;
    detector->htrs_icoordinates2pixel[4][4] = 16;
    detector->htrs_icoordinates2pixel[4][5] = 33;

    detector->htrs_icoordinates2pixel[5][1] = 27;
    detector->htrs_icoordinates2pixel[5][2] = 13;
    detector->htrs_icoordinates2pixel[5][3] = 14;
    detector->htrs_icoordinates2pixel[5][4] = 15;
    detector->htrs_icoordinates2pixel[5][5] = 32;

    detector->htrs_icoordinates2pixel[6][1] = 28;
    detector->htrs_icoordinates2pixel[6][2] = 29;
    detector->htrs_icoordinates2pixel[6][3] = 30;
    detector->htrs_icoordinates2pixel[6][4] = 31;

    for(xi=0 ; xi<detector->width; xi++) {
      for(yi=0 ; yi<detector->width; yi++) {
	if (detector->htrs_icoordinates2pixel[xi][yi] > 0) {
	  // In the program the pixel index starts at 0, not at 1 !!
	  // So shift all given pixel indices!
	  detector->htrs_icoordinates2pixel[xi][yi]--;
	  detector->htrs_pixel2icoordinates
	    [detector->htrs_icoordinates2pixel[xi][yi]].x = xi;   
	  detector->htrs_pixel2icoordinates
	    [detector->htrs_icoordinates2pixel[xi][yi]].y = yi;		 
	}
      }
    }

    // Now the 2 different numbering schemes can be easily converted 
    // among each other.
    


    // Calculate the centers of the hexagonal HTRS pixels.
    centers = (struct Point2d*)malloc(HTRS_N_PIXELS * sizeof(struct Point2d));
    if (centers == NULL) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error: Not enough memory available for HTRS initialization!\n");
      HD_ERROR_THROW(msg, status);
      break;
    }
    
    for(xi=0; xi<detector->width; xi++) {
      for(yi=0 ; yi<detector->width; yi++) {
	pixel = detector->htrs_icoordinates2pixel[xi][yi];
	if (pixel != INVALID_PIXEL) {
	  centers[pixel].x = (xi-3)                *1.5*detector->a;
	  centers[pixel].y = (yi-3+((xi+1)%2)*0.5) *2.0*detector->h;
	}
      }
    }


    int l[3];
    // Get memory and clear the array which is used to find the pixel
    // for given coordinates on the detector.
    detector->htrs_lines2pixel = 
      (int***) malloc(2*detector->width * sizeof(int**));
    for (l[0]=0; l[0]<2*detector->width; l[0]++) {
      detector->htrs_lines2pixel[l[0]] = 
	(int**) malloc(2*detector->width * sizeof(int*));
      for (l[1]=0; l[1]<2*detector->width; l[1]++) {
	detector->htrs_lines2pixel[l[0]][l[1]] = 
	  (int*) malloc(2*detector->width * sizeof(int));
	for (l[2]=0; l[2]<2*detector->width; l[2]++) {
	  detector->htrs_lines2pixel[l[0]][l[1]][l[2]] = INVALID_PIXEL;
	}
      }
    }

    int direction;
    struct Point2d point;
    const double sin30 = sin(M_PI/6.);
    const double cos30 = cos(M_PI/6.);

    for(pixel=0; pixel<HTRS_N_PIXELS; pixel++) {
      // For each pixel choose 6 points located around the center and
      // determine the line indices which define this pixel section.
    
      for(direction=0; direction<6; direction++) {
	point = centers[pixel];
	
	switch(direction) {
	case 0: // above center
	  point.y += detector->h/2;
	  break;
	case 1: // upper right section
	  point.x += detector->h/2*cos30;
	  point.y += detector->h/2*sin30;
	  break;
	case 2: // lower right section
	  point.x += detector->h/2*cos30;
	  point.y -= detector->h/2*sin30;
	  break;
	case 3: // below center
	  point.y -= detector->h/2;
	  break;
	case 4: // lower left section
	  point.x -= detector->h/2*cos30;
	  point.y -= detector->h/2*sin30;
	  break;
	case 5: // upper left section
	  point.x -= detector->h/2*cos30;
	  point.y += detector->h/2*sin30;
	  break;
	}

	htrs_get_lines(point, detector, l);
	detector->htrs_lines2pixel[l[0]][l[1]][l[2]] = pixel;
      }
    }

  } while (0);  // END of error handling loop


  // --- clean up ---

  free(centers);

  return(status);
}




/*
// None-Gaussian Charge clouds:
// Calculate the carge splitting according to a concept proposed by Konrad Dennerl.
int get_pixel_square(
		     Detector* detector, 
		     struct Point2d position, 
		     int* x, int* y, 
		     double* fraction
		     )
{
  int npixels = 0;  // number of affected pixels

  // The following array entries are used to transform between 
  // different array indices.
  int xe[4] = {1,0,-1,0};   
  int ye[4] = {0,1,0,-1};

  // -- Determine the affected pixels:

  // Calculate pixel indices (integer) of central affected pixel:
  x[0] = (int)(position.x/detector->pixelwidth + (double)(detector->width/2) +1.)-1;
  y[0] = (int)(position.y/detector->pixelwidth + (double)(detector->width/2) +1.)-1;
  
  // Calculate the distances from the event to the borders of the 
  // surrounding pixel (in [m]):
  double distances[4] = { 
    // distance to right pixel edge
    (x[0]-detector->offset+1)*detector->pixelwidth - position.x,
    // distance to upper edge
    (y[0]-detector->offset+1)*detector->pixelwidth - position.y,
    // distance to left pixel edge
    position.x - (x[0]-detector->offset)*detector->pixelwidth,
    // distance to lower edge
    position.y - (y[0]-detector->offset)*detector->pixelwidth
  };

  int mindist = min_dist(distances, 4);
  x[1] = x[0] + xe[mindist];
  y[1] = y[0] + ye[mindist];

  double minimum = distances[mindist];
  // search for the next to minimum distance to an edge
  distances[mindist] = -1.;
  int secmindist = min_dist(distances, 4);
  distances[mindist] = minimum;

  x[2] = x[0] + xe[secmindist];
  y[2] = y[0] + ye[secmindist];
  x[3] = x[1] + xe[secmindist];
  y[3] = y[1] + ye[secmindist];

  // --- Now we know the affected pixels and can determine the charge fractions.

  fraction[0] = exp(-pow((sqrt(pow(detector->pixelwidth/2-distances[mindist],2.)+
			       pow(detector->pixelwidth/2-distances[secmindist],2.)))
			 /(0.35*75.e-6),2.));
  fraction[1] = exp(-pow((sqrt(pow(detector->pixelwidth/2+distances[mindist],2.)+
			       pow(detector->pixelwidth/2-distances[secmindist],2.)))
			 /(0.35*75.e-6),2.));
  fraction[2] = exp(-pow((sqrt(pow(detector->pixelwidth/2-distances[mindist],2.)+
			       pow(detector->pixelwidth/2+distances[secmindist],2.)))
			 /(0.35*75.e-6),2.));
  fraction[3] = exp(-pow((sqrt(pow(detector->pixelwidth/2+distances[mindist],2.)+
			       pow(detector->pixelwidth/2+distances[secmindist],2.)))
			 /(0.35*75.e-6),2.));
  double sum = fraction[0]+fraction[1]+fraction[2]+fraction[3];
  fraction[0] /= sum;
  fraction[1] /= sum;
  fraction[2] /= sum;
  fraction[3] /= sum;


  npixels = 4;


  // --- Check whether all pixels lie inside the detector:
  int count;
  for(count=0; count<npixels; count++) {
    if ((x[count]<0) || (x[count]>=detector->width) ||
	(y[count]<0) || (y[count]>=detector->width)) {
      x[count] = INVALID_PIXEL;
      y[count] = INVALID_PIXEL;
    }
  }

  // Return the number of affected pixel, and 0, if the position lies outside
  // the detector array.
  return(npixels);
}
*/



///////////////////////////////////////////////////////
void htrs_free_Detector(Detector* detector)
{
  free(detector->htrs_pixel2icoordinates);

  int count;
  for(count=0; count<detector->width; count++) {
    free(detector->htrs_icoordinates2pixel[count]);
  }
  free(detector->htrs_icoordinates2pixel);

  int count2;
  for(count=0; count<detector->width; count++) {
    for(count2=0; count2<detector->width; count2++) {
      free(detector->htrs_lines2pixel[count][count2]);
    }
    free(detector->htrs_lines2pixel[count]);
  }
  free(detector->htrs_lines2pixel);

}





