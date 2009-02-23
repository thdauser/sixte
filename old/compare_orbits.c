#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "vector.h"
#include "fitsio.h"
#include "fits_ctlg.h"

//////////////////////////////////////////////////////////////////////////////////////////
// This programs reads two different orbit files (FITS files) and compares the positions,
// i.e. calculates the distances of the satellites at each point of time.
//
// Usage:
// compare_orbits.c <scale1> <orbitfile1> <scale2> <orbitfile2a> <orbitfile2b> <orbitfile2c> ...
// scale: 1.0 if position in FITS file has unit [km], 0.001 if [m]
//
// Output:
// time, distance between the individual positions [km]
////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
  //const long double mu = 398601.3;  // G*Msolar; units: (km^3/s^2)

  int filecounter;
  fitsfile *fptr;
  int fitsstatus=0;
  int hdunum, hdutype;
  long nfrows, narows, fitscounter, arraycounter;
  struct orbit_entry oposition;
  struct orbit_entry *orbit_catalog;    // orbit catalog
  double scale1, scale2;
  


  // check if there are enough parameters
  if (argc < 5) {
    printf("You must specify an orbit file!\n");
    return EXIT_FAILURE;
  }

  // read the scaling factors for the data in both FITS files
  scale1 = atof(argv[1]);
  scale2 = atof(argv[3]);
  
  // read the reference data from the first FITS file
  // open 1st FITS file
  if (fits_open_file(&fptr, argv[2], READONLY, &fitsstatus)) {
    printf("Error: could not open orbit file '%s'!\n", argv[2]);
  } else {
    // after opening the orbit file, get the number of the current HDU
    if (fits_get_hdu_num(fptr, &hdunum) == 1) {
      // this is the primary array
      // try to move to the first extension and see if it is a table
      fits_movabs_hdu(fptr, 2, &hdutype, &fitsstatus);
    } else {
      // get the HDU type
      fits_get_hdu_type(fptr, &hdutype, &fitsstatus);
    }

    // if we have an image HDU, we cannot read any data, so print an error message
    if (hdutype == IMAGE_HDU) {
      printf("Error: extension in source file '%s' is not a table (default: second extension)!\n", argv[filecounter]);
    } else {
      // get the number of rows in the FITS table
      fits_get_num_rows(fptr, &narows, &fitsstatus);
	
      // get memory for the orbit catalog:
      orbit_catalog = (struct orbit_entry *) malloc(narows * sizeof(struct orbit_entry));
      if (!orbit_catalog) {
	printf("Error: Not enough memory available to store the orbit catalog!\n");
	return EXIT_FAILURE;
      }

      // perform a loop over the timesteps in the FITS file
      // (in the RXTE orbit files, there is an entry for each 60s)
      for(arraycounter = 0; arraycounter < narows; arraycounter++) {

	// read position and velocity from the FITS file:
	get_orbtbl_row(fptr, arraycounter, &oposition.time, &oposition.r, &oposition.v, &fitsstatus);
	// rescale from to [km]:
	oposition.r.x = oposition.r.x*scale1;
	oposition.r.y = oposition.r.y*scale1;
	oposition.r.z = oposition.r.z*scale1;    
	oposition.v.x = oposition.v.x*scale1;
	oposition.v.y = oposition.v.y*scale1;
	oposition.v.z = oposition.v.z*scale1;

	orbit_catalog[arraycounter] = oposition;
      }
      
      // close the 1st FITS file
      fits_close_file(fptr, &fitsstatus);
    }
  }



  
  // open the 2nd FITS file to calculate the position differences
  for(arraycounter=0,filecounter=4; filecounter < argc; filecounter++) {
    // open FITS file to read the orbit positions
    if (fits_open_file(&fptr, argv[filecounter], READONLY, &fitsstatus)) {
      printf("Error: could not open orbit file '%s'!\n", argv[filecounter]);
    } else {
      // after opening the orbit file, get the number of the current HDU
      if (fits_get_hdu_num(fptr, &hdunum) == 1) {
	// this is the primary array
	// try to move to the first extension and see if it is a table
	fits_movabs_hdu(fptr, 2, &hdutype, &fitsstatus);
      } else {
	// get the HDU type
	fits_get_hdu_type(fptr, &hdutype, &fitsstatus);
      }

      // if we have an image HDU, we cannot read any data, so print an error message
      if (hdutype == IMAGE_HDU) {
	printf("Error: extension in source file '%s' is not a table (default: second extension)!\n", argv[filecounter]);
      } else {
	// get the number of rows in the FITS table
	fits_get_num_rows(fptr, &nfrows, &fitsstatus);
	

	// perform a loop over the timesteps in the FITS file
	// (in the RXTE orbit files, there is an entry for each 60s)
	for(fitscounter = 0; fitscounter < nfrows; fitscounter++) {

	  // read position and velocity from the FITS file:
	  get_orbtbl_row(fptr, fitscounter, &oposition.time, &oposition.r, &oposition.v, &fitsstatus);
	  // rescale from to [km]:
	  oposition.r.x = oposition.r.x*scale2;
	  oposition.r.y = oposition.r.y*scale2;
	  oposition.r.z = oposition.r.z*scale2;    
	  oposition.v.x = oposition.v.x*scale2;
	  oposition.v.y = oposition.v.y*scale2;
	  oposition.v.z = oposition.v.z*scale2;    

	  // find the same point of time in the reference orbit catalog
	  while ( (fabs(orbit_catalog[arraycounter].time-oposition.time)>0.1) && (orbit_catalog[arraycounter].time < oposition.time) && (arraycounter < narows-1)) {
	    arraycounter++;
	  }

	  while ( (fabs(orbit_catalog[arraycounter].time-oposition.time)>0.1) && (orbit_catalog[arraycounter].time > oposition.time) && (arraycounter > 0)) {
	    arraycounter--;
	  }

	  if (fabs(orbit_catalog[arraycounter].time-oposition.time)>0.1) {
	    printf("#Error: no reference time available! (file %s)\n", argv[filecounter]);
	    return EXIT_FAILURE;
	  }

          // calculate and print  the distance between the position in the actual FITS file and in the orbit catalog:
          printf("%lf\t%lf\n", oposition.time, sqrt(pow(orbit_catalog[arraycounter].r.x-oposition.r.x,2.)+pow(orbit_catalog[arraycounter].r.y-oposition.r.y,2.)+pow(orbit_catalog[arraycounter].r.z-oposition.r.z,2.)));
	}

	// close the FITS file
	fits_close_file(fptr, &fitsstatus);
      }
    }
  } // end of loop over all file arguments

  // print any error that have occurred on FITS access
  if (fitsstatus) {
    fits_report_error(stderr, fitsstatus);
  }
  return(fitsstatus);

}

