/** 
 * This file is the header file of 'orbatt.c'
 * contains all definitions for orbit and attitude handling.
 */

#ifndef ORBATT_H
#define ORBATT_H 1

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fitsio.h"
#include "headas.h"
#include "headas_error.h"

#include "fits_attitude.h"
#include "fits_ctlg.h"
#include "vector.h"
#include "telescope.h"


// maximum number of entries in the orbit-/attitude-catalog (satellite_catalog)
#define MAX_NORBIT_ENTRIES 500000   
// time intervals [s] to update the preselected source catalog
#define ORBIT_UPDATE_TIME 5400.     



/** 
 * Structure that contains the direction of the telescope and information about the motion of the satellite:
struct telescope_data {
  double time;
  struct vector r;       // position vector of the satellite
  struct vector v;       // velocity vector of the satellite
  struct vector nz;      // direction of the telescope axis (in simplest case same direction as r)
  struct vector nx, ny;  // nx usually points in the direction of motion, but is perpendicular to nz,
                         // ny is perpendicular to nx and nz
};
*/



// This function allocates memory for the satellite catalog (containing orbit and attitude information), 
// reads this data from the corresponding FITS files and stores it in the catalog
int get_satellite_catalog(struct Telescope **, long *nentries, double t0, 
			  double timespan, const char orbit_filename[], 
			  const char attitude_filename[]);



#endif

