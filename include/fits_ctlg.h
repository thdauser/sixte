#ifndef FITS_CTLG_H
#define FITS_CTLG_H

#include "fitsio.h"
#include "vector.h"

// number of fields in the rnd_source table (= number of columns)
#define N_SOURCE_FIELDS 6    
// number of fields in the orbit table (t, \vec{r}, \vec{v})
#define N_ORBIT_FIELDS 7     



/////////////////////
// SOURCES:

// define a structure which contains the necessary data for each source
struct source {
  Vector r;  // position of the source
  double flux;      // energy flux
};

// creates the necessary parameters to generate the table in the source FITS file
void create_srctbl_parameters(char *ftype[N_SOURCE_FIELDS], 
			      char *fform[N_SOURCE_FIELDS], 
			      char *funit[N_SOURCE_FIELDS]);

// writes a row of data into the source FITS file
void add_srctbl_row(fitsfile *fptr, long row, const double rasc, 
		    const double dec, const float countrate, int *status);

// reads a row of data from the source FITS file
int get_srctbl_row(fitsfile *fptr, long row, const int columns[], 
		   double *rasc, double *dec, float *countrate,int *status);



///////////////////
// ORBITS:


// define a structure to store the orbit data (time, position, velocity)
struct orbit_entry {
  double time;      // time
  Vector r;  // position
  Vector v;  // velocity
};


// writes a row of orbit data into the orbit FITS file
void add_orbtbl_row(fitsfile *fptr, long row, double time, 
		    Vector r, Vector v, int *status);

// creates the necessary parameters to generate the table in the FITS file
void create_orbtbl_parameter(char *ftype[N_ORBIT_FIELDS], 
			     char *fform[N_ORBIT_FIELDS], 
			     char *funit[N_ORBIT_FIELDS]);

// reads a row of data from the orbit FITS file
int get_orbtbl_row(fitsfile *fptr, long row, double *time, 
		   Vector *r, Vector *v, int *status);



#endif

