/** 
 * This fiel belongs to 'measurement.c' and
 * contains all source code for orbit and attitude handling.
 */

#include "sixt.h"
#include "orbatt.h"


// This function allocates memory for the satellite catalog 
// (containing orbit and attitude information), reads this data 
// from the corresponding FITS files and stores it in the catalog.
int get_satellite_catalog(
			  struct Telescope **sat_catalog,  // pointer to the array of orbit-/attitude-entries
			  long *nentries,      // number of entries in the satellite catalog (return value)
			  double t0,           // requested start time for the catalog array
			  double timespan,     // requested end time for the catalog array
			  const char orbit_filename[], 
			  const char attitude_filename[]
			  )
{
  fitsfile *orbit_fptr;  // FITS file pointer to orbit file
  fitsfile *att_fptr;    // -"-               to attitude file
  int hdunum, hdutype;   // HDU number and type for FITS access
  long nrows;            // number of lines in the FITS file
  long entry=0;          // At the end this variable contains the number of entries
                         // in the satellite catalog.
    
  char msg[MAXMSG];
  int status=0;                 // error handling variable

  do {    // beginning of error handling loop

    // get memory for the satellite catalog:
    *sat_catalog = (struct Telescope *) 
      malloc(MAX_NORBIT_ENTRIES * sizeof(struct Telescope));
    if (*sat_catalog==NULL) {
      status = EXIT_FAILURE;
      sprintf(msg, "Not enough memory available to store the satellite catalog!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }


    // Read-in the orbit and attitude data from the FITS files and store 
    // it in satellite catalog array.
    // First open the orbit file:
    if (fits_open_table(&orbit_fptr, orbit_filename, READONLY, &status)) break;

    // After opening the orbit file, get the number of the current HDU.
    if (fits_get_hdu_num(orbit_fptr, &hdunum) == 1) {
      // This is the primary array, so try to move to the first extension 
      // and see if it is a table
      if (fits_movabs_hdu(orbit_fptr, 2, &hdutype, &status)) break;
    } else {
      // get the HDU type
      if (fits_get_hdu_type(orbit_fptr, &hdutype, &status)) break;
    }

    // image HDU results in an error message
    if (hdutype==IMAGE_HDU) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: FITS extension in orbit file '%s' is not a table "
	      "but an image (HDU number: %d)\n", orbit_filename, hdunum);
      HD_ERROR_THROW(msg,status);
      break;
    }

    // get the number of rows in the orbit file (number of given points of time)
    fits_get_num_rows(orbit_fptr, &nrows, &status);


    
    long counter;
    double time;                // buffer for the time
    struct vector sat_r, sat_v; // buffer for position and velocity of the satellite
    
    // Subsequently read lines from orbit file:
    for(counter=0, entry=0; (counter<nrows)&&(!status) ; counter++) {
      if (get_orbtbl_row(orbit_fptr, counter, &time, &sat_r, &sat_v, &status)) 
	break;

      if (time >= t0) {
	// check, if we are already at the end of the orbit file:
	if (counter>=nrows) {
	  status=EXIT_FAILURE;
	  sprintf(msg, "Error: Not enough orbit data available for the "
		  "specified period!");
	  HD_ERROR_THROW(msg,status);
	  break;
	}

	// store the satellite orbit data in the array
	(*sat_catalog)[entry].time = time;
	(*sat_catalog)[entry].r = sat_r;
	(*sat_catalog)[entry++].v = sat_v;
      }
	  
      // end of reading process (end of specified period)
      if (time>t0+timespan) {
	break;
      }
    }  // end of the orbit readout loop
   

    // Then open attitude file:
    if (fits_open_table(&att_fptr, attitude_filename, READONLY, &status)) break;

    // After opening the attitude file, get the number of the current HDU
    if (fits_get_hdu_num(att_fptr, &hdunum) == 1) {
      // This is the primary array, so try to move to the first extension 
      // and see if it is a table
      if (fits_movabs_hdu(att_fptr, 2, &hdutype, &status)) break;
    } else {
      // get the HDU type
      if (fits_get_hdu_type(att_fptr, &hdutype, &status)) break;
    }

    // image HDU results in an error message
    if (hdutype==IMAGE_HDU) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: FITS extension in attitude file '%s' is not a table "
	      "but an image (HDU number: %d)\n", attitude_filename, hdunum);
      HD_ERROR_THROW(msg,status);
      break;
    }

    // get the number of rows in the orbit file (number of given points of time)
    fits_get_num_rows(att_fptr, &nrows, &status);


    char valtime[19];
    double view_ra, view_dec, rollangle;  // angles specifying the direction of the telescope (attitude data)  
    struct vector nx, ny;                 // unit vectors spanning the satellite's intrinsic coordinate system
    // read lines from attitude file subsequently
    for (counter=0, entry=0; (counter<nrows)&&(!status) ; counter++) {
      if (get_atttbl_row(att_fptr, counter, valtime, &time, &view_ra, &view_dec, &rollangle, &status)) break;

      if (time >= t0) {
	// check, if we are already at the end of the attitude file:
	if (counter>=nrows) {
	  status=EXIT_FAILURE;
	  sprintf(msg, "Error: Not enough attitude data available for the "
		  "specified period!");
	  HD_ERROR_THROW(msg,status);
	  break;
	}

	// store the satellite orbit data in the array
	if ((*sat_catalog)[entry].time != time) {
	  status=EXIT_FAILURE;
	  sprintf(msg, "Error: time values not equivalent for orbit and "
		  "attitude data: %lf vs. %lf!", (*sat_catalog)[entry].time, time);
	  HD_ERROR_THROW(msg,status);
	  break;
	}
	  
	// calculate and store attitude data:
	// rescale from degrees to radians:
	view_ra = view_ra*M_PI/180.;
	view_dec = view_dec*M_PI/180.;
	rollangle = rollangle*M_PI/180.;

	// telescope direction nz:
	(*sat_catalog)[entry].nz = unit_vector(view_ra, view_dec);

	// nx:
	nx = normalize_vector((*sat_catalog)[entry].v);
	double scp = scalar_product(nx, (*sat_catalog)[entry].nz);
	nx.x -= scp*(*sat_catalog)[entry].nz.x;
	nx.y -= scp*(*sat_catalog)[entry].nz.y;
	nx.z -= scp*(*sat_catalog)[entry].nz.z;	
	nx = normalize_vector(nx);

	// ny:
	ny = normalize_vector(vector_product((*sat_catalog)[entry].nz, nx));

	(*sat_catalog)[entry].nx.x = cos(rollangle)*nx.x + sin(rollangle)*ny.x;
	(*sat_catalog)[entry].nx.y = cos(rollangle)*nx.y + sin(rollangle)*ny.y;
	(*sat_catalog)[entry].nx.z = cos(rollangle)*nx.z + sin(rollangle)*ny.z;
	(*sat_catalog)[entry].ny.x = -sin(rollangle)*nx.x + cos(rollangle)*ny.x;
	(*sat_catalog)[entry].ny.y = -sin(rollangle)*nx.y + cos(rollangle)*ny.y;
	(*sat_catalog)[entry].ny.z = -sin(rollangle)*nx.z + cos(rollangle)*ny.z;

	entry++;
      }
	  
      // end of reading process (end of specified period)
      if (time>t0+timespan) {
	break;
      }
    }  // End of the orbit readout loop

  } while (0); // End of error handling loop



  // --- clean up ---

  // close FITS files
  if (orbit_fptr) fits_close_file(orbit_fptr, &status);
  if (att_fptr) fits_close_file(att_fptr, &status);

  // return number of stored data lines and error status variable
  *nentries = entry;
  return(status);
}


