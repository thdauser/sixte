#include <malloc.h>
#include "fitsio.h"

/** Data structure containing the mirror vignetting function. */
typedef struct {
  int nenergies; /**< Number of energy bins. */
  int ntheta; /**< Number of off-axis angles. */
  int nphi; /**< Number of azimuth angles. */

  float* energ_lo; /**< Minimum energy of bin in [keV]. */
  float* energ_hi; /**< Maximum energy of bin in [keV]. */
  float* theta; /**< Off-axis angle in [rad]. */
  float* phi; /**< Azimuth angle in [rad]. */
  float*** vignet; /**< Vignetting data. Array[energy, theta, phi] */
} Vignetting;



/////////////////////////////////////////////////////////////
/** Constructor of the Vignetting data structure. 
 * Loads the vignetting function from a given FITS file. 
 * The format of the FITS file is defined by 
 * OGIP Memo CAL/GEN/92-021. */
Vignetting* get_Vignetting(char* filename, int* status) {
  Vignetting* vignetting=NULL;
  fitsfile* fptr=NULL;
  float* data_buffer=NULL;

  char msg[MAXMSG]; 


  do {

    // Allocate memory for the Vignetting data STRUCTURE:
    vignetting = (Vignetting*)malloc(sizeof(Vignetting));
    if (NULL==vignetting) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: could not allocate memory for storing "
	      "the vignetting data!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }


    // Open the FITS file for reading the vignetting function.
    if (fits_open_table(&fptr, filename, READONLY, status)) break;
    
    // Determine the column numbers of the individual columns.
    int column_energ_lo=0, column_energ_hi=0, column_theta=0;
    int column_phi=0, column_vignet=0;
    if(fits_get_colnum(fptr, CASEINSEN, "ENERG_LO", &column_energ_lo, status)) break;
    if(fits_get_colnum(fptr, CASEINSEN, "ENERG_HI", &column_energ_hi, status)) break;
    if(fits_get_colnum(fptr, CASEINSEN, "THETA", &column_theta, status)) break;
    if(fits_get_colnum(fptr, CASEINSEN, "PHI", &column_phi, status)) break;
    if(fits_get_colnum(fptr, CASEINSEN, "VIGNET", &column_vignet, status)) break;
    // Determine the format of the multi-dimensional array from the corresponding
    // header keyword.
    char tdim_name[MAXMSG], tdim_value[MAXMSG], comment[MAXMSG];
    sprintf(tdim_name, "TDIM%d", column_vignet);
    if (fits_read_key(fptr, TSTRING, tdim_name, tdim_value, comment, status)) break;
    sscanf(tdim_value, "(%d,%d,%d)", &vignetting->nenergies, 
	   &vignetting->ntheta, &vignetting->nphi);

    // Allocate memory for the Vignetting data:
    vignetting->energ_lo = (float*)malloc(vignetting->nenergies*sizeof(float));
    vignetting->energ_hi = (float*)malloc(vignetting->nenergies*sizeof(float));
    vignetting->theta    = (float*)malloc(vignetting->ntheta   *sizeof(float));
    vignetting->phi      = (float*)malloc(vignetting->nphi     *sizeof(float));
    if ((NULL==vignetting->energ_lo) || (NULL==vignetting->energ_hi) ||
	(NULL==vignetting->theta) || (NULL==vignetting->phi)) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: could not allocate memory for storing "
	      "the vignetting data!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    } else {
      int count;
      for(count=0; count<vignetting->nenergies; count++) {
	vignetting->energ_lo[count] = 0.;
	vignetting->energ_hi[count] = 0.;
      }
      for(count=0; count<vignetting->ntheta; count++) {
	vignetting->theta[count] = 0.;
      }
      for(count=0; count<vignetting->nphi; count++) {
	vignetting->phi[count] = 0.;
      }
    }

    vignetting->vignet = (float***)malloc(vignetting->nenergies*sizeof(float**));
    if (NULL!=vignetting->vignet) {
      int count;
      for(count=0; count<vignetting->nenergies; count++) {
	vignetting->vignet[count] = (float**)malloc(vignetting->ntheta*sizeof(float*));
	if (NULL!=vignetting->vignet[count]) {
	  int count2;
	  for(count2=0; count2<vignetting->ntheta; count2++) {
	    vignetting->vignet[count][count2] = 
	      (float*)malloc(vignetting->nphi*sizeof(float));
	    if (NULL!=vignetting->vignet[count][count2]) {
	      int count3;
	      for(count3=0; count3<vignetting->nphi; count3++) {
		vignetting->vignet[count][count2][count3] = 0.;
	      }
	    } else {
	      *status=EXIT_FAILURE;
	      sprintf(msg, "Error: could not allocate memory for storing "
		      "the vignetting data!\n");
	      HD_ERROR_THROW(msg, *status);
	      break;
	    }
	  }
	} else {
	  *status=EXIT_FAILURE;
	  sprintf(msg, "Error: could not allocate memory for storing "
		  "the vignetting data!\n");
	  HD_ERROR_THROW(msg, *status);
	  break;
	}
      }
    } else {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: could not allocate memory for storing "
	      "the vignetting data!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }
    if (EXIT_SUCCESS!=*status) break;

    data_buffer = (float*)malloc(vignetting->nenergies*
				 vignetting->ntheta*
				 vignetting->nphi*sizeof(float));
    if (NULL==data_buffer) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: could not allocate memory for storing "
	      "the vignetting data!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    } else {
      int count;
      for(count=0; 
	  count<vignetting->nenergies*vignetting->ntheta*vignetting->nphi; count++) {
	data_buffer[count]=0.;
      }
    }
    // Now all memory is allocated successfully!

    
    // READ the data from the FITS table.
    int anynul=0;
    fits_read_col(fptr, TFLOAT, column_energ_lo, 1, 1, vignetting->nenergies, 
		  vignetting->energ_lo, vignetting->energ_lo, &anynul, status);
    fits_read_col(fptr, TFLOAT, column_energ_hi, 1, 1, vignetting->nenergies, 
		  vignetting->energ_hi, vignetting->energ_hi, &anynul, status);
    fits_read_col(fptr, TFLOAT, column_theta, 1, 1, vignetting->ntheta, 
		  vignetting->theta, vignetting->theta, &anynul, status);
    fits_read_col(fptr, TFLOAT, column_phi, 1, 1, vignetting->nphi, 
		  vignetting->phi, vignetting->phi, &anynul, status);
    fits_read_col(fptr, TFLOAT, column_vignet, 1, 1, 
		  vignetting->nenergies*vignetting->ntheta*vignetting->nphi, 
		  data_buffer, data_buffer, &anynul, status);

    // Transfer the data from the data buffer to the Vignetting data structure:
    if (EXIT_SUCCESS!=*status) break;
    int count1, count2, count3; 
    for(count1=0; count1<vignetting->nenergies; count1++) {
      for(count2=0; count2<vignetting->ntheta; count2++) {
	for(count3=0; count3<vignetting->nphi; count3++) {
	  vignetting->vignet[count1][count2][count3] = 
	    data_buffer[count1+
			count2*vignetting->nenergies+
			count3*vignetting->nenergies*vignetting->ntheta];
	}
      }
    }

  } while(0); // END of Error handling loop

  // --- Clean up ---

  if (NULL!=data_buffer) free(data_buffer);
  
  if (NULL!=fptr) fits_close_file(fptr, status);

  if (EXIT_SUCCESS!=*status) vignetting=NULL;
  return(vignetting);
}


