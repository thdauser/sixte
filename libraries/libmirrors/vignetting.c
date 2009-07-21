#include "vignetting.h"


/////////////////////////////////////////////////////////////
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

    // Determine the minimum and maximum available energy:
    int count; 
    vignetting->Emin=-1.;
    vignetting->Emax=-1.;
    for (count=0; count<vignetting->nenergies; count++) {
      if ((vignetting->energ_lo[count] < vignetting->Emin) || (vignetting->Emin<0.)) {
	vignetting->Emin = vignetting->energ_lo[count];
      }
      if (vignetting->energ_hi[count] > vignetting->Emax) {
	vignetting->Emax = vignetting->energ_hi[count];
      }
    }						
    // Scale from [deg] -> [rad]:
    for (count=0; count<vignetting->ntheta; count++) {
      vignetting->theta[count] *= M_PI/180.;
    }
    for (count=0; count<vignetting->nphi; count++) {
      vignetting->phi[count] *= M_PI/180.;
    }

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



//////////////////////////////////////////////
void free_Vignetting(Vignetting* vi) {
  if (NULL!=vi) {
    if (NULL!=vi->energ_lo) free(vi->energ_lo);
    if (NULL!=vi->energ_hi) free(vi->energ_hi);
    if (NULL!=vi->theta)    free(vi->theta);
    if (NULL!=vi->phi)      free(vi->phi);
    
    if (NULL!=vi->vignet) {
      int count1, count2;
      for (count1=0; count1<vi->nenergies; count1++) {
	if (NULL!=vi->vignet[count1]) {
	  for (count2=0; count2<vi->ntheta; count2++) {
	    if (NULL!=vi->vignet[count1][count2]) {
	      free(vi->vignet[count1][count2]);
	    }
	  }
	  free(vi->vignet[count1]);
	}
      }
      free(vi->vignet);
    }
  }
  free(vi);
}




///////////////////////////////////////////////////////////////////////////////
float get_Vignetting_Factor(Vignetting* vi, float energy, float theta, float phi) {
  float factor=0.;

  // At the moment this routine can only handle the case with phi = 0.
  if (phi!=0.) {
    HD_ERROR_THROW("Error: vignetting can only be determined for phi=0!\n", EXIT_FAILURE);
    return(0.);
  }

  if ((energy<vi->Emin) || (energy>vi->Emax)) {
    factor = -1.; // Energy is out of range!
  } else {
    // Find the right energy bin:
    int count;
    for(count=0; count<vi->nenergies; count++) {
      if ((energy>vi->energ_lo[count]) && (energy<=vi->energ_hi[count])) {
	// Find the best fitting theta:
	int count2;
	float dtheta_min=-1.;
	for(count2=0; count2<vi->ntheta; count2++) {
	  if ((fabs(theta - vi->theta[count2]) < dtheta_min) || (dtheta_min<0.)) {
	    dtheta_min = fabs(theta - vi->theta[count2]);
	    factor = vi->vignet[count][count2][0];
	  }
	} // Loop to find the best theta.
	break;
      }
    } // Loop to find the right energy bin.
    assert (count<vi->nenergies);
  }

  return(factor);
}

