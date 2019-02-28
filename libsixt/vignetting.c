/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
*/

#include "vignetting.h"


Vignetting* newVignetting(const char* const filename, int* const status) 
{
  Vignetting* vignetting=NULL;
  fitsfile* fptr=NULL;
  float* data_buffer=NULL;
  int count1, count2, count3; 

  do {

    // Allocate memory for the Vignetting data STRUCTURE:
    vignetting=(Vignetting*)malloc(sizeof(Vignetting));
    if (NULL==vignetting) {
      *status=EXIT_FAILURE;
      SIXT_ERROR("could not allocate memory for storing the vignetting data");
      break;
    }


    // Open the FITS file for reading the vignetting function.
    headas_chat(5, "open Vignetting FITS file '%s' ...\n", filename);
    if (fits_open_table(&fptr, filename, READONLY, status)) break;
    
    // Determine the column numbers of the individual columns.
    int column_energ_lo=0, column_energ_hi=0, column_theta=0;
    int column_energy=0;
    int column_phi=0, column_vignet=0;

    if(fits_get_colnum(fptr, CASEINSEN, "ENERGY", &column_energy, status)){
        if(fits_get_colnum(fptr, CASEINSEN, "ENERG_LO", &column_energ_lo, status)) break;
        if(fits_get_colnum(fptr, CASEINSEN, "ENERG_HI", &column_energ_hi, status)) break;
    }


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
    vignetting->energy  =(float*)malloc(vignetting->nenergies*sizeof(float));
    vignetting->theta   =(float*)malloc(vignetting->ntheta   *sizeof(float));
    vignetting->phi     =(float*)malloc(vignetting->nphi     *sizeof(float));
    if ((NULL==vignetting->energy)  || (NULL==vignetting->theta) || (NULL==vignetting->phi)) {
    	*status=EXIT_FAILURE;
    	SIXT_ERROR("could not allocate memory for storing the vignetting data");
    	break;
    } else {
    	for(count1=0; count1<vignetting->nenergies; count1++) {
    		vignetting->energy[count1]=0.;
    	}
    	for(count1=0; count1<vignetting->ntheta; count1++) {
    		vignetting->theta[count1]=0.;
    	}
    	for(count1=0; count1<vignetting->nphi; count1++) {
    		vignetting->phi[count1]=0.;
    	}
    }

    vignetting->vignet=(float***)malloc(vignetting->nenergies*sizeof(float**));
    if (NULL!=vignetting->vignet) {
    	for(count1=0; count1<vignetting->nenergies; count1++) {
    		vignetting->vignet[count1]=(float**)malloc(vignetting->ntheta*sizeof(float*));
    		if (NULL!=vignetting->vignet[count1]) {
    			for(count2=0; count2<vignetting->ntheta; count2++) {
    				vignetting->vignet[count1][count2]=
    						(float*)malloc(vignetting->nphi*sizeof(float));
    				if (NULL!=vignetting->vignet[count1][count2]) {
    					for(count3=0; count3<vignetting->nphi; count3++) {
    						vignetting->vignet[count1][count2][count3]=0.;
    					}
    				} else {
    					*status=EXIT_FAILURE;
    					SIXT_ERROR("could not allocate memory for storing the vignetting data");
    					break;
    				}
    			}
    		} else {
    			*status=EXIT_FAILURE;
    			SIXT_ERROR("could not allocate memory for storing the vignetting data");
    			break;
    		}
    	}
    } else {
    	*status=EXIT_FAILURE;
    	SIXT_ERROR("could not allocate memory for storing the vignetting data");
    	break;
    }
    CHECK_STATUS_BREAK(*status);

    data_buffer=(float*)malloc(vignetting->nenergies*
    		vignetting->ntheta*
			vignetting->nphi*sizeof(float));
    if (NULL==data_buffer) {
    	*status=EXIT_FAILURE;
    	SIXT_ERROR("could not allocate memory for storing the vignetting data");
    	break;
    } else {
    	for(count1=0;
    			count1<vignetting->nenergies*vignetting->ntheta*vignetting->nphi;
    			count1++) {
    		data_buffer[count1]=0.;
    	}
    }
    // Now all memory is allocated successfully!

    // READ the data from the FITS table.
    int anynul=0;
    if (column_energy){
        fits_read_col(fptr, TFLOAT, column_energ_lo, 1, 1, vignetting->nenergies,
        		vignetting->energy, vignetting->energy, &anynul, status);
    } else {
    	// if only ENERG_LO and ENERG_HI is given (old format), we need to calculate ENERGY
        float* energ_lo=(float*)malloc(vignetting->nenergies*sizeof(float));
        float* energ_hi=(float*)malloc(vignetting->nenergies*sizeof(float));

        if (energ_lo==NULL || energ_hi==NULL){
        	*status=EXIT_FAILURE;
        	SIXT_ERROR("could not allocate memory for storing the vignetting data");
        	break;
        }

        fits_read_col(fptr, TFLOAT, column_energ_lo, 1, 1, vignetting->nenergies,
        		energ_lo, energ_lo, &anynul, status);
        fits_read_col(fptr, TFLOAT, column_energ_hi, 1, 1, vignetting->nenergies,
        		energ_hi, energ_hi, &anynul, status);

        for (int ii=0; ii<vignetting->nenergies; ii++){
        	vignetting->energy[ii] = 0.5*(energ_lo[ii]+energ_hi[ii]);
        }
        free(energ_lo);
        free(energ_hi);
    }


    fits_read_col(fptr, TFLOAT, column_theta, 1, 1, vignetting->ntheta, 
    		vignetting->theta, vignetting->theta, &anynul, status);
    fits_read_col(fptr, TFLOAT, column_phi, 1, 1, vignetting->nphi, 
    		vignetting->phi, vignetting->phi, &anynul, status);
    fits_read_col(fptr, TFLOAT, column_vignet, 1, 1, 
    		vignetting->nenergies*vignetting->ntheta*vignetting->nphi,
			data_buffer, data_buffer, &anynul, status);

    // Determine the minimum and maximum available energy:
    vignetting->Emin=vignetting->energy[0];
    vignetting->Emax=vignetting->energy[vignetting->nenergies-1];

    // Scale from [deg] -> [rad]:
    for (count1=0; count1<vignetting->ntheta; count1++) {
    	vignetting->theta[count1]*=M_PI/180.;
    }
    for (count1=0; count1<vignetting->nphi; count1++) {
    	vignetting->phi[count1]*=M_PI/180.;
    }

    // Plot debug information about available energies, off-axis
    // angles, and azimuthal angles.
    headas_chat(5, "Vignetting - available energies:\n");
    for (count1=0; count1<vignetting->nenergies; count1++) {
      headas_chat(5, " %.3lf keV  \n",
          vignetting->energy[count1]);
    }
    headas_chat(5, "Vignetting - available off-axis angles:\n");
    for (count1=0; count1<vignetting->ntheta; count1++) {
      headas_chat(5, " %.3lf arc min\n", vignetting->theta[count1]/M_PI*180.*60.);
    }
    headas_chat(5, "Vignetting - available azimuthal angles:\n");
    for (count1=0; count1<vignetting->nphi; count1++) {
      headas_chat(5, " %.3lf deg\n", vignetting->phi[count1]/M_PI*180.);
    }

    // Transfer the data from the data buffer to the Vignetting data structure:
    if (EXIT_SUCCESS!=*status) break;
    for(count1=0; count1<vignetting->nenergies; count1++) {
      for(count2=0; count2<vignetting->ntheta; count2++) {
	for(count3=0; count3<vignetting->nphi; count3++) {
	  vignetting->vignet[count1][count2][count3] = 
	    data_buffer[count1+
			count2*vignetting->nenergies+
			count3*vignetting->nenergies*vignetting->ntheta];

	  // Output of the vignetting value for this particular parameters.
	  headas_chat(5, "Vignetting: %.2lf%% for "
		      "%.1lf keV, %.4lf arc min, %.4lf deg, \n",
		      vignetting->vignet[count1][count2][count3]*100., 
		      vignetting->energy[count1],
		      vignetting->theta[count2]/M_PI*180.*60.,
		      vignetting->phi[count3]/M_PI*180.);
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


void destroyVignetting(Vignetting** const vi) {
  if (NULL!=*vi) {
    if (NULL!=(*vi)->energy) free((*vi)->energy);
    if (NULL!=(*vi)->theta)    free((*vi)->theta);
    if (NULL!=(*vi)->phi)      free((*vi)->phi);
    
    if (NULL!=(*vi)->vignet) {
      int count1, count2;
      for (count1=0; count1<(*vi)->nenergies; count1++) {
	if (NULL!=(*vi)->vignet[count1]) {
	  for (count2=0; count2<(*vi)->ntheta; count2++) {
	    if (NULL!=(*vi)->vignet[count1][count2]) {
	      free((*vi)->vignet[count1][count2]);
	    }
	  }
	  free((*vi)->vignet[count1]);
	}
      }
      free((*vi)->vignet);
    }
    free((*vi));
    *vi=NULL;
  }
}


// At the moment SIXTE can only handle the cases without phi dependece phi = 0.
static float interpol_vign_theta(const float theta, const float* arr_theta, float** arr_vign, const int ntheta){


	// Check if the required angle is larger than the biggest in the
	// vignetting data.
	if (theta>=arr_theta[ntheta-1]) {
		return(arr_vign[ntheta-1][0]);
	}

	// Find the two values in the vignetting data surrounding
	// the required angle.
	int jj;
	for(jj=1; jj<ntheta; jj++) {
		if (arr_theta[jj]>=theta) {
			// Interpolate between both values.
			return(arr_vign[jj-1][0]+
					(arr_vign[jj][0]-arr_vign[jj-1][0])*
					(theta-arr_theta[jj-1])/(arr_theta[jj]-arr_theta[jj-1]));
		}
	}
	return 1.;
}

float get_Vignetting_Factor(const Vignetting* const vi, const float energy, 
			    const float theta, const float phi) 
{
  // Check if any vignetting is specified. 
  // If not, return a default value of 1.
  if (NULL==vi) return(1.);

  (void)phi;
  
  /*
  // At the moment this routine can only handle the case with phi = 0.
  if (phi!=0.) {
    SIXT_ERROR("vignetting can only be determined for phi=0");
    return(0.);
  }
  */

  /* if we are above or below the defined energies, we return the vignetting
   * value at the edge of the grid  */
  if (energy<=vi->Emin){
	  return interpol_vign_theta(theta, vi->theta, vi->vignet[0], vi->ntheta);
  }
  if (energy>=vi->Emax){
	  return interpol_vign_theta(theta, vi->theta, vi->vignet[vi->nenergies-1], vi->ntheta);
  }

  // Find the right energy bin.
  int ind = binary_search_float(energy,vi->energy,vi->nenergies );

  double vign_val[2];
  // assume linear interpolation
  vign_val[0] = interpol_vign_theta(theta, vi->theta, vi->vignet[ind], vi->ntheta);
  vign_val[1] = interpol_vign_theta(theta, vi->theta, vi->vignet[ind+1], vi->ntheta);


  double ifac = (energy-vi->energy[ind]) / ( vi->energy[ind+1] - vi->energy[ind] );

  return interp_lin_1d(ifac, vign_val[0], vign_val[1]);

}

