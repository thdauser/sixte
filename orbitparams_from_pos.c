#include <stdio.h>
#include <math.h>
#include "vector.h"
#include "fitsio.h"
#include "fits_ctlg.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// This programs reads the positions and velocities of a satellite from an orbit file (FITS file) and  //
// calculates the Keplerian orbital parameters for the corresponding orbit at each point of time.      //
//                                                                                                     //
// Usage:                                                                                              //
// orbitparams_from_pos <scale> <orbitfile>                                                            //
// scale: "1.0" if orbit positions have unit [km], "0.001" if unit is [m]
// orbitfile: FITS file containing the orbit to be analysed
//
// Output:
// time, position, Keplerian orbital elements
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char *argv[]) 
{
  const long double mu = 398601.3;  // G*Msolar; units: (km^3/s^2)

  double a,p,e,i,Omega,omega,M;
  struct vector ex, ey, ez, er, hv, qv, n, b;
  double h, q;
  double theta,E;

  int filecounter;
  double time;
  struct vector r, v;
  double scale;
  fitsfile *fptr;
  int fitsstatus = 0;
  int hdunum, hdutype;
  long nrows, counter;


  // check if there is a FITS file specified
  if (argc < 3) {
    printf("You must specify an orbit file!\n");
    return(-1);
  }

  scale = atof(argv[1]);
  
  for(filecounter=2; filecounter < argc; filecounter++) {
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
	printf("%s\n", argv[filecounter]);
	
	// get the number of rows in the FITS table
	fits_get_num_rows(fptr, &nrows, &fitsstatus);
	

	// perform a loop over the timesteps in the FITS file
	// (in the RXTE orbit files, there is an entry for each 60s)
	for(counter = 0; counter < nrows; counter++) {

	  // read position and velocity from the FITS file:
	  get_orbtbl_row(fptr, counter, &time, &r, &v, &fitsstatus);
	  // rescale from [m] to [km]:
	  r.x = r.x*scale;
	  r.y = r.y*scale;
	  r.z = r.z*scale;    
	  v.x = v.x*scale;
	  v.y = v.y*scale;
	  v.z = v.z*scale;    
      
	  /*	  if((filecounter == 2) && (counter == 0)) {
	    time_offset = time;
	  }
	  time -= time_offset;*/
    

	  // calculate orbital elements from this data (according to Steiner p. 86)
	  // step 1 + 3:
	  er = normalize_vector(r);
      
	  hv = vector_product(r, v);
	  ez = normalize_vector(hv);

	  b = vector_product(v, hv);
	  qv.x = b.x - mu*er.x;
	  qv.y = b.y - mu*er.y;
	  qv.z = b.z - mu*er.z;
	  ex = normalize_vector(qv);

	  ey = vector_product(ez, ex);
	  
	  h = sqrt(scalar_product(hv,hv));
	  q = sqrt(scalar_product(qv,qv));

	  n.x = -ez.y;
	  n.y = ez.x;
	  n.z = 0.;
	  n = normalize_vector(n);

	  // step 2:
	  p = pow(h,2.)/mu;
	  e = q/mu;
	  a = p/(1.-pow(e,2.));

	  // step 4:
	  i = acos(ez.z);
	  
	  // step 5:
	  if ((Omega = atan2(n.y, n.x))<0.0) {
	    Omega += 2.*M_PI;
	  }
      
	  // step 6:
	  if((omega = atan2(-scalar_product(n,ey),scalar_product(n,ex)))<0.0) {
	    omega += 2.*M_PI;
	  }

	  // step 7:
	  theta = atan2(scalar_product(er,ey), scalar_product(er,ex));
	  E = 2.*atan(sqrt((1-e)/(1+e))*tan(theta/2.));
	  M = E-e*sin(E);                 // Kepler equation
	  if (M<0.) M += 2.*M_PI;



	  // print the resulting orbital elements:
	  printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", time, r.x, r.y, r.z, a, e, i*180./M_PI, Omega*180./M_PI, omega*180./M_PI, M*180./M_PI);
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

