/** 
 * This file contains routines/functions for photon handling.
 */
#include "photon.h"


///////////////////////////////////////////////////////
int create_lightcurve(
		      struct source_cat_entry *src,// source catalog entry
		      double time,           // start time for lightcurve
     		      gsl_rng *gsl_random_g  // pointer to GSL random number generator
		      )
{
  int count;
  int status = EXIT_SUCCESS;    // error status
  char msg[MAXMSG];             // error message output buffer


  // At first light curve creation for this source, get memory for lightcurve.
  if (src->lightcurve == NULL) {
    src->lightcurve = (struct lightcurve_entry *) 
      malloc((N_LIGHTCURVE_BINS+1) * sizeof(struct lightcurve_entry));
    if (!(src->lightcurve)) {
      status = EXIT_FAILURE;
      sprintf(msg, "Not enough memory available to store the lightcurve!\n");
      HD_ERROR_THROW(msg,status);
    }
  }


  // Generate the light curve:
#ifdef CONSTANT_LIGHTCURVE // Create constant light curve with average count rate.
  for (count=0; count<N_LIGHTCURVE_BINS; count++) {
    src->lightcurve[count].t = time + count * LIGHTCURVE_BINWIDTH;
    src->lightcurve[count].rate = src->rate;
  }
  src->lightcurve[N_LIGHTCURVE_BINS].t = 
    time + N_LIGHTCURVE_BINS*LIGHTCURVE_BINWIDTH;

#else  // Create a light curve with power law PSD and read noise.

  // Repeat the light curve creation until we have proper data 
  // (without zero-rates).
  int check;
  do { 
    check=0;

    // First create frequencies.
    double freq[N_LIGHTCURVE_BINS];
    for (count=0; count<N_LIGHTCURVE_BINS; count++) {
      freq[count] = (double)count;
    }
    // Create PSD.
    double psd[N_LIGHTCURVE_BINS];
    for (count=0; count<N_LIGHTCURVE_BINS; count++) {
      // create red noise: ~1/f² (Uttley, McHardy 2001)
      psd[count] = pow(1./freq[count], 2.);  
    }
    // Create Fourier components.
    double fcomp[N_LIGHTCURVE_BINS];
    fcomp[0] = 0.;
    for (count=1; count<N_LIGHTCURVE_BINS/2; count++) {
      REAL(fcomp, count) = gsl_ran_ugaussian(gsl_random_g)*sqrt(psd[count]);
      IMAG(fcomp, count) = gsl_ran_ugaussian(gsl_random_g)*sqrt(psd[count]);
    }
    fcomp[N_LIGHTCURVE_BINS/2] = 
      gsl_ran_ugaussian(gsl_random_g)*sqrt(psd[N_LIGHTCURVE_BINS/2]); // TODO 

    // Perform Fourier transformation.
    gsl_fft_halfcomplex_radix2_inverse(fcomp, 1, N_LIGHTCURVE_BINS);
  
    // Calculate mean and variance
    double mean=0., variance=0.;
    for (count=0; count<N_LIGHTCURVE_BINS; count++) {
      mean += fcomp[count];
      variance += pow(fcomp[count], 2.); 
    }
    mean = mean/(double)N_LIGHTCURVE_BINS;
    variance = variance/(double)N_LIGHTCURVE_BINS;
    variance -= pow(mean,2.);     // var = <x^2>-<x>^2

    // desired properties of the light curve
    const float mean_rate = src->rate;
    const float sigma = mean_rate/3.; // (mean_rate-10.)/3.;

    // normalize and copy FFT array to light curve array
    if (status != EXIT_FAILURE) {
      for (count=0; count<N_LIGHTCURVE_BINS; count++) {
	src->lightcurve[count].t = time + count * LIGHTCURVE_BINWIDTH;
	src->lightcurve[count].rate = 
	  (fcomp[count]-mean) *sigma/sqrt(variance) + mean_rate;  // TODO

	// Avoid negative or close2zero count rates (no physical meaning):
	if (src->lightcurve[count].rate < 0.01*mean) { check = -1; }
      }

      src->lightcurve[N_LIGHTCURVE_BINS].t = 
	time + N_LIGHTCURVE_BINS * LIGHTCURVE_BINWIDTH;
    }

  } while (check==-1);  // Repeat the light curve creation until we have proper data.

#endif // END of #ifndef CONSTANT_LIGHTCURVE

  /*
  // plot light curve for debugging
  double avg=0.;
  FILE *pfile=NULL;
  pfile = fopen("lightcurve.dat", "w+");
  for (count=0; count < N_LIGHTCURVE_BINS; count++) {
    avg += src->lightcurve[count].rate/N_LIGHTCURVE_BINS;
    fprintf(pfile, "%lf %lf\n", src->lightcurve[count].t, src->lightcurve[count].rate);
  }
  printf("avg: %lf\n", avg);
  fclose(pfile);
  */

  return(status);
}



///////////////////////////////////////////////////////////////////////////////////
// Create a randomly chosen photon energy according to the spectrum of the 
// specified source.
// The used spectrum is given in [PHA channels].
// The function returns the photon energy in [keV].
float photon_energy(
		    struct source_cat_entry src,     // source data
		    // Detector information like Nchannels and ebounds
		    Detector* detector
		    )
{
  // get a random PHA value according to given PHA distribution
  double rand = get_random_number();
  long upper = detector->Nchannels-1, lower=0, mid;
  
  // determine the energy of the photon
  while (upper-lower>1) {
    mid = (int)((lower+upper)/2);
    if (src.spectrum->data[mid] < rand) {
      lower = mid;
    } else {
      upper = mid;
    }
  }
    
  if (src.spectrum->data[lower] < rand) {
    lower = upper;
  }

  // return energy chosen randomly from the determined PHA bin
  return(detector->ebounds.row[lower].E_min + 
	 (detector->ebounds.row[lower].E_max-detector->ebounds.row[lower].E_min)*
	 get_random_number());
}





///////////////////////////////////////////////////////////////////////////////////
int create_photons(
		   struct source_cat_entry* src,  // source data
		   // current time and time interval for photon creation
		   double time, double dt,       
		   struct Photon_Entry** pl,      // time ordered photon list
		   // Detector information (Nchannels, ebounds)
		   Detector* detector,     
		   gsl_rng *gsl_random_g
		   )
{
  int bin=0;                         // light curve bin counter
  // Flag: first photon in the generation process for this source, i.e.,
  // in the run of this function. This is needed, because the first photon
  // can be inserted at the beginning of the photon list.
  int first_photon=1;              
  // Alternative pointer to photon list, that can be moved along the list,   
  // without loosing the first entry.
  struct Photon_Entry* pe = NULL;  
  int status=EXIT_SUCCESS;


  // Set pe to point to the same value as pl 
  // (pl is the pointer to the first entry in the entire photon list.).
  // So pe can be updated and modified without loosing the first entry 
  // of the photon list.
  // The pointer pe is updated by insert_photon, so the scanning of the 
  // photon list does not have to start at the first 
  // entry for each individual photon.
  // For one and the same source, the photons are created in temporal order.

  // Copy the pointer to the photon list.
  pe = *pl;
  
  // if there is no photon time stored so far
  if (src->t_last_photon < 0.) {
    src->t_last_photon = time;
  }


  // create photons and insert them in the given time ordered list
  while (src->t_last_photon < time+dt) {
    struct Photon new_photon;        // buffer for new photon
    new_photon.ra = src->ra;
    new_photon.dec = src->dec; 
    new_photon.direction = src->r;   // REMOVE

    // Create the energy of the new photon
    new_photon.energy = photon_energy(*src, detector);

    // Determine the current count rate of the light curve.
    // If the source has no light curve, or the assigned light curve is 
    // too old, create a new one.
    if (src->lightcurve == NULL) {
      if ((status=create_lightcurve(src, time, gsl_random_g)) != EXIT_SUCCESS) break;
    } else if (src->lightcurve[N_LIGHTCURVE_BINS].t < src->t_last_photon) {
      if ((status=create_lightcurve(src, time, gsl_random_g)) != EXIT_SUCCESS) break;
    }

    // Find light curve time bin corresponding to the current time.
    for ( ; bin<N_LIGHTCURVE_BINS; bin++) {
      if (src->lightcurve[bin+1].t > src->t_last_photon) {
	break;
      }
    }

    // calculate arrival time depending on former photon creation
    src->t_last_photon += rndexp(1./(src->lightcurve[bin].rate));
    //printf("%lf %lf\n", src->t_last_photon, src->lightcurve[bin].rate);
    //    printf("%d %f\n", bin, src->lightcurve[bin].rate);
    new_photon.time = src->t_last_photon;


    if (first_photon == 1) {
      // The first photon for the list might be inserted at the 
      // beginning of the time-ordered photon list.
      // Therefore this case requires a special treatment.
      if ((status=insert_photon(pl, new_photon)) != EXIT_SUCCESS) break;
      pe = *pl;
      first_photon = 0;
    } else {
      // If this is not the first photon for this source in the current 
      // run of this function, give the variable pointer 'pe' to the 
      // insertion procedure. 'pe' might point to later events in the 
      // time-ordered photon list, whereas the beginning of the list
      // is stored in '*pl'. The use of 'pe' can improve the program performance.

      // Insert photon to the global photon list:
      if ((status=insert_photon(&pe, new_photon)) != EXIT_SUCCESS) break;  
    }
  }

  return(status);
}



////////////////////////////////////////////////////////////////
// Clears the photon list.
void clear_photon_list(struct Photon_Entry **pe) {
  if ((*pe) != NULL) {
    // this is not the last entry in the list, so call routine recursively
    clear_photon_list(&((*pe)->next_entry));

    // free memory and reset pointer to NULL
    free(*pe);
    *pe = NULL;  
  }
}



////////////////////////////////////////////////////////////////
// Inserts a new photon into the time-ordered photon list.
// The return value is the value of the error status.
int insert_photon(struct Photon_Entry **pe, struct Photon ph) {
  char msg[MAXMSG];

  while (((*pe)!=NULL) && (ph.time > (*pe)->photon.time)) {
    pe = &((*pe)->next_entry);
  }

  // create a new Photon_Entry and insert it at the current position
  struct Photon_Entry *buffer=NULL;
  buffer = (struct Photon_Entry *) malloc(sizeof(struct Photon_Entry));

  if (buffer==NULL) {
    sprintf(msg, "Error: Could not allocate memory for photon entry!\n");
    HD_ERROR_THROW(msg,EXIT_FAILURE);
    return(EXIT_FAILURE);
  } else {
    buffer->photon = ph;
    buffer->next_entry = *pe;
    *pe = buffer;
    return(EXIT_SUCCESS);
  }
}




///////////////////////////////////////////////////////////////////
int create_photonlist_file(
			   fitsfile **fptr,
			   char filename[],
			   int *status
			   )
{
  char *ftype[N_PHOTON_FIELDS];
  char *fform[N_PHOTON_FIELDS];
  char *funit[N_PHOTON_FIELDS];
  int counter;

  char msg[MAXMSG];  // error output buffer

  do { // Beginning of ERROR handling loop

    // Create a new FITS file:
    if (fits_create_file(fptr, filename, status)) break;

    // To create a FITS table, the format of the individual columns has to 
    // be specified.
    for(counter=0; counter<N_PHOTON_FIELDS; counter++) {
      // Allocate memory
      ftype[counter] = (char *) malloc(8 * sizeof(char));
      fform[counter] = (char *) malloc(4 * sizeof(char));
      funit[counter] = (char *) malloc(20 * sizeof(char));

      // Check if all memory was allocated successfully:
      if ((!ftype[counter]) || (!fform[counter]) || (!funit[counter])) {
	*status = EXIT_FAILURE;
	sprintf(msg, "Error: no memory allocation for FITS table parameters "
		"failed (photon list)!\n");
	HD_ERROR_THROW(msg, *status);
      }
    }

    // If an error has occurred during memory allocation, 
    // skip the following part.
    if (*status != EXIT_SUCCESS) break;

    // Set the field types of the table in the FITS file.
    // 1. time
    strcpy(ftype[0], "TIME");
    strcpy(fform[0], "D");
    strcpy(funit[0], "s");

    // 2. energy
    strcpy(ftype[1], "ENERGY");
    strcpy(fform[1], "E");
    strcpy(funit[1], "keV");

    // 3. right ascension
    strcpy(ftype[2], "RA");
    strcpy(fform[2], "D");
    strcpy(funit[2], "decimal degrees");

    // 4. declination
    strcpy(ftype[3], "DEC");
    strcpy(fform[3], "D");
    strcpy(funit[3], "decimal degrees");

    // create the table
    if (fits_create_tbl(*fptr, BINARY_TBL, 0, N_PHOTON_FIELDS, 
			ftype, fform, funit, "PHOTONLIST", status)) break;
    

    // write descriptory data into the header of the FITS file
    if (fits_write_key(*fptr, TSTRING, "COMMENT", "PHOTONLIST",
		       "", status)) break;

    // If desired by the user, print all program parameters to HISTORY of 
    // FITS file (HDU number 1).
    HDpar_stamp(*fptr, 2, status);
    
  } while (0);  // END of ERROR handling loop


  //----------------
  // clean up
  for (counter=0; counter<N_PHOTON_FIELDS; counter++) {
    if (ftype[counter]) free(ftype[counter]);
    if (fform[counter]) free(fform[counter]);
    if (funit[counter]) free(funit[counter]);
  }

  return(*status);
}



///////////////////////////////////////////////////////////////////
int create_impactlist_file(
			   fitsfile **fptr,
			   char filename[],
			   int *status
			   )
{
  char *ftype[N_IMPACT_FIELDS];
  char *fform[N_IMPACT_FIELDS];
  char *funit[N_IMPACT_FIELDS];
  int counter;

  char msg[MAXMSG];  // error output buffer

  do { // Beginning of ERROR handling loop

    // Create a new FITS file:
    if (fits_create_file(fptr, filename, status)) break;

    // To create a FITS table, the format of the individual columns has to 
    // be specified.
    for(counter=0; counter<N_PHOTON_FIELDS; counter++) {
      // Allocate memory
      ftype[counter] = (char *) malloc(8 * sizeof(char));
      fform[counter] = (char *) malloc(4 * sizeof(char));
      funit[counter] = (char *) malloc(20 * sizeof(char));

      // Check if all memory was allocated successfully:
      if ((!ftype[counter]) || (!fform[counter]) || (!funit[counter])) {
	*status = EXIT_FAILURE;
	sprintf(msg, "Error: no memory allocation for FITS table parameters "
		"failed (impact list)!\n");
	HD_ERROR_THROW(msg, *status);
      }
    }

    // If an error has occurred during memory allocation, 
    // skip the following part.
    if (*status != EXIT_SUCCESS) break;

    // Set the field types of the table in the FITS file.
    // 1. time
    strcpy(ftype[0], "TIME");
    strcpy(fform[0], "D");
    strcpy(funit[0], "s");

    // 2. energy
    strcpy(ftype[1], "ENERGY");
    strcpy(fform[1], "E");
    strcpy(funit[1], "keV");

    // 3. right ascension
    strcpy(ftype[2], "X");
    strcpy(fform[2], "D");
    strcpy(funit[2], "m");

    // 4. declination
    strcpy(ftype[3], "y");
    strcpy(fform[3], "D");
    strcpy(funit[3], "m");

    // create the table
    if (fits_create_tbl(*fptr, BINARY_TBL, 0, N_PHOTON_FIELDS, 
			ftype, fform, funit, "IMPACTLIST", status)) break;
    

    // write descriptory data into the header of the FITS file
    if (fits_write_key(*fptr, TSTRING, "COMMENT", "IMPACTLIST", "", status)) break;

    // If desired by the user, print all program parameters to HISTORY of 
    // FITS file (HDU number 1).
    HDpar_stamp(*fptr, 2, status);
    
  } while (0);  // END of ERROR handling loop


  //----------------
  // clean up
  for (counter=0; counter<N_PHOTON_FIELDS; counter++) {
    if (ftype[counter]) free(ftype[counter]);
    if (fform[counter]) free(fform[counter]);
    if (funit[counter]) free(funit[counter]);
  }

  return(*status);
}


