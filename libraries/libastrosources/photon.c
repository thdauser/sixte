#include "photon.h"


///////////////////////////////////////////////////////
int create_lightcurve(
		      PointSource* ps, 
		      double t0,             // start time for lightcurve
     		      gsl_rng *gsl_random_g  // pointer to GSL random number generator
		      )
{
  int count;

  int status = EXIT_SUCCESS;    // error status
  char msg[MAXMSG];             // error message output buffer


  // At first light curve creation for this source, get memory for lightcurve.
  if (ps->lightcurve == NULL) {
    ps->lightcurve = (struct lightcurve_entry *) 
      malloc((N_LIGHTCURVE_BINS+1) * sizeof(struct lightcurve_entry));
    if (!(ps->lightcurve)) {
      status = EXIT_FAILURE;
      sprintf(msg, "Not enough memory available to store the lightcurve!\n");
      HD_ERROR_THROW(msg,status);
    }
  }


  // Generate the light curve:
#ifdef CONSTANT_LIGHTCURVE // Create constant light curve with average count rate.
  for (count=0; count<N_LIGHTCURVE_BINS; count++) {
    ps->lightcurve[count].t = t0 + count * LIGHTCURVE_BINWIDTH;
    ps->lightcurve[count].rate = ps->rate;
  }
  ps->lightcurve[N_LIGHTCURVE_BINS].t = 
    t0 + N_LIGHTCURVE_BINS*LIGHTCURVE_BINWIDTH;

#else  // Create a light curve with power-law PSD and read noise.

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
    const float mean_rate = ps->rate;
    const float sigma = mean_rate/3.; // (mean_rate-10.)/3.;

    // normalize and copy FFT array to light curve array
    if (status != EXIT_FAILURE) {
      for (count=0; count<N_LIGHTCURVE_BINS; count++) {
	ps->lightcurve[count].t = t0 + count * LIGHTCURVE_BINWIDTH;
	ps->lightcurve[count].rate = 
	  (fcomp[count]-mean) *sigma/sqrt(variance) + mean_rate;  // TODO

	// Avoid negative or close2zero count rates (no physical meaning):
	if (ps->lightcurve[count].rate < 0.01*mean) { check = -1; }
      }

      ps->lightcurve[N_LIGHTCURVE_BINS].t = 
	t0 + N_LIGHTCURVE_BINS * LIGHTCURVE_BINWIDTH;
    }

  } while (check==-1);  // Repeat the light curve creation until we have proper data.

#endif // END of #ifndef CONSTANT_LIGHTCURVE

  return(status);
}



///////////////////////////////////////////////////////////////////////////////////
// Create a randomly chosen photon energy according to the spectrum of the 
// specified source.
// The used spectrum is given in [PHA channels].
// The function returns the photon energy in [keV].
float photon_energy(Spectrum* pha_spectrum, struct RMF* rmf)
{
  // Get a random PHA channel according to the given PHA distribution.
  float rand = (float)get_random_number();
  long upper = pha_spectrum->NumberChannels-1, lower=0, mid;
  
  if(rand > pha_spectrum->rate[pha_spectrum->NumberChannels-1]) {
    printf("PHA sum: %f < RAND: %f\n", 
	   pha_spectrum->rate[pha_spectrum->NumberChannels-1], rand);
  }

  // Determine the energy of the photon.
  while (upper-lower>1) {
    mid = (long)((lower+upper)/2);
    if (pha_spectrum->rate[mid] < rand) {
      lower = mid;
    } else {
      upper = mid;
    }
  }
    
  if (pha_spectrum->rate[lower] < rand) {
    lower = upper;
  }

  // Return an energy chosen randomly out of the determined PHA bin:
  return(rmf->ChannelLowEnergy[lower] + 
	 get_random_number()*(rmf->ChannelHighEnergy[lower]-
			      rmf->ChannelLowEnergy[lower]));
}



///////////////////////////////////////////////////////////////////////////////////
int create_photons(
		   PointSource* ps /**< Source data. */,
		   double time /**< Current time. */, 
		   double dt /**< Time interval for photon generation. */,       
		   /** Address of pointer to time-ordered photon list.*/
		   struct PhotonOrderedListEntry** list_first,
		   struct RMF* rmf,     
		   gsl_rng *gsl_random_g
		   )
{
  int bin=0; // light curve bin counter

  // Second pointer to photon list, that can be moved along the list,   
  // without loosing the first entry.
  struct PhotonOrderedListEntry* list_current = *list_first;

  int status=EXIT_SUCCESS;

  // If there is no photon time stored so far, set the current time.
  if (ps->t_last_photon < 0.) {
    ps->t_last_photon = time;
  }


  // Create photons and insert them into the given time-ordered list:
  while (ps->t_last_photon < time+dt) {
    Photon new_photon; // buffer for new photon
    new_photon.ra = ps->ra;
    new_photon.dec = ps->dec; 
    new_photon.direction = unit_vector(ps->ra, ps->dec); // REMOVE

    // Create the energy of the new photon
    new_photon.energy = photon_energy(ps->spectrum, rmf);

    // Determine the current count rate of the light curve.
    // If the source has no light curve, or the assigned light curve is 
    // too old, create a new one.
    if (ps->lightcurve == NULL) {
      if ((status=create_lightcurve(ps, time, gsl_random_g)) != EXIT_SUCCESS) break;
    } else if (ps->lightcurve[N_LIGHTCURVE_BINS].t < ps->t_last_photon) {
      if ((status=create_lightcurve(ps, time, gsl_random_g)) != EXIT_SUCCESS) break;
    }

    // Find light curve time bin corresponding to the current time.
    for ( ; bin<N_LIGHTCURVE_BINS; bin++) {
      if (ps->lightcurve[bin+1].t > ps->t_last_photon) {
	break;
      }
    }

    // Calculate arrival time depending on previous photon.
    ps->t_last_photon += rndexp(1./(ps->lightcurve[bin].rate));
    new_photon.time = ps->t_last_photon;


    // Insert photon to the global photon list:
    if ((status=insert_Photon2TimeOrderedList(list_first, &list_current, &new_photon)) 
	!= EXIT_SUCCESS) break;  
    
  } // END of loop 'while(t_last_photon < time+dt)'

  return(status);
}




////////////////////////////////////////////////////////////////
void clear_PhotonList(struct PhotonOrderedListEntry** pole) {
  if ((*pole) != NULL) {
    // This is not the last entry in the list, so call routine recursively.
    clear_PhotonList(&((*pole)->next));

    // Free memory and reset pointer to NULL.
    free(*pole);
    *pole = NULL;  
  }
}



////////////////////////////////////////////////////////////////
int insert_Photon2TimeOrderedList(struct PhotonOrderedListEntry** first,
				  struct PhotonOrderedListEntry** current, 
				  Photon* ph) 
{
  // The iterator is shifted over the time-ordered list to find the right entry.
  struct PhotonOrderedListEntry** iterator = current;

  // Find the right entry where to insert the new photon.
  while ((NULL!=*iterator) && ((*iterator)->photon.time < ph->time)) {
    iterator = &((*iterator)->next);
  }
  // Now '*iterator' points to the entry before which the new photon has to be inserted.
  // I.e. 'iterator' is equivalent to '&[fromer_entry]->next'. 
  // Therefore changeing '*iterator' means redirecting the chain.
  // '*iterator' has to be redirected to the new entry.
    
  // Create a new PhotonOrderedListEntry and insert it before '**iterator'.
  struct PhotonOrderedListEntry* new_entry = 
    (struct PhotonOrderedListEntry*)malloc(sizeof(struct PhotonOrderedListEntry));
  if (NULL==new_entry) {
    HD_ERROR_THROW("Error: Could not allocate memory for new photon entry!\n", EXIT_FAILURE);
    return(EXIT_FAILURE);
  }

  // Set the values of the new entry:
  new_entry->photon = *ph;
  new_entry->next = *iterator;
  if (*iterator == *first) {
    // If the new photon was inserted as first entry of the list, the pointer to this 
    // first entry has to be redirected.
    *first = new_entry;
  }
  *iterator = new_entry;

  // The pointer '*current' should point to the new entry in the time-ordered list.
  *current = new_entry;
  
  return(EXIT_SUCCESS);
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
    strcpy(ftype[3], "Y");
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




/////////////////////////////////////////////////////////////////////
int insert_Photon2BinaryTree(struct PhotonBinaryTreeEntry** ptr, Photon* ph)
{
  int status = EXIT_SUCCESS;

  if (NULL==*ptr) {
    // Reached an end of the tree. So create a new entry and insert it here.
    *ptr = (struct PhotonBinaryTreeEntry*)malloc(sizeof(struct PhotonBinaryTreeEntry));
    if (NULL==*ptr) { // Check if memory allocation was successfull
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for binary tree failed!\n", status);
      return(status);
    }

    (*ptr)->photon = *ph;
    (*ptr)->sptr = NULL;
    (*ptr)->gptr = NULL;

  } else {
    // We have to go deeper into the tree. So decide whether the new photon
    // is earlier or later than the current entry.
    if ((*ptr)->photon.time > ph->time) {
      status=insert_Photon2BinaryTree(&(*ptr)->sptr, ph);
    } else {
      status=insert_Photon2BinaryTree(&(*ptr)->gptr, ph);
    }
  }

  return(status);
}



////////////////////////////////////////////////////////////////////
int CreateOrderedPhotonList(struct PhotonBinaryTreeEntry** tree_ptr,
			    struct PhotonOrderedListEntry** list_first,
			    struct PhotonOrderedListEntry** list_current)
{
  int status=EXIT_SUCCESS;

  // Check if the current tree entry exists.
  if (NULL != *tree_ptr) {

    // Add the entries before the current one to the time-ordered list.
    status = CreateOrderedPhotonList(&(*tree_ptr)->sptr, list_first, list_current);
    if (EXIT_SUCCESS != status) return(status);

    // Insert the current entry into the time-ordered list at the right position.
    status = insert_Photon2TimeOrderedList(list_first, list_current, &(*tree_ptr)->photon);
    if (EXIT_SUCCESS != status) return(status);
    
    // Add the entries after the current one to the time-ordered list.
    status = CreateOrderedPhotonList(&(*tree_ptr)->gptr, list_first, list_current);
    if (EXIT_SUCCESS != status) return(status);

    // Delete the binary tree entry as it is not required any more.
    free(*tree_ptr);
    *tree_ptr=NULL;

  } // END if current tree entry exists.

  return(status);
}
						       

