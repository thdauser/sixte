#include "rmf.h"


struct RMF* loadRMF(const char* const filename, int* const status) 
{
  fitsfile* fptr=NULL;

  struct RMF* rmf = (struct RMF*)malloc(sizeof(struct RMF));
  if (NULL==rmf) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: could not allocate memory for RMF!\n", *status);
    return(rmf);
  }

  // Load the RMF from the FITS file using the HEAdas RMF access routines
  // (part of libhdsp).
  fits_open_file(&fptr, filename, READONLY, status);
  if (*status!=EXIT_SUCCESS) return(rmf);
  
  // Read the 'SPECRESP MATRIX' or 'MATRIX' extension:
  if ((*status=ReadRMFMatrix(fptr, 0, rmf))!=EXIT_SUCCESS) return(rmf);

  // Print some information:
  headas_chat(5, "RMF loaded with %ld energy bins and %ld channels\n",
	      rmf->NumberEnergyBins, rmf->NumberChannels);

#ifdef NORMALIZE_RMF
  // Normalize the RMF:
  headas_printf("### Warning: RSP/RMF is explicitly normalized!                 ###\n");
  headas_printf("### ARF contributions must be contained in the input spectrum! ###\n");
  NormalizeRMF(rmf);
#else
  // Check if the RSP file contains matrix rows with a sum of more than 1.
  // In that case the RSP probably also contains the mirror ARF, what should 
  // normally not be the case for this simulation.
  long chancount, bincount;
  double maxsum = 0.;
  for (bincount=0; bincount<rmf->NumberEnergyBins; bincount++) {
    double sum = 0.;
    for (chancount=0; chancount<rmf->NumberChannels; chancount++) {
      sum += ReturnRMFElement(rmf, chancount, bincount);
    }
    if (sum > maxsum) {
      maxsum = sum;
    }
  }
  if (fabs(maxsum-1.)>1.e-3) {
    headas_printf("### Warning: RSP probably contains ARF "
		  "(row-sum = %lf)! ###\n", maxsum);
  }
#endif

  // Read the EBOUNDS extension:
  if ((*status=ReadRMFEbounds(fptr, 0, rmf))!=EXIT_SUCCESS) return(rmf);

  // Close the open FITS file.
  fits_close_file(fptr, status);
  
  return(rmf);
}



void freeRMF(struct RMF* const rmf) 
{
  if (NULL!=rmf) {
    free(rmf);
  }
}



long getEBOUNDSChannel(const float energy, const struct RMF* const rmf)
{
  // Check if the charge is outside the range of the energy bins defined
  // in the EBOUNDS table. In that case the return value of this function is '-1'.
  if (rmf->ChannelLowEnergy[0] > energy) {
    return(0); // TODO
  } else if (rmf->ChannelHighEnergy[rmf->NumberChannels-1] < energy) {
    return(rmf->NumberChannels - 1 + rmf->FirstChannel);
  }
  
  // Perform a binary search to obtain the detector PHA channel 
  // that corresponds to the given detector charge.
  long min, max, mid;
  min = 0;
  max = rmf->NumberChannels-1;
  while (max > min) {
    mid = (min+max)/2;
    if (rmf->ChannelHighEnergy[mid] < energy) {
      min = mid+1;
    } else {
      max = mid;
    }
  }
  
  // Return the PHA channel.
  return(min + rmf->FirstChannel);
}



float getEBOUNDSEnergy(long channel, const struct RMF* const rmf, 
		       const int boundary)
{
  // Subtract the channel offset (EBOUNDS may either start at 0 or at 1).
  channel -= rmf->FirstChannel;
  if ((channel < 0) || (channel >= rmf->NumberChannels)) {
    return(-1.);
  }

  // Return the randomized mean / lower / upper energy value corresponding 
  // to the specified PHA channel according to the EBOUNDS table.
  if (0==boundary) {
    return(rmf->ChannelLowEnergy[channel] +
	   sixt_get_random_number()*(rmf->ChannelHighEnergy[channel]-
				     rmf->ChannelLowEnergy[channel]));
  } else if (-1==boundary) {
    return(rmf->ChannelLowEnergy[channel]);
  } else {
    return(rmf->ChannelHighEnergy[channel]);
  }
}


