/* Routine to return the channel for a photon of given input energy */

#include <stdio.h>
#include "heasp.h"
#include "headas_rand.h"

void ReturnChannel(struct RMF *rmf, float energy, int NumberPhoton, long *channel)
{
  long upper, lower, middle, energybin;
  int i, j, igrp, ichan;
  float *sumresponse;
  float *RandomNumber;

  /* initialize the output array to -1s in the event that either the input energy is
     outside the response range or that the response does not sum to unity and events
     can all off the end of the channels. */

  for (i=0; i<NumberPhoton; i++) channel[i] = -1;

  lower = 0;
  upper = rmf->NumberEnergyBins-1;

  /* trap the case of the energy being outside the response range */

  if ( energy < rmf->LowEnergy[lower] || energy > rmf->HighEnergy[upper] ) return;

  /* find the energy bin associated with the input energy - assumes the energies are in increasing order */

  while ( upper - lower > 1 ) {
    middle = (upper + lower)/2;
    if ( energy < rmf->HighEnergy[middle] ) {
      upper = middle;
    } else {
      lower = middle;
    }
  }
  if ( energy > rmf->HighEnergy[lower] ) {
    energybin = upper;
  } else {
    energybin = lower;
  }

  /* generate an array of size channel each element of which is the integrated response up to and
     including that channel */

  sumresponse = (float *) malloc(rmf->NumberChannels*sizeof(float));
  for (i=0; i<rmf->NumberChannels; i++) sumresponse[i] = 0.;
  for (i=0; i<rmf->NumberGroups[energybin]; i++) {
    igrp = i + rmf->FirstGroup[energybin];
    for (j=0; j<rmf->NumberChannelGroups[igrp]; j++) {
      ichan = j + rmf->FirstChannelGroup[igrp];
      sumresponse[ichan] = rmf->Matrix[j+rmf->FirstElement[igrp]];
    }
  }
  for (i=1; i<rmf->NumberChannels; i++) sumresponse[i] += sumresponse[i-1];

  /* generate random numbers between 0 and 1 */

  RandomNumber = (float *) malloc(NumberPhoton*sizeof(float));
  for (i=0; i<NumberPhoton; i++) RandomNumber[i] = (float) HDmtDrand();

  /* loop round the photons */

  for (i=0; i<NumberPhoton; i++) {

  /* find the array element containing this random number. note that we do
     not assume that the total response sums to 1 - if the random number
     exceeds the total response then we assume that the event fell off the
     end of the channel array and return a -1 */

    lower = 0;
    upper = rmf->NumberChannels-1;
    if ( RandomNumber[i] <= sumresponse[upper] ) {
      while ( upper - lower > 1 ) {
	middle = (upper + lower)/2;
	if ( RandomNumber[i] < sumresponse[middle] ) {
	  upper = middle;
	} else {
	  lower = middle;
	}
      }
      if ( RandomNumber[i] > sumresponse[lower] ) {
	channel[i] = upper;
      } else {
	channel[i] = lower;
      }

  /* correct the channel number for the first channel number in use in the response matrix */

      channel[i] += rmf->FirstChannel;

    }

  }

  /* memory tidy-up */

  free(sumresponse);
  free(RandomNumber);

  return;
}
