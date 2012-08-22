/* General RMF routines
   RMFCompressLine      Compress a line of an RMF
   RMFExpandLine        Expand the RMF for a line
   ReturnRMFElement     Return the RMF element for a channel and energy bin
   ReturnRMFchanElement Same but from the transposed structure
   NormalizeRMF         Ensure that RMF sums to 1.0 for each energy bin
*/


#include <string.h>
#include <stdio.h>
#include "heasp.h"

/************************ prototypes ************************/

void RMFCompressLine(float *Vector, long VectorLength, float Threshold, 
		     long CurrentGroup, long CurrentElement,
		     long *GroupStart, long *GroupSize, long *GroupFirstElement,
		     long *NumberGroups, float *MatrixElements, 
		     long *NumberElements);

void RMFExpandLine(long *GroupStart, long *GroupSize, long *GroupFirstElement, 
		   long NumberGroups, long FirstGroup, float *MatrixElements, 
		   long VectorLength, float *Vector);

float ReturnRMFElement(struct RMF *rmf, long channel, long energybin);

float ReturnRMFchanElement(struct RMFchan *rmfchan, long channel, long energybin);

void NormalizeRMF(struct RMF *rmf);

/************************ RMFCompressLine ************************/

/* Compress a vector and store only those elements above the threshold value
   Arguments :
     float *Vector            i: Input vector
     long VectorLength        i: Size of input vector
     float Threshold          i: Store elements above this threshold value
     long CurrentGroup        i: Current group
     long CurrentElement      i: Current element
     long *GroupStart         r: Start position in vector of each group
                                 NB all pointers to output should be for
				 first element of the original array
     long *GroupSize          r: Size of each group
     long *GroupFirstElement  r: Start position in MatrixElements of group
     long *NumberGroups       r: Number of groups created from the vector
     float *MatrixElements    r: Stored elements
     long *NumberElements     r: Number of stored elements from this vector
*/


void RMFCompressLine(float *Vector, long VectorLength, float Threshold, 
		     long CurrentGroup, long CurrentElement,
                     long *GroupStart, long *GroupSize, long *GroupFirstElement,
		     long *NumberGroups, float *MatrixElements, 
		     long *NumberElements)
{

  long i;
  long igroup=CurrentGroup;
  long ielement=CurrentElement;
  int counting=0;

  /* loop round input vector elements */

  for (i=0; i<VectorLength; i++) {

    /* if this element is above the threshold and we aren't in a group
       then start one. If we are in a group then increment the size */

    if ( Vector[i] > Threshold ) {

      ielement++;
      MatrixElements[ielement] = Vector[i];

      if ( !counting ) {
	igroup++;
	GroupStart[igroup] = i;
	GroupFirstElement[igroup] = ielement;
      }
      counting++;

    } else {

      /* If this element is not above the threshold then close out the group */

      if ( counting ) {
	GroupSize[igroup] = counting+1;
	counting = 0;
      }      

    }

    /* end loop over vector entries */

  }

  /* if we are still in a group then close it out */

  if ( counting ) GroupSize[igroup] = counting+1;

  /* return the number of groups and number of elements created from this
     vector */

  *NumberGroups = igroup - CurrentGroup;
  *NumberElements = ielement - CurrentElement;
	
  return;
}

/************************ RMFExpandLine ************************/

/* Expand a line in the response from the compressed format to a vector
   Arguments :
     long *GroupStart         i: Start position in vector of each group
                                 NB all pointers to input arrays should be for
				 first element of the original array
     long *GroupSize          i: Size of each group
     long *GroupFirstElement  i: Start position in MatrixElements of group
     long NumberGroups        i: Number of groups created from the vector
     long FirstGroup          i: First group in this line
     float *MatrixElements    r: Stored elements
     long *NumberElements     r: Number of stored elements from this vector
     long VectorLength        i: Size of output vector
     float *Vector            r: Output vector
*/


void RMFExpandLine(long *GroupStart, long *GroupSize, long *GroupFirstElement,
		     long NumberGroups, long FirstGroup, float *MatrixElements, 
		     long VectorLength, float *Vector)
{

  long i, j, igroup, ivec, ielt;

  /* initialize the output array */

  for (i=0; i<VectorLength; i++) Vector[i] = 0.0;

  /* loop round input groups */

  for (i=0; i<NumberGroups; i++) {

    igroup = i + FirstGroup;
    ivec = GroupStart[igroup];
    ielt = GroupFirstElement[igroup];

    /* loop round elements in this group - adding them to the output vector */

    for (j=0; j<GroupSize[igroup]; j++) Vector[ivec+j] += MatrixElements[ielt+j];

  }
	
  return;
}

/************************ ReturnRMFElement ************************/

/* Return the RMF element for a channel and energy bin
   Arguments :
     struct *rmf              i: RMF structure
     long channel             i: channel number 
     long energybin           i: energy bin number 
     float ReturnRMFElement   r: Matrix value for this channel and energy
*/


float ReturnRMFElement(struct RMF *rmf, long channel, long energybin)
{
  int i;

  if ( channel < rmf->FirstChannel || channel >= rmf->FirstChannel+rmf->NumberChannels || 
       energybin < 0 || energybin >= rmf->NumberEnergyBins ) return(0.);

  /* loop round groups for this energy bin */

  for(i=rmf->FirstGroup[energybin];i<rmf->FirstGroup[energybin]+rmf->NumberGroups[energybin];i++) {

    if(channel >= rmf->FirstChannelGroup[i] && 
       channel < rmf->FirstChannelGroup[i]+rmf->NumberChannelGroups[i]) {

      return(rmf->Matrix[rmf->FirstElement[i]+channel-rmf->FirstChannelGroup[i]]);

    }

  }

  return(0.);

}

/************************ ReturnRMFchanElement ************************/

/* Return the RMF element for a channel and energy bin using the transposed
   structure.
   Arguments :
     struct *rmfchan             i: RMFchan structure
     long channel                i: channel number 
     long energybin              i: energy bin number 
     float ReturnRMFchanElement  r: Matrix value for this channel and energy
*/

float ReturnRMFchanElement(struct RMFchan *rmfchan, long channel, long energybin)
{
  int i;

  channel -= rmfchan->FirstChannel;

  if ( channel < 0 || channel >= rmfchan->NumberChannels || 
       energybin < 0 || energybin >= rmfchan->NumberEnergyBins ) return(0.);

  /* loop round groups for this channel */

  for(i=rmfchan->FirstGroup[channel];i<rmfchan->FirstGroup[channel]+rmfchan->NumberGroups[channel];i++) {

    if(energybin >= rmfchan->FirstEnergyGroup[i] && 
       energybin < rmfchan->FirstEnergyGroup[i]+rmfchan->NumberEnergyGroups[i]) {

      return(rmfchan->Matrix[rmfchan->FirstElement[i]+energybin-rmfchan->FirstEnergyGroup[i]]);

    }

  }

  return(0.);

}

/************************ ReturnRMFElement ************************/

/* Normalize the RMF so it sums to 1.0 for each energy bin
   Arguments :
     struct *rmf             i/r: RMF structure
*/

void NormalizeRMF(struct RMF *rmf)
{
  int ie, i, j, igrp;
  float sumresp;

  /* Loop over energies */

  for (ie=0; ie<rmf->NumberEnergyBins; ie++) {

    /* sum up the response in this energy */

    sumresp = 0.0;

    for (i=0; i<rmf->NumberGroups[ie]; i++) {
       igrp = i + rmf->FirstGroup[ie];
       for (j=0; j<rmf->NumberChannelGroups[igrp]; j++) {
         sumresp += rmf->Matrix[j+rmf->FirstElement[igrp]];
       }
    }

    /* divide through by the summed response */

    for (i=0; i<rmf->NumberGroups[ie]; i++) {
       igrp = i + rmf->FirstGroup[ie];
       for (j=0; j<rmf->NumberChannelGroups[igrp]; j++) {
         rmf->Matrix[j+rmf->FirstElement[igrp]] /= sumresp;
       }
    }

  }

  return;
}
