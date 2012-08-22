/* Function to transpose a matrix with energy rows to channel rows */

#include <string.h>
#include <stdio.h>
#include "heasp.h"

void TransposeRMF(struct RMF *rmf, struct RMFchan *rmfchan)
{
  int i, j, igroup, ielement, counting;
  long *TempFirstEnergyGroup;
  long *TempNumberEnergyGroups;
  long *TempFirstElement;
  long *TempOrderGroup;
  long *CurrentGroup;
  long *GroupOffset;
  float *ChannelVector;

  /* Copy all the standard information which is in common between
     the two structures */

  rmfchan->NumberChannels = rmf->NumberChannels;
  rmfchan->NumberEnergyBins = rmf->NumberEnergyBins;
  rmfchan->NumberTotalElements = 2*rmf->NumberTotalElements;
  rmfchan->FirstChannel = rmf->FirstChannel;
  rmfchan->isOrder = rmf->isOrder;

  rmfchan->LowEnergy = (float *) malloc(rmfchan->NumberEnergyBins*sizeof(float));
  rmfchan->HighEnergy = (float *) malloc(rmfchan->NumberEnergyBins*sizeof(float));

  for (i=0; i<rmf->NumberEnergyBins; i++) {
    rmfchan->LowEnergy[i] = rmf->LowEnergy[i];
    rmfchan->HighEnergy[i] = rmf->HighEnergy[i];
  }

  rmfchan->ChannelLowEnergy = (float *) malloc(rmfchan->NumberChannels*sizeof(float));
  rmfchan->ChannelHighEnergy = (float *) malloc(rmfchan->NumberChannels*sizeof(float));

  for (i=0; i<rmf->NumberChannels; i++) {
    rmfchan->ChannelLowEnergy[i] = rmf->ChannelLowEnergy[i];
    rmfchan->ChannelHighEnergy[i] = rmf->ChannelHighEnergy[i];
  }

  rmfchan->AreaScaling = rmf->AreaScaling;
  rmfchan->ResponseThreshold = rmf->ResponseThreshold;

  strcpy(rmfchan->ChannelType, rmf->ChannelType);
  strcpy(rmfchan->RMFVersion, rmf->RMFVersion);
  strcpy(rmfchan->EBDVersion, rmf->EBDVersion);
  strcpy(rmfchan->Telescope, rmf->Telescope);
  strcpy(rmfchan->Instrument, rmf->Instrument);
  strcpy(rmfchan->Detector, rmf->Detector);
  strcpy(rmfchan->Filter, rmf->Filter);
  strcpy(rmfchan->RMFType, rmf->RMFType);
  strcpy(rmfchan->RMFExtensionName, rmf->RMFExtensionName);
  strcpy(rmfchan->EBDExtensionName, rmf->EBDExtensionName);

  /* Grab the memory for the rmfchan arrays - create temporary arrays for those of
     size NumberTotalGroups since we can't assume that this will be the same as for
     the rmf structure. Set size of these arrays to twice that of the rmf structure. */

  rmfchan->NumberGroups = (long *) malloc(rmfchan->NumberChannels*sizeof(long));
  rmfchan->FirstGroup   = (long *) malloc(rmfchan->NumberChannels*sizeof(long));

  TempFirstEnergyGroup  = (long *) malloc(10*rmf->NumberTotalGroups*sizeof(long));
  TempNumberEnergyGroups= (long *) malloc(10*rmf->NumberTotalGroups*sizeof(long));
  TempFirstElement      = (long *) malloc(10*rmf->NumberTotalGroups*sizeof(long));
  TempOrderGroup        = (long *) malloc(10*rmf->NumberTotalGroups*sizeof(long));

  rmfchan->Matrix       = (float *) malloc(rmfchan->NumberTotalElements*sizeof(float));

  /* Grab memory for and initialize temporary arrays for the channel */

  CurrentGroup  = (long *) malloc(rmfchan->NumberEnergyBins*sizeof(long));
  GroupOffset   = (long *) malloc(rmfchan->NumberEnergyBins*sizeof(long));
  ChannelVector = (float *) malloc(rmfchan->NumberEnergyBins*sizeof(float));

  /* Initialize counters for groups and matrix elements */

  igroup = 0;
  ielement = -1;

  /* Loop round channels */

  for (i=0; i<rmfchan->NumberChannels; i++) {

    for (j=0; j<rmfchan->NumberEnergyBins; j++ ) ChannelVector[j] = 0.0;

    /* Special case for first channel to initialize the temporary arrays */

    if ( i == 0 ) {

      for (j=0; j<rmfchan->NumberEnergyBins; j++) {

	if (rmf->NumberGroups[j] > 0 && 
            rmf->FirstChannelGroup[rmf->FirstGroup[j]] == rmf->FirstChannel) {
	  CurrentGroup[j] = rmf->FirstGroup[j];
	  GroupOffset[j] = 0;
	} else {
	  CurrentGroup[j] = -1;
	  GroupOffset[j] = -1;
	}

      }

    } else {

      /* Calculate the new values of the temporary arrays */

      for (j=0; j<rmfchan->NumberEnergyBins; j++) {

	/* If we haven't started a group yet for this energy */

	if ( CurrentGroup[j] < 0 ) {
	  if (rmf->NumberGroups[j] > 0 && 
              rmf->FirstChannelGroup[rmf->FirstGroup[j]] == (i+rmf->FirstChannel)) {
	    CurrentGroup[j] = rmf->FirstGroup[j];
	    GroupOffset[j] = 0;
	  }
		  
	} else {

	  /* We have started a group. There are two possibilities for the last channel at this
             energy - either it was part of a group in which case GroupOffset >= 0 and CurrentGroup
             gives the group number or it was not part of a group in which case GroupOffset < 0
	     and CurrentGroup gives the next group number */

	  if ( GroupOffset[j] >= 0 ) {

	    /* Check whether this channel takes us off the end of the current response energy group -
	       if so increment CurrentGroup and set GroupOffset to -1, if not just increment GroupOffset */

	    if ( GroupOffset[j] == rmf->NumberChannelGroups[CurrentGroup[j]]-1 ) {
	      GroupOffset[j] = -1;
	      CurrentGroup[j]++;
	    } else {
	      GroupOffset[j]++;
	    }

	  } else {

	    /* We are not in a group so check whether this channel is the start of the next response
	       energy group - if so set GroupOffset to 0 otherwise do nothing */

	    if ( i+rmf->FirstChannel == rmf->FirstChannelGroup[CurrentGroup[j]] ) GroupOffset[j] = 0;

	  }

	}

	/* end loop over energy bins */

      }

    }

    /* Generate the response vector for this channel */

    for ( j=0; j<rmfchan->NumberEnergyBins; j++ ) {

      if ( GroupOffset[j] >= 0 ) {
	ChannelVector[j] = rmf->Matrix[rmf->FirstElement[CurrentGroup[j]]+GroupOffset[j]];
      }

    }

    /* Convert to the compressed form of the response */

    rmfchan->NumberGroups[i] = 0;
    counting = 0;
    for ( j=0; j<rmf->NumberEnergyBins; j++) {

      if (ChannelVector[j] > rmf->ResponseThreshold) {

	ielement++;
	rmfchan->Matrix[ielement] = ChannelVector[j];
	if (counting == 0) {
	  TempFirstEnergyGroup[igroup] = j;
	  TempFirstElement[igroup] = ielement;
	  rmfchan->NumberGroups[i]++;
	  if (rmfchan->NumberGroups[i]==1) rmfchan->FirstGroup[i] = igroup;
	}
	counting++;

      } else {

	if (counting > 0) {
	  TempNumberEnergyGroups[igroup] = counting;
	  igroup++;
	  counting = 0;
	}

      }
    
    }

    /* if we are still in a group then close it out */

    if (counting > 0) {
      TempNumberEnergyGroups[igroup] = counting+1;
      igroup++;
      counting = 0;
    }

    /* end of loop over channels */

  }

  /* The total number of groups created is now igroup + 1 so we can now create the permanent
     arrays and get rid of the temporarys */

  rmfchan->NumberTotalGroups = igroup;
  rmfchan->NumberTotalElements = ielement + 1;

  printf("%ld %ld %ld %ld\n", rmf->NumberTotalGroups, rmfchan->NumberTotalGroups, rmf->NumberTotalElements, rmfchan->NumberTotalElements);

  rmfchan->FirstEnergyGroup = (long *) malloc(rmfchan->NumberTotalGroups*sizeof(long));
  rmfchan->NumberEnergyGroups = (long *) malloc(rmfchan->NumberTotalGroups*sizeof(long));
  rmfchan->FirstElement = (long *) malloc(rmfchan->NumberTotalGroups*sizeof(long));
  rmfchan->OrderGroup = (long *) malloc(rmfchan->NumberTotalGroups*sizeof(long));

  for (i=0; i<rmfchan->NumberTotalGroups; i++) {
    rmfchan->FirstEnergyGroup[i] = TempFirstEnergyGroup[i];
    rmfchan->NumberEnergyGroups[i] = TempNumberEnergyGroups[i];
    rmfchan->FirstElement[i] = TempFirstElement[i];
    rmfchan->OrderGroup[i] = TempOrderGroup[i];
  }

  free(TempFirstEnergyGroup);
  free(TempNumberEnergyGroups);
  free(TempFirstElement);
  free(TempOrderGroup);

  free(CurrentGroup);
  free(GroupOffset);
  free(ChannelVector);

  return;
}
