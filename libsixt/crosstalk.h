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


   Copyright 2016 Thomas Dauser, FAU
*/

#ifndef CROSSTALK_H
#define CROSSTALK_H 1

#include "sixt.h"
#include "advdet.h"
#include "pixelimpact.h"

#define CROSSTALK_ID_NONE 0
#define CROSSTALK_ID_ALL 1
#define CROSSTALK_ID_THERM 2
#define CROSSTALK_ID_ELEC 3
#define CROSSTALK_ID_IMOD 4

////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////

typedef struct {
	int* pixid;
	int* chan;
	double* freq;
	int len;
} channel_list;



////////////////////////////////////////////////////////////////////////
// Function Declarations.
////////////////////////////////////////////////////////////////////////

/** Destructor of the ARF library structure */
void init_crosstalk(AdvDet* det, int* const status);

/** Load the Channel-Frequency List */
channel_list* load_channel_list(char* fname, int* status);

/** Load the Channel-Frequency List */
ReadoutChannels* get_readout_channels(AdvDet* det, int* status);

/** free the channel list */
void free_channel_list(channel_list** chans);

/** Compute influence of the crosstalk event on an impact using the timedependence table */
int computeCrosstalkInfluence(AdvDet* det,PixImpact* impact,PixImpact* crosstalk,double* influence);

#endif /* CROSSTALK_H */


