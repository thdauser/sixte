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


   Copyright 2014 Jelle de Plaa, SRON, Thorsten Brand, FAU
*/

#ifndef TESDATASTREAM_H
#define TESDATASTREAM_H 1

#include "sixt.h"
#include "tesproftemplates.h"
#include "pixelimpact.h"
#include "pixelimpactfile.h"
#include "tesnoisespectrum.h"
#include <stdint.h>

#define TESFITSMAXPIX 40

/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Structure containing the data stream */
typedef struct{
  
  /** Number of pixels. */
  int Npix;
  
  /** Number of time steps. */
  long Ntime;
  
  /** Time stamp. */
  double *time;
  
  /** Signal array for all pixels. [time][pixID] */
  uint16_t **adc_value;
  
}TESDataStream;

/** Structure containing the calorimeter pixel properties. */
typedef struct{
  
  /** Number of pixels. */
  int Npix;
  
  /** Array containing the indices of the pixels. */
  int *pixID;
  
  /** Array containing the indices of the version IDs of the pixels. */
  int *versionID;
  
  /** Templates. */
  TESProfiles *templates;
  
}TESPulseProperties;

/** Structure containing the information for a FITS table. */
typedef struct{
  
  /** Name of the stream. */
  char name[9];
  
  /** Number of pixels in the struct. Maximum TESFITSMAXPIX. */
  int Npix;
  
  /** Number of time bins in the struct. */
  long Ntime;
  
  /** Array containing the IDs of the pixels. */
  int *pixID;
  
  /** Time array. */
  double *time;
  
  /** Signal array for all pixels. [pixID][time] */
  uint16_t **adc_value;
  
}TESFitsStream;

/** Linked list containing active pulses */
typedef struct node{
  
  /** Time array */
  double *time;
  
  /** Pulse array containing scaled pulse */
  double *adcpulse;
  
  /** Number of time steps */ 
  long Nt;
  
  /** Current position (index) in the pulse */
  long count;  
  
  /** Next node */
  struct node * next;

}EvtNode;

/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////

/** Function which initializes the data structure TESDataStream. */
TESDataStream* newTESDataStream(int* const status);

/** Function which allocates memory for the TESDataStream structure */
void allocateTESDataStream(TESDataStream* stream, 
			   long Nt, 
			   int Npix, 
			   int* const status);

/** Destructor for the data structure TESDataStream. */
void destroyTESDataStream(TESDataStream* stream);

/** Function which initializes the data structure TESPulseProperties. */
TESPulseProperties* newTESPulseProperties(int* const status);

/** Destructor for the data structure TESPulseProperties. */
void destroyTESPulseProperties(TESPulseProperties* prop);

/** Function which initializes the data structure TESFitsStream. */
TESFitsStream* newTESFitsStream(int* const status);

/** Function which allocates memory for the TESDataStream structure */
void allocateTESFitsStream(TESFitsStream* stream, 
			   long Nt, 
			   int Npix, 
			   int* const status);

/** Destructor for the data structure TESFitsStream. */
void destroyTESFitsStream(TESFitsStream* stream);

/** Creates a new TESFitsStream file. */
void createTESFitsStreamFile(fitsfile **fptr, 
			     char *filename,
			     char* const telescop,
			     char* const instrume,
			     char* const filter,
			     char* const ancrfile,
			     char* const respfile,
			     char* const xmlfile,
			     char* const impactlist,
			     const double mjdref,
			     const double timezero,
			     const double tstart,
			     const double tstop,
			     const char clobber, 
			     int* const status);

/** Function to write a TESFitsStream into a FITS Extension. */
void writeTESFitsStream(fitsfile *fptr, 
			TESFitsStream *stream,
			double tstart,
			double tstop,
			double timeres,
			int* const status);

/** Main engine generating TES data stream */
void getTESDataStream(TESDataStream* TESData, 
		      PixImpFile* PixFile, 
		      TESProfiles* TESProf,
		      AdvDet* det, 
		      double tstart, 
		      double tstop,
		      int Ndetpix,
		      int Nactive,
		      int* activearray,
		      int* const status);

/** Add an event to the node list */
int addEventToNode(EvtNode** ActPulses, 
		   TESProfiles* Pulses, 
                   PixImpact* impact,
		   int pixno,
		   int versionID,
		   int EnID,
		   int* const status);

/** Remove an event from the node list */
void removeEventFromNode(EvtNode** ActPulses, 
			int* pixel);

/** Destroy array of linked lists */
void destroyEventNodes(EvtNode** ActPulses);

/** Initialize array of linked lists */
EvtNode** newEventNodes(int *NPixel, int* const status);

/** Checks if the pixID is in the list of active pixels */
int checkPixIfActive(int pixID, int Npix, int* activearray);

#endif /* TESDATASTREAM_H */