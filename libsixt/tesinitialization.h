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


   Copyright 2014 Philippe Peille, IRAP
*/

#ifndef TESINITIALIZATION_H
#define TESINITIALIZATION_H 1

#include "sixt.h"
#include "pixelimpact.h"
#include "advdet.h"
#include "pixelimpactfile.h"
#include "tesproftemplates.h"
#include "tesnoisespectrum.h"
#include "tesdatastream.h"
#include "testriggerfile.h"

////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

/** Structure containing all the possible paramaters of the tools
    simulating the TES data chain*/
typedef struct {

  char PixImpList[MAXFILENAME];
  char XMLFile[MAXFILENAME];
  char streamname[MAXFILENAME];
  char tesTriggerFile[MAXFILENAME];
  char TesEventFile[MAXFILENAME];
  
  char activePixels[9];
  int Nactive;
  int Npix;
  int nlo;
  int nhi;
  int triggerSize;
  int preBufferSize;
  
  double tstart;
  double tstop;
  
  char writeStreamFile;
  char clobber;
  char history;
  char Reconstruct;
  char WriteRecordFile;
  char check_times;

  unsigned long int seed;

  double tstart_stream;
  
} TESGeneralParameters;

/** Structure containing all the pointers and values necessary to run a
    simulation of the TES data stream optionnaly wrapping
    the generation of triggers*/
typedef struct {
  /** Structure around the impact file */
  PixImpFile* impfile;

  /** TES profiles to load */
  TESProfiles* profiles;

  /** Advanced detector structure */
  AdvDet *det;

  /** Record file */
  TesTriggerFile* record_file;

  /** Event File */
  TesEventFile* event_file;
  
  /** Array of active pixels */
  int *activearray;
  
  /** Array of event numbers */
  long *Nevts;
  
  //Keywords
  char telescop[MAXMSG];
  char instrume[MAXMSG];
  char filter[MAXMSG];
  char ancrfile[MAXMSG];
  char respfile[MAXMSG];
  double mjdref;
  double timezero;
  double tstart;
  double tstop;
  

} TESInitStruct;

////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////

/** Initializes the different variables necessary fo the simulations. Depending
    on the tool calling this function, not all the variables are set. */
void tesinitialization(TESInitStruct* const init,TESGeneralParameters* const par, int* const status);

/** Constructor. Returns a pointer to an empty TESInitStruct data
    structure. */
TESInitStruct* newInitStruct(int* const status);

/** Destructor. */
void freeTESInitStruct(TESInitStruct** const init, int* const status);

#endif /* TESINITIALIZATION_H */
