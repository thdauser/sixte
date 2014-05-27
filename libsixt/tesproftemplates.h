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


   Copyright 2014 Thorsten Brand, FAU
*/

#ifndef TESPROFTEMPLFILE_H
#define TESPROFTEMPLFILE_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Structure containing the calorimeter profile templates of one extension. */
typedef struct{

  /** Number of time steps. */
  long Nt;
  
  /** Number of Energy steps. */
  int NE;
  
  /** time array. */
  double *time;
  
  /** Energy array. */
  double *energy;
  
  /** Signal array. */
  double **adc_value;  
  
}TESProfilesEntries;

/** Structure containing the calorimeter profile templates of one file. */
typedef struct{
  
  /** Number of versions (pixels etc.). */
  int Nv;
  
  /** version string array. */
  char **version;
  
  /** Array containing the profiles. */
  TESProfilesEntries *profiles;
  
}TESProfiles;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////

/** Function which initializes the data structure TESProfilesEntries. */
void newTESProfilesEntries(TESProfilesEntries *prof);

/** Function which initializes the data structure TESProfiles. */
TESProfiles* newTESProfiles(int* const status);

/** Destructor for the data structure TESProfilesEntries. */
void destroyTESProfilesEntries(TESProfilesEntries* prof);

/** Destructor for the data structure TESProfiles. */
void destroyTESProfiles(TESProfiles* prof);

/** Function which returns 1 if a version is already loaded, otherwise 0*/
int testTESProfilesExist(TESProfiles *prof, char *version);

/** Function which reads one FITS table with pulse profile templates. */
void readTESProfiles(char *filename, 
			     char *version, 
			     TESProfiles *prof, 
			     int* const status);

/** Function which looks for a specific version and returns the index 
    or -1 if the version is not yet present in the template collection */
int findTESProfileVersionIndex(TESProfiles* prof,
			       char *version);

/** Function which returns the index of the best-matching energy
    index in the proper version array. */
int findTESProfileEnergyIndex(TESProfiles* prof, 
			      int version, 
			      double energy);

#endif /* TESPROFTEMPLFILE_H */