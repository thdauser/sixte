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


   Copyright 2015 Philippe Peille, IRAP
   Copyright 2016-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef TESTRIGGER_H
#define TESTRIGGER_H 1

#include "testriggerfile.h"
#include "tesinitialization.h"

/** Save pixels, NES/NET and monoen keywords to the given FITS file */
void saveTriggerKeywords(fitsfile* fptr,int firstpix,int lastpix,int numberpix,float monoen,
		int* const numberSimulated,int* const numberTrigger,int* const status);

/** Copy pixels, NES/NET and monoen keywords from one file to another */
void copyTriggerKeywords(fitsfile* fptr,fitsfile* fptr2,int* const status);

/** Writes the ADC curves in the TES trigger format*/
void triggerWithImpact(TESDataStream* const stream,TESGeneralParameters * par,
		TESInitStruct* init,float monoen,ReconstructInit* reconstruct_init,int event_list_size,
		const char identify,int* const status);

#endif /* TESTRIGGER_H */
