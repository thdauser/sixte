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
*/

#ifndef TESRECONSTRUCTION_H
#define TESRECONSTRUCTION_H 1

#include "optimalfilters.h"
#include "testriggerfile.h"
#include "teseventlist.h"
#include "tesproftemplates.h"

#define TOOLSUB tesreconstruction_main
#include "headas_main.c"

struct Parameters {
	//File containing the optimal filter
	char OptimalFilterFile[MAXFILENAME];

	//File to reconstruct
	char RecordFile[MAXFILENAME];

	//Ouput event list
	char TesEventFile[MAXFILENAME];

	//File containing the pulse template
	char PulseTemplateFile[MAXFILENAME];

	//Pulse Length
	int PulseLength;

	//Threshold level
	double Threshold;

	//Calibration factor
	double Calfac;

	//Default size of the event list
	int EventListSize;

	//Boolean to choose whether to erase an already existing event list
	char clobber;

	//Boolean to choose to save the run parameters in the output file
	char history;

};

int getpar(struct Parameters* const par);

#endif /* TESRECONSTRUCTION_H */
