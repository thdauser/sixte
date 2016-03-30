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

#ifndef TESWITHRMF_H
#define TESWITHRMF_H 1

#include "sixt.h"
#include "pixelimpactfile.h"
#include "advdet.h"
#include "crosstalk.h"

#define INVGRADE -1
#define PILEUP -2
#define CROSSTALK -3
#define WRONGE -4


typedef struct{
	double next,current,previous;
}gradingTimeStruct;

typedef struct{
	PixImpact *impact;
	PixImpact * crosstalk;
	gradingTimeStruct * times;
	long row;
	int n_active_crosstalk;
	double crosstalk_energy;
	int nb_crosstalk_influence;
}GradeProxy;

typedef struct {
	gradingTimeStruct *times;
	long row;
	double totalenergy;
}pixGrade;

typedef struct{
	int* pixelHits;
	int num_pix;
}imodProxy;

/** given grade1 and grade 2, make a decision about the high/mid/los res events **/
int makeGrading(long grade1,long grade2,AdvPix* pixel);

/** calculate the grading in samples from the a given impact, and its previous and next impact **/
void calcGradingTimes(double sample_length, gradingTimeStruct pnt,long *grade1, long *grade2, int* status);

/** writes the grading to an existing piximpact file **/
void writeGrading2PixImpactFile(AdvDet *det,PixImpFile *piximpacfile,int *status);

/** Process the impacts contained in the piximpacts file with the RMF method */
void processImpactsWithRMF(AdvDet* det,PixImpFile* piximpacfile,TesEventFile* event_file,int* const status);

/** Processes the impacts, including crosstalk and RMF energy randomization **/
void impactsToEvents(AdvDet *det,PixImpFile *piximpactfile,TesEventFile* event_file,int save_crosstalk,int* const status);

/** Apply matrix cross talk: create new events on concerned pixels if corresponding energy is above the detection threshold, affect previous event otherwise */
void applyMatrixCrossTalk(MatrixCrossTalk* cross_talk,GradeProxy* grade_proxys,const double sample_length,
		PixImpact* impact,AdvDet* det,TesEventFile* event_file,int save_crosstalk,int* const status);

/** Same as "applyMatrixCrosstalk", but for the more complicated intermodulation crosstalk */
void applyIntermodCrossTalk(IntermodulationCrossTalk* cross_talk,GradeProxy* grade_proxys,const double sample_length,
		PixImpact* impact, AdvDet* det,TesEventFile* event_file,int save_crosstalk,int* const status);

/** Processes a crosstalk event using addCrosstalkEvent or processGradedEvent depending on whether it is above threshold or not */
void processCrosstalkEvent(GradeProxy* grade_proxy,const double sample_length,PixImpact* impact,AdvDet* det,TesEventFile* event_file,int save_crosstalk,int* const status);

/** Processes a graded event : update grading proxy and save previous event */
void processGradedEvent(GradeProxy* grade_proxy,const double sample_length,PixImpact* next_impact,AdvDet* det,TesEventFile* event_file,int is_crosstalk,int* const status);

/** Adds a below threshold crosstalk event to the grading proxy */
void addCrosstalkEvent(GradeProxy* grade_proxy,const double sample_length,PixImpact* impact,AdvDet* det,TesEventFile* event_file,int save_crosstalk,int* const status);

#endif /* TESWITHRMF_H */
