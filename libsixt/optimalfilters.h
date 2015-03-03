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


#ifndef OPTIMALFILTERS_H
#define OPTIMALFILTERS_H 1

#include "sixt.h"
#include "tesrecord.h"
#include "tesproftemplates.h"
#include "teseventlist.h"

////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

typedef struct {
  /** Duration of the optimal filter. */
  int filter_duration;

  /** Array containing the filter data. */
  double* filter;

} OptimalFilter;

typedef struct {
  /** Number of optimal filters in the structure. */
  int nfilters;

  /** Oth order coeff in pulse height to keV conversion */
  double ph_b;

  /** 1th order coeff in pulse height to keV conversion */
  double ph_a;

  /** Array containing all the optimal filters. */
  OptimalFilter* optimal_filters;

} OptimalFilterCollection;

/** Structure containing all the pointers and values to run the reconstruction stage */
typedef struct {
	/** OptimalFilterCollection structure */
	OptimalFilterCollection* opt_filter_collection;

	/** Pulse template */
	double * pulse_template;

	/** Derivated pulse template */
	double * derivated_template;

	/** Pulse template height */
	double pulse_template_height;

	/** Pulse length */
	int pulse_length;

	/** Calibration factor */
	double calfac;

	/** Trigger threshold */
	double threshold;

	/** Minimal distance before using OFs after a mireconstruction */
	int normal_exclusion;

	/** Minimal distance before reconstructing any event after a mireconstruction */
	int derivate_exclusion;

	/** Saturation level of the ADC curves */
	double saturation_value;


} ReconstructInit;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** OptimalFilterCollection constructor. Returns a pointer to an empty OptimalFilterCollection data
    structure. */
OptimalFilterCollection* newOptimalFilterCollection(int* const status);

/** OptimalFilterCollection Destructor. */
void freeOptimalFilterCollection(OptimalFilterCollection** const opt_filter_collection);

/** Allocates memory for an OptimalFilterCollection structure. */
void allocateOptimalFilterCollection(OptimalFilterCollection* opt_filter_collection,int nfilters,int* const status);


/** OptimalFilter constructor. Returns a pointer to an empty OptimalFilter data structure. */
OptimalFilter* newOptimalFilter(int* const status);

/** Allocates memory for an OptimalFilter structure. */
void allocateOptimalFilter(OptimalFilter* opt_filter,int filter_duration,int* const status);

/** OptimalFilter destructor. */
void freeOptimalFilter(OptimalFilter* opt_filter);



/** Create and retireve an OptimalFilterCollection from a file. */
OptimalFilterCollection*  getOptimalFilterCollection(const char* const filename, int nfilters, int* const status);


/** Derivate a stream */
void derivate_stream(double * data_stream,double * derivated_stream,int stream_length);

/** Remove a pulse from a data stream */
void subtractPulse(double * data_stream,int pulse_time,double * pulse_template,double factor,int pulse_length,int stream_length);

/** Filter the pulse and return the energy */
double filterPulse(double * data_stream,int pulse_time,double * filter,int filter_length);

/** Trigger on the pulses and updates the event_list accordingly */
int triggerEvents(TesRecord* record,TesEventList* event_list,double * derivated_pulse,int derivated_pulse_length,
		double threshold,double pulse_template_height,double saturation_value,int derivate_exclusion,
		int normal_exclusion,int* const status);

/** Computes the energy of the detected pulses and save the result in the event list */
void computeEnergy(TesRecord* record,TesEventList* event_list,OptimalFilterCollection* opt_filter_collection,double * pulse_template,
		int pulse_length,double calfac,const char identify,int* const status);

/** Wrapper around the whole pulse reconstruction */
void reconstructRecord(TesRecord* record,TesEventList* event_list,ReconstructInit* reconstruct_init,const char identify,int* const status);


/** Constructor. Returns a pointer to an empty ReconstructInit data
    structure. */
ReconstructInit* newReconstructInit(int* const status);

/** Destructor. */
void freeReconstructInit(ReconstructInit* reconstruct_init);

/** Initializes the different variables necessary for the reconstruction */
void initializeReconstruction(ReconstructInit* reconstruct_init,char* const optimal_filter_file,int pulse_length,
		char* const pulse_template_file,double threshold,double calfac,int normal_exclusion,int derivate_exclusion,
		double saturation_value,int* const status);

#endif /* OPTIMALFILTERS_H */
