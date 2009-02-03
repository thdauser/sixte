#ifndef SCAN_TES_EVENTS_H
#define SCAN_TES_EVENTS_H 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <limits.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "global_constants.h"
#include "event_list.h"


#define TOOLSUB scan_tes_events_main
#include "headas_main.c"

// Number of photons stored in the event list for each individual pixel.
#define TES_EVENT_LIST_ENTRIES 5

int scan_tes_events_getpar(char inputfile[], double* time_unit,
			   int* units_before_pulse, int* units_after_pulse); 


#endif
