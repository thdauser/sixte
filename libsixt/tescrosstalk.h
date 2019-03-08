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


   Copyright 2016 Christian Kirsch, FAU
   Copyright 2017-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef TESCROSSTALK_H
#define TESCROSSTALK_H 1

#include "sixt.h"
#include "advdet.h"
#include "crosstalk.h"
#include <gsl/gsl_complex_math.h>

////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
// Function Declarations.
////////////////////////////////////////////////////////////////////////

/** Asign resonance frequencies to advdet pixels */
void get_resonance_frequencies(AdvDet* det, int* status);

/** Alloc an empty FDM system */
FDMSystem* newFDMSystem(int num_pixels, int* status);

/** Initialize this channel's FDM system */
void init_FDMSystem(Channel* chan, double L_Common, double C_Common, double TTR, int* status);

/** Solve this channel's FDM system for the current channel state and output to the pixels */
void solve_FDM_old(Channel *chan);
void solve_FDM(Channel *chan);



#endif /* TESCROSSTALK_H */
