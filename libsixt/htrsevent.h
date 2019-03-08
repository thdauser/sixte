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


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef HTRSEVENT_H
#define HTRSEVENT_H 1


/** HTRS-specific event. */
typedef struct {
  /** Detection time. */
  double time;

  /** PHA channel. */
  long pha;

  /** Photon energy [keV]. */
  float energy;

  /** Number of the pixel, where the event is detected. Numbering
      starts at 0. */
  int pixel;

  /** Event grade. This gives information about the energy and time
      resolution of the event. There are the following event grades:

      -  0: event is measured with the slow shaper

      -  1: event is measured with the fast shaper

      -  2: event is mixed up with another event (pile-up) because they
            cannot be distinguished even with the fast shaper

      - -1: lost during reset

  */
  int grade;


  /** Exact impact position of the photon on the detector. */
  double x, y;

} HTRSEvent;


#endif /* HTRSEVENT_H */
