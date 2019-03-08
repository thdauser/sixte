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

#ifndef GAUSSIANCHARGECLOUD_H
#define GAUSSIANCHARGECLOUD_H 1

typedef struct {

  /** Sigma value for Gaussian shape charge clouds (given in [m]).
      This quantity is used to calculate size/extension of the
      Gaussian shape charge cloud. */
  double ccsigma;
  /** Size of the charge cloud (given in [m]). Quantity to estimate
      the extension of a Gaussian shape charge cloud. The value is
      defined to be three times ccsigma. It is assumed that
      approximately all charge is with in a radius of this size (3
      sigma) around the center of the charge cloud. */
  double ccsize;

} GaussianChargeCloud;

#endif /* GAUSSIANCHARGECLOUD_H */
