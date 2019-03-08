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

#ifndef KDTREEELEMENT_H
#define KDTREEELEMENT_H 1

#include "sixt.h"
#include "check_fov.h"
#include "simput.h"
#include "source.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Element of a KDTree (multidimensional binary tree). */
struct structKDTreeElement {
  Source* src;
  struct structKDTreeElement* left;
  struct structKDTreeElement* right;
};
typedef struct structKDTreeElement KDTreeElement;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
KDTreeElement* newKDTreeElement(int* const status);

/** Destructor. */
void freeKDTreeElement(KDTreeElement** el);

/** Build up the KDTree from the given list of Sources. */
KDTreeElement* buildKDTree2(Source* const list,
			    const long nelements,
			    const int depth,
			    int* const status);

/** Perform a range search on the given kdTree, i.e., return all X-ray
    sources lying within a certain radius around the reference
    point. This region is defined by the minimum cosine value for the
    scalar product of the source direction and the reference
    vector. The function returns a time-ordered list of newly
    generated photons. */
LinkedPhoListElement* KDTreeRangeSearch(KDTreeElement* const node,
					const int depth,
					const Vector* const ref,
					const double min_align,
					const double t0, const double t1,
					const double mjdref,
					SimputCtlg* const simputcat,
					int* const status);


#endif /* KDTREEELEMENT_H */
