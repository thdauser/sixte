#include "lad.h"


////////////////////////////////////////////////////////////////////
// Program Code.
////////////////////////////////////////////////////////////////////


LAD* newLAD(int* const status) 
{
  // Allocate memory.
  LAD* lad=(LAD*)malloc(sizeof(LAD));
  if (NULL==lad) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for LAD failed");
    return(lad);
  }

  // Initialize all pointers with NULL and set initial values.
  lad->panel  =NULL;
  lad->npanels=0;
  lad->fov_diameter=0.;
  lad->arf         =NULL;
  lad->arf_filename=NULL;
  lad->rmf         =NULL;
  lad->rmf_filename=NULL;
  lad->temperature =0.;
  lad->efield      =0.;
  lad->mobility    =0.;
  lad->deadtime    =0.;
  lad->threshold_readout_lo_keV=NULL;
  lad->threshold_readout_up_keV=NULL;
  lad->filename    =NULL;
  lad->filepath    =NULL;

  return(lad);
}


void freeLAD(LAD** const lad)
{
  if (NULL!=*lad) {
    if (NULL!=(*lad)->panel) {
      long ii;
      for (ii=0; ii<(*lad)->npanels; ii++) {
	freeLADPanel(&((*lad)->panel[ii]));
      }
      free((*lad)->panel);
    }
    if (NULL!=(*lad)->arf) {
      freeARF((*lad)->arf);
    }
    if (NULL!=(*lad)->arf_filename) {
      free((*lad)->arf_filename);
    }
    if (NULL!=(*lad)->rmf) {
      freeRMF((*lad)->rmf);
    }
    if (NULL!=(*lad)->rmf_filename) {
      free((*lad)->rmf_filename);
    }
    if (NULL!=(*lad)->threshold_readout_lo_keV) {
      free((*lad)->threshold_readout_lo_keV);
    }
    if (NULL!=(*lad)->threshold_readout_up_keV) {
      free((*lad)->threshold_readout_up_keV);
    }
    if (NULL!=(*lad)->filename) {
      free((*lad)->filename);
    }
    if (NULL!=(*lad)->filepath) {
      free((*lad)->filepath);
    }

    free(*lad);
    *lad=NULL;
  }
}


LADPanel* newLADPanel(int* const status) 
{
  // Allocate memory.
  LADPanel* panel=(LADPanel*)malloc(sizeof(LADPanel));
  if (NULL==panel) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for LADPanel failed");
    return(panel);
  }

  // Initialize all pointers with NULL and set initial values.
  panel->module  =NULL;
  panel->nmodules=0;
  panel->id      =0;
  panel->nx      =0;
  panel->ny      =0;
  panel->xdim    =0.;
  panel->ydim    =0.;

  return(panel);
}


void freeLADPanel(LADPanel** const panel)
{
  if (NULL!=*panel) {
    if (NULL!=(*panel)->module) {
      long ii;
      for (ii=0; ii<(*panel)->nmodules; ii++) {
	freeLADModule(&((*panel)->module[ii]));
      }
      free((*panel)->module);
    }

    free(*panel);
    *panel=NULL;
  }
}


LADModule* newLADModule(int* const status) 
{
  // Allocate memory.
  LADModule* module=(LADModule*)malloc(sizeof(LADModule));
  if (NULL==module) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for LADModule failed");
    return(module);
  }

  // Initialize all pointers with NULL and set initial values.
  module->element  =NULL;
  module->nelements=0;
  module->id       =0;
  module->nx       =0;
  module->ny       =0;
  module->xdim     =0.;
  module->ydim     =0.;

  return(module);
}


void freeLADModule(LADModule** const module)
{
  if (NULL!=*module) {
    if (NULL!=(*module)->element) {
      long ii;
      for (ii=0; ii<(*module)->nelements; ii++) {
	freeLADElement(&((*module)->element[ii]));
      }
      free((*module)->element);
    }

    free(*module);
    *module=NULL;
  }
}


LADElement* newLADElement(int* const status) 
{
  // Allocate memory.
  LADElement* element=(LADElement*)malloc(sizeof(LADElement));
  if (NULL==element) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for LADElement failed");
    return(element);
  }

  // Initialize all pointers with NULL and set initial values.
  element->id       =0;
  element->xdim     =0.;
  element->ydim     =0.;
  element->xborder  =0.;
  element->yborder  =0.;
  element->nanodes  =0;

  return(element);
}


void freeLADElement(LADElement** const element)
{
  if (NULL!=*element) {
    free(*element);
    *element=NULL;
  }
}


void LADCollimatorHoleIdx(const struct Point2d position,
			  long* col, long* row)
{
  // Distance between holes in the collimator.
  const double pitch = 28.e-6; // [m]
  // Diameter of holes in the collimator.
  const double diameter = 25.e-6; // [m]
  // Square radius of the holes in the collimator.
  double radius2 = pow(diameter*0.5, 2.);

  // Distance between 2 subsequent rows.
  double h = pitch*0.5/tan(M_PI/6.);

  // Determine the row.
  long rowidx = (long)(position.y/h);
  // Determine the column.
  long colidx = (long)(position.x/(pitch*0.5));

  // Center indices of the circles that have to be taken into account.
  long c[2], r[2];
  
  // Check if row is an even or an odd number.
  if (rowidx % 2 == 1) {
    // Odd (as first row).
    
    // Check if column is an even or an odd number.
    if (colidx % 2 == 1) {
      c[0] = colidx;
      r[0] = rowidx;
      c[1] = colidx+1;
      r[1] = rowidx+1;
    } else {
      c[0] = colidx+1;
      r[0] = rowidx;
      c[1] = colidx;
      r[1] = rowidx+1;
    }
    // END of column is an even or odd number.

  } else {
    // Row is even number (as second row).

    // Check if column is an even or an odd number.
    if (colidx % 2 == 1) {
      c[0] = colidx+1;
      r[0] = rowidx;
      c[1] = colidx;
      r[1] = rowidx+1;
    } else {
      c[0] = colidx;
      r[0] = rowidx;
      c[1] = colidx+1;
      r[1] = rowidx+1;
    }
    // END of column is an even or odd number.

  }
  // END of row is odd or even.

  // Calculate the distances to the 2 circles that might contain 
  // the specified position.
  int ii;
  for (ii=0; ii<2; ii++) {
    struct Point2d center;
    center.x = c[ii]*pitch/2.;
    center.y = r[ii]*h;
    // Check if the specified position lies within the regarded circle (hole).
    if (pow(position.x-center.x, 2.) + pow(position.y-center.y, 2.) < radius2) {
      *col = c[ii];
      *row = r[ii];
      return;  
    }
  }

  // Position is not open.
  *col = -1;
  *row = -1;
  return;
}

