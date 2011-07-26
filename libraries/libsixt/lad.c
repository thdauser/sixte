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

  return(element);
}


void freeLADElement(LADElement** const element)
{
  if (NULL!=*element) {
    free(*element);
    *element=NULL;
  }
}

