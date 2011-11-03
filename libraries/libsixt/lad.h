#ifndef LAD_H
#define LAD_H 1

#include "sixt.h"
#include "arf.h"
#include "point.h"
#include "rmf.h"
#include "xmlbuffer.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Element of the LAD detector for LOFT. */
typedef struct {
  
  /** Unique identifier within the parent module. */
  long id;

  /** Width of the element in x-direction (along the drift as in
      Campana, 2010) [m]. */
  float xdim;

  /** Width of the element in y-direction (parallel to the anodes as
      in Campana, 2010) [m]. */
  float ydim;

  /** Border of the element at the top and at the bottom of the anode
      lines [m]. */
  float xborder;

  /** Border of the element at both sides of the drift strips [m]. */
  float yborder;

  /** Number of anodes. Must be an even number, since each element
      contains the two rows of anodes at the top and at the bottom
      edge. */
  long nanodes;

} LADElement;


/** Module of the LAD detector for LOFT. */
typedef struct {

  /** Unique identifier within the parent panel. */
  long id;

  /** Number of child elements in x-direction. */
  long nx;

  /** Number of child elements in y-direction. */
  long ny;

  /** Number of elements. */
  long nelements;

  /** Width of the module in x-direction [m]. */
  float xdim;

  /** Width of the module in y-direction [m]. */
  float ydim;
  
  /** Array containing pointers to the individual elements. */
  LADElement** element;

} LADModule;


/** Panel of the LAD detector for LOFT. */
typedef struct {

  /** Unique identifier within the LAD. */
  long id;

  /** Number of child modules in x-direction. */
  long nx;

  /** Number of child modules in y-direction. */
  long ny;

  /** Number of modules. */
  long nmodules;
  
  /** Width of the panel in x-direction [m]. */
  float xdim;

  /** Width of the panel in y-direction [m]. */
  float ydim;

  /** Array containing pointers to the individual modules. */
  LADModule** module;

} LADPanel;


/** LAD detector for LOFT. */
typedef struct {

  /** Number of LAD panels. */
  long npanels;

  /** Array containing pointers to the individual panels. */
  LADPanel** panel;

  /** Diameter of the FoV defined by the collimator [rad]. In the XML
      file the diameter is given in [deg], but it is converted to
      [rad] when parsing the XML file. */
  float fov_diameter;

  /** ARF containing the on-axis effective area of the instrument. */
  char* arf_filename;
  struct ARF* arf;

  /** RMF of the detector. */
  char* rmf_filename;
  struct RMF* rmf;

  /** Temperature of the SDDs [K]. */
  float temperature;
  
  /** Strength of the electric drift field [V]. */
  float efield;

  /** Mobility [m**2/V/s]. */
  float mobility; 

  /** Lower read-out threshold [keV]. */
  float* threshold_readout_lo_keV;

  /** Upper read-out threshold [keV]. */
  float* threshold_readout_up_keV;

  /** File name (without path contributions) of the FITS file
      containing the XML LAD definition. */
  char* filename;

  /** Path to the FITS file containing the XML detector definition. */
  char* filepath;

} LAD;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Load the LAD detector definition from the XML file and create a
    new LAD data structure based on these data. */
LAD* getLADfromXML(const char* const filename, 
		   int* const status);

/** Constructor. Allocates memory for a new LAD data structure. */
LAD* newLAD(int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the LAD data structure to NULL. */
void freeLAD(LAD** const lad);


/** Constructor. Allocates memory for a new LADPanel data structure. */
LADPanel* newLADPanel(int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the LADPanel data structure to NULL. */
void freeLADPanel(LADPanel** const panel);


/** Constructor. Allocates memory for a new LADModule data structure. */
LADModule* newLADModule(int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the LADModule data structure to NULL. */
void freeLADModule(LADModule** const module);


/** Constructor. Allocates memory for a new LADElement data structure. */
LADElement* newLADElement(int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the LADElement data structure to NULL. */
void freeLADElement(LADElement** const element);

/** Determine the column and row indices of the hole, in which the
    specified position on the collimator lies. If the position is not
    inside a hole but on the absorbing material, the return values of
    the column and row indices are set to -1. */
void LADCollimatorHoleIdx(const struct Point2d position,
			  long* col, long* row);


#endif
