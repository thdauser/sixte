#ifndef ATTITUDECATALOG2_H
#define ATTITUDECATALOG2_H 1

#include "sixt.h"
#include "vector.h"
#include "telescope.h"
#include "attitudefile.h"


//Type Declarations

typedef struct {
  double time; 

  Vector nx;   //telescope motion direction
  Vector nz;   //telescope pointing direction

  float roll_angle; //in [rad]
}AttEntry;

typedef struct {
  long nentries;        //number of Attitude Entry-elements
  AttEntry* entry;      //AttitudeEntry-elements
  long current_entry;   //number of current entry
  int alignment;        //value 0: rollangle determined with respect
                        //         to equatorial plane.
                        //value 1: rollangle with respect to telescope
                        //         motion direction
}AttCatalog;


//Function Declarations

AttCatalog* getAttCatalog(int* const status);

AttEntry initializeAttitudeEntry ();

AttCatalog* loadAttCatalog(const char* filename, int* const status);

Vector GetTelescopeNz (AttCatalog* ac, const double time, 
		       int* const status);

void freeAttCatalog(AttCatalog** const ac);

void getTelAxes(AttCatalog* const ac,
		Vector* const nx,
		Vector* const ny, 
		Vector* const nz, 
		const double time,
		int* const status);

float GetRollAngle(AttCatalog* const ac,
		   const double time,
		   int* const status);


#endif /*ATTITUDECATALOG2_H*/
