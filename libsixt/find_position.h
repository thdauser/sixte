#ifndef FIND_POSITION_H
#define FIND_POSITION_H 1

#include "sixt.h"
#include "sourceimage.h"
#include "squarepixels.h"

////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////
typedef struct {
  int posX, posY;           //position of source in pixel-coordinates
  double posRA, posDEC;     //position of source in ra,dec
  double pixval;            //value of source-pixel
  /* double errorRA, errorDEC;*/ //error of found ra,dec compared to originally given value (source.simput)
} PixPosition;

typedef struct {
  PixPosition** entry; //pointer to current entry of PixPosition
  int entryCount;     //count for all PixPosition-elements (all sources)
} PixPositionList;

/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////
PixPosition* getPixPosition(); //returns PixPosition structure

PixPositionList* getPixPositionList(); //returns PixPositionList structure

//returns value of current brightest source
double findBrightestPix(int threshold, SourceImage* sky_pixels, double pixval, PixPositionList* ppl,
			float delta, double telRA, double telDEC);

//returns RA/DEC of current brightest source (identified as one pixel)  TODO:more than one pix?!
double getPosRa(int pix, SourceImage* sky_pixels, float delta, double telRA);
double getPosDec(int pix, SourceImage* sky_pixels, float delta, double telDEC);

//returns '1' for pix that is still 20% (arbitrary) above mean of all left pix, '2' else
int getThresholdForSources(double pix, PixPositionList* ppl, SourceImage* sky_pixels);

//returns mean value of all left pix (all except already identified sources)
double getMeanValue(PixPositionList* ppl, SourceImage* sky_pixels);

//saves PixPositionList to a FITS-file
void savePositionList(PixPositionList* ppl, char* filename, int* status);

void FreePixPositionList(PixPositionList* ppl);

#endif /* FIND_POSITION_H */

