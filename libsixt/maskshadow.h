#ifndef MASKSHADOW_H
#define MASKSHADOW_H 1

#include "sixt.h"
#include "find_position.h"
#include "vector.h"
#include "reconstruction.h"
#include "find_position.h"
#include "sourceimage.h"
#include "fft_array.h"
#include "eventarray.h"

////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////
typedef struct {
  double** map;    /**whole map of re-pixeled mask derived from ReconArray (values betw 0...1)*/ 
  double** shadow; /**part of re-pixeled mask corresponding to current source */ 
} MaskShadow;

/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////

MaskShadow* getEmptyMaskShadowElement(int* const status);

MaskShadow* getMaskShadowElement(int const Size1, int const Size2, int* const status);

void getMaskRepix(ReconArray* recon, MaskShadow* ms);

void getMaskShadow(MaskShadow* ms, PixPositionList* ppl, SourceImage* sky_pixels, SquarePixels* det_pix,
		    ReconArray* r, const Vector* const nx, const Vector* const ny, double const distance);
void getMaskShadow2(MaskShadow* ms, struct wcsprm* wcs2, PixPositionList* ppl, SourceImage* sky_pixels, SquarePixels* det_pix, ReconArray* r, int* const status);

void testImageShadow(MaskShadow* ms, SquarePixels* det_pix, int* status);

void testImageMap(MaskShadow* ms, ReconArray* r, int* status);

double getNormalization1(MaskShadow* ms, ReadEvent* ea, SquarePixels* det_pix, int const xdiff, int const ydiff);
double getNormalization2(MaskShadow* ms, ReadEvent* ea, SquarePixels* det_pix, int const xdiff, int const ydiff);

void testImageEventArray(ReadEvent* ea, SquarePixels* det_pix, int const xdiff, int const ydiff, char* filename, int* status);

void FreeMaskShadow(MaskShadow* ms,int const Size1);

#endif /* MASKSHADOW_H */
