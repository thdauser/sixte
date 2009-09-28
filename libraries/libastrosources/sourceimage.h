#ifndef SOURCEIMAGE_H
#define SOURCEIMAGE_H 1

#include "sixt.h"
#include "spectrum.h"


/** One pixel in the extended SourceImage. */
struct SourceImagePixel {
  float rate; /**< Photon rate in this pixel */
  double t_next_photon; /**< Time of next photon emitted from this Pixel. */
};


/** Object containing the extended source image from one FITS image extension. */
typedef struct {
  struct SourceImagePixel** pixel;

  int naxis1, naxis2;    /**< Width of the image [pixel]. */
  double cdelt1, cdelt2; /**< Width of one pixel [rad]. */
  double crpix1, crpix2; /**< [pixel] */
  double crval1, crval2; /**< [rad] */

  double minra, maxra;   /**< Maximum right ascension covered by the image [rad]. */
  double mindec, maxdec; /**< Maximum declination covered by the image [rad]. */

  /** SpectrumStore containing the source spectra used for this image of the sky. */
  SpectrumStore spectrumstore;
} SourceImage;


/** Catalog of all extended SourceImage Objects. */
typedef struct {
  /** Total number of ClusterImage elements in the catalog. */
  int nimages;
  /** Catalog of ClusterImage elements. */
  SourceImage** images; /* nimages */
} SourceImageCatalog;


//////////////////////////////////////////////////////////////////


/** Constructor for SourceImage objects. */
SourceImage* get_SourceImage();

/** Constructor for SourceImage objects. 
 * Reads a SourceImage from a FITS file and stores it in the 
 * corresponding data structure. */
SourceImage* get_SourceImage_fromFile(char* filename, int* status);

/** Destructor for SourceImage objects. */
void free_SourceImage(SourceImage* si);

/** Constructor for the SourceImageCatalog. */
SourceImageCatalog* get_SourceImageCatalog();

/** Destructor for the SourceImageCatalog. */
void free_SourceImageCatalog(SourceImageCatalog* sic);


#endif /* SOURCEIMAGE_H */

