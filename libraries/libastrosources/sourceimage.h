#ifndef SOURCEIMAGE_H
#define SOURCEIMAGE_H 1

#include "sixt.h"
#include "sixt_random.h"
#include "spectrum.h"


/** One pixel in the extended SourceImage. */
/*
struct SourceImagePixel {
  float rate; //< Photon rate in this pixel
};
*/

/** Object containing the extended source image from one FITS image extension. */
typedef struct {
  /** Pixel probability distribution. */
  double** pixel;

  /** If this flag is set (1), the SourceImage pixels contain a
      normalized probability distribution. If not (0), the pixels
      contain a real image. */
  int accumulated;

  /** Total photon rate of the entire SourceImage. This is the sum of
      the pixels read from the FITS file. In the computer memory the
      image is normalized to 1. */
  float total_rate;

  /** Time of the last photon emitted from this SourceImage. */
  double t_last_photon;

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


struct SourceImageParameters {
  /** Dimensions of the source image. */
  int naxis1, naxis2;
  /** Width of the image pixels [rad]. */
  double cdelt1, cdelt2;
};


//////////////////////////////////////////////////////////////////


/** Constructor for SourceImage objects returning the bare data
    structure without any memory allocated for the pixels. */
SourceImage* get_SourceImage();

/** Constructor for SourceImage objects generating an empty source
    image with allocated memory for the requested image dimensions. */
SourceImage* getEmptySourceImage(struct SourceImageParameters* sip, int* status);

/** Constructor for SourceImage objects. Reads a SourceImage from a
    FITS file and stores it in the corresponding data structure. */
SourceImage* get_SourceImage_fromFile(char* filename, int* status);

/** Constructor for SourceImage objects. Reads a SourceImage from a
    FITS file and stores it in the corresponding data structure. */
SourceImage* get_SourceImage_fromHDU(fitsfile* fptr, int* status);

/** Destructor for SourceImage objects. */
void free_SourceImage(SourceImage* si);

/** Store the SourceImage in a FITS file. */
void saveSourceImage(SourceImage* si, char* filename, int* status);

/** Constructor for the SourceImageCatalog. */
SourceImageCatalog* get_SourceImageCatalog();

/** Destructor for the SourceImageCatalog. */
void free_SourceImageCatalog(SourceImageCatalog* sic);

/** Add new SourceImage object to the SourceImageCatalog. The new
    image is loaded from the FITS HDU designated by fptr and added to
    the SourceImageCatalog. The return value is the error status. */
int addSourceImage2Catalog(SourceImageCatalog* sic, fitsfile* fptr);

/** Determine a random pixel according to the probability distribution
    given by the SourceImage. */
void getRandomSourceImagePixel(SourceImage* si, int* x, int* y);


#endif /* SOURCEIMAGE_H */

