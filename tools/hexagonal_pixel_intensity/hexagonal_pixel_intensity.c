#include "sixt.h"
#include "htrsdetector.h"
#include "impact.h"
#include "impactlistfile.h"


int main(int argc, char* argv[])
{
  HexagonalPixels hexagonalPixels;
  ImpactListFile impactlistfile;
  Impact impact;
  int pixel; // Pixel hit by an incident photon.

  long ntotal_photons; // Total number of photons.
  long nphotons[37];   // Number of photons per pixel.

  int status = EXIT_SUCCESS;


  if (argc < 2) {
    printf("Error: have to specify impact list!\n");
    return(EXIT_FAILURE);
  }

  do { // Error handling loop 

    // Initialization of HexagonalPixels data structure.
    struct HexagonalPixelsParameters hpparameters = {
      .npixels = 37,
      .pixelwidth = 4.e-3  // in [m]
    };
    if(EXIT_SUCCESS!=(status=initHexagonalPixels(&hexagonalPixels, &hpparameters))) break;

    // Open the impact list FITS file.
    status = openImpactListFile(&impactlistfile, argv[1], READONLY);
    if (EXIT_SUCCESS!=status) break;

    ntotal_photons=impactlistfile.nrows;

    // Loop over all impacts in the FITS file.
    while ((EXIT_SUCCESS==status)&&(0==ImpactListFile_EOF(&impactlistfile))) {

      status=getNextImpactListFileRow(&impactlistfile, &impact);
      if (EXIT_SUCCESS!=status) break;

      // Determine the pixel that is hit by the photon.
      getHexagonalPixel(&hexagonalPixels, impact.position, &pixel);

      // Record data for statistics.
      assert(pixel<37);
      if (INVALID_PIXEL!=pixel) {
	nphotons[pixel]++;
      }

    } // END of scanning the impact list.

    // Release memory.
    cleanupHexagonalPixels(&hexagonalPixels);


    // Determine the statistics and print them.
    long ndetected=0;
    double mean=0., mean2;
    for (pixel=0; pixel<37; pixel++) {
      ndetected += nphotons[pixel];
      mean2 += pow((double)nphotons[pixel], 2.)/37.;
    }
    mean = ndetected*1./37.;

    printf("%lf %ld %lf %lf\n",
	   hpparameters.pixelwidth, ntotal_photons,
	   (double)ndetected/ntotal_photons, // fraction of detected photons
	   sqrt(mean2-pow(mean,2.)) // rms
	   );

  } while (0); // END of error handling loop


  // --- Clean up ---

  return(status);
}
