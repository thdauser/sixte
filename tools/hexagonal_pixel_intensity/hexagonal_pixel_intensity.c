#include "sixt.h"
#include "htrsdetector.h"
#include "impact.h"
#include "impactlistfile.h"


int main(int argc, char* argv[])
{
  ImpactListFile impactlistfile;
  Impact impact;
  int pixel; // Pixel hit by an incident photon.

  int count_pixelwidth;                // Counter for the different pixel widths.
  const int n_pixelwidths=2;          // Number of different pixelwidths to simulate.
  const double min_pixelwidth=0.5e-3;  // Minimum pixel width ([m]).
  const double step_pixelwidth=0.5e-3; // Increment step for the pixel width ([m]).

  HexagonalPixels hexagonalPixels[n_pixelwidths];

  long ntotal_photons;                // Total number of photons.
  long nphotons[n_pixelwidths][37];   // Number of photons per pixel.

  int status = EXIT_SUCCESS;


  if (argc < 2) {
    printf("Error: You have to specify an impact list!\n");
    return(EXIT_FAILURE);
  }

  do { // Error handling loop 

    // Loop over the different pixel widths.
    for (count_pixelwidth=0; count_pixelwidth<n_pixelwidths; count_pixelwidth++) {
      // Init buffers for statistical data.
      for (pixel=0; pixel<37; pixel++) {
	nphotons[count_pixelwidth][pixel] = 0;
      }

      // Init data structures for different pixel widths.
      struct HexagonalPixelsParameters hpparameters = {
	.npixels = 37,
	.pixelwidth = min_pixelwidth + count_pixelwidth*step_pixelwidth  // [m]
      };

      // Initialization of HexagonalPixels data structure.
      if(EXIT_SUCCESS!=(status=initHexagonalPixels(&hexagonalPixels[count_pixelwidth], 
						   &hpparameters))) break;
    } // END of loop over different pixel widths.
    
    
    // Open the impact list FITS file.
    status = openImpactListFile(&impactlistfile, argv[1], READONLY);
    if (EXIT_SUCCESS!=status) break;

    ntotal_photons=impactlistfile.nrows;

    // Loop over all impacts in the FITS file.
    while ((EXIT_SUCCESS==status)&&(0==ImpactListFile_EOF(&impactlistfile))) {

      status=getNextImpactListFileRow(&impactlistfile, &impact);
      if (EXIT_SUCCESS!=status) break;

      // Loop over the different pixel widths.
      for (count_pixelwidth=0; count_pixelwidth<n_pixelwidths; count_pixelwidth++) {
	// Determine the pixel that is hit by the photon.
	getHexagonalPixel(&hexagonalPixels[count_pixelwidth], impact.position, &pixel);

	// Record data for statistics.
	assert(pixel<37);
	if (INVALID_PIXEL!=pixel) {
	  nphotons[count_pixelwidth][pixel]++;
	}
      } // END of loop over different pixel widths.

    } // END of scanning the impact list.


    // Loop over the different pixel widths.
    for (count_pixelwidth=0; count_pixelwidth<n_pixelwidths; count_pixelwidth++) {
      // Determine the statistics and print them.
      long ndetected=0;
      double mean=0., mean2;
      for (pixel=0; pixel<37; pixel++) {
	ndetected += nphotons[count_pixelwidth][pixel];
	mean2 += pow((double)nphotons[count_pixelwidth][pixel], 2.)/37.;
      }
      mean = ndetected*1./37.;
      
      printf("%lf %ld %ld %lf %lf\n",
	     min_pixelwidth+count_pixelwidth*step_pixelwidth, // Pixel width
	     ntotal_photons,                   // Total number of photons
	     ndetected,                        // Number of detected photons
	     ndetected*1./ntotal_photons,      // Fraction of detected photons
	     sqrt(mean2-pow(mean,2.))/mean     // rms/mean
	     );

    } // END of loop over different pixel widths.

  } while (0); // END of error handling loop.


  // --- Clean up ---

  // Close the impact list file.
  status += closeImpactListFile(&impactlistfile);

  // Loop over the different pixel widths.
  for (count_pixelwidth=0; count_pixelwidth<n_pixelwidths; count_pixelwidth++) {
    // Release memory.
    cleanupHexagonalPixels(&hexagonalPixels[count_pixelwidth]);
  } // END of loop over different pixel widths.

  return(status);
}
