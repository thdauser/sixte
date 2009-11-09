#include "sixt.h"
#include "hexagonalpixels.h"
#include "arcpixels.h"
#include "impact.h"
#include "impactlistfile.h"


int main(int argc, char* argv[])
{
  ImpactListFile impactlistfile;
  Impact impact;
  int pixel; // Pixel hit by an incident photon.

#ifdef HTRS_HEXPIXELS
  const int n_pixels=37;
  int count_pixelwidth;                  // Counter for the different pixel widths.
  const double min_pixelwidth=0.5e-3;    // Minimum pixel width ([m]).
  const double step_pixelwidth=0.125e-3; // Increment step for the pixel width ([m]).
  const int n_pixelwidths=37;            // Number of different pixelwidths to simulate.
  HexagonalPixels hexagonalPixels[n_pixelwidths];
  long hex_nphotons[n_pixelwidths][n_pixels];  // Number of photons per hexagonal pixel.
#endif
#ifdef HTRS_ARCPIXELS
  const int n_pixels=37;
  ArcPixels arcPixels;
  long arc_nphotons[n_pixels]; // Number of photons per arc pixel.
  int ring;
#endif

  long ntotal_photons;   // Total number of photons.

  int status = EXIT_SUCCESS;

  if (argc < 2) {
    printf("Error: You have to specify an impact list!\n");
    return(EXIT_FAILURE);
  }

  do { // Error handling loop.

#ifdef HTRS_HEXPIXELS
    // Loop over the different pixel widths.
    for (count_pixelwidth=0; count_pixelwidth<n_pixelwidths; count_pixelwidth++) {
      // Init buffers for statistical data.
      for (pixel=0; pixel<n_pixels; pixel++) {
	hex_nphotons[count_pixelwidth][pixel] = 0;
      }

      // Init data structures for different pixel widths.
      struct HexagonalPixelsParameters hpparameters = {
	.npixels = n_pixels,
	.pixelwidth = min_pixelwidth + count_pixelwidth*step_pixelwidth  // [m]
      };

      // Initialization of HexagonalPixels data structure.
      status=initHexagonalPixels(&hexagonalPixels[count_pixelwidth], 
				 &hpparameters);
      if (EXIT_SUCCESS!=status) break;
    } // END of loop over different pixel widths.
    if (EXIT_SUCCESS!=status) break;
#endif

#ifdef HTRS_ARCPIXELS
    // Setup the configuration of the ArcPixels structure.
    for (pixel=0; pixel<n_pixels; pixel++) {
      arc_nphotons[pixel] = 0;
    }

    // Configuration with 31 pixels.
    int npixels[4] = { 1, 6, 12, 12 };
    double radii[4] = { 2.26e-3, 5.5e-3, 8.85e-3, 12.0e-3 };
    double offset_angles[4] = { 0., 0., M_PI/12, 0. };
    /*// Configuration with 37 pixels.
    int npixels[4] = { 1, 7, 14, 15 };
    double radii[4] = { 2.3e-3, 5.44e-3, 8.77e-3, 12.0e-3 };
    double offset_angles[4] = { 0., 0., 0., 0. };*/
    /*// Configuration with 50 pixels.
    int npixels[4] = { 1, 9, 25, 25 };
    double radii[4] = { 2.1e-3, 5.47e-3, 8.79e-3, 12.0e-3 };
    double offset_angles[4] = { 0., 0., 0., 0. };*/
    struct ArcPixelsParameters apparameters = {
      .nrings = 4, 
      .npixels = npixels,
      .radius = radii,
      .offset_angle = offset_angles
    };
    // Initialization of ArcPixels data structure.
    if(EXIT_SUCCESS!=(status=initArcPixels(&arcPixels, &apparameters))) break;
    // END setup of ArcPixels.
#endif
    

    // Open the impact list FITS file.
    status = openImpactListFile(&impactlistfile, argv[1], READONLY);
    if (EXIT_SUCCESS!=status) break;

    ntotal_photons=impactlistfile.nrows;

    // Loop over all impacts in the FITS file.
    int pixel;
    while ((EXIT_SUCCESS==status)&&(0==ImpactListFile_EOF(&impactlistfile))) {

      status=getNextImpactListFileRow(&impactlistfile, &impact);
      if (EXIT_SUCCESS!=status) break;

#ifdef HTRS_HEXPIXELS
      // Loop over the different pixel widths.
      for (count_pixelwidth=0; count_pixelwidth<n_pixelwidths; count_pixelwidth++) {
	// Determine the hexagonal pixel that is hit by the photon.
	getHexagonalPixel(&hexagonalPixels[count_pixelwidth], impact.position, &pixel);

	// Record data for statistics.
	assert(pixel<n_pixels);
	if (INVALID_PIXEL!=pixel) {
	  hex_nphotons[count_pixelwidth][pixel]++;
	}
      } // END of loop over different pixel widths.
#endif

#ifdef HTRS_ARCPIXELS
      // Determine the ArcPixel that is hit by the photon.
      getArcPixel(&arcPixels, impact.position, &ring, &pixel);
      // Determine the absolute pixel index (numbering for all pixels of this
      // detector) from the given ring and internal pixel number in this ring.
      arc_nphotons[getArcPixelIndex(&arcPixels, ring, pixel)]++;
#endif
    } // END of scanning the impact list.


    // Calculate and print statistics.
    long ndetected=0;
    double mean=0.;
    double mean2=0.;

#ifdef HTRS_HEXPIXELS
    printf("# Hexagonal Pixel Structure (different pixel sizes):\n");
    // Loop over the different pixel sizes.
    for (count_pixelwidth=0; count_pixelwidth<n_pixelwidths; count_pixelwidth++) {
      // Determine the statistics and print them.
      ndetected=0;
      mean2=0.;
      for (pixel=0; pixel<n_pixels; pixel++) {
	ndetected += hex_nphotons[count_pixelwidth][pixel];
	mean2 += pow((double)hex_nphotons[count_pixelwidth][pixel], 2.)/n_pixels;
      }
      mean = ndetected*1./n_pixels;

      printf("# pixel width, N_photons, N_detected_photons, "
	     "fraction of detected photons, sigma/average\n");
      printf("%lf %ld %ld %.10lf %lf\t pixel:",
	     min_pixelwidth+count_pixelwidth*step_pixelwidth, // Pixel width
	     ntotal_photons,               // Total number of photons
	     ndetected,                    // Number of detected photons
	     ndetected*1./ntotal_photons,  // Fraction of detected photons
	     sqrt(mean2-pow(mean,2.))/mean // rms/mean
	     );
      for (pixel=0; pixel<n_pixels; pixel++) {
	printf(" %lf,", (double)hex_nphotons[count_pixelwidth][pixel]/ntotal_photons);
      }
      printf("\n");
    } // END of loop over different pixel sizes.
#endif

#ifdef HTRS_ARCPIXELS
    printf("# ArcPixel Structure:\n");
    // Determine the statistics and print them.
    ndetected=0;
    mean2=0.;
    for (pixel=0; pixel<n_pixels; pixel++) {
      ndetected += arc_nphotons[pixel];
      mean2 += pow((double)arc_nphotons[pixel], 2.)/n_pixels;
    }
    mean = ndetected*1./n_pixels;
      
    printf("# N_photons, N_detected_photons, "
	   "fraction of detected photons, sigma/average\n");
    printf("%ld %ld %.10lf %lf\n",
	   ntotal_photons,               // Total number of photons
	   ndetected,                    // Number of detected photons
	   ndetected*1./ntotal_photons,  // Fraction of detected photons
	   sqrt(mean2-pow(mean,2.))/mean // rms/mean
	   );
    printf("# radii: %lf %lf %lf %lf\n# pixel:", radii[0], radii[1], radii[2], radii[3]);
    printf("\n %ld, ", ntotal_photons); // RM
    for (pixel=0; pixel<n_pixels; pixel++) {
      //      printf(" %lf,", (double)arc_nphotons[pixel]/ntotal_photons);
      printf(" %ld,", arc_nphotons[pixel]); // RM
    }
    printf("\n");
#endif

  } while (0); // END of error handling loop.


  // --- Clean up ---

  // Close the impact list file.
  status += closeImpactListFile(&impactlistfile);

#ifdef HTRS_HEXPIXELS
  // Loop over the different pixel widths.
  for (count_pixelwidth=0; count_pixelwidth<n_pixelwidths; count_pixelwidth++) {
    // Release memory.
    cleanupHexagonalPixels(&hexagonalPixels[count_pixelwidth]);
  } // END of loop over different pixel widths.
#endif

#ifdef HTRS_ARCPIXELS
  cleanupArcPixels(&arcPixels);
#endif

  return(status);
}

