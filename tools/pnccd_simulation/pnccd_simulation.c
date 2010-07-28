#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "pnccd_simulation.h"

////////////////////////////
/** Main procedure */
int pnccd_simulation_main() {

	struct Parameters parameters; // containing all program parameters read by PIL

	// Detector data structure (containing the pixel array, width, ...)

	// Background structure

	ImpactListFile impactlistfile;

	int status=EXIT_SUCCESS; // Error status

	// Register HEATOOL:
	set_toolname("pnccd_simulation");
	set_toolversion("0.01");


	return(status);
}

int getpar(struct Parameters* parameters) {

	int status=EXIT_SUCCESS; // error status
