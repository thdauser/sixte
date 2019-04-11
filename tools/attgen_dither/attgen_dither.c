/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "sixt.h"
#include "attitude.h"
#include "simput.h"

#define TOOLSUB attgen_dither_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  /** Attitude file. */
  char* fname_attitude;

  double amplitude;
  double ra, dec;


  double tstart;
  double expos;
  int nbins;  // dt will be calculated

  double mjdref;

  int clobber;
};


int attgen_dither_getpar(struct Parameters *parameters);


int attgen_dither_main() {
  // Program parameters.
  struct Parameters par;


  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("attgen_dither");
  set_toolversion("0.01");

  Attitude* ac = NULL;

  do { // Beginning of the ERROR handling loop.

	  // --- Initialization ---

	  // Read the program parameters using PIL library.
	  if ((status=attgen_dither_getpar(&par))) break;

	  ac = get_default_attitude_lissajous(par.amplitude, par.ra, par.dec,
	  		par.tstart, par.tstart+par.expos, par.mjdref, &status);

	  assert(ac!=NULL);

	  // now we need to write the attitude file
	  write_attitude(ac, par.fname_attitude, par.mjdref, par.clobber, &status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  freeAttitude(&ac);

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int attgen_dither_getpar(struct Parameters *par) {

  int status=EXIT_SUCCESS; // Error status

  query_simput_parameter_file_name("Attitude", &par->fname_attitude, &status);

  query_simput_parameter_double("Amplitude", &par->amplitude, &status);

  // only load srcRA and srcDec if Simput is not given
  query_simput_parameter_double("SrcRA",&par->ra,&status);
  query_simput_parameter_double("SrcDec",&par->dec,&status);

  query_simput_parameter_double("Exposure", &par->expos, &status);
  query_simput_parameter_double("TSTART", &par->tstart, &status);
  query_simput_parameter_int("nbins", &par->nbins, &status);
  query_simput_parameter_bool("clobber", &par->clobber, &status);

  return(status);
}
