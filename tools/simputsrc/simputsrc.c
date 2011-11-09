#include "simputsrc.h"


int simputsrc_main() 
{
  // Filename constants.
  const char CMDFILE[] = "isis.tmp";
  const char SPECFILE[] = "spec%d.tmp";

  // Program parameters.
  struct Parameters par;

  // Temporary files for ISIS interaction.
  FILE* cmdfile=NULL;
  fitsfile* specfile=NULL;

  // Buffers for spectral components.
  float* flux=NULL;

  // SIMPUT data structures.
  SimputMissionIndepSpec* simputspec=NULL;
  SimputSource* src=NULL;
  SimputCatalog* cat=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("simputsrc");
  set_toolversion("0.03");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----
    
    // Read the parameters using PIL.
    status=simputsrc_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // ---- END of Initialization ----


    // ---- Main Part ----

    // Create a new SIMPUT catalog.
    remove(par.Simput);
    cat=openSimputCatalog(par.Simput, READWRITE, &status);
    CHECK_STATUS_BREAK(status);

    // Insert a point-like source.
    float totalFlux = par.plFlux+par.bbFlux+par.flFlux+par.rflFlux;
    char src_name[MAXMSG];
    strcpy(src_name, par.Src_Name);
    strtoupper(src_name);
    if ((strlen(src_name)>0)&&(strcmp(src_name, "NONE"))) {
      // Copy again the name, this time without upper-case conversion.
      strcpy(src_name, par.Src_Name);
    } else {
      strcpy(src_name, "");
    }
    
    // Get a new source entry.
    src=getSimputSourceV(1, src_name, 0., 0., 0., 1., 
			 par.Emin, par.Emax, totalFlux, 
			 "[SPECTRUM,1]", "", "", &status);
    CHECK_STATUS_BREAK(status);
    appendSimputSource(cat, src, &status);
    CHECK_STATUS_BREAK(status);

    // END of creating the catalog.


    // Create the spectrum.

    // Open the ISIS command file.
    // TODO Get a random temporary name instead of using a constant.
    cmdfile=fopen(CMDFILE,"w+");
    CHECK_NULL_BREAK(cmdfile, status, "opening temporary file failed");

    // Write the header.
    fprintf(cmdfile, "require(\"isisscripts\");\n");
    fprintf(cmdfile, "()=xspec_abund(\"wilm\");\n");
    fprintf(cmdfile, "use_localmodel(\"relline\");\n");
    
    // Define the energy grid.
    fprintf(cmdfile, "variable grid=[0.05:100.0:0.01];\n");
    fprintf(cmdfile, "variable lo=grid[[0:length(grid)-2]];\n");
    fprintf(cmdfile, "variable hi=grid[[1:length(grid)-1]];\n");
    fprintf(cmdfile, "variable flux;\n");
    fprintf(cmdfile, "variable spec;\n");

    // Loop over the different components of the spectral model.
    int ii;
    for (ii=0; ii<4; ii++) {

      // Define the spectral model and set the parameters.
      switch(ii) {
      case 0:
	fprintf(cmdfile, "fit_fun(\"phabs(1)*powerlaw(1)\");\n");
	fprintf(cmdfile, "set_par(\"powerlaw(1).PhoIndex\", %e);\n", 
		par.plPhoIndex);
	break;
      case 1:
	fprintf(cmdfile, "fit_fun(\"phabs(1)*bbody(1)\");\n");
	fprintf(cmdfile, "set_par(\"bbody(1).kT\", %e);\n", par.bbkT);	
	break;
      case 2:
	fprintf(cmdfile, "fit_fun(\"phabs(1)*egauss(1)\");\n");
	fprintf(cmdfile, "set_par(\"egauss(1).center\", 6.4);\n");	
	fprintf(cmdfile, "set_par(\"egauss(1).sigma\", %e);\n", par.flSigma);	
	break;
      case 3:
	fprintf(cmdfile, "fit_fun(\"phabs(1)*relline(1)\");\n");
	fprintf(cmdfile, "set_par(\"relline(1).lineE\", 6.4);\n");	
	fprintf(cmdfile, "set_par(\"relline(1).a\", %f);\n", par.rflSpin);	
	break;
      default:
	status=EXIT_FAILURE;
	break;
      }
      CHECK_STATUS_BREAK(status);

      // Absorption is the same for all spectral components.
      fprintf(cmdfile, "set_par(\"phabs(1).nH\", %e);\n", par.nH);

      // Evaluate the spectral model and store the data in a temporary 
      // FITS file.
      fprintf(cmdfile, "flux=eval_fun_keV(lo, hi)/(hi-lo);\n");
      fprintf(cmdfile, "spec=struct{ENERGY=0.5*(lo+hi), FLUX=flux};\n");
      char command[MAXMSG];
      sprintf(command, 
	      "fits_write_binary_table(\"%s\",\"SPECTRUM\", spec);\n", SPECFILE);
      fprintf(cmdfile, command, ii);
    } 
    CHECK_STATUS_BREAK(status);
    // END of loop over the different spectral components.

    fprintf(cmdfile, "exit;\n");

    // End of writing the ISIS file.
    fclose(cmdfile);
    cmdfile=NULL;

    // Construct the shell command to run ISIS.
    char command[MAXMSG];
    strcpy(command, "/data/system/software/local/bin/isis ");
    strcat(command, CMDFILE);
      
    // Run ISIS.
    status=system(command);
    CHECK_STATUS_BREAK(status);

    // Add the spectra and append the total spectrum to the SIMPUT file.
    simputspec=getSimputMissionIndepSpec(&status);
    CHECK_STATUS_BREAK(status);

    // Loop over the different components of the spectral model.
    for (ii=0; ii<4; ii++) {

      // Read the spectrum.
      char filename[MAXFILENAME];
      sprintf(filename, SPECFILE, ii);
      fits_open_table(&specfile, filename, READONLY, &status);
      CHECK_STATUS_BREAK(status);

      // Read the data.
      long nrows=0;
      int anynull;
      if (0==ii) {
	// Determine the number of rows.
	fits_get_num_rows(specfile, &nrows, &status);
	CHECK_STATUS_BREAK(status);

	// Allocate memory.
	simputspec->nentries=nrows;
	simputspec->energy=(float*)malloc(nrows*sizeof(float));
	CHECK_NULL_BREAK(simputspec->energy, status, "memory allocation failed");
	simputspec->pflux=(float*)malloc(nrows*sizeof(float));
	CHECK_NULL_BREAK(simputspec->energy, status, "memory allocation failed");

	// Read the energy column.
	fits_read_col(specfile, TFLOAT, 1, 1, 1, nrows, 0, simputspec->energy, 
		      &anynull, &status);
	CHECK_STATUS_BREAK(status);
      }

      flux=(float*)malloc(nrows*sizeof(float));
      CHECK_NULL_BREAK(flux, status, "memory allocation failed");

      fits_read_col(specfile, TFLOAT, 2, 1, 1, nrows, 0, flux, 
		    &anynull, &status);
      CHECK_STATUS_BREAK(status);

      fits_close_file(specfile, &status);
      specfile=NULL;
      CHECK_STATUS_BREAK(status);

      // Normalize the spectrum.
      // Determine the required flux in the reference band.
      float shouldflux=0.;
      switch(ii) {
      case 0:
	shouldflux = par.plFlux;
	break;
      case 1:
	shouldflux = par.bbFlux;
	break;
      case 2:
	shouldflux = par.flFlux;
	break;
      case 3:
	shouldflux = par.rflFlux;
	break;
      default:
	status = EXIT_FAILURE;
	break;
      }
      CHECK_STATUS_BREAK(status);

      float factor;
      if (shouldflux==0.) {
	factor=0.;
      } else {
	// Determine the current flux in the reference band.
	// Conversion factor from [keV] -> [erg].
	const float keV2erg = 1.602e-9;
	// Flux in reference energy band.
	float isflux=0.; 
	long jj;
	for (jj=0; jj<nrows; jj++) {
	  float binmin, binmax;
	  if (0==jj) {
	    binmin=simputspec->energy[0];
	  } else {
	    binmin=0.5*(simputspec->energy[jj]+simputspec->energy[jj-1]);
	  }
	  if (nrows-1==jj) {
	    binmax=simputspec->energy[nrows-1];
	  } else {
	    binmax=0.5*(simputspec->energy[jj]+simputspec->energy[jj+1]);
	  }
	  if ((binmax>par.Emin)&&(binmin<par.Emax)) {
	    float min = MAX(binmin, par.Emin);
	    float max = MIN(binmax, par.Emax);
	    assert(max>min);
	    isflux += (max-min)*flux[jj]*simputspec->energy[jj];
	  }
	}
	// Convert units of 'flux' from [keV/s/cm^2] -> [erg/s/cm^2].
	isflux *= keV2erg;
	
	factor = shouldflux/isflux;
      }

      long jj;
      for (jj=0; jj<nrows; jj++) {
	// Normalize.
	flux[jj] *= factor;

	// Add the component to the total spectrum.
	simputspec->pflux[jj] += flux[jj];
      }

    }
    CHECK_STATUS_BREAK(status);
    // END of loop over the different spectral components.

    long jj;
    for (jj=0; jj<simputspec->nentries; jj++) {
      // Check if the flux has a physically reasonable value.
      if ((simputspec->pflux[jj]<0.)||(simputspec->pflux[jj]>1.e12)) {
	SIXT_ERROR("flux out of boundaries");
	status=EXIT_FAILURE;
      }
    }
    saveSimputMissionIndepSpec(simputspec, par.Simput, "SPECTRUM", 1, &status);
    CHECK_STATUS_BREAK(status);
    // END of creating the spectrum.

    // ---- END of Main Part ----

  } while(0); // END of error handling loop.

  // Close the temporary files.
  if (NULL!=cmdfile) {
    fclose(cmdfile);
    cmdfile=NULL;
  }
  if (NULL!=specfile) {
    fits_close_file(specfile, &status);
    specfile=NULL;
  }
  // Remove the temporary files.
  remove(CMDFILE);
  int ii;
  for (ii=0; ii<4; ii++) {
    char filename[MAXFILENAME];
    sprintf(filename, SPECFILE, ii);
    remove(filename);
  }

  // Release memory.
  freeSimputMissionIndepSpec(&simputspec);
  freeSimputSource(&src);
  freeSimputCatalog(&cat, &status);

  if (NULL!=flux) {
    free(flux);
    flux=NULL;
  }

  if (status==EXIT_SUCCESS) headas_chat(3, "finished successfully!\n\n");
  return(status);
}



int simputsrc_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_float("plPhoIndex", &par->plPhoIndex);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the plPhoIndex parameter failed");
    return(status);
  }

  status=ape_trad_query_float("plFlux", &par->plFlux);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the plFlux parameter failed");
    return(status);
  }

  status=ape_trad_query_float("bbkT", &par->bbkT);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the bbkT parameter failed");
    return(status);
  }

  status=ape_trad_query_float("bbFlux", &par->bbFlux);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the bbFlux parameter failed");
    return(status);
  }

  status=ape_trad_query_float("flSigma", &par->flSigma);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the flSigma parameter failed");
    return(status);
  }

  status=ape_trad_query_float("flFlux", &par->flFlux);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the flFlux parameter failed");
    return(status);
  }

  status=ape_trad_query_float("rflSpin", &par->rflSpin);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the rflSpin parameter failed");
    return(status);
  }

  status=ape_trad_query_float("rflFlux", &par->rflFlux);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the rflFlux parameter failed");
    return(status);
  }

  status=ape_trad_query_float("nH", &par->nH);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the nH parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Emin", &par->Emin);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the Emin parameter failed");
    return(status);
  }

  status=ape_trad_query_float("Emax", &par->Emax);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the Emax parameter failed");
    return(status);
  }

  status=ape_trad_query_file_name("Simput", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("reading the name of the output SIMPUT catalog file failed");
    return(status);
  } 
  strcpy(par->Simput, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }

  return(status);
}


