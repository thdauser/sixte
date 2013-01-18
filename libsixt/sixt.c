#include "sixt.h"
#include "headas_rand.h"


#ifdef USE_RCL
// Use the Remeis random number server.
#include <rcl.h>
#endif


double sixt_get_random_number(int* const status)
{
  // Return a random value out of the interval [0,1).

#ifdef USE_RCL

  // Use the Remeis random number server running on leo.
  // The library for accessing the server is maintained
  // by Fritz-Walter Schwarm.
  int rcl_status=0;
  double rand=rcl_rand_net_get_double(NULL, NULL, &rcl_status);

  if(RCL_RANDOM_SUCCESS!=rcl_status) {
    SIXT_ERROR("failed getting random number from RCL");
    *status=EXIT_FAILURE;
  }

  return(rand);  
#else

  // Use the HEAdas random number generator.
  return(HDmtDrand());

  // Status variable is not needed.
  (void)(*status);
#endif
}


void sixt_init_rng(const int seed, int* const status)
{
  // Initialize HEAdasS random number generator.
  // Note that this has to be done in any case, even
  // if the RCL random number server is used, because
  // the HEAdas routines (like heasp) rely on HDmtDrand().
  HDmtInit(seed);

#ifdef USE_RCL

  // Call the RCL random number generator specifying the
  // server and method.
  int rcl_status=0;
  rcl_rand_net_get_double("draco", "rand", &rcl_status);

  if(RCL_RANDOM_SUCCESS!=rcl_status) {
    SIXT_ERROR("failed getting random number from RCL");
    *status=EXIT_FAILURE;
  }

#else 

  // The status variable is not used.
  (void)(*status);

#endif
}


void sixt_destroy_rng()
{
  // Release HEADAS random number generator:
  HDmtFree();
}


void sixt_get_gauss_random_numbers(double* const x, 
				   double* const y, 
				   int* const status)
{
  double sqrt_2rho=sqrt(-log(sixt_get_random_number(status))*2.);
  CHECK_STATUS_VOID(*status);
  double phi=sixt_get_random_number(status)*2.*M_PI;
  CHECK_STATUS_VOID(*status);

  *x=sqrt_2rho * cos(phi);
  *y=sqrt_2rho * sin(phi);
}


double rndexp(const double avgdist, int* const status)
{
  assert(avgdist>0.);

  double rand=sixt_get_random_number(status);
  CHECK_STATUS_RET(*status, 0.);
  if (rand<1.e-15) {
    rand=1.e-15;
  }

  return(-log(rand)*avgdist);
}


void strtoupper(char* const string) 
{
  int count=0;
  while (string[count]!='\0') {
    string[count]=toupper(string[count]);
    count++;
  }
}


void sixt_error(const char* const func, const char* const msg)
{
  // Use the HEADAS error output routine.
  char output[MAXMSG];
  sprintf(output, "Error in %s: %s!\n", func, msg);
  HD_ERROR_THROW(output, EXIT_FAILURE);
}


void sixt_warning(const char* const msg)
{
  // Print the formatted output message.
  headas_chat(1, "### Warning: %s!\n", msg);
}


void sixt_get_XMLFile(char* const filename,
		      const char* const xmlfile,
		      const char* const mission,
		      const char* const instrument,
		      const char* const mode,
		      int* const status)
{
  // Convert the user input to capital letters.
  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];
  strcpy(Mission, mission);
  strcpy(Instrument, instrument);
  strcpy(Mode, mode);
  strtoupper(Mission);
  strtoupper(Instrument);
  strtoupper(Mode);

  // Check the available missions, instruments, and modes.
  char XMLFile[MAXFILENAME];
  strcpy(XMLFile, xmlfile);
  strtoupper(XMLFile);
  if (0==strcmp(XMLFile, "NONE")) {
    // Determine the base directory containing the XML
    // definition files.
    strcpy(filename, SIXT_DATA_PATH);
    strcat(filename, "/instruments");

    // Determine the XML filename according to the selected
    // mission, instrument, and mode.
    if (0==strcmp(Mission, "SRG")) {
      strcat(filename, "/srg");
      if (0==strcmp(Instrument, "EROSITA")) {
	strcat(filename, "/erosita.xml");
      } else {
	*status=EXIT_FAILURE;
	SIXT_ERROR("selected instrument is not supported");
	return;
      }

    } else if (0==strcmp(Mission, "IXO")) {
      strcat(filename, "/ixo");
      if (0==strcmp(Instrument, "WFI")) {
	strcat(filename, "/wfi");
	if (0==strcmp(Mode, "FULLFRAME")) {
	  strcat(filename, "/fullframe.xml");
	} else {
	  *status=EXIT_FAILURE;
	  SIXT_ERROR("selected mode is not supported");
	  return;
	}
      } else {
	*status=EXIT_FAILURE;
	SIXT_ERROR("selected instrument is not supported");
	return;
      }

    } else if (0==strcmp(Mission, "ATHENA")) {
      strcat(filename, "/athena");
      if (0==strcmp(Instrument, "WFI")) {
	strcat(filename, "/wfi");
	if (0==strcmp(Mode, "FULLFRAME")) {
	  strcat(filename, "/fullframe.xml");
	} else {
	  *status=EXIT_FAILURE;
	  SIXT_ERROR("selected mode is not supported");
	  return;
	}
      } else {
	*status=EXIT_FAILURE;
	SIXT_ERROR("selected instrument is not supported");
	return;
      }

    } else if (0==strcmp(Mission, "GRAVITAS")) {
      strcat(filename, "/gravitas");
      if (0==strcmp(Instrument, "HIFI")) {
	strcat(filename, "/hifi.xml");
      } else {
	*status=EXIT_FAILURE;
	SIXT_ERROR("selected instrument is not supported");
	return;
      }
      
    } else if (0==strcmp(Mission, "NUSTAR")) {
      strcat(filename, "/nustar");
      if (0==strcmp(Instrument, "NUSTAR")) {
	strcat(filename, "/nustar.xml");
      } else {
	*status=EXIT_FAILURE;
	SIXT_ERROR("selected instrument is not supported");
	return;
      }

    } else {
      *status=EXIT_FAILURE;
      SIXT_ERROR("selected mission is not supported");
      return;
    }
    
  } else {
    // The XML filename has been given explicitly.
    strcpy(filename, xmlfile);
  }
}


void sixt_get_LADXMLFile(char* const filename,
			 const char* const xmlfile)
{
  // Check the available missions, instruments, and modes.
  char XMLFile[MAXFILENAME];
  strcpy(XMLFile, xmlfile);
  strtoupper(XMLFile);
  if (0==strcmp(XMLFile, "NONE")) {
    // Set default LAD XML file.
    strcpy(filename, SIXT_DATA_PATH);
    strcat(filename, "/instruments/loft/lad.xml");
    
  } else {
    // The XML filename has been given explicitly.
    strcpy(filename, xmlfile);
  }
}

