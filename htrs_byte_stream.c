#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>


#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#define TOOLSUB htrs_byte_stream_main
#include "headas_main.c"



//////////////////////////////////
int htrs_byte_stream_main()
{
  int status=EXIT_SUCCESS;


  // HEATOOLs: register program
  set_toolname("htrs_byte_stream");
  set_toolversion("0.01");



  return(status);
}

