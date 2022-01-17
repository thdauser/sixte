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


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

//
// sixt_main.c: modification of headas_main.c that preserves
// argc/argv as sixt_argc and sixt_argv for those SIXT tools that
// need argc/argv access
//

// Otherwise the tool works identically to headas_main


/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
#include "headas.h"
/******************************************************************************/

int sixt_argc;
char **sixt_argv;

/* C/C++ compatibility. */
#ifdef __cplusplus
extern "C" {
#endif

  /****************************************************************************
   * Function declarations.                                                   *
   ****************************************************************************/
  int headas_init(int argc, char **argv);
  int headas_close(int errNum);
  /****************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
}
#endif

// TODO: why do we need this??
#ifndef TOOLSUB
#error: TOOLSUB is not defined
#endif

int TOOLSUB(void);

/* Universal main function. */
int main(int argc, char **argv) {
  int status;

  sixt_argc=argc;
  sixt_argv=argv;
    
  printf("SIXTE version %s\n",PACKAGE_VERSION);

  status = headas_init(argc, argv);

  if(0 == status) status = TOOLSUB();

  status = headas_close(status);

  return status;
}
