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

#include "radec2xy.h"


int radec2xy_main() {
    // Containing all programm parameters read by PIL
    struct Parameters par;

    // Input event file.
    EventFile* elf = NULL;

    // Error status.
    int status = EXIT_SUCCESS;

    // Register HEATOOL:
    set_toolname("radec2xy");
    set_toolversion("0.10");

    do { // Beginning of the ERROR handling loop (will at most be run once).

        // --- Initialization ---
        headas_chat(3, "initialization ...\n");

        // Read parameters using PIL library:
        radec2xy_getpar(&par,&status);
        CHECK_STATUS_BREAK(status);

        // Open the input event file.
        elf = openEventFile(par.EvtFile, READWRITE, &status);
        CHECK_STATUS_BREAK(status);

        addXY2eventfile(elf, &par.RefRA, &par.RefDec, par.Projection, &status);
        CHECK_STATUS_BREAK(status);

        // Append a check sum to the header of the event extension.
        int hdutype = 0;
        fits_movabs_hdu(elf->fptr, 2, &hdutype, &status);
        fits_write_chksum(elf->fptr, &status);
        CHECK_STATUS_BREAK(status);

    } while (0); // END of the error handling loop.

    // --- Cleaning up ---
    headas_chat(3, "cleaning up ...\n");

    // Close the files.
    freeEventFile(&elf, &status);

    if (EXIT_SUCCESS == status) {
        headas_chat(3, "finished successfully!\n\n");
        return (EXIT_SUCCESS);
    } else {
        return (EXIT_FAILURE);
    }
}


void radec2xy_getpar(struct Parameters* const par, int* const status){
    query_simput_parameter_file_name_buffer("EvtFile", par->EvtFile, MAXFILENAME, status);
    query_simput_parameter_string_buffer("Projection", par->Projection, MAXFILENAME, status);
    query_simput_parameter_float("RefRa", &par->RefRA, status);
    query_simput_parameter_float("RefDec", &par->RefDec, status);
}
