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

 Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
 Erlangen-Nuernberg
 */

#include "radec2xylib.h"


// Convert world coordinates to image coordinates X and Y.
ImgPos radec2xy( Event* const event, struct wcsprm* const wcs, int* const status ) {

    double world[2] = {
            event->ra * 180. / M_PI,
            event->dec * 180. / M_PI
    };

    ImgPos pos;
    double imgcrd[2], pixcrd[2];
    double phi, theta;
    int wcsstatus = 0;

    wcss2p(wcs, 1, 2, world, &phi, &theta, imgcrd, pixcrd, &wcsstatus);
    if (0 != wcsstatus) {
        char msg[MAXMSG];
        sprintf(msg,
                "WCS coordinate conversion failed (RA=%lf, Dec=%lf, error code %d)",
                world[0], world[1], wcsstatus);

        char *projection, *wcstype;
        wcstype = strdup(wcs->ctype[0]);
        if( wcstype != NULL ) {
            while ((projection = strsep(&wcstype, "-")) != NULL);
            if (projection != NULL && strcmp(projection, "AIT") != 0) {
                char tmpmsg[MAXMSG];
                sprintf(tmpmsg, "\n(---> You might want to consider the AIT projection type "
                                "instead of %s)", projection);
                strcat(msg, tmpmsg);
            }
        }
        SIXT_ERROR(msg);
        *status = EXIT_FAILURE;
    }

    pos.x = (long) pixcrd[0];
    if (pixcrd[0] < 0.){
        pos.x--;
    }
    pos.y = (long) pixcrd[1];
    if (pixcrd[1] < 0.){
        pos.y--;
    }
    return pos;
}

// WCS data structure used for projection.
struct wcsprm getRadec2xyWCS( float* RefRA, float* RefDec, char* Projection, int* const status ) {
    /** RefRa: Right ascension of reference point [deg].
        RefDec: Declination of reference point [deg].
        Projection type (usually SIN).
    */
    struct wcsprm wcs = {.flag=-1};

    // Check whether an appropriate WCS projection has been selected.
    if (strlen(Projection) != 3) {
        SIXT_ERROR("invalid WCS projection type");
        *status = EXIT_FAILURE;
        return wcs;
    }

    if (0 != wcsini(1, 2, &wcs)) {
        SIXT_ERROR("initalization of WCS data structure failed");
        *status = EXIT_FAILURE;
        return wcs;
    }

    wcs.crpix[0] = 0.0;
    wcs.crpix[1] = 0.0;
    wcs.crval[0] = (double) *RefRA;
    wcs.crval[1] = (double) *RefDec;
    wcs.cdelt[0] = -0.05 / 3600.;
    wcs.cdelt[1] = 0.05 / 3600.;
    strcpy(wcs.cunit[0], "deg");
    strcpy(wcs.cunit[1], "deg");
    strcpy(wcs.ctype[0], "RA---");
    strcat(wcs.ctype[0], Projection);
    strcpy(wcs.ctype[1], "DEC--");
    strcat(wcs.ctype[1], Projection);

    return wcs;
}


void addXY2eventfile(EventFile* const evtfile, float* RefRA, float* RefDec, char* Projection, int* const status) {

    headas_chat(3, "add XY coordinates to event file ...\n");

    // Check if the input file contains recombined event patterns.
    fits_movnam_hdu(evtfile->fptr, BINARY_TBL, "EVENTS", 0, status);
    CHECK_STATUS_VOID(*status);

    char evtype[MAXMSG], comment[MAXMSG];
    fits_read_key(evtfile->fptr, TSTRING, "EVTYPE", evtype, comment, status);
    if (EXIT_SUCCESS != *status) {
        SIXT_ERROR("could not read FITS keyword 'EVTYPE'");
        return;
    }
    strtoupper(evtype);
    if (0 != strcmp(evtype, "PATTERN")) {
        *status = EXIT_FAILURE;
        char msg[MAXMSG];
        sprintf(msg, "event type of input file is '%s' (must be 'PATTERN')", evtype);
        SIXT_ERROR(msg);
        return;
    }

    // Determine WCS
    struct wcsprm wcs = getRadec2xyWCS(RefRA, RefDec, Projection, status);
    CHECK_STATUS_VOID(*status)

    // Add 'X', 'Y' column to evtfile if necessary
    int cx;
    fits_get_colnum(evtfile->fptr, CASEINSEN, "X", &cx, status);
    if( *status == COL_NOT_FOUND ){
        fits_clear_errmsg();
        *status=EXIT_SUCCESS;
        cx = evtfile->cra + 1;
        addCol2EventFile(evtfile, &cx, "X", "J", "", status);
        CHECK_STATUS_VOID(*status);
    }

    int cy;
    fits_get_colnum(evtfile->fptr, CASEINSEN, "Y", &cy, status);
    if( *status == COL_NOT_FOUND ){
        fits_clear_errmsg();
        *status=EXIT_SUCCESS;
        cy = evtfile->cra + 2;
        addCol2EventFile(evtfile, &cy, "Y", "J", "", status);
        CHECK_STATUS_VOID(*status);
    }

    // Update the WCS keywords in the output file.
    char keyword[MAXMSG];
    sprintf(keyword, "TCTYP%d", cx);
    fits_update_key(evtfile->fptr, TSTRING, keyword, wcs.ctype[0],
                    "projection type", status);
    sprintf(keyword, "TCTYP%d", cy);
    fits_update_key(evtfile->fptr, TSTRING, keyword, wcs.ctype[1],
                    "projection type", status);
    sprintf(keyword, "TCRVL%d", cx);
    fits_update_key(evtfile->fptr, TDOUBLE, keyword, &wcs.crval[0],
                    "reference value", status);
    sprintf(keyword, "TCRVL%d", cy);
    fits_update_key(evtfile->fptr, TDOUBLE, keyword, &wcs.crval[1],
                    "reference value", status);
    sprintf(keyword, "TCRPX%d", cx);
    fits_update_key(evtfile->fptr, TFLOAT, keyword, &wcs.crpix[0],
                    "reference point", status);
    sprintf(keyword, "TCRPX%d", cy);
    fits_update_key(evtfile->fptr, TFLOAT, keyword, &wcs.crpix[1],
                    "reference point", status);
    sprintf(keyword, "TCDLT%d", cx);
    fits_update_key(evtfile->fptr, TDOUBLE, keyword, &wcs.cdelt[0],
                    "pixel increment", status);
    sprintf(keyword, "TCDLT%d", cy);
    fits_update_key(evtfile->fptr, TDOUBLE, keyword, &wcs.cdelt[1],
                    "pixel increment", status);
    sprintf(keyword, "TCUNI%d", cx);
    fits_update_key(evtfile->fptr, TSTRING, keyword, wcs.cunit[0],
                    "axis units", status);
    sprintf(keyword, "TCUNI%d", cy);
    fits_update_key(evtfile->fptr, TSTRING, keyword, wcs.cunit[1],
                    "axis units",status);
    CHECK_STATUS_VOID(*status);

    fits_update_key(evtfile->fptr, TSTRING, "REFXCTYP", wcs.ctype[0],
                    "projection type", status);
    fits_update_key(evtfile->fptr, TSTRING, "REFYCTYP", wcs.ctype[1],
                    "projection type", status);
    fits_update_key(evtfile->fptr, TSTRING, "REFXCUNI", wcs.cunit[0],
                    "axis units", status);
    fits_update_key(evtfile->fptr, TSTRING, "REFYCUNI", wcs.cunit[1],
                    "axis units", status);
    fits_update_key(evtfile->fptr, TFLOAT, "REFXCRPX", &wcs.crpix[0],
                    "reference value", status);
    fits_update_key(evtfile->fptr, TFLOAT, "REFYCRPX", &wcs.crpix[1],
                    "reference value", status);
    fits_update_key(evtfile->fptr, TDOUBLE, "REFXCRVL", &wcs.crval[0],
                    "reference value", status);
    fits_update_key(evtfile->fptr, TDOUBLE, "REFYCRVL", &wcs.crval[1],
                    "reference value", status);
    fits_update_key(evtfile->fptr, TDOUBLE, "REFXCDLT", &wcs.cdelt[0],
                    "pixel increment", status);
    fits_update_key(evtfile->fptr, TDOUBLE, "REFYCDLT", &wcs.cdelt[1],
                    "pixel increment", status);
    CHECK_STATUS_VOID(*status);


    // Loop over all events in the input list.
    Event* event = getEvent(status);
    CHECK_STATUS_VOID(*status);
    long ii;
    for (ii = 1; ii <= evtfile->nrows; ii++) {
        // Get event
        getEventFromFile(evtfile, (int)ii, event, status);
        CHECK_STATUS_BREAK(*status);

        ImgPos pos = radec2xy( event, &wcs, status );
        CHECK_STATUS_BREAK(*status);

        // Save changes to eventfile
        fits_write_col(evtfile->fptr, TLONG, cx, ii, 1, 1, &pos.x, status);
        fits_write_col(evtfile->fptr, TLONG, cy, ii, 1, 1, &pos.y, status);
        CHECK_STATUS_BREAK(*status);

    }
    // Release memory.
    freeEvent(&event);
    wcsfree(&wcs);

    if (*status == EXIT_SUCCESS) {
        headas_chat(3, " ... adding X, Y coordinates successful!\n");
    } else {
        char msg[MAXMSG];
        sprintf(msg,"*** ERROR occurred while adding X, Y coordinates to eventfile");
        SIXT_ERROR(msg);
        *status = EXIT_FAILURE;
    }
}