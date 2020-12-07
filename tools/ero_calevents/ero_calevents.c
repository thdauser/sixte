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

#include "ero_calevents.h"


int ero_calevents_main() {
    // Containing all programm parameters read by PIL
    struct Parameters par;

    // Input event file.
    EventFile *elf = NULL;

    // File pointer to the output eROSITA event file.
    fitsfile *fptr = NULL;

    // WCS data structure used for projection.
    struct wcsprm wcs = {.flag=-1};
    // String buffer for FITS header.
    char *headerstr = NULL;

    GTI *gti = NULL;
    Attitude *ac = NULL;

    // Error status.
    int status = EXIT_SUCCESS;


    // Register HEATOOL:
    set_toolname("ero_calevents");
    set_toolversion("0.20");


    do { // Beginning of the ERROR handling loop (will at most be run once).

        // --- Initialization ---

        headas_chat(3, "initialization ...\n");

        // Read parameters using PIL library:
        if ((status = getpar(&par))) break;

        // Check whether an appropriate WCS projection has been selected.
        if (strlen(par.Projection) != 3) {
            SIXT_ERROR("invalid WCS projection type");
            status = EXIT_FAILURE;
            break;
        }


        // Open the input event file.
        elf = openEventFile(par.EvtFile, READONLY, &status);
        CHECK_STATUS_BREAK(status);

        // Check if the input file contains recombined event patterns.
        char evtype[MAXMSG], comment[MAXMSG];
        fits_read_key(elf->fptr, TSTRING, "EVTYPE", evtype, comment, &status);
        if (EXIT_SUCCESS != status) {
            SIXT_ERROR("could not read FITS keyword 'EVTYPE'");
            break;
        }
        strtoupper(evtype);
        if (0 != strcmp(evtype, "PATTERN")) {
            status = EXIT_FAILURE;
            char msg[MAXMSG];
            sprintf(msg, "event type of input file is '%s' (must be 'PATTERN')", evtype);
            SIXT_ERROR(msg);
            break;
        }

        // Read keywords from the input file.
        double mjdref = 0.0;
        fits_read_key(elf->fptr, TDOUBLE, "MJDREF", &mjdref, comment, &status);
        if (EXIT_SUCCESS != status) {
            char msg[MAXMSG];
            sprintf(msg, "could not read FITS keyword 'MJDREF' from input "
                         "event list '%s'", par.EvtFile);
            SIXT_ERROR(msg);
            break;
        }

        double timezero = 0.0;
        fits_write_errmark();
        int status2 = EXIT_SUCCESS;
        fits_read_key(elf->fptr, TDOUBLE, "TIMEZERO", &timezero, comment, &status2);
        fits_clear_errmark();
        if (EXIT_SUCCESS != status2) {
            timezero = 0.;
        }

        char date_obs[MAXMSG];
        fits_read_key(elf->fptr, TSTRING, "DATE-OBS", date_obs, comment, &status);
        if (EXIT_SUCCESS != status) {
            char msg[MAXMSG];
            sprintf(msg, "could not read FITS keyword 'DATE-OBS' from input "
                         "event list '%s'", par.EvtFile);
            SIXT_ERROR(msg);
            break;
        }

        char time_obs[MAXMSG];
        fits_read_key(elf->fptr, TSTRING, "TIME-OBS", time_obs, comment, &status);
        if (EXIT_SUCCESS != status) {
            char msg[MAXMSG];
            sprintf(msg, "could not read FITS keyword 'TIME-OBS' from input "
                         "event list '%s'", par.EvtFile);
            SIXT_ERROR(msg);
            break;
        }

        char date_end[MAXMSG];
        fits_read_key(elf->fptr, TSTRING, "DATE-END", date_end, comment, &status);
        if (EXIT_SUCCESS != status) {
            char msg[MAXMSG];
            sprintf(msg, "could not read FITS keyword 'DATE-END' from input "
                         "event list '%s'", par.EvtFile);
            SIXT_ERROR(msg);
            break;
        }

        char time_end[MAXMSG];
        fits_read_key(elf->fptr, TSTRING, "TIME-END", time_end, comment, &status);
        if (EXIT_SUCCESS != status) {
            char msg[MAXMSG];
            sprintf(msg, "could not read FITS keyword 'TIME-END' from input "
                         "event list '%s'", par.EvtFile);
            SIXT_ERROR(msg);
            break;
        }

        double tstart = 0.0;
        fits_read_key(elf->fptr, TDOUBLE, "TSTART", &tstart, comment, &status);
        if (EXIT_SUCCESS != status) {
            char msg[MAXMSG];
            sprintf(msg, "could not read FITS keyword 'TSTART' from input "
                         "event list '%s'", par.EvtFile);
            SIXT_ERROR(msg);
            break;
        }

        double tstop = 0.0;
        fits_read_key(elf->fptr, TDOUBLE, "TSTOP", &tstop, comment, &status);
        if (EXIT_SUCCESS != status) {
            char msg[MAXMSG];
            sprintf(msg, "could not read FITS keyword 'TSTOP' from input "
                         "event list '%s'", par.EvtFile);
            SIXT_ERROR(msg);
            break;
        }

        // Verify values of MJDREF and TIMEZERO.
        verifyMJDREF(eromjdref, mjdref, "in event file", &status);
        CHECK_STATUS_BREAK(status);
        verifyTIMEZERO(timezero, &status);
        CHECK_STATUS_BREAK(status);

        // Determine the file creation date for the header.
        char creation_date[MAXMSG];
        int timeref;
        fits_get_system_time(creation_date, &timeref, &status);
        CHECK_STATUS_BREAK(status);

        // Check if the output file already exists.
        int exists;
        fits_file_exists(par.eroEvtFile, &exists, &status);
        CHECK_STATUS_BREAK(status);
        if (0 != exists) {
            if (0 != par.clobber) {
                // Delete the file.
                remove(par.eroEvtFile);
            } else {
                // Throw an error.
                char msg[MAXMSG];
                sprintf(msg, "file '%s' already exists", par.eroEvtFile);
                SIXT_ERROR(msg);
                status = EXIT_FAILURE;
                break;
            }
        }

        // Read arf from header key to determine filter
        char arffile[MAXMSG];
        fits_read_key(elf->fptr, TSTRING, "ANCRFILE", arffile, comment, &status);
        if (EXIT_SUCCESS != status) {
            char msg[MAXMSG];
            sprintf(msg, "could not read FITS keyword 'ANCRFILE' from input "
                         "event list '%s'", par.EvtFile);
            SIXT_ERROR(msg);
            break;
        }

        // Read the FILTER key word (from event file)
        char filter[MAXMSG];
        fits_read_key(elf->fptr, TSTRING, "FILTER", filter, comment, &status);
        if (EXIT_SUCCESS != status) {
            char msg[MAXMSG];
            sprintf(msg, "could not read FITS keyword 'FILTER' from input "
                         "event list '%s'", par.EvtFile);
            SIXT_ERROR(msg);
            break;
        }

        // Create and open a new FITS file.
        headas_chat(3, "create new eROSITA event list file '%s' ...\n",
                    par.eroEvtFile);
        fits_create_file(&fptr, par.eroEvtFile, &status);
        CHECK_STATUS_BREAK(status);

        // Create the event table.
        char *ttype[] = {"TIME", "RA", "DEC", "X", "Y", "ENERGY",
                         "EV_WEIGHT", "RAWX", "RAWY", "SUBX", "SUBY",
                         "PHA", "PAT_TYP", "PAT_INF", "TM_NR",
                         "FLAG", "FRAME", "RECORDTIME", "FRAMETIME","PI"};
        char *tform[] = {"D", "D", "D", "D", "D", "E",
                         "E", "I", "I", "D", "D",
                         "I", "I", "B", "B",
                         "J", "J","D","D","E"};
        char *tunit[] = {"", "", "", "", "", "eV",
                         "", "", "", "", "",
                         "adu", "", "", "",
                         "", "", "s", "s","eV"};
        fits_create_tbl(fptr, BINARY_TBL, 0, 20, ttype, tform, tunit,
                        "EVENTS", &status);
        if (EXIT_SUCCESS != status) {
            char msg[MAXMSG];
            sprintf(msg, "could not create binary table for events "
                         "in file '%s'", par.eroEvtFile);
            SIXT_ERROR(msg);
            break;
        }

        // Insert header keywords.
        char hduclass[MAXMSG] = "OGIP";
        fits_update_key(fptr, TSTRING, "HDUCLASS", hduclass, "", &status);
        char hduclas1[MAXMSG] = "EVENTS";
        fits_update_key(fptr, TSTRING, "HDUCLAS1", hduclas1, "", &status);
        CHECK_STATUS_BREAK(status);

        // Insert the standard eROSITA header keywords.
        sixt_add_fits_erostdkeywords(fptr, 1, filter, creation_date, date_obs, time_obs,
                                     date_end, time_end, tstart, tstop,
                                     mjdref, timezero, par.CCDNr, &status);
        CHECK_STATUS_BREAK(status);
        sixt_add_fits_erostdkeywords(fptr, 2, filter, creation_date, date_obs, time_obs,
                                     date_end, time_end, tstart, tstop,
                                     mjdref, timezero, par.CCDNr, &status);
        CHECK_STATUS_BREAK(status);
	
        // Determine the column numbers.
        int ctime, cra, cdec, cx, cy, cenergy, cev_weight, crawx, crawy,
	  csubx, csuby, cpha, cpat_typ, cpat_inf, cccdnr, cflag,
	  cframe, crecordtime, cframetime, cpi;
        fits_get_colnum(fptr, CASEINSEN, "TIME", &ctime, &status);
        fits_get_colnum(fptr, CASEINSEN, "RA", &cra, &status);
        fits_get_colnum(fptr, CASEINSEN, "DEC", &cdec, &status);
        fits_get_colnum(fptr, CASEINSEN, "X", &cx, &status);
        fits_get_colnum(fptr, CASEINSEN, "Y", &cy, &status);
        fits_get_colnum(fptr, CASEINSEN, "ENERGY", &cenergy, &status);
        fits_get_colnum(fptr, CASEINSEN, "EV_WEIGHT", &cev_weight, &status);
        fits_get_colnum(fptr, CASEINSEN, "RAWX", &crawx, &status);
        fits_get_colnum(fptr, CASEINSEN, "RAWY", &crawy, &status);
        fits_get_colnum(fptr, CASEINSEN, "SUBX", &csubx, &status);
        fits_get_colnum(fptr, CASEINSEN, "SUBY", &csuby, &status);
        fits_get_colnum(fptr, CASEINSEN, "PHA", &cpha, &status);
        fits_get_colnum(fptr, CASEINSEN, "PAT_TYP", &cpat_typ, &status);
        fits_get_colnum(fptr, CASEINSEN, "PAT_INF", &cpat_inf, &status);
        fits_get_colnum(fptr, CASEINSEN, "TM_NR", &cccdnr, &status);
        fits_get_colnum(fptr, CASEINSEN, "FLAG", &cflag, &status);
        fits_get_colnum(fptr, CASEINSEN, "FRAME", &cframe, &status);
	fits_get_colnum(fptr, CASEINSEN, "RECORDTIME", &crecordtime, &status);
	fits_get_colnum(fptr, CASEINSEN, "FRAMETIME", &cframetime, &status);
        fits_get_colnum(fptr, CASEINSEN, "PI", &cpi, &status);
        CHECK_STATUS_BREAK(status);

        // Set the TLMIN and TLMAX keywords.
        // For the PHA column.
        char keyword[MAXMSG];
        int tlmin_pha = 0, tlmax_pha = 4095;
        sprintf(keyword, "TLMIN%d", cpha);
        fits_update_key(fptr, TINT, keyword, &tlmin_pha, "", &status);
        sprintf(keyword, "TLMAX%d", cpha);
        fits_update_key(fptr, TINT, keyword, &tlmax_pha, "", &status);
        CHECK_STATUS_BREAK(status);

        // For the ENERGY column.
        float tlmin_energy = 0.0, tlmax_energy = 20.48;
        sprintf(keyword, "TLMIN%d", cenergy);
        fits_update_key(fptr, TFLOAT, keyword, &tlmin_energy, "", &status);
        sprintf(keyword, "TLMAX%d", cenergy);
        fits_update_key(fptr, TFLOAT, keyword, &tlmax_energy, "", &status);
        CHECK_STATUS_BREAK(status);

        // For the X and Y column.
        long tlmin_x = -12960000, tlmax_x = 12960000;
        long tlmin_y = -6480000, tlmax_y = 6480000;
        sprintf(keyword, "TLMIN%d", cx);
        fits_update_key(fptr, TLONG, keyword, &tlmin_x, "", &status);
        fits_update_key(fptr, TLONG, "REFXLMIN", &tlmin_x, "", &status);
        sprintf(keyword, "TLMAX%d", cx);
        fits_update_key(fptr, TLONG, keyword, &tlmax_x, "", &status);
        fits_update_key(fptr, TLONG, "REFXLMAX", &tlmax_x, "", &status);
        sprintf(keyword, "TLMIN%d", cy);
        fits_update_key(fptr, TLONG, keyword, &tlmin_y, "", &status);
        fits_update_key(fptr, TLONG, "REFYLMIN", &tlmin_y, "", &status);
        sprintf(keyword, "TLMAX%d", cy);
        fits_update_key(fptr, TLONG, keyword, &tlmax_y, "", &status);
        fits_update_key(fptr, TLONG, "REFYLMAX", &tlmax_y, "", &status);
        CHECK_STATUS_BREAK(status);

        // Set up the WCS data structure.
        struct wcsprm wcs = getRadec2xyWCS( &par.RefRA, &par.RefDec, par.Projection, &status );
        CHECK_STATUS_BREAK(status);

        // Update the WCS keywords in the output file.
        sprintf(keyword, "TCTYP%d", cx);
        fits_update_key(fptr, TSTRING, keyword, wcs.ctype[0],
                        "projection type", &status);
        sprintf(keyword, "TCTYP%d", cy);
        fits_update_key(fptr, TSTRING, keyword, wcs.ctype[1],
                        "projection type", &status);
        sprintf(keyword, "TCRVL%d", cx);
        fits_update_key(fptr, TDOUBLE, keyword, &wcs.crval[0],
                        "reference value", &status);
        sprintf(keyword, "TCRVL%d", cy);
        fits_update_key(fptr, TDOUBLE, keyword, &wcs.crval[1],
                        "reference value", &status);
        sprintf(keyword, "TCRPX%d", cx);
        fits_update_key(fptr, TFLOAT, keyword, &wcs.crpix[0],
                        "reference point", &status);
        sprintf(keyword, "TCRPX%d", cy);
        fits_update_key(fptr, TFLOAT, keyword, &wcs.crpix[1],
                        "reference point", &status);
        sprintf(keyword, "TCDLT%d", cx);
        fits_update_key(fptr, TDOUBLE, keyword, &wcs.cdelt[0],
                        "pixel increment", &status);
        sprintf(keyword, "TCDLT%d", cy);
        fits_update_key(fptr, TDOUBLE, keyword, &wcs.cdelt[1],
                        "pixel increment", &status);
        sprintf(keyword, "TCUNI%d", cx);
        fits_update_key(fptr, TSTRING, keyword, wcs.cunit[0],
                        "axis units", &status);
        sprintf(keyword, "TCUNI%d", cy);
        fits_update_key(fptr, TSTRING, keyword, wcs.cunit[1],
                        "axis units", &status);
        CHECK_STATUS_BREAK(status);

        fits_update_key(fptr, TSTRING, "REFXCTYP", wcs.ctype[0],
                        "projection type", &status);
        fits_update_key(fptr, TSTRING, "REFYCTYP", wcs.ctype[1],
                        "projection type", &status);
        fits_update_key(fptr, TSTRING, "REFXCUNI", wcs.cunit[0],
                        "axis units", &status);
        fits_update_key(fptr, TSTRING, "REFYCUNI", wcs.cunit[1],
                        "axis units", &status);
        fits_update_key(fptr, TFLOAT, "REFXCRPX", &wcs.crpix[0],
                        "reference value", &status);
        fits_update_key(fptr, TFLOAT, "REFYCRPX", &wcs.crpix[1],
                        "reference value", &status);
        fits_update_key(fptr, TDOUBLE, "REFXCRVL", &wcs.crval[0],
                        "reference value", &status);
        fits_update_key(fptr, TDOUBLE, "REFYCRVL", &wcs.crval[1],
                        "reference value", &status);
        fits_update_key(fptr, TDOUBLE, "REFXCDLT", &wcs.cdelt[0],
                        "pixel increment", &status);
        fits_update_key(fptr, TDOUBLE, "REFYCDLT", &wcs.cdelt[1],
                        "pixel increment", &status);
        CHECK_STATUS_BREAK(status);

        // Set the TZERO keywords for the columns SUBX and SUBY. Note that
        // the TZERO values also have to be set with the routine
        // fits_set_tscale(). Otherwise CFITSIO will access the raw values
        // in the file.
        double tzero_subx_suby = -0.843333333333333;
        double tscal_subx_suby = 6.66666666666667e-3;
        sprintf(keyword, "TZERO%d", csubx);
        fits_update_key(fptr, TDOUBLE, keyword, &tzero_subx_suby, "", &status);
        sprintf(keyword, "TSCAL%d", csubx);
        fits_update_key(fptr, TDOUBLE, keyword, &tscal_subx_suby, "", &status);
        fits_set_tscale(fptr, csubx, tscal_subx_suby, tzero_subx_suby, &status);
        sprintf(keyword, "TZERO%d", csuby);
        fits_update_key(fptr, TDOUBLE, keyword, &tzero_subx_suby, "", &status);
        sprintf(keyword, "TSCAL%d", csuby);
        fits_update_key(fptr, TDOUBLE, keyword, &tscal_subx_suby, "", &status);
        fits_set_tscale(fptr, csuby, tscal_subx_suby, tzero_subx_suby, &status);
        CHECK_STATUS_BREAK(status);

        // Set the TZERO and TSCAL keywords for the columns RA and DEC.
        // Note that both values also have to be set with the routine
        // fits_set_tscale(). Otherwise CFITSIO will access the raw values
        // in the file.
        double tzero_ra_dec = 0.0, tscal_ra_dec = 1.e-6;
        sprintf(keyword, "TZERO%d", cra);
        fits_update_key(fptr, TDOUBLE, keyword, &tzero_ra_dec, "", &status);
        sprintf(keyword, "TSCAL%d", cra);
        fits_update_key(fptr, TDOUBLE, keyword, &tscal_ra_dec, "", &status);
        fits_set_tscale(fptr, cra, tscal_ra_dec, tzero_ra_dec, &status);
        sprintf(keyword, "TZERO%d", cdec);
        fits_update_key(fptr, TDOUBLE, keyword, &tzero_ra_dec, "", &status);
        sprintf(keyword, "TSCAL%d", cdec);
        fits_update_key(fptr, TDOUBLE, keyword, &tscal_ra_dec, "", &status);
        fits_set_tscale(fptr, cdec, tscal_ra_dec, tzero_ra_dec, &status);
        CHECK_STATUS_BREAK(status);

        // --- END of initialization ---

        // --- Beginning of copy events ---

        headas_chat(3, "copy events ...\n");

        // Actual minimum and maximum values of X and Y.
        long refxdmin, refxdmax, refydmin, refydmax;
        double ra_min, ra_max, dec_min, dec_max;

        // Loop over all events in the FITS file.
        long input_row, output_row = 0;
        for (input_row = 0; input_row < elf->nrows; input_row++) {

            // Read the next event from the input file.
            Event event;
            getEventFromFile(elf, input_row + 1, &event, &status);
            CHECK_STATUS_BREAK(status);

            // Determine the event data based on the event information.
            eroCalEvent ev;

            // Time and frame.
            ev.time = event.time;
            ev.frame = event.frame;

            ev.ra = event.ra * 180. / M_PI;
            if (event.ra < 0.) {
                SIXT_WARNING("value for right ascension <0.0deg");
            }
            ev.dec = event.dec * 180. / M_PI;

            // Determine the minimum and maximum values of RA and Dec in [rad].
            if (0 == input_row) {
                ra_min = event.ra;
                ra_max = event.ra;
                dec_min = event.dec;
                dec_max = event.dec;
            }
            if (event.ra < ra_min) {
                ra_min = event.ra;
            }
            if (event.ra > ra_max) {
                ra_max = event.ra;
            }
            if (event.dec < dec_min) {
                dec_min = event.dec;
            }
            if (event.dec > dec_max) {
                dec_max = event.dec;
            }

            // Convert world coordinates to image coordinates X and Y.
            ImgPos pos = radec2xy( &event, &wcs, &status );
            CHECK_STATUS_BREAK(status);
            ev.x = pos.x;
            ev.y = pos.y;

            // Determine the actual minimum and maximum values of X and Y.
            if (0 == input_row) {
                refxdmin = ev.x;
                refxdmax = ev.x;
                refydmin = ev.y;
                refydmax = ev.y;
            }
            if (ev.x < refxdmin) {
                refxdmin = ev.x;
            }
            if (ev.x > refxdmax) {
                refxdmax = ev.x;
            }
            if (ev.y < refydmin) {
                refydmin = ev.y;
            }
            if (ev.y > refydmax) {
                refydmax = ev.y;
            }

            // TODO In the current implementation the value of FLAG is set
            // by default. This needs to be changed later.
            ev.flag = 0xC00001C0;

            // TODO Inverse vignetting correction factor is not used.
            ev.ev_weight = 1.0; // Invers vignetting correction factor.

            // CCD number.
            ev.ccdnr = par.CCDNr;

	    // RECORDTIME (TDOUBLE, time in seconds): Time of the
	    // telemetry record, produced by ITC, responsible for
	    // (next to the 8bit FrameCounter) getting the data in the
	    // correct order. A bright source has many events and one
	    // CCD frame can be distributed over multiple telemetry
	    // frames. The recordtime can (but needn't)
	    // increase. Reference point is 1.1.2000.
	    ev.recordtime = ev.time+(mjdref-eromjdref)*86400.;

	    // FRAMETIME (TDOUBLE, time in seconds): Time stamp of the
	    // CCD Frames, reference point is 1.1.2000, produced by CE
	    ev.frametime = ev.time+(mjdref-eromjdref)*86400.;
	    
            // Loop over all split partners contributing to the event.
            int ii;
            for (ii = 0; ii < 9; ii++) {

                // Only regard split events with a non-vanishing contribution.
                if (event.signals[ii] <= 0.0) continue;

                // Raw pixel coordinates.

                /** eROSITA standard is that the CCD is viewed from the bottom
                 *  (towards the sky), which is the opposite of the Sixte standard
                 *  (from the mirrors onto the CCD). Need to flip the y-axis to
                 *  take this change into account
                 */
                int ibuffer = 384;  // hard coded here, bad style but ero CCDs will not change
                ev.rawx = ibuffer - 1 - event.rawx;  // rawx/y defined from 0 to ibuffer

                ev.rawx += ii % 3;

                ev.rawy = event.rawy + ii / 3;

                // TODO Sub-pixel resolution is not implemented.
                ev.subx = 0;
                ev.suby = 0;

                // Detected channel.
                ev.pha = event.phas[ii];
		
                // Calibrated and recombined amplitude in [eV].
                // The amplitude is positive for the main event only. For
                // split partners it is negative.
                if (4 == ii) {
                    ev.energy = event.signal * 1000.;
		    // Copy the raw energy grid into a new column named
		    // PI. The eSASS reads in the PI column (apperently in
		    // units of eV and not columns). This hack. however,
		    // bypasses all detector corrections. A proper PHA2PI
		    // correction has to be done!
		    ev.pi = event.signal * 1000.;
                } else {
                    ev.energy = -event.signal * 1000.;
		    ev.pi = -event.signal * 1000.;
                }
		
                // Event type.
                if (event.type >= 0) {
                    ev.pat_typ = event.npixels;
                } else {
                    // Invalid events.
                    ev.pat_typ = 0;
                }

                // Event type and alignment.
                if (event.type >= 1) {
                    int pixelnr = (ii + 1) - ((ii / 3) - 1) * 6;
                    ev.pat_inf = event.type * 10 + pixelnr;
                } else {
                    ev.pat_inf = 0;
                }

                // Store the event in the output file.
                output_row++;
                fits_write_col(fptr, TDOUBLE, ctime, output_row, 1, 1, &ev.time, &status);
                fits_write_col(fptr, TLONG, cframe, output_row, 1, 1, &ev.frame, &status);
                fits_write_col(fptr, TLONG, cpha, output_row, 1, 1, &ev.pha, &status);
                fits_write_col(fptr, TFLOAT, cenergy, output_row, 1, 1, &ev.energy, &status);
                fits_write_col(fptr, TINT, crawx, output_row, 1, 1, &ev.rawx, &status);
                fits_write_col(fptr, TINT, crawy, output_row, 1, 1, &ev.rawy, &status);
                fits_write_col(fptr, TDOUBLE, cra, output_row, 1, 1, &ev.ra, &status);
                fits_write_col(fptr, TDOUBLE, cdec, output_row, 1, 1, &ev.dec, &status);
                fits_write_col(fptr, TLONG, cx, output_row, 1, 1, &ev.x, &status);
                fits_write_col(fptr, TLONG, cy, output_row, 1, 1, &ev.y, &status);
                fits_write_col(fptr, TINT, csubx, output_row, 1, 1, &ev.subx, &status);
                fits_write_col(fptr, TINT, csuby, output_row, 1, 1, &ev.suby, &status);
                fits_write_col(fptr, TLONG, cflag, output_row, 1, 1, &ev.flag, &status);
                fits_write_col(fptr, TUINT, cpat_typ, output_row, 1, 1, &ev.pat_typ, &status);
                fits_write_col(fptr, TBYTE, cpat_inf, output_row, 1, 1, &ev.pat_inf, &status);
                fits_write_col(fptr, TFLOAT, cev_weight, output_row, 1, 1, &ev.ev_weight, &status);
                fits_write_col(fptr, TINT, cccdnr, output_row, 1, 1, &ev.ccdnr, &status);
                fits_write_col(fptr, TDOUBLE, crecordtime, output_row, 1, 1, &ev.recordtime, &status);
                fits_write_col(fptr, TDOUBLE, cframetime, output_row, 1, 1, &ev.frametime, &status);
                fits_write_col(fptr, TFLOAT, cpi, output_row, 1, 1, &ev.pi, &status);
                CHECK_STATUS_BREAK(status);
            }
            CHECK_STATUS_BREAK(status);
            // End of loop over all split partners.
        }
        CHECK_STATUS_BREAK(status);
        // END of loop over all events in the FITS file.

        // Set the RA_MIN, RA_MAX, DEC_MIN, DEC_MAX keywords (in [deg]).
        ra_min *= 180. / M_PI;
        ra_max *= 180. / M_PI;
        dec_min *= 180. / M_PI;
        dec_max *= 180. / M_PI;
        fits_update_key(fptr, TDOUBLE, "RA_MIN", &ra_min, "", &status);
        fits_update_key(fptr, TDOUBLE, "RA_MAX", &ra_max, "", &status);
        fits_update_key(fptr, TDOUBLE, "DEC_MIN", &dec_min, "", &status);
        fits_update_key(fptr, TDOUBLE, "DEC_MAX", &dec_max, "", &status);
        CHECK_STATUS_BREAK(status);

	int pat_sel=15;
	int flagsel=0;
        fits_update_key(fptr, TLONG, "FLAGSEL", &flagsel, "Flag selection", &status);
        fits_update_key(fptr, TINT, "PAT_SEL", &pat_sel, "Pattern selection (1=s,3=s+d,7=s+d+t,15=all)", &status);
        CHECK_STATUS_BREAK(status);

        // Set the number of unique events to the number of entries in the table.
        // long uniq_evt;
        // fits_get_num_rows(fptr, &uniq_evt, &status);
        // CHECK_STATUS_BREAK(status);
        // fits_update_key(fptr, TLONG, "UNIQ_EVT", &uniq_evt,
        //	    "Number of unique events inside", &status);
        // CHECK_STATUS_BREAK(status);

        // Set the REF?DMIN/MAX keywords.
        fits_update_key(fptr, TLONG, "REFXDMIN", &refxdmin, "", &status);
        fits_update_key(fptr, TLONG, "REFXDMAX", &refxdmax, "", &status);
        fits_update_key(fptr, TLONG, "REFYDMIN", &refydmin, "", &status);
        fits_update_key(fptr, TLONG, "REFYDMAX", &refydmax, "", &status);
        CHECK_STATUS_BREAK(status);

        // Determine the relative search threshold for split partners.
        fits_write_errmark();
        float spltthr;
        int opt_status = EXIT_SUCCESS;
        fits_read_key(elf->fptr, TFLOAT, "SPLTTHR", &spltthr, comment, &opt_status);
        if (EXIT_SUCCESS == opt_status) {
            fits_update_key(fptr, TFLOAT, "SPLTTHR", &spltthr,
                            "Relative search level for split events", &status);
            CHECK_STATUS_BREAK(status);
        }
        opt_status = EXIT_SUCCESS;
        fits_clear_errmark();

        // --- End of copy events ---

        // --- Begin of append GTI extension ---

        headas_chat(3, "append GTI extension ...\n");

        // Load the GTI extension from the input file.
        gti = loadGTI(par.EvtFile, &status);
        CHECK_STATUS_BREAK(status);

        // Make sure that the MJDREF of the GTI extension agrees with
        // the value in the input event file.
        verifyMJDREF(mjdref, gti->mjdref, "in GTI file", &status);
        CHECK_STATUS_BREAK(status);

        // Store the GTI extension in the output file.
        char gti_extname[MAXMSG];
        sprintf(gti_extname, "GTI%d", par.CCDNr);
        saveGTIExt(fptr, gti_extname, gti, &status);
        CHECK_STATUS_BREAK(status);

        // --- End of append GTI extension ---

        // --- Beginning of append DEADCOR extension ---

        headas_chat(3, "append DEADCOR extension ...\n");

        // Create the DEADCOR table.
        char deadcor_extname[MAXMSG];
        sprintf(deadcor_extname, "DEADCOR%d", par.CCDNr);
        char *deadcor_ttype[] = {"TIME", "DEADC"};
        char *deadcor_tform[] = {"D", "E"};
        char *deadcor_tunit[] = {"", ""};
        fits_create_tbl(fptr, BINARY_TBL, 0, 2,
                        deadcor_ttype, deadcor_tform, deadcor_tunit,
                        deadcor_extname, &status);
        if (EXIT_SUCCESS != status) {
            SIXT_ERROR("could not create binary table for DEADCOR extension");
            break;
        }

        // Insert header keywords.
        fits_update_key(fptr, TSTRING, "HDUCLASS", "OGIP", "", &status);
        fits_update_key(fptr, TSTRING, "HDUCLAS1", "TEMPORALDATA", "", &status);
        fits_update_key(fptr, TSTRING, "HDUCLAS2", "TSI", "", &status);
        CHECK_STATUS_BREAK(status);

        // Determine the individual column numbers.
        int cdeadcor_time, cdeadc;
        fits_get_colnum(fptr, CASEINSEN, "TIME", &cdeadcor_time, &status);
        fits_get_colnum(fptr, CASEINSEN, "DEADC", &cdeadc, &status);
        CHECK_STATUS_BREAK(status);

        // Store the data in the table.
        double dbuffer[2] = {tstart, tstop};
        fits_write_col(fptr, TDOUBLE, cdeadcor_time, 1, 1, 2, dbuffer, &status);
        float fbuffer[2] = {1., 1.};
        fits_write_col(fptr, TFLOAT, cdeadc, 1, 1, 2, fbuffer, &status);
        CHECK_STATUS_BREAK(status);

        // --- End of append DEADCOR extension ---

        // --- Beginning of append BADPIX extension ---

        headas_chat(3, "append BADPIX extension ...\n");

        // Create the BADPIX table.
        char badpix_extname[MAXMSG];
        sprintf(badpix_extname, "BADPIX%d", par.CCDNr);
        char *badpix_ttype[] = {"RAWX", "RAWY", "YEXTENT", "TYPE", "BADFLAG",
                                "TIMEMIN", "TIMEMAX", "PHAMIN", "PHAMAX", "PHAMED"};
        char *badpix_tform[] = {"I", "I", "I", "I", "I",
                                "D", "D", "I", "I", "E"};
        char *badpix_tunit[] = {"", "", "", "", "",
                                "", "", "", "", ""};
        fits_create_tbl(fptr, BINARY_TBL, 0, 10,
                        badpix_ttype, badpix_tform, badpix_tunit,
                        badpix_extname, &status);
        if (EXIT_SUCCESS != status) {
            SIXT_ERROR("could not create binary table for BADPIX extension");
            break;
        }

        // Insert header keywords.
        fits_update_key(fptr, TSTRING, "HDUCLASS", "OGIP", "", &status);
        fits_update_key(fptr, TSTRING, "HDUCLAS1", "BADPIX", "", &status);
        fits_update_key(fptr, TSTRING, "HDUCLAS2", "STANDARD", "", &status);
        CHECK_STATUS_BREAK(status);

        // Determine the individual column numbers.
        int cbadpix_rawx, cbadpix_rawy, cbadpix_yextent, cbadpix_type,
                cbadflag, ctimemin, ctimemax, cphamin, cphamax, cphamed;
        fits_get_colnum(fptr, CASEINSEN, "RAWX", &cbadpix_rawx, &status);
        fits_get_colnum(fptr, CASEINSEN, "RAWY", &cbadpix_rawy, &status);
        fits_get_colnum(fptr, CASEINSEN, "YEXTENT", &cbadpix_yextent, &status);
        fits_get_colnum(fptr, CASEINSEN, "TYPE", &cbadpix_type, &status);
        fits_get_colnum(fptr, CASEINSEN, "BADFLAG", &cbadflag, &status);
        fits_get_colnum(fptr, CASEINSEN, "TIMEMIN", &ctimemin, &status);
        fits_get_colnum(fptr, CASEINSEN, "TIMEMAX", &ctimemax, &status);
        fits_get_colnum(fptr, CASEINSEN, "PHAMIN", &cphamin, &status);
        fits_get_colnum(fptr, CASEINSEN, "PHAMAX", &cphamax, &status);
        fits_get_colnum(fptr, CASEINSEN, "PHAMED", &cphamed, &status);
        CHECK_STATUS_BREAK(status);

        // --- End of append BADPIX extension ---

        // --- Beginning of append CORRATT extension ---

        headas_chat(3, "append CORRATT extension ...\n");

        // Set up the Attitude.
        if (par.Attitude == NULL) {
            // Set up a simple pointing attitude.
            ac = getPointingAttitude(mjdref, tstart, tstop,
                                     par.RA * M_PI / 180., par.Dec * M_PI / 180., par.rollangle * M_PI / 180., &status);
            CHECK_STATUS_BREAK(status);

        } else {
            // Load the attitude from the given file.
            ac = loadAttitude(par.Attitude, &status);
            CHECK_STATUS_BREAK(status);

            // Check if the required time interval for the simulation
            // is a subset of the period covered by the attitude file.
            checkAttitudeTimeCoverage(ac, mjdref, tstart, tstop,
                                      &status);
            CHECK_STATUS_BREAK(status);
        }
        // END of setting up the attitude.


        if (NULL != ac) {
            // Create the CORRATT table.
            char corratt_extname[MAXMSG];
            sprintf(corratt_extname, "CORRATT%d", par.CCDNr);
            char *corratt_ttype[] = {"TIME", "RA", "DEC", "ROLL"};
            char *corratt_tform[] = {"D", "D", "D", "D"};
            char *corratt_tunit[] = {"", "deg", "deg", "deg"};
            fits_create_tbl(fptr, BINARY_TBL, 0, 4,
                            corratt_ttype, corratt_tform, corratt_tunit,
                            corratt_extname, &status);
            if (EXIT_SUCCESS != status) {
                SIXT_ERROR("could not create binary table for CORRATT extension");
                break;
            }

            // Insert header keywords.
            fits_update_key(fptr, TSTRING, "HDUCLASS", "OGIP", "", &status);
            fits_update_key(fptr, TSTRING, "HDUCLAS1", "TEMPORALDATA", "", &status);
            fits_update_key(fptr, TSTRING, "HDUCLAS2", "ASPECT", "", &status);
            CHECK_STATUS_BREAK(status);

            // Determine the individual column numbers.
            int ccorratt_time, ccorratt_ra, ccorratt_dec, croll;
            fits_get_colnum(fptr, CASEINSEN, "TIME", &ccorratt_time, &status);
            fits_get_colnum(fptr, CASEINSEN, "RA", &ccorratt_ra, &status);
            fits_get_colnum(fptr, CASEINSEN, "DEC", &ccorratt_dec, &status);
            fits_get_colnum(fptr, CASEINSEN, "ROLL", &croll, &status);
            CHECK_STATUS_BREAK(status);

            // Determine the rotation of the CCD from the keyword in the event file.
            float ccdrotation;
            fits_read_key(elf->fptr, TFLOAT, "CCDROTA", &ccdrotation, comment, &status);
            if (EXIT_SUCCESS != status) {
                SIXT_ERROR("failed reading keyword CCDROTA in input file");
                break;
            }

            // Insert the data.
            // Number of rows in the output attitude extension.
            long nrows = 0;
            // Loop over all intervals in the GTI collection.
            int gtibin = 0;
            do {
                // Currently regarded interval.
                double t0, t1;

                // Determine the currently regarded interval.
                if (NULL == gti) {
                    t0 = tstart;
                    t1 = tstop;
                } else {
                    t0 = gti->start[gtibin];
                    t1 = gti->stop[gtibin];
                }

                // Note that the attitude is stored in steps of 1s
                // according to the official event file format definition.
                double currtime;
                for (currtime = t0; currtime <= t1; currtime += 1.0) {
                    Vector nx, ny, nz;
                    getTelescopeAxes(ac, &nx, &ny, &nz, currtime, &status);
                    CHECK_STATUS_BREAK(status);

                    double ra, dec;
                    calculate_ra_dec(nz, &ra, &dec);
                    ra *= 180. / M_PI;
                    dec *= 180. / M_PI;

                    // Determine the roll angle.
                    Vector x1, x2;
                    Vector z = {0.0, 0.0, 1.0};
                    x2 = normalize_vector(vector_product(nz, z));
                    x1 = vector_product(x2, nz);

                    float rollangle =
                            atan2(scalar_product(&nx, &x2), scalar_product(&nx, &x1)) * 180. / M_PI;

                    // Apply the rotation angle of the CCD.
                    rollangle += ccdrotation;

                    // Store the data in the file.
                    nrows++;
                    fits_write_col(fptr, TDOUBLE, ccorratt_time, nrows, 1, 1, &currtime, &status);
                    fits_write_col(fptr, TDOUBLE, ccorratt_ra, nrows, 1, 1, &ra, &status);
                    fits_write_col(fptr, TDOUBLE, ccorratt_dec, nrows, 1, 1, &dec, &status);
                    fits_write_col(fptr, TFLOAT, croll, nrows, 1, 1, &rollangle, &status);
                    CHECK_STATUS_BREAK(status);
                }
                CHECK_STATUS_BREAK(status);

                // Proceed to the next GTI interval.
                if (NULL != gti) {
                    gtibin++;
                    if (gtibin >= gti->ngti) break;
                }

            } while (NULL != gti);
            CHECK_STATUS_BREAK(status);
            // End of loop over the individual GTI intervals.
        }

        // --- End of append CORRATT extension ---

        // Append a check sum to the header of the event extension.
        int hdutype = 0;
        fits_movabs_hdu(fptr, 2, &hdutype, &status);
        fits_write_chksum(fptr, &status);
        CHECK_STATUS_BREAK(status);

    } while (0); // END of the error handling loop.


    // --- Cleaning up ---
    headas_chat(3, "cleaning up ...\n");

    // Close the files.
    freeEventFile(&elf, &status);
    if (NULL != fptr) fits_close_file(fptr, &status);

    // Release memory.
    wcsfree(&wcs);
    if (NULL != headerstr) free(headerstr);
    freeGTI(&gti);
    freeAttitude(&ac);

    if (EXIT_SUCCESS == status) {
        headas_chat(3, "finished successfully!\n\n");
        return (EXIT_SUCCESS);
    } else {
        return (EXIT_FAILURE);
    }
}


int getpar(struct Parameters *const par) {
    // String input buffer.
    char *sbuffer = NULL;

    // Error status.
    int status = EXIT_SUCCESS;

    // check if any obsolete keywords are given
    sixt_check_obsolete_keyword(&status);
    CHECK_STATUS_RET(status, EXIT_FAILURE);

    status = ape_trad_query_file_name("EvtFile", &sbuffer);
    if (EXIT_SUCCESS != status) {
        SIXT_ERROR("failed reading the name of the input pattern list");
        return (status);
    }
    strcpy(par->EvtFile, sbuffer);
    free(sbuffer);

    status = ape_trad_query_file_name("eroEvtFile", &sbuffer);
    if (EXIT_SUCCESS != status) {
        SIXT_ERROR("failed reading the name of the output event list");
        return (status);
    }
    strcpy(par->eroEvtFile, sbuffer);
    free(sbuffer);

    status = ape_trad_query_int("CCDNr", &par->CCDNr);
    if (EXIT_SUCCESS != status) {
        SIXT_ERROR("failed reading the CCDNr parameter");
        return (status);
    }

    status = ape_trad_query_string("Projection", &sbuffer);
    if (EXIT_SUCCESS != status) {
        SIXT_ERROR("failed reading the name of the projection type");
        return (status);
    }
    strcpy(par->Projection, sbuffer);
    free(sbuffer);

    status = ape_trad_query_float("RefRA", &par->RefRA);
    if (EXIT_SUCCESS != status) {
        SIXT_ERROR("failed reading RefRA");
        return (status);
    }

    status = ape_trad_query_float("RefDec", &par->RefDec);
    if (EXIT_SUCCESS != status) {
        SIXT_ERROR("failed reading RefDEC");
        return (status);
    }

    query_simput_parameter_file_name("Attitude", &(par->Attitude), &status);

    // only load RA,Dec if Attitude is not given
    if (par->Attitude) {
        // set to default values
        par->RA = 0.0;
        par->Dec = 0.0;
        par->rollangle = 0.0;
        headas_chat(3, "using Attitude File: %s \n", par->Attitude);
    } else {
        query_simput_parameter_float("RA", &(par->RA), &status);
        query_simput_parameter_float("Dec", &(par->Dec), &status);
        query_simput_parameter_float("rollangle", &(par->rollangle), &status);
    }


    status = ape_trad_query_bool("clobber", &par->clobber);
    if (EXIT_SUCCESS != status) {
        SIXT_ERROR("failed reading the clobber parameter");
        return (status);
    }

    return (status);
}
