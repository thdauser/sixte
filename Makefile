##/////////////////////////////////////////////////////////////////////////////
## General Information:
## file:             Makefile
## project:          eROSITA - NRTA - simulation
## author:           Christian Schmid
## date:             12/2007
##/////////////////////////////////////////////////////////////////////////////
## Description:
## This is the Makefile of the NRTA-simulation for the eROSITA mission.
## It builds several runnables (some actual software and some testing software)
## and librarys with commonly used functions/routines/structures/constants.
##/////////////////////////////////////////////////////////////////////////////



##============================================================================
##  definitions and flags
##============================================================================

#   GNU C compiler:
CC=gcc

# "-O3" means best optimization of the code
#OPT=-O3

# "-g" creates debug information
# "-Wall" creates warnings, e.g. for unused variables
DEBUG=-g -W -Wall -O0
# -Wstrict-prototypes -Wmissing-prototypes -pedantic



#   PROGRAM NAMES:
# names of all programms:
# names of all own headas applications / ftools
OWN_HEATOOLS= 	conv_rosat2fits conv_psf2fits create_rnd_sctlg create_orbit \
		create_attitude create_rmf create_spectrum measurement plot_psf \
		plot_det_images plot_eventlist create_psf \
		scan_tes_events byte_stream
OWN_LIBRARIES= 	libfitsformats.a
ALL= 		$(OWN_HEATOOLS) $(OWN_LIBRARIES) \
		sources_in_fov orbitparams_from_pos compare_orbits \
		calcJ2perturbations create_htrs_psf htrs_pixel \
		test_tle_output test_fov test_distrndsources test_simulation \
		test_lightcurve test_fabs



#   Include Paths:
# HEADAS libraries:
LHPATH=-L$(FTOOLS)/lib
IHPATH=-I$(FTOOLS)/include
# own libraries:
#LOPATH = -L.
#IOPATH = -I.



#   Libraries:
#   CFITSIO-library
LIBFITS=-lcfitsio_3.12 -lm -lnsl
#   PIL-library (Parameter Interface Library)
LIBPIL=-lape_2.4.0 -lreadline
#   HEADAS libraries for creating HEATOOLS
LIBHEATOOLS=-lhdinit_2.5 -lhdutils_2.5 -lhdio_2.5 $(LIBFITS) $(LIBPIL)
#   PNG-library
LIBPNG=-lpng
#   GNU-scientific library
LIBGSL=-lgsl -lgslcblas


# Set the standard compiler/linker flags
CFLAGS=$(DEBUG) $(OPT) $(IHPATH)
LFLAGS=$(DEBUG) $(OPT) $(LHPATH)


##============================================================================
##  generating runables
##============================================================================

all: $(ALL)

sim_headas: $(OWN_HEATOOLS) 

conv_rosat2fits: conv_rosat2fits.o strftcpy.o
	$(CC) $(LFLAGS) -o conv_rosat2fits $^ $(LIBHEATOOLS)

conv_psf2fits: conv_psf2fits.o psf.o vector.o random.o photon.o detector.o \
	strftcpy.o event_list.o fits_pha.o
	$(CC) $(LFLAGS) -o conv_psf2fits $^ $(LIBHEATOOLS) $(LIBGSL)

create_htrs_psf: create_htrs_psf.o psf.o detector.o event_list.o fits_pha.o \
	photon.o random.o vector.o
	$(CC) $(LFLAGS) -o create_htrs_psf $^ $(LIBHEATOOLS) $(LIBGSL)

create_psf: create_psf.o psf.o random.o vector.o
	$(CC) $(LFLAGS) -o create_psf $^ $(LIBGSL) $(LIBHEATOOLS)

create_rnd_sctlg: create_rnd_sctlg.o fits_ctlg.o create_rnd_source.o strftcpy.o \
	random.o
	$(CC) $(LFLAGS) -o create_rnd_sctlg $^ $(LIBHEATOOLS)

create_orbit: create_orbit.o vector.o fits_ctlg.o
	$(CC) $(LFLAGS) -o create_orbit create_orbit.o vector.o fits_ctlg.o \
	$(LHPATH) $(LIBHEATOOLS)

create_attitude: create_attitude.o vector.o fits_ctlg.o fits_attitude.o
	$(CC) $(LFLAGS) -o create_attitude $^ $(LIBHEATOOLS)

create_rmf: create_rmf.o fits_pha.o
	$(CC) $(LFLAGS) -o create_rmf $^ $(LIBHEATOOLS)

create_spectrum: create_spectrum.o detector.o photon.o random.o fits_pha.o \
	event_list.o vector.o
	$(CC) $(LFLAGS) -o create_spectrum $^ $(LHPATH) $(LIBHEATOOLS) $(LIBGSL)

compare_orbits: compare_orbits.o vector.o fits_ctlg.o
	$(CC) $(LFLAGS) -o compare_orbits $^ $(LIBFITS)

calcJ2perturbations: calcJ2perturbations.o
	$(CC) $(LFLAGS) -o calcJ2perturbations $^ -lm

byte_stream: byte_stream.o libfitsformats.a
	$(CC) $(LFLAGS) -o byte_stream $^ $(LIBHEATOOLS)

htrs_pixel: htrs_pixel.o
	$(CC) $(LFLAGS) -o htrs_pixel $^ -lm

measurement: measurement.o vector.o psf.o detector.o sources.o photon.o orbatt.o \
	random.o spectrum.o check_fov.o fits_ctlg.o imglib.o strftcpy.o \
	measurement_array.o event_list.o fits_attitude.o fits_pha.o
	$(CC) $(LFLAGS) -o measurement $^ $(LIBPNG) $(LIBGSL) $(LIBHEATOOLS)

orbitparams_from_pos: orbitparams_from_pos.o vector.o fits_ctlg.o
	$(CC) $(LFLAGS) -o orbitparams_from_pos $^ $(LIBFITS)

plot_eventlist: plot_eventlist.o imglib.o event_list.o
	$(CC) $(LFLAGS) -o plot_eventlist $^ $(LIBHEATOOLS) $(LIBPNG)
# -lm

plot_det_images: plot_det_images.o imglib.o event_list.o
	$(CC) $(LFLAGS) -o plot_det_images $^ $(LIBHEATOOLS) $(LIBPNG)

#plot_spectrum: plot_spectrum.o detector.o random.o photon.o event_list.o \
	fits_pha.o
#	$(CC) $(LFLAGS) -o plot_spectrum $^ $(LIBHEATOOLS) $(LIBPNG) $(LIBGSL) -lm

plot_psf: plot_psf.o imglib.o strftcpy.o
	$(CC) $(LFLAGS) -o plot_psf $^ $(LHPATH) $(LIBHEATOOLS) $(LIBPNG)

random: random.o
	$(CC) $(LFLAGS) -o random $^ $(LIBHEATOOLS)

scan_tes_events: scan_tes_events.o event_list.o 
	$(CC) $(LFLAGS) -o scan_tes_events $^ $(LIBHEATOOLS) $(LIBGSL)

sources_in_fov: sources_in_fov.o check_fov.o fits_ctlg.o vector.o
	$(CC) $(LFLAGS) -o sources_in_fov $^ $(LIBFITS) $(LIBPIL)

test_tle_output: test_tle_output.o tle.o strftcpy.o
	$(CC) $(LFLAGS) -o test_tle_output $^ -lm

test_simulation: test_simulation.o vector.o
	$(CC) $(LFLAGS) -o test_simulation $^ -lm $(LHPATH) $(LIBPIL)

test_distrndsources: test_distrndsources.o create_rnd_source.o random.o
	$(CC) $(LFLAGS) -o test_distrndsources $^ $(LIBHEATOOLS) -lm

test_fov: test_fov.o check_fov.o vector.o
	$(CC) $(LFLAGS) -o test_fov $^ -lm

test_lightcurve: test_lightcurve.o vector.o fits_ctlg.o random.o strftcpy.o \
	fits_pha.o sources.o photon.o detector.o event_list.o
	$(CC) $(LFLAGS) -o test_lightcurve $^ $(LIBHEATOOLS) $(LIBGSL) -lm

test_fabs: test_fabs.o
	$(CC) $(LFLAGS) -o test_fabs $^ -lm


##============================================================================
##  generating libraries
##============================================================================


libfitsformats.a: event_list.o
	ar -r libfitsformats.a $^



##============================================================================
##  generating object files
##============================================================================

DEPENDFILE = make.dependencies
SRC = conv_rosat2fits.c conv_psf2fits.c create_psf.c event_list.c \
	fits_attitude.c \
	fits_pha.c create_rnd_source.c create_rnd_sctlg.c create_rmf.c \
	create_spectrum.c fits_ctlg.c sources_in_fov.c vector.c check_fov.c tle.c \
	create_orbit.c create_attitude.c orbitparams_from_pos.c compare_orbits.c \
	calcJ2perturbations.c measurement.c psf.c detector.c sources.c photon.c \
	orbatt.c spectrum.c split.c random.c strftcpy.c plot_eventlist.c \
	byte_stream.c \
	plot_det_images.c measurement_array.c scan_tes_events.c htrs_pixel.c \
	imglib.c plot_psf.c test_tle_output.c test_simulation.c create_htrs_psf.c \
	test_distrndsources.c test_fov.c test_lightcurve.c test_fabs.c

dep: $(SRC)
	$(CC) $(CFLAGS) -MM $(IHPATH) $(SRC) > $(DEPENDFILE)

-include $(DEPENDFILE)


##============================================================================
##  cleaning
##============================================================================

.PHONY: clean

clean:
	rm -f $(ALL) $(DEPENDFILE) *.o *.a *~


