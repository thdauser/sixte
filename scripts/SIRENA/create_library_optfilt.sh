#! /bin/bash

# % This file is part of SIXTE/SIRENA software.
# %
# % SIXTE is free software: you can redistribute it and/or modify it
# % under the terms of the GNU General Public License as published by
# % the Free Software Foundation, either version 3 of the License, or
# % any later version.
#
# %  SIXTE is distributed in the hope that it will be useful,
# %  but WITHOUT ANY WARRANTY; without even the implied warranty of
# %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# %  GNU General Public License for more details.
#
# %  For a copy of the GNU General Public License see
# %  <http://www.gnu.org/licenses/>.
#
# %  Copyright 2020:  SIRENA scripts have been developed by the
# %  INSTITUTO DE FISICA DE  CANTABRIA (CSIC-UC) with funding 
# %  from the Spanish Ministry of Science and  Innovation (MICIN) 
# %  under project RTI2018-096686-B-C21
#
# ==================================================
#
#
# SIRENA procedure to reconstruct TESSIM or XIFUSIM-simulated data
# https://sirena.readthedocs.io/
#
# 1. Create Noise spectrum (create_noise_spectrum.sh)
# 2. Simulate pulses for Pulse template (simulate_pulses_template.sh)
# 3. Create library with OPTIMAL FILTER (create_library_optfilt.sh)
# 4. reconstruct xifusim/tessim-simulated data example (example_reconstruction.sh)
#
#
#
#
# ====================================================
# This script runs SIRENA utilities to create library with OPTIMAL FILTER (step 3) )
# ====================================================
#

# 3.1) define simulation parameters
#----------------------------------
largeFilter=8192 # largest filter to be created (samples); power of 2 smaller filters will also be created
tempeV=6000. # (energy of optimal filter (eV)

# FILENAMES
#---------------
noise_spec="mynoise.spec"
temp_sim="myTemp.fits"
libfile="myLib.fits"
xmlfile="${XIFUSIM}/share/xifusim/instruments/1pix_lpa2.5a_fll_SIRENAintegration.xml"

# Use noise spectrum and pulses for template and create a library of optimal filters
#echo "teslib Recordfile=${temp_sim} TesEventFile=myLibEvents.fits monoenergy=${tempeV} LibraryFile=${libfile} largeFilter=${largeFilter} NoiseFile=${noise_spec} XMLFile=${xmlfile}"
teslib Recordfile=${temp_sim} \
                                    TesEventFile=myLibEvents.fits\
                                    monoenergy=${tempeV} \
                                    LibraryFile=${libfile} \
                                    largeFilter=${largeFilter} \
                                    NoiseFile=${noise_spec} \
                                    XMLFile=${xmlfile} \
                                    preBuffer=n
 
                            
                                
                            
