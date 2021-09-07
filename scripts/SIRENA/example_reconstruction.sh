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
# ==================================================
# This script runs SIRENA utilities to reconstruct xifusim-simulated data  (step 4) )
# ==================================================
#

# SIXTE XML files *** see note in README file ***
xmldirSX=${XIFUSIM}/share/xifusim/instruments/
xmlfileSX=${xmldirSX}/1pix_lpa2.5a_fll.xml # any file with grading information (pixel info is not relevant here)

tempeV=6000. # (energy of optimal filter (eV)

data_sim="myData.fits"
evt_recons="myEvents.fits"
libfile="myLib.fits"


# Reconstruct data using optimal filter (for different reconstruction methods have a look to SIRENA documentation)
# https://sirena.readthedocs.io/
tesreconstruction Recordfile=${data_sim} \
                                    TesEventFile=${evt_recons} \
                                    monoenergy=${tempeV} \
                                    LibraryFile=${libfile} \
                                    opmode=1 \
                                    XMLFile=${xmlfileSX}

                            
                                
                            
