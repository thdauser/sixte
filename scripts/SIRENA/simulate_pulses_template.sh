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
#==================================================
# This script runs SIRENA utilities to simulate pulses to build a template (step 2) )
# ==================================================
#
#
# 2.1)  define simulation parameters
#----------------------------------

samprate=156250  # sampling frequency
pulse_length=8192 # length for reconstruction (samples)
pulse_distance=20000 # separation between consecutive isolated pulses (samples) >> pulse_length
recordSeparation=20000 # distance (samples ) between consecutive records
preBuffer=1000 # free samples before pulse triggering (samples)
TsizeTC=`python -c "print (${preBuffer}+${pulse_distance} + ${recordSeparation}+1000)"` # simulated double-pulse record length (samples)
TsizeSM=`python -c "print (${preBuffer} + ${pulse_length} + 1000)"` # trigger (isolated pulse) record length 
npulses=20000 # number of simulated pulses to build the template
time=`python -c "print (${npulses} / 2. * ${TsizeTC}/ ${samprate})"` #  simulation time (s)
tempKeV=6. # (energy for template and then for optimal filter (keV)

# FILENAMES
#---------------
temp_imp="myTemp.piximpact"
temp_sim="myTemp.fits"

# 2.2) Create piximpact list to simulate isolated pulses (tesconstpileup must be used if jitter wants to be simulated, i.e. offsets in the arrival time of photons)
#-----------------------------------------------------------------
xmldirSX=${XIFUSIM}/share/xifusim/instruments
xmlfileSX=${xmldirSX}/1pix_lpa2.5a_fll.xml
# 
tesconstpileup PixImpList=${temp_imp} \
                                XMLFile=${xmlfileSX} \
                                timezero=3.E-7 \
                                tstop=${time} \
                                offset=-1 \
                                energy=${tempKeV} \
                                pulseDistance=${pulse_distance} \
                                TriggerSize=${TsizeTC} \
                                sample_freq=${samprate}\
                                clobber=yes

# 2.3) Simulate isolated pulses (same conditions as those used to simulate data)
#----------------------------------------------------------------------------
# (if data have been simulated with TESSIM)
pixfile='file:${SIXTE}/share/sixte/instruments/athena-xifu/newpix_LPA75um.fits'
triggerType='diff:3:100:${pulse_length}'
tessim PixID=1 \
                PixImpList=${temp_imp} \
                Streamfile=${temp_sim} \
                tstart=0. \
                tstop=${time} \
                triggertype=${triggerType}\
                triggerSize=${TsizeSM} \
                prebuffer=${preBuffer}\
                PixType=${pixfile} \
                acbias=yes

# (if data have been simulated with XIFUSIM)
# define  XML files
xmldirXF=${XIFUSIM}/share/xifusim/instruments
xmlfile=${xmldirXF}/1pix_lpa2.5a_fll.xml

xifusim PixImpList=${temp_imp} \
                                Streamfile=${temp_sim} \
                                tstop=${time} \
                                acbias=no \
                                XMLfilename=${xmlfile} \
                                trig_reclength=${TsizeSM} \
                                trig_n_pre=${preBuffer} \
                                trig_n_suppress=${pulse_length} \
                                simnoise=y \
                                clobber=yes
                                
                                
                                
                            
