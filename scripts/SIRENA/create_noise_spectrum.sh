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
#==============================================
# This script runs SIRENA utilities to create a NOISE SPECTRUM (step 1) )
# ==============================================
#

# 1.1) define simulation parameters
#----------------------------------
# define  XML files 
xmldirXF=${XIFUSIM}/share/xifusim/instruments
# create XML XF file for noise with <Trigger  model=TriggerNoise>
xmlfile=${xmldirXF}/8pix_nobbfb.xml
xmlfilenoise=${xmldirXF}/8pix_noise_nobbfb.xml

if [ ! -e ${xmlfilenoise} ]; then
    cp  ${xmlfile} ${xmlfilenoise}
    sed -i '/TriggerDiff/c\<Trigger model="TriggerNoise" filename="pars_8.fits" hduname="TrigNoise"/>' ${xmlfilenoise}
fi

samprate=156250  # sampling frequency
pulse_length=8192 # length for reconstruction (samples)
rlength=8200 # record length (samples) >= pulse_length
preBuffer=1000 # Number of samples before the impact in the trigger
nintervals=2000 #noise intervals for the noise spectrum
time=`python -c "print(${nintervals} * ${pulse_length} / ${samprate})"`  # noise simulation time (s)
acbias="yes"

# FILENAMES
#---------------
noise_imp="mynoise.piximpact"
noise_sim="mynoise_sim.fits"
noise_spec="mynoise.spec"

# 1.2) Create piximpact list to simulate noise streams
#-----------------------------------------------
tesconstpileup PixImpList=${noise_imp} \
                                XMLFile=${xmlfile}\
                                tstop=${time}\
                                offset=0\
                                energy=0.\
                                pulseDistance=1\
                                TriggerSize=${rlength}\
                                sample_freq=${samprate}\
                                preBufferSize=0 \
                                clobber=yes

# 1.3) Simulate noise stream (same conditions as those used to simulate data)
#-----------------------------------------------------------------------
# (if data have been simulated with TESSIM)
pixfile='file:${SIXTE}/share/sixte/instruments/athena-xifu/newpix_LPA75um.fits'
tessim PixID=1 \
                PixImpList=${noise_imp} \
                Streamfile=${noise_sim} \
                tstart=0. \
                tstop=${time} \
                triggertype=noise \
                triggersize=${rlength} \
                prebuffer=0 \
                PixType=${pixfile} \
                acbias=${acbias}

# (if data have been simulated with XIFUSIM)

xifusim PixImpList=${noise_imp}\
                                Streamfile=${noise_sim}\
                                tstop=${time}\
                                acbias=${acbias}\
                                XMLfilename=${xmlfilenoise}\
                                trig_reclength=${rlength}\
                                trig_n_pre=0 \
                                trig_thresh=0.\
                                simnoise=y
                                
# 1.4) Build noise spectrum
#--------------------------
gennoisespec inFile=${noise_sim} \
                            outFile=${noise_spec} \
                            nintervals=${nintervals} \
                            pulse_length=${pulse_length} \
                            intervalMinSamples=${pulse_length}\

                            
                            
