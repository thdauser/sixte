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

samprate=156250  # sampling frequency
rlength=8200 # record length (samples) >= High resolution length
nintervals=1000 # default noise intervals for the noise spectrum
time=`python -c "print(${nintervals} * ${rlength} / ${samprate})"`  # noise simulation time (s)

# FILENAMES
#---------------
noise_imp="mynoise.piximpact"
noise_sim="mynoise_sim.fits"
noise_spec="mynoise.spec"

# 1.2) Create piximpact list to simulate noise streams (both tesconstpileup or tesgenimpacts)
#----------------------------------------------------------------------------------------
#echo "tesgenimpacts PixImpList = ${noise_imp} opmode = const tstart = 0 tstop = ${time} EConst = 0. dtau = 1"
tesgenimpacts \
        PixImpList = ${noise_imp} \
        opmode = const \
        tstart = 0 \
        tstop = ${time} \
        EConst = 0. \
        dtau = 1
        
# 1.3) Simulate noise stream (same conditions as those used to simulate data)
#-----------------------------------------------------------------------
# 1.3.1) (if your data have been simulated with TESSIM)
#pixfile='file:${SIXTE}/share/sixte/instruments/athena-xifu/newpix_LPA75um.fits'
pixfile='file:mypixel_configuration.fits'
######################################################################################################################
# From SIXTE manual: Input parameters of the model are the heat capacity of the sensor, the properties of the thermal
# link, the operating point resistance and temperature, and depending on whether the TES is AC- or DC-biased, the ﬁlter
# inductance and parasitic resistance, or the shunt resistance and the eﬀective inductance of the readout circuit. These
# input parameters can either be given on the command line using the Parameter Interface Library (PIL) syntax, or can
# be read from a FITS ﬁle ('mypixel_configuration.fits') with a TESDATASTREAM extension whose header contains the
# parameters. Such a ﬁle can be obtained by running tessim with the propertiesonly=yes option.
######################################################################################################################
#echo "tessim PixID=1 PixImpList=${noise_imp} Streamfile=${noise_sim} tstart=0. tstop=${time} triggertype=noise triggersize=${rlength} prebuffer=0 PixType=${pixfile} acbias=yes"
tessim PixID=1 \
                PixImpList=${noise_imp} \
                Streamfile=${noise_sim} \
                tstart=0. \
                tstop=${time} \
                triggertype=noise \
                triggersize=${rlength} \
                prebuffer=0 \
                PixType=${pixfile} \
                acbias=yes

# 1.3.2) (if your data have been simulated with XIFUSIM)
# define  XML files 
xmldirXF=${XIFUSIM}/share/xifusim/instruments
# create XML XF file for noise with <Trigger  model=TriggerNoise>
xmlfile=${xmldirXF}/1pix_lpa2.5a_fll.xml
xmlfilenoise=${xmldirXF}/1pix_lpa2.5a_fll_noise.xml

# create a noise-trigger XML file if not already present
if [ ! -e ${xmlfilenoise} ]; then
    cp  ${xmlfile} ${xmlfilenoise}
    sed -i '/TriggerDiff/c\<Trigger model="TriggerNoise" filename="pars_8.fits" hduname="TrigNoise"/>' ${xmlfilenoise}
fi


#echo "xifusim PixImpList=${noise_imp} Streamfile=${noise_sim} tstop=${time} acbias=no XMLfilename=${xmlfilenoise} trig_reclength=${rlength} trig_n_pre=0 trig_thresh=0. simnoise=y"
xifusim PixImpList=${noise_imp}\
                                Streamfile=${noise_sim}\
                                tstop=${time}\
                                acbias=no\
                                XMLfilename=${xmlfilenoise}\
                                trig_reclength=${rlength}\
                                trig_n_pre=0 \
                                trig_thresh=0.\
                                simnoise=y \
                                writeEP=0	# EP is not going to run in XIFUSIM
                                
# 1.4) Build noise spectrum
#--------------------------
gennoisespec inFile=${noise_sim} \
                            outFile=${noise_spec} 
                            
            

                            
                            
