version [2.6.2]
  - new tool: sixte_arfgen
    * ARF generator to calculate off-axis ARFs for any location and attitude.
  - updates internal handling of RMFs

version [2.6.1]
  - ero_calevents: Make SIXTE output consistent with eSASS' evtool/srctool
  - makespec
    * Add GTIFile option to produce multiple lightcurves from one eventfile
    * Add RMF consistency check
  - Add CMake support
  - Add exposure_map end-to-end test
  - Set clobber=yes as default
  - Multiple bug fixes

version [2.6.0]
  - changes input parameters of exposure_map tool
    * the exposure_map tool now uses WCS keyrecords to specify the grid on
      which the exposure map is calculated instead of an RA/Dec
      interval (previous input). See simulator manual (v1.3.10, p. 60), for
      more details.
    * for backward compatibility the old interface still works, but will not be
      actively supported anymore in the future.
  - adds new tool radec2xy
    * calculates backprojected X,Y instrument coordinates from RA/Dec sky
      coordinates and adds these coordinates to given event file.

version [2.5.12]
  - fix azimuthal angle range in get_psf_pos
  - fix attgen_dither tool (MJDREF, ALIGNMEN keywords not set in 2.5.11)

version [2.5.11]
  - adds warning if no photons were imaged
  - fixes crash upon writing INSTRUME/TELESCOP keys to GTI
    extension
  - fixes use of wrong RMF for PHA background energies if using
    an Advanced Detector XML file
  - updates erosim
    * erosim can now be given a comma separated list of filenames
      for each Simput parameter
  - updates SIRENA
    * 0-padding (8192-length filters filled in with 0's) and
      filters with preBuffer

version [2.5.10]
 - adds new tool attgen_dither which creates a Lissajous attitude
 - adds roll angle parameter for pointing observations via command line
 - adds missing std header keys to conform HEASARC
 - adjusts parameter input for runsixt. If an xml file is given,
   the parameters "Mission", "Instrument" and "Mode" are not queried anymore
 - SIRENA
   * fixes several bugs reagarding writing/reading HEADER
   * testriggerfile now compatible with RECORDS and
     TESRECORDS HDU
   * improves handling is no maximum is found for lags
   * improves list of input files handling

version [2.5.9]
 - fixes bug in last version regarding pha2pi correction
 - updates SIRENA to conform new file format of xifusim

version [2.5.8]
 - update & improve pha2pi correction handling
   * detector xml now may include a PIRMF or SPECARF
     entry to specify the pha2pi corrected RMF and the
     calibrated ARF
   * adjust makespec to automatically use PIRMF and
     SPECARF if available

version [2.5.7]
 - updates copyright and license notices
 - fixes bug in reading older vignetting files
 - changes vignetting handling (energies are interpolated linearly)
 - fixes bug in xmm timing mode

version [2.5.6]
 - ero_calevents: CORRATT Attitude changes reverted, given now
   in the Sixte system (in accordance with eSASS)
 - tessim: discard empty records at end ofsimulation
 - fix pha2pi tool crashing without message

version [2.5.5]
 - fixes bug in lib installation for Mac

version [2.5.4]
 - adds compatibitlity with new, April 2019 eRO Attitude
 - fixes Event Mode for NuSTAR simulations

version [2.5.3]
 - ero_calevents takes coord. trans. into account for attitude
   (CORRATT extension correctly written in eRO system now)
 - ero_calevents changes RAWX orientatation to eRO system
   (needs XML file update to srg-1.6.0)

version [2.5.2]
  - Pha2Pi correction automatically works with MC RMF (provided
    information is given in the XML file)
  - phogen tool warns if pseudo RNG is used (just meant for testing)
  - fixes minor bug in "psfgen" for uneven number of pixels

version [2.5.1]
  - fixes minor input reading "ero_vis"
  - updates display of e2e test setup
