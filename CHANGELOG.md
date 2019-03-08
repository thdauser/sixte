version [2.5.6]
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
