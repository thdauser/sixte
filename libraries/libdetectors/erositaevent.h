#ifndef EROSITAEVENT_H
#define EROSITAEVENT_H 1


/** eROSITA specific event. Contains the elements specifying an
    eROSITA event, i.e., the typcial column entries in the FITS event
    file. */
typedef struct {
  double time;
  long pha;
  float energy; /**< Event energy [eV] */
  /** RAWX and RAWY  pixel coordinates, starting at 0. */
  int xi, yi; 
  long frame;

  /** Number of the CCD. For eROSITA the CCDs are numbered from 1 to
      7. */
  char ccdnr;

  /** Information about split events. */
  int pat_typ, pat_inf;


  /** Right ascension and declination [degree]. */
  double ra, dec; 
  /** Integer sky pixel coordinates. */
  long sky_xi, sky_yi; 

} eROSITAEvent;


#endif /* EROSITAEVENT_H */

