#ifndef EROSITAEVENT_H
#define EROSITAEVENT_H 1


/** eROSITA specific event. */
typedef struct {
  double time;
  long pha;
  float energy; /**< Event energy [eV] */
  int xi, yi; /**< Pixel coordinates starting at 0. */
  long frame;

  /** Number of the CCD.
   * For eROSITA the CCDs are numbered from 1 to 7. */
  char ccdnr;

  double ra, dec; /**< Right ascension and declination [degree]. */
  long sky_xi, sky_yi; /**< Integer sky pixel coordinates. */
} eROSITAEvent;


#endif /* EROSITAEVENT_H */

