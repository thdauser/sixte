#ifndef WFIDETECTOR_H
#define WFIDETECTOR_H 1

#include "sixt.h"
#include "genericdetector.h"
#include "squarepixels.h"
#include "eventfile.h"
#include "wfieventfile.h"
#include "impact.h"


/** Model for the WFI detector on IXO.  This data structure contains
    the data required for the simulation of the WFI.  It inherits some
    properties of the GenericDetector and SquarePixels data
    structures.  The WFI-specific properties are mainly related to the
    particular readout mode with 1 or 2 readout lines. (The WFI
    prototype used for comparison of simulation data with measurements
    from the laboratory has only 1 readout line.)  The WFIDetector
    data structure can be initialized by calling the initWFIDetector()
    function with a WFIDetectorParameters data structure containing
    the desired setup.  A new photon Impact can be added to the
    WFIDetector array by the function addImpact2WFIDetector().
    Finally after the simulation when the data structure is not
    required any more, the cleanupWFIDetector() routine should be
    called to release allocated memory and close open file
    connections. */
typedef struct {

  /** Generic Detector properties like, e.g., the detector response. */
  GenericDetector generic;
  /** Array of square pixels. */
  SquarePixels pixels;


  /** Number of readout directions. In the current implementation
      this can be either 1 or 2. */
  int readout_directions;

  /** Time required to read out one line of the WIF detector (given in
      [s]). This is the time interval between the readout of two
      subsequent detector lines. */
  double line_readout_time;

  /** Time required to clear a line of pixels on the detector. (given
      in [s]). In the current implementation the clear time is a part
      of the readout time. */
  double line_clear_time;

  /** Current readout time (given in [s]). This is used as a counter
      for the current readout time of the WFIDetector, i.e., the point
      of time, when the last line readout has been performed. */
  double readout_time;

  /** Current readout lines of the WFI detector. Gives the indices of
      the lines that have been read out recently. */
  int readout_lines[2];

  /** Number of the current frame. The currently active readout frame
      that will be assigned to the events in the event file. */
  long frame; 

  /** Output event list. The events read out from the detector array
      are written to this event file that must have the WFI-specific
      format. */
  WFIEventFile eventlist;

} WFIDetector;


/** Parameters of the WFIDetector model.  This data structure contains
    the parameters for setting up the WFIDetector data structure.  It
    is used as input for the initWFIDetector() routine.  Besides some
    generic detector parameters especially the readout mode and the
    event file for data output have to be defined.  For documentation
    of the inidividual parameters see WFIDetector. */
struct WFIDetectorParameters {
  struct GenericDetectorParameters generic;
  struct SquarePixelsParameters pixels;

  double line_readout_time, line_clear_time;
  int readout_directions;
  char* eventlist_filename;
  char* eventlist_template;

  double t0;
};


////////////////////////////////////////////////////////


/** Set up the configuration of a WFIDetector. The routine is
    responsible to set up the initial the WFIDetector configuration
    which is given in the WFIDetectorParameters data structure. It
    has to take care of allocating the required memory for the pixel
    array and to create an event file for the output of the measured
    data. For some of these tasks it simply calls the init routines
    of the underlying data structures. */
int initWFIDetector(WFIDetector*, struct WFIDetectorParameters*);

/** Clean up the WFIDetector data structure. This routine should be
    called if the WFIDetector data structure is not required any more.
    It takes care of releasing allocated memory and closes open file
    connections. If applicable it calls clean-up routines of
    underlying data structures. */
int cleanupWFIDetector(WFIDetector*);

/** Check out whether the WFI detector needs to be read out.
    Regarding the given time and the current readout time of the
    WFIDetector this function checks, whether and which lines have to
    be read out. In case it is necessary it calls the
    readoutLinesWFIDetector() function, which actually performs the
    readout. */
int checkReadoutWFIDetector(WFIDetector*, double time);

/** Read out the charge from the currently active readout lines (one
    or two) of the WFI detector. This function implements the
    different readout modes of the WFIDetector model. According to
    the chosen number of readout lines it performs the readout of the
    currently active detector lines. The events read from the
    detector pixels are stored in the output even file. */
inline int readoutLinesWFIDetector(WFIDetector*);

/** Add a photon impact to the WFIDetector pixel array. This is the
    standard routine to be called for the simulation of the
    WFIDetector. For a new photon incident on the detector this
    routine determines the resulting generated charge from the
    detector response, sums the charge to the affected detector pixels
    (one or several according to the selected split model), and checks
    whether the detector has to be read out according to the chosen
    readout mode. In the latter case the right routines for the
    detector readout are called and the events are stored in the
    output event file. */
int addImpact2WFIDetector(WFIDetector*, Impact*);

/** Check the detector sensitivity in the given detector line at the
    specified time. This routine checks, whether the WFIDetector is
    sensitive in line 'y' at the given time. According to the
    currently implemented readout strategy the WFIDetector is
    insensitive to new photons in an interval immediately after the
    line readout. This interval is called 'line_clear_time'. This
    function checks whether the designated line was just read out and
    the time lies within the clear time.  It returns 1 if the detector
    is sensitive in the line 'y' at the specified time. Otherwise the
    return value is '0'. */
inline int WFIDetectorIsSensitive(int y, WFIDetector*, double time);


#endif /* WFIDETECTOR_H */
