/******************************************************************************
 *   File name:                                                               *
 *                                                                            *
 * Description:                                                               *
 *                                                                            *
 *    Language: C or C++                                                      *
 *                                                                            *
 *      Author: Bryan Irby, for HEASARC/GSFC/NASA                             *
 *              Wrapper to "Mersenne Twister" (mt) routines provided by       *
 *              Bob Wiegand                                                   *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/
#ifndef HEADAS_RAND_H
#define HEADAS_RAND_H

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
/******************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
extern "C" {
#endif

  /****************************************************************************
   * Constants.                                                               *
   ****************************************************************************/
  /****************************************************************************/

  /****************************************************************************
   * Type declarations/definitions.                                           *
   ****************************************************************************/
	struct mt_state_t;
	typedef struct mt_state_t HDmt_state;

  /****************************************************************************/

  /****************************************************************************
   * Global variable forward declarations.                                    *
   ****************************************************************************/
  /****************************************************************************/

  /****************************************************************************
   * Function declarations.                                                   *
   ****************************************************************************/

	/* HDmt_srand returns a state which must be freed by the caller
	 * using HDmt_destroy_state */
	HDmt_state *HDmt_srand (unsigned long int s);
	void HDmt_destroy_state (HDmt_state *state);
	unsigned long int HDmt_rand (HDmt_state *state);
	double HDmt_drand (HDmt_state *state);

        /* Versions of routines which do not have the HDmt_state
	   pointer as an argument */
        unsigned long int HDmtRand ();
        double HDmtDrand ();
        void HDmtInit (unsigned long int s);
        void HDmtFree ();

  /****************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
}
#endif

#endif

/******************************************************************************
 * $Log: headas_rand.h,v $
 * Revision 1.2  2007/04/19 22:23:12  kaa
 * Added prototypes for new HDmt* routines.
 *
 * Revision 1.1  2005/03/23 17:24:31  irby
 * Add routines for temporary file and pseudo-random number generation.
 *
 ******************************************************************************/
