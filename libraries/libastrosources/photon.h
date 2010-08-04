#ifndef PHOTON_H
#define PHOTON_H 1

#include "sixt.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif

#include "vector.h"
#include "sixt_random.h"
#include "pointsources.h"
#include "extendedsources.h"


////////////////////////////////////////////////////////////////////////
// Definitions.
////////////////////////////////////////////////////////////////////////


// Counter for the number of entirely generated photons.
long photon_counter;


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Contains all information about a single photon in the sky. */
typedef struct {
  double time;  /**< Real time, when the photon is falling on the detector (in [s]). */
  float energy; /**< Photon energy in [keV]. */

  /** Right ascension and declination of photon position [rad]. */
  double ra, dec; 

  // REMOVE (obsolete)
  Vector direction; // direction from which the photon originates 
                    // (source direction)
} Photon;


// Structure containing a photon and a pointer to the next photon in the 
// time-ordered photon list.
struct PhotonOrderedListEntry {
  Photon photon; 
  struct PhotonOrderedListEntry *next;  // pointer to the next entry
};


/** Entry in the binary tree that stores the generated photons. */
struct PhotonBinaryTreeEntry {
  Photon photon; /**< Photon data. */

  struct PhotonBinaryTreeEntry* sptr; /**< Pointer to entry with smaller time value. */
  struct PhotonBinaryTreeEntry* gptr; /**< Pointer to entry with greater time value. */
};


//////////////////////////////////////////////////////////////////////////
//   Function declarations
//////////////////////////////////////////////////////////////////////////


/** Create photons for a particular PointSource. Take into account
    the particular rate specified by the light curve and add the
    photons to the time ordered photon list. The return value is the
    value of the error status variable. */
int create_PointSourcePhotons(PointSource* ps, 
			      double time, double dt,
			      struct PhotonOrderedListEntry** pl, 
			      struct RMF*);

/** Create photons for a particular ExtendedSource. Take into account
    the particular rate specified by the light curve and add the
    photons to the time ordered photon list. The return value is the
    value of the error status variable. */
int create_ExtendedSourcePhotons(ExtendedSource* es, 
				 double time, double dt,
				 struct PhotonOrderedListEntry** pl, 
				 struct RMF*);

/** Create random photon energy. The routine randomly determines a
    photon energy according to the given source spectrum.  The
    spectrum must be in [PHA channels].  The function returns the
    photon energy in [keV]. */
float photon_energy(Spectrum* spectrum, struct RMF* rmf);

/** Inserts a new photon into the time-ordered photon list.
 * If the list does not exist so far, the routine creates a new list.
 * The return value is the error status. */
int insert_Photon2TimeOrderedList
(/** Address of the pointer to the absolutely first entry in the time-ordered list. 
  * The pointer might be  NULL, if the list is empty. The routine might modify the pointer. */
 struct PhotonOrderedListEntry** first,
 /** Address of the pointer to the entry in the time-ordered list where the insert 
  * routine should start searching. That may not be the absolutely first entry of 
  * the list. Might be NULL, if it points to the end of the list. The routine 
  * might modify the pointer. */
 struct PhotonOrderedListEntry** current,
 /** Data of the photon that should be inserted. */
 Photon* ph);

/** Clear the time-ordered photon list. */
void clear_PhotonList(/** Address of the pointer to the first entry of the list. 
		       * Might be NULL. */
		      struct PhotonOrderedListEntry** pole);

/** Insert a new photon to a binary tree. The routine checks, whether
    the binary tree already exists. If not (the NULL==*first_entry),
    it uses the new photon as the first entry in a new tree and
    changes *first_entry respectively. For that reason the address of
    the pointer to the first PhotonBinaryTreeEntry has been given as
    parameter, such that the function can redirect the pointer to the
    first entry. If the tree already exists, the appropriate position
    is searched, where the new photon can be inserted. */
int insert_Photon2BinaryTree(/** Address of the pointer to the first
			         entry of the binary tree. Can be
			         NULL. */
			     struct PhotonBinaryTreeEntry** first_entry,
			     Photon* photon/**< Data of the photon that should be inserted. */ );

/** Creates a time-ordered photon list from a given binary tree. The
    return value is the error status. The routine deletes the binary
    tree after the readout. */
int CreateOrderedPhotonList
(struct PhotonBinaryTreeEntry** tree_ptr /**< Pointer to the binary tree. 
					  * Will be reset to NULL by this routine. 
					  * (*tree_ptr) might be NULL. */,
 struct PhotonOrderedListEntry** list_first /**< Address of the pointer to the pointer to 
					     * the absolutely 
					     * first entry in the time-ordered list. 
					     * (*list_ptr) might be NULL, if the list is 
					     * empty. */,
 struct PhotonOrderedListEntry** list_current /**< Address of the pointer to the current entry
					       * in the time-ordered photon list.
					       * Might be NULL, if it points to the end of the
					       * list. */ );


#endif  /* PHOTON_H */

