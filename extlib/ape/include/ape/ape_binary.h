/** \file ape_binary.h
    \brief Declaration of facilities to support binaries using ape, including especially pget, pset etc..
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_ape_binary_h
#define ape_ape_binary_h

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ApeBinary {
  int (*func)(int, char **);
  const char * name;
  const char * usage;
  int min_argc;
  int max_argc;
} ApeBinary;

/** \brief Run a binary top level function, using the given arguments.
    \param argc Number of command line arguments, including executable name.
    \param argv Command line arguments, including argv[0] == executable name.
    \param binary Structure containing all that is needed to run an ape binary.
*/
int ape_binary_run(int argc, char ** argv, ApeBinary * binary);

/** \brief Perform pget action, for retrieving individual parameter values.
    \param argc Number of command line arguments, not including executable name.
    \param argv Command line arguments, not including argv[0] == executable name.
*/
int ape_binary_pget(int argc, char ** argv);

/** \brief Perform plist action, for displaying all values in a parameter file.
    \param argc Number of command line arguments, not including executable name.
    \param argv Command line arguments, not including argv[0] == executable name.
*/
int ape_binary_plist(int argc, char ** argv);

/** \brief Perform pquery action, for prompting for an individual parameter as needed.
    \param argc Number of command line arguments, not including executable name.
    \param argv Command line arguments, not including argv[0] == executable name.
*/
int ape_binary_pquery(int argc, char ** argv);

/** \brief Perform pset action, for setting individual parameter values.
    \param argc Number of command line arguments, *not* including name of executable.
    \param argv Command line arguments, *not* including name of executable.
*/
int ape_binary_pset(int argc, char ** argv);

/** \brief Perform punlearn action, for reverting local parameter file to system settings.
    \param argc Number of command line arguments, *not* including name of executable.
    \param argv Command line arguments, *not* including name of executable.
*/
int ape_binary_punlearn(int argc, char ** argv);

/** \brief Perform pcheck action, for assessing validity of a parameter file format.
    \param argc Number of command line arguments, *not* including name of executable.
    \param argv Command line arguments, *not* including name of executable.
*/
int ape_binary_pcheck(int argc, char ** argv);

/** \brief Run some benchmarks.
    \param argc Number of command line arguments, *not* including name of executable.
    \param argv Command line arguments, *not* including name of executable.
*/
int ape_binary_performance(int argc, char ** argv);

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: ape_binary.h,v $
 * Revision 1.6  2007/07/26 16:16:19  peachey
 * Add ape_binary_performance, for measuring Ape's speed.
 *
 * Revision 1.5  2006/08/22 20:34:14  peachey
 * Move body of pcheck binary to a function ape_binary_pcheck.
 *
 * Revision 1.4  2006/05/23 19:31:11  peachey
 * Add punlearn facilities through ape_binary_punlearn function.
 *
 * Revision 1.3  2006/05/23 16:23:35  peachey
 * Add pset facility through ape_binary_pset function.
 *
 * Revision 1.2  2006/05/22 17:36:16  peachey
 * Add pquery functionality.
 *
 * Revision 1.1  2006/05/19 17:42:22  peachey
 * Add ape_binary module.
 *
*/
