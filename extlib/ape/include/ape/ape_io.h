/** \file ape_io.h
    \brief Declaration of input/output facilities.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_ape_io_h
#define ape_ape_io_h

#include "ape/ape_list.h"
#include "ape/ape_par.h"

#ifdef __cplusplus
extern "C" {
#endif

struct ApeParFile;
typedef struct ApeParFile ApeParFile;

/** \brief Create an exact duplicate of an ApeParFile. 
    \param orig_file Pointer to the original file.
    \param clone_file Pointer to the output cloned file pointer.
*/
int ape_io_clone_file(ApeParFile * orig_file, ApeParFile ** clone_file);

/** \brief Free up resources associated with parameter file abstraction. Does not delete parameter
    file on disk.
    \param par_file The file structure to destroy.
*/
void ape_io_destroy_file(ApeParFile * par_file);

/** \brief Find the named parameter file.
    \param par_name The name of the parameter to find.
    \param par_file The parameter file to search.
    \param par Output iterator pointing to found parameter.
*/
int ape_io_find_par(const char * par_name, ApeParFile * par_file, ApeListIterator * par_itor);

/** \brief Get container of parameters from the given file.
    \param par_file The parameter file whose parameter container to get.
    \param par_cont Pointer to output container pointer.
*/
int ape_io_get_par_cont(ApeParFile * par_file, ApeList ** par_cont);

/** \brief Get the name of the given parameter file structure.
    \param par_file The parameter file structure.
    \param file_name The name of the parameter file.
*/
int ape_io_get_file_name(const ApeParFile * par_file, const char ** file_name);

/** \brief Set the name of the given parameter file object.
    \param par_file The parameter file whose name is being set.
    \param file_name Pointer to the input new file name string.
*/
int ape_io_set_file_name(ApeParFile * par_file, const char * file_name);

/** \brief Get the default mode, i.e. the value of the mode parameter for the given file.
    \param par_file The file.
    \param mode The value of the mode.
*/
int ape_io_get_default_mode(ApeParFile * par_file, char * mode);

/** \brief Merge values from system and local parameter files to produce a single merged parameter file.
    \param sys_par_file The input system parameter file.
    \param loc_par_file The input local parameter file.
    \param merged_par_file The output (merged) parameter file.
*/
int ape_io_merge_par_files(ApeParFile * sys_par_file, ApeParFile * loc_par_file, ApeParFile ** merged_par_file);

/** \brief Restore to their previous state any unlearned parameters in a parameter file.
    \param current The current set of parameters, i.e. the ones to be unlearned.
    \param previous The previous set of parameters, whose values are restored to the current set.
*/
int ape_io_revert_unlearned(ApeParFile * current, ApeParFile * previous);

/** \brief Get PFILES variable setting currently in use by ape. Note this may or may not be the same as
    the environment variable PFILES, depending on whether the client has called ape_io_set_pfiles.
    The client should free the returned string.
    \param pfiles Pointer to output variable holding the current value for PFILES.
*/
int ape_io_get_pfiles(char ** pfiles);

/** \brief Set PFILES variable. Note this does not alter the environment variable, just the internal
    setting used by ape to locate parameter files.
    \param pfiles The new value for PFILES, which will be copied. If pfiles is 0, ape *will* use the
    value of the environment variable PFILES.
*/
int ape_io_set_pfiles(const char * pfiles);

/** \brief Get a list whose elements can be cast to const char * and interpreted as the left-to-right
    contents of the local parameter file path. (This is the left path in the PFILES variable.)
    \param loc_path Pointer to list pointer containing the directories in the path. This is read-only
    memory; client is NOT responsible.
*/
int ape_io_get_loc_path(ApeList ** loc_path);

/** \brief Get a list whose elements can be cast to const char * and interpreted as the left-to-right
    contents of the system parameter file path. (This is the right path in the PFILES variable.)
    \param sys_path Pointer to list pointer containing the directories in the path. This is read-only
    memory; client is NOT responsible.
*/
int ape_io_get_sys_path(ApeList ** sys_path);

/** \brief Read the named file and interpret its contents to create a parameter container.
    \param file_name The name of the input parameter file.
    \param par_file The output parameter file object.
*/
int ape_io_read(const char * file_name, ApeParFile ** par_file);

/** \brief Read the named file and interpret its contents to create a parameter container.
    \param file_name The file name only of the input parameter file.
    \param path Path to search for the file.
    \param file The parameter file object to read.
*/
int ape_io_read_file_path(const char * file_name, ApeList * path, ApeParFile ** par_file);

/** \brief Write a parameter container to the named file.
    \param par_file The parameter file to write.
    \param force_write Flag indicating whether to override read-only status of file and force the file to be written.
*/
int ape_io_write(ApeParFile * par_file, char force_write);

/** \brief Interpret the given command line arguments in the context of the given parameter file object,
    and apply the arguments to modify the parameter file object..
    \param par_file The par file to be modified.
    \param argc The number of arguments (not including the executable/parameter file name.)
    \param argv The arguments (not including the executable/parameter file name.)
*/
int ape_io_apply_command_line(ApeParFile * par_file, int argc, char ** argv);

/** \brief Check contents of parameter file and report any errors in the file format.
    \param par_file The file to check.
    \param check_value Flag specifying whether the parameter value itself should be checked.
*/
int ape_io_check_file_format(ApeParFile * par_file, char check_value);

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: ape_io.h,v $
 * Revision 1.20  2006/05/23 16:26:56  peachey
 * Add read-only flag to par file structure. Add force_write arguement to
 * ape_io_write for overriding read-only status.
 *
 * Revision 1.19  2006/05/19 17:29:32  peachey
 * Update TODO which was done.
 *
 * Revision 1.18  2006/05/18 03:15:43  peachey
 * Add ape_io_set_file_name to allow the file name used by a
 * file object to be modified.
 *
 * Revision 1.17  2006/05/17 02:19:48  peachey
 * Add ape_io_get_par_cont, for getting the container of parameters from a file.
 *
 * Revision 1.16  2006/05/16 19:46:16  peachey
 * Add ape_io_get_pfiles, and automatically perform default
 * initialization if ape_io_get_loc_path or ape_io_get_sys_path are called
 * and PFILES is not defined.
 *
 * Revision 1.15  2006/05/12 17:22:53  peachey
 * Add ape_io_revert_unlearned, to forget unlearned parameters when requested.
 *
 * Revision 1.14  2006/05/12 03:30:17  peachey
 * Add function ape_io_get_default_mode, for getting mode parameter value from file.
 *
 * Revision 1.13  2006/05/10 16:48:49  peachey
 * Generalize format tests to allow them to be run with or without
 * value checking, to let them be used both in contexts where the whole parameter
 * needs to be checked (checking user input) and contexts in which a bad value
 * is not a show-stopper, because the user or client may yet correct the problem.
 *
 * Revision 1.12  2006/05/09 19:31:42  peachey
 * Add ape_io_check_file_format, for checking the whole parameter file.
 *
 * Revision 1.11  2006/05/02 23:48:28  peachey
 * Add ape_io_apply_command_line for handling command line parameters,
 * directly imported from pil's PIL_allow_spaces_cmndarg2value.
 *
 * Revision 1.10  2006/05/01 19:11:58  peachey
 * Correct typo.
 *
 * Revision 1.9  2006/04/24 13:28:43  peachey
 * Add ape_io_clone_file and the beginnings of ape_io_merge_par_files.
 *
 * Revision 1.8  2006/04/21 14:52:48  peachey
 * Change signature of ape_io_find_par so that it returns an iterator
 * instead of the parameter itself.
 *
 * Revision 1.7  2006/04/21 01:30:57  peachey
 * Add and test ape_io_find_par, to locate a parameter by name.
 *
 * Revision 1.6  2006/04/19 02:30:09  peachey
 * Refactor ape_io_write to use ApeParFile instead of ApeList.
 *
 * Revision 1.5  2006/04/19 01:33:02  peachey
 * Change signature of ape_io_read so that it returns an ApeParFile rather than
 * an ApeList of parameters.
 *
 * Revision 1.4  2006/04/15 03:06:25  peachey
 * Add ape_io_destroy_file for cleaning up ApeParFiles.
 *
 * Revision 1.3  2006/04/15 00:47:06  peachey
 * Add more support for dealing with parameter files. Specifically, a
 * struct ApeParFile, which contains a file name and a container of
 * parameters, and two new functions:
 * o ape_io_get_file_name: retrieves the file name from an ApeParFile.
 * o ape_io_read_file_path: searches a path for a file. If it is found,
 *   open it and load its parameters.
 *
 * Revision 1.2  2006/04/14 03:43:21  peachey
 * Add facilities for handling PFILES: ape_io_set_pfiles, ape_io_get_loc_path,
 * ape_io_get_sys_path.
 *
 * Revision 1.1.1.1  2006/04/05 13:45:19  peachey
 * Initial import of All-purpose Parameter Environment (APE).
 *
*/
