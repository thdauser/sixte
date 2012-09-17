/** \file ape_test.h
    \brief Declaration of unit test functions.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_ape_test_h
#define ape_ape_test_h

#ifdef __cplusplus
extern "C" {
#endif

/** \brief Perform unit tests of i/o facilities. */
void ape_io_test(void);

/** \brief Perform unit tests of list facility. */
void ape_list_test(void);

/** \brief Perform unit tests of message handling facility. */
void ape_msg_test(void);

/** \brief Perform unit tests of parameter facility. */
void ape_par_test(void);

/** \brief Perform unit tests of parameter group facility. */
void ape_par_group_test(void);

/** \brief Perform unit tests of session interface. */
void ape_session_test(void);

/** \brief Perform unit tests of traditional interface. */
void ape_trad_test(void);

/** \brief Perform unit tests of binary support interface. */
void ape_binary_test(void);

/** \brief Perform unit tests of general internal utilities. */
void ape_util_test(void);

/** \brief Perform unit tests of general internal utilities. */
void ape_pil_test(void);

/** \brief Perform all unit tests. */
int ape_test(int argc, char ** argv);

/** \brief Compare actual outcome of a test producing a double and a status to its expected
    outcome, using the supplied hint to customize the message.
    \param hint Test-specific hint.
    \param result Actual result of a test.
    \param expected_result Expected result of a test.
    \param status Actual status code returned by test.
    \param expected_status Expected status code returned by test.
*/
void ape_test_cmp_double(const char * hint, double result, double expected_result, int status, int expected_status);

/** \brief Compare actual outcome of a test producing a integral value and a status to its expected
    outcome, using the supplied hint to customize the message.
    \param hint Test-specific hint.
    \param result Actual result of a test.
    \param expected_result Expected result of a test.
    \param status Actual status code returned by test.
    \param expected_status Expected status code returned by test.
*/
void ape_test_cmp_long(const char * hint, long result, long expected_result, int status, int expected_status);

/** \brief Compare actual outcome of a test producing a pointer value and a status to its expected
    outcome, using the supplied hint to customize the message.
    \param hint Test-specific hint.
    \param result Actual result of a test.
    \param expected_result Expected result of a test.
    \param status Actual status code returned by test.
    \param expected_status Expected status code returned by test.
*/
void ape_test_cmp_ptr(const char * hint, void * result, void * expected_result, int status, int expected_status);

/** \brief Compare actual outcome of a test producing an unsigned integral value and a status to its expected
    outcome, using the supplied hint to customize the message.
    \param hint Test-specific hint.
    \param result Actual result of a test.
    \param expected_result Expected result of a test.
    \param status Actual status code returned by test.
    \param expected_status Expected status code returned by test.
*/
void ape_test_cmp_ulong(const char * hint, unsigned long result, unsigned long expected_result,
  int status, int expected_status);

/** \brief Compare actual outcome of a test producing a string and a status to its expected outcome,
    using the supplied hint to customize the message.
    \param hint Test-specific hint.
    \param result Actual result of a test.
    \param expected_result Expected result of a test.
    \param status Actual status code returned by test.
    \param expected_status Expected status code returned by test.
*/
void ape_test_cmp_string(const char * hint, const char * result, const char * expected_result, int status, int expected_status);

/** \brief Compare actual outcome of a test producing an array of strings and a status to its expected outcome,
    using the supplied hint to customize the message.
    \param hint Test-specific hint.
    \param result Actual result of a test.
    \param expected_result Expected result of a test.
    \param status Actual status code returned by test.
    \param expected_status Expected status code returned by test.
*/
void ape_test_cmp_string_array(const char * hint, char ** result, const char ** expected_result,
  int status, int expected_status);

/** \brief Report failure of a test. Sets global status flag. */
void ape_test_failed(const char * fmt, ...);

/** \brief Set status of unit test.
    \param status The status of the test. If global status is already set, nothing happens.
*/
void ape_test_set_status(int status);

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: ape_test.h,v $
 * Revision 1.11  2010/06/02 19:22:09  peachey
 * Add initial (partial) implementation of ape_session module.
 *
 * Revision 1.10  2007/10/09 16:45:01  peachey
 * Use command line arguments in test binary to name the log file.
 * If no log file, do not redirect output for the convenience of debugging.
 *
 * Revision 1.9  2006/06/06 13:28:34  peachey
 * Add ape_pil_test.
 *
 * Revision 1.8  2006/05/19 17:43:46  peachey
 * Add ape_binary_test.
 *
 * Revision 1.7  2006/04/26 14:31:09  peachey
 * Add ape_test_cmp_string_array for comparing arrays of strings.
 *
 * Revision 1.6  2006/04/20 03:43:23  peachey
 * Add ape_test_cmp_ptr for comparing pointers in tests.
 *
 * Revision 1.5  2006/04/14 14:15:22  peachey
 * Add ape_test_cmp_ulong for comparing unsigned results.
 *
 * Revision 1.4  2006/04/13 18:42:24  peachey
 * Add ape_test_cmp_* family of functions to facilitate smaller and more
 * consistent tests.
 *
 * Revision 1.3  2006/04/12 14:19:30  peachey
 * Test traditional interface.
 *
 * Revision 1.2  2006/04/10 21:12:17  peachey
 * Add test for ape_util module.
 *
 * Revision 1.1.1.1  2006/04/05 13:45:19  peachey
 * Initial import of All-purpose Parameter Environment (APE).
 *
*/
