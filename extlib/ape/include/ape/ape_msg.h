/** \file ape_msg.h
    \brief Declaration of message facilities.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_ape_msg_h
#define ape_ape_msg_h

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** \brief Display a debugging message. Only has effect if debugging mode enabled.
    \param fmt Format for output, see documentation for printf for details.
    \param ... Additional (variable) arguments giving the items to display.
*/
void ape_msg_debug(const char * fmt, ...);

/** \brief Enable/disable debugging mode, which results in many more messages.
    \param enable Flag indicating whether to enable (non-0) or disable (0) debugging messages.
*/
void ape_msg_debug_enable(char enable);

/** \brief Get current debug mode.
*/
char ape_msg_get_debug_mode(void);

/** \brief Display an error message.
    \param fmt Format for output, see documentation for printf for details.
    \param ... Additional (variable) arguments giving the items to display.
*/
void ape_msg_error(const char * fmt, ...);

/** \brief Display a (suppressible) informational message.
    \param chatter Importance level for the message, the lower the more important.
    \param fmt Format for output, see documentation for printf for details.
    \param ... Additional (variable) arguments giving the items to display.
*/
void ape_msg_info(unsigned int chatter, const char * fmt, ...);

/** \brief Display an unsuppressible output message.
    \param fmt Format for output, see documentation for printf for details.
    \param ... Additional (variable) arguments giving the items to display.
*/
void ape_msg_out(const char * fmt, ...);

/** \brief Redirect output destined for stderr to the given stream.
    \param new_stream Stream to use for output.
*/
void ape_msg_set_err_stream(FILE * new_stream);

/** \brief Redirect output destined for stdout to the given stream.
    \param new_stream Stream to use for output.
*/
void ape_msg_set_out_stream(FILE * new_stream);

/** \brief Set the function which handles output from ape_msg_out.
    \param func The new handler function.
*/
void ape_msg_set_out_handler(void (*func)(const char *));

/** \brief Display a (suppressible) warning message.
    \param chatter Importance level for the message, the lower the more important.
    \param fmt Format for output, see documentation for printf for details.
    \param ... Additional (variable) arguments giving the items to display.
*/
void ape_msg_warn(unsigned int chatter, const char * fmt, ...);

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: ape_msg.h,v $
 * Revision 1.6  2006/06/16 01:18:48  peachey
 * Add ape_msg_get_debug_mode, for getting current debug mode.
 *
 * Revision 1.5  2006/06/03 02:11:37  peachey
 * Allow redirection of output stream using functions.
 *
 * Revision 1.4  2006/05/19 17:30:33  peachey
 * Add ape_msg_set_out_stream, for redirecting stdout the way ape already
 * redirects stderr.
 *
 * Revision 1.3  2006/05/18 13:58:40  peachey
 * Add ape_msg_out, for output.
 *
 * Revision 1.2  2006/04/12 17:56:33  peachey
 * Add ERROR banner to ape_test_failed, and add ape_msg_info function.
 *
 * Revision 1.1.1.1  2006/04/05 13:45:19  peachey
 * Initial import of All-purpose Parameter Environment (APE).
 *
*/
