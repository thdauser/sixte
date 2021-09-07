/***********************************************************************
   This file is part of SIXTE/SIRENA software.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

   Copyright 2014:  TASKSSIRENA has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01, 
   ESP2013-48637-C2-1-P, ESP2014-53672-C3-1-P and RTI2018-096686-B-C21.

***********************************************************************
*                      LOG
*
*  File:       log.h
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef LOG_H
#define LOG_H

#undef log_test
#undef log_trace
#undef log_debug
#undef log_info
#undef log_warning
#undef log_error
#undef log_fatal

#define LOG
#define TEST_LOG

#ifndef TEST_LOG

#define log_test(fmt, ...) ((void)0)

#else

#define log_test tlog::test

namespace tlog 
{
  void test(const char* format, ...);
}

#endif

#ifndef LOG

#define log_trace(fmt, ...) ((void)0)
#define log_debug(fmt, ...) ((void)0)
#define log_info(fmt, ...) ((void)0)
#define log_warning(fmt, ...) ((void)0)
#define log_error(fmt, ...) ((void)0)
#define log_fatal(fmt, ...) ((void)0)

#else

#define log_trace slog::trace
#define log_debug slog::debug
#define log_info slog::info
#define log_warning slog::warning
#define log_error slog::error
#define log_fatal slog::fatal

namespace slog
{
  enum level
  {
    TRACE = 0,
    DEBUG = 1,
    INFO = 2,
    WARNING = 3,
    ERR = 4,
    FATAL = 5,
  };
  void trace(const char* format, ...);
  void debug(const char* format, ...);
  void info(const char* format, ...);
  void warning(const char* format, ...);
  void error(const char* format, ...);
  void fatal(const char* format, ...);
}

#endif

#endif
