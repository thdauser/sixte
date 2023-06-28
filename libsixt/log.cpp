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

   Copyright 2023:  LOG has been developed by the INSTITUTO DE FISICA DE
   CANTABRIA (CSIC-UC) with funding under different projects:
   ESP2006-13608-C02-01, AYA2012-39767-C02-01, ESP2013-48637-C2-1-P,
   ESP2014-53672-C3-1-P, RTI2018-096686-B-C21 and PID2021-122955OB-C41.

***********************************************************************
*                      LOG
*
*  File:       log.cpp
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#include "log.h"

#ifdef LOG

#include <cstdarg>
#include <stdio.h>
#include <string>
#include <sstream>

#include <thread>
#include <mutex>

namespace tlog
{
  static std::mutex test_mutex;
  static std::string test_header("%llu (thread id: %s: ");

  static unsigned long long test_get_timestamp()
  {
    return static_cast<unsigned long long> (std::clock()/(CLOCKS_PER_SEC/1000));
  }

  static std::string test_get_thread_id()
  {
    std::thread::id id (std::this_thread::get_id());
    std::stringstream ss;
    ss << id;
    return (ss.str());
  }

  static void test_log(const char* fmt, va_list args)
  {
    std::lock_guard<std::mutex> guard(test_mutex);
    printf(test_header.c_str(),
            test_get_timestamp(),test_get_thread_id().c_str());
    std::string format(fmt);
    format.append("\n");
    vprintf(format.c_str(), args);
  }

  void test(const char* format, ...)
  {
    va_list args;
    va_start(args, format);
    test_log(format, args);
    va_end(args);
  }
}

namespace slog
{
  static const level default_level = level::FATAL;
  //static const level default_level = level::TRACE;
  static std::string header("%llu (thread id: %s [%7.7s]: ");

  static std::mutex level_mutex;
  static std::mutex write_mutex;

  static level get_default_level()
  {
    std::lock_guard<std::mutex> guard(level_mutex);
    return default_level;
  }
  
  static unsigned long long get_timestamp()
  {
    return static_cast<unsigned long long> (std::clock()/(CLOCKS_PER_SEC/1000));
  }

  static std::string get_thread_id()
  {
    std::thread::id id (std::this_thread::get_id());
    std::stringstream ss;
    ss << id;
    return (ss.str());
  }

  static const char* levels[] = 
    {
      "  TRACE",
      "  DEBUG",
      "   INFO",
      "WARNING",
      "  ERROR",
      "  FATAL"
    };

  static const char* level_to_string(level lvl)
  {
    if (static_cast<std::size_t>(lvl) < sizeof(levels) / sizeof(*levels))
      return levels[lvl];
    return "";
    //return static_cast<const char*>(lvl);
  }

  static void log(level lvl, const char* fmt, va_list args)
  {
    std::lock_guard<std::mutex> guard(write_mutex);
    if (lvl >= get_default_level()){
      printf(header.c_str(),
             get_timestamp(),get_thread_id().c_str(),level_to_string(lvl));
      std::string format(fmt);
      format.append("\n");
      vprintf(format.c_str(), args);
    }
  }

  void trace(const char* format, ...)
  {
    va_list args;
    va_start(args, format);
    log(level::TRACE, format, args);
    va_end(args);
  }
  void debug(const char* format, ...)
  {
    va_list args;
    va_start(args, format);
    log(level::DEBUG, format, args);
    va_end(args);
  }
  void info(const char* format, ...)
  {
    va_list args;
    va_start(args, format);
    log(level::INFO, format, args);
    va_end(args);
  }
  void warning(const char* format, ...)
  {
    va_list args;
    va_start(args, format);
    log(level::WARNING, format, args);
    va_end(args);
  }
  void error(const char* format, ...)
  {
    va_list args;
    va_start(args, format);
    log(level::ERR, format, args);
    va_end(args);
  }
  void fatal(const char* format, ...)
  {
    va_list args;
    va_start(args, format);
    log(level::FATAL, format, args);
    va_end(args);
  }
}

#endif
