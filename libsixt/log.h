
#ifndef LOG_H
#define LOG_H

#undef log_trace
#undef log_debug
#undef log_info
#undef log_warning
#undef log_error
#undef log_fatal

#define LOG

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
