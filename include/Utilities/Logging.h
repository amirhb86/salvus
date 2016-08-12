#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

// The global logging state. Define at the top of `main()` only once
#define INIT_LOGGING_STATE() GlobalLoggerState GLOBAL_LOGGER_STATE;

struct GlobalLoggerState;
extern GlobalLoggerState GLOBAL_LOGGER_STATE;

enum class LogProc {
  ROOTONLY, ALLPROCS
    };

enum class LogLevel {
  STD, DEBUG, VERBOSE, ERROR
    };

enum class LogWhere {
  STDOUT, FILE, STDOUT_AND_FILE  
    };

class Logger {
protected:
  std::ostringstream os;
  LogLevel level;
private:
  Logger(const Logger &);
  Logger &operator=(const Logger &);

public:
  Logger();
  virtual ~Logger();
  std::ostringstream& Get(LogLevel level_) { level = level_; return os; }
};

#define LOG() Logger().Get(LogLevel::STD)
#define DEBUG() Logger().Get(LogLevel::DEBUG)
#define VERBOSE() Logger().Get(LogLevel::VERBOSE)
#define ERROR() Logger().Get(LogLevel::ERROR)

struct GlobalLoggerState {
  LogProc proc = LogProc::ROOTONLY;
  LogLevel level = LogLevel::STD;
  std::ofstream output_file;
  LogWhere log_where = LogWhere::STDOUT;
};

