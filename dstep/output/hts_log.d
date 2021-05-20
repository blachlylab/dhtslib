/// \file htslib/hts_log.h
/// Configuration of log levels.
/* The MIT License
Copyright (C) 2017 Genome Research Ltd.

Author: Anders Kaplan

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

extern (C):

/// Log levels.
enum htsLogLevel
{
    HTS_LOG_OFF = 0, ///< All logging disabled.
    HTS_LOG_ERROR = 1, ///< Logging of errors only.
    HTS_LOG_WARNING = 3, ///< Logging of errors and warnings.
    HTS_LOG_INFO = 4, ///< Logging of errors, warnings, and normal but significant events.
    HTS_LOG_DEBUG = 5, ///< Logging of all except the most detailed debug events.
    HTS_LOG_TRACE = 6 ///< All logging enabled.
}

/// Sets the selected log level.
void hts_set_log_level(htsLogLevel level);

/// Gets the selected log level.
htsLogLevel hts_get_log_level();

/// Selected log level.
/*!
 * One of the HTS_LOG_* values. The default is HTS_LOG_WARNING.
 * \note Avoid direct use of this variable. Use hts_set_log_level and hts_get_log_level instead.
 */
extern __gshared int hts_verbose;

/*! Logs an event.
* \param severity      Severity of the event:
*                      - HTS_LOG_ERROR means that something went wrong so that a task could not be completed.
*                      - HTS_LOG_WARNING means that something unexpected happened, but that execution can continue, perhaps in a degraded mode.
*                      - HTS_LOG_INFO means that something normal but significant happened.
*                      - HTS_LOG_DEBUG means that something normal and insignificant happened.
*                      - HTS_LOG_TRACE means that something happened that might be of interest when troubleshooting.
* \param context       Context where the event occurred. Typically set to "__func__".
* \param format        Format string with placeholders, like printf.
*/
void hts_log(
    htsLogLevel severity,
    const(char)* context,
    const(char)* format,
    ...);

/*! Logs an event with severity HTS_LOG_ERROR and default context. Parameters: format, ... */

/*! Logs an event with severity HTS_LOG_WARNING and default context. Parameters: format, ... */

/*! Logs an event with severity HTS_LOG_INFO and default context. Parameters: format, ... */

/*! Logs an event with severity HTS_LOG_DEBUG and default context. Parameters: format, ... */

/*! Logs an event with severity HTS_LOG_TRACE and default context. Parameters: format, ... */

// #ifndef HTS_LOG_H
