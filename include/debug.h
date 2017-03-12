/*

    debug.h provides logging facilities for mpcgs.

    Copyright (C) {2017}  {Philip Davis}

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef MPCGS_DEBUG_H
#define MPCGS_DEBUG_H

#include<stdlib.h>
#include<stdio.h>

enum mpcgs_log_level_t {
    MPCGS_LOG_NONE,
    MPCGS_LOG_ERR,
    MPCGS_LOG_WARN,
    MPCGS_LOG_DEBUG,
    MPCGS_LOG_HIDEBUG
};

static enum mpcgs_log_level_t log_threshold;
static enum mpcgs_log_level_t err_threshold;
static FILE *mpcgs_logfile;
static FILE *mpcgs_errfile;

void mpcgs_set_log_threshold(const enum mpcgs_log_level_t log_level)
{

    log_threshold = log_level;

}

void mpcgs_set_err_threshold(const enum mpcgs_log_level_t err_level)
{

   err_threshold = err_level;

}

void mpcgs_set_log_file(FILE *log_file)
{

    mpcgs_logfile = log_file;

}

void mpcgs_set_err_file(FILE *err_file)
{

    mpcgs_errfile = err_file;

}

void mpcgs_log_init()
{

    log_threshold = MPCGS_LOG_NONE;
    err_threshold = MPCGS_LOG_ERR;
    mpcgs_logfile = stdout;
    mpcgs_errfile = stderr;

}

#define mpcgs_log(log_level, ...) {\
    if(log_level <= log_threshold) {\
        fprintf(mpcgs_logfile, __VA_ARGS__);\
        fflush(mpcgs_logfile);\
    }\
}

#define mpcgs_err(log_level, ...) {\
    if(log_level <= err_threshold) {\
        fprintf(mpcgs_errfile, __VA_ARGS__);\
        fflush(mpcgs_errfile);\
    }\
}

#define log_err(...) mpcgs_log(MPCGS_LOG_ERR, __VA_ARGS__)
#define log_warn(...) mpcgs_log(MPCGS_LOG_WARN, __VA_ARGS__)
#define log_debug(...) mpcgs_log(MPCGS_LOG_DEBUG, __VA_ARGS__)
#define log_hidebug(...) mpcgs_log(MPCGS_LOG_HIDEBUG, __VA_ARGS__)
#define err_err(...) mpcgs_err(MPCGS_LOG_ERR, __VA_ARGS__)
#define err_warn(...) mpcgs_err(MPCGS_LOG_WARN, __VA_ARGS__)
#define err_debug(...) mpcgs_debug(MPCGS_LOG_DEBUG, __VA_ARGS__)
#define err_hidebug(...) mpcgs_hidebug(MPCGS_LOG_HIDEBUG, __VA_ARGS__)

//TODO Adapt CUDA debugging facilities.

#endif /* MPCGS_DEBUG_H */
