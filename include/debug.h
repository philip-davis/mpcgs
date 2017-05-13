/*

    debug.h provides logging facilities for mpcgs.

    Copyright (C) 2017  Philip Davis

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

#include <stdio.h>
#include <stdlib.h>

enum mpcgs_log_level_t
{
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

#define mpcgs_log(log_level, ...)                                              \
    {                                                                          \
        if (log_level <= mpcgs_get_log_threshold()) {                          \
            fprintf(mpcgs_get_log_file(), __VA_ARGS__);                        \
            fflush(mpcgs_logfile);                                             \
        }                                                                      \
    }

#define mpcgs_err(log_level, ...)                                              \
    {                                                                          \
        if (log_level <= mpcgs_get_err_threshold()) {                          \
            fprintf(mpcgs_get_err_file(), __VA_ARGS__);                        \
            fflush(mpcgs_errfile);                                             \
        }                                                                      \
    }

#define log_err(...) mpcgs_log(MPCGS_LOG_ERR, __VA_ARGS__)
#define log_warn(...) mpcgs_log(MPCGS_LOG_WARN, __VA_ARGS__)
#define log_debug(...) mpcgs_log(MPCGS_LOG_DEBUG, __VA_ARGS__)
#define log_hidebug(...) mpcgs_log(MPCGS_LOG_HIDEBUG, __VA_ARGS__)
#define err_err(...) mpcgs_err(MPCGS_LOG_ERR, __VA_ARGS__)
#define err_warn(...) mpcgs_err(MPCGS_LOG_WARN, __VA_ARGS__)
#define err_debug(...) mpcgs_err(MPCGS_LOG_DEBUG, __VA_ARGS__)
#define err_hidebug(...) mpcgs_err(MPCGS_LOG_HIDEBUG, __VA_ARGS__)

#define debug_var_decl(errname)                                                \
    int err, param_ok;                                                         \
    const char *err_name = errname;                                            \
    char *err_str

#define param_chk(param_ok, target)                                            \
    {                                                                          \
        if (!(param_ok)) {                                                     \
            err = EINVAL;                                                      \
            err_str = "passed invalid parameters";                             \
            goto target;                                                       \
        }                                                                      \
    }

#define alloc_chk(ptr, target)                                                 \
    {                                                                          \
        if (!ptr) {                                                            \
            err_str = strerror(ENOMEM);                                        \
            err = ENOMEM;                                                      \
            goto target;                                                       \
        }                                                                      \
    }

#define debug_err_out() err_out(err_name, err_str, -err)

void mpcgs_set_log_threshold(const enum mpcgs_log_level_t log_level);
void mpcgs_set_err_threshold(const enum mpcgs_log_level_t err_level);
void mpcgs_set_log_file(FILE *log_file);
void mpcgs_set_err_file(FILE *err_file);
int mpcgs_get_log_threshold();
int mpcgs_get_err_threshold();
FILE *mpcgs_get_log_file();
FILE *mpcgs_get_err_file();
void mpcgs_log_init();
void err_out(const char *err_name, const char *err_str, int rc);

// TODO Adapt CUDA debugging facilities.

#endif /* MPCGS_DEBUG_H */
