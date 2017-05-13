/*
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

#include <stdio.h>
#include <stdlib.h>

#include "debug.h"

extern enum mpcgs_log_level_t log_threshold;
extern enum mpcgs_log_level_t err_threshold;
extern FILE *mpcgs_logfile;
extern FILE *mpcgs_errfile;

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

inline int mpcgs_get_log_threshold()
{

    return (log_threshold);
}

inline int mpcgs_get_err_threshold()
{

    return (err_threshold);
}

FILE *mpcgs_get_log_file()
{

    return (mpcgs_logfile);
}

FILE *mpcgs_get_err_file()
{

    return (mpcgs_errfile);
}

void mpcgs_log_init()
{

    log_threshold = MPCGS_LOG_HIDEBUG;
    err_threshold = MPCGS_LOG_ERR;
    mpcgs_logfile = stdout;
    mpcgs_errfile = stderr;
}

void err_out(const char *err_name, const char *err_str, int rc)
{

    err_err("Error %s: %s.\n", err_name, err_str);
    exit(rc);
}
