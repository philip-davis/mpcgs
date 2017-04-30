/*
    mpcgs - multiple-proposal coalescent genealogy sampler

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

#include<argp.h>
#include<errno.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "debug.h"
#include "mpcgs.h"
#include "phylip.h"
#include "tree.h"

#define MPCGS_FLAG_NBURN 0x01
#define MPCGS_FLAG_NCHAIN 0x02
#define MPCGS_FLAG_GDATFILE 0x04
#define MPCGS_FLAG_NITER 0x08
#define MPCGS_FLAG_THETA 0x10

//Argument parsing
const char *argp_program_version = "mpcgs 0.2.0a";
static char doc[] = "mpcgs -- a parallel coalescent genealogy sampler.";

static struct argp_option options[] = {
    {"verbose", 'v', "LEVEL",  OPTION_ARG_OPTIONAL, "Set verbosity LEVEL" },
    {"quiet",   'q', 0, 0, "Produce minimal output" },
    {"gpu",     'g', 0, OPTION_HIDDEN, "use GPU for theta estimation"},
	{"seed",	's', "INTEGER", OPTION_HIDDEN, "random number generator seed"},
    {"mpi",     'm', 0, OPTION_HIDDEN, "distribute processing using MPI"},
    {"threads", 'N', "INTEGER", OPTION_HIDDEN, "run INTEGER parallel threads"},
    {"numiter", 'n', "INTEGER", 0, "run INTEGER estimation iterations"},
    {"chainlen", 'c', "INTEGER", 0, "generated INTEGER samples per iteration"},
    {"burninlen", 'b', "INTEGER", 0, "length of burn-in phase"},
    {"theta",   't', "VALUE", 0, "Set initial driving theta to VALUE" },
    {"input",   'i', "FILE", 0, "Sequence data contained in FILE (phylib format)"},
    { 0 }
};


static error_t parse_opt(int key, char *arg, struct argp_state *state)
{

    struct mpcgs_opt_t *arguments = state->input;
    static unsigned int optflags = 0;

    switch(key) {
        case 'b':
            arguments->nburn = atoi(arg);
            break;
        case 'c':
            arguments->nchain = atoi(arg);
            break;
        case 'i':
            arguments->gdatfile = arg;
            break;
        case 'n':
            arguments->niter = atoi(arg);
            break;
        case 's':
        	arguments->seed = atol(arg);
        	break;
        case 't':
            arguments->init_theta = atof(arg);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }

    

    return 0;

}

static struct argp argp = {options, parse_opt, 0, doc};


///end argument parsing

int main(int argc, char *argv[])
{

    struct mpcgs_opt_t args = {0};
    int err;

    mpcgs_log_init();
    mpcgs_set_log_threshold(MPCGS_LOG_HIDEBUG);
	mpcgs_set_err_threshold(MPCGS_LOG_HIDEBUG);

    err = argp_parse(&argp, argc, argv, 0, 0, &args);
    if(err) {
        err_out("parsing arguments", strerror(err), -err);
    }

    mpcgs_estimate(&args);

    return 0;

}
