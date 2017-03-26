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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "debug.h"
#include "mpcgs.h"
#include "phylip.h"

void mpcgs_estimate(struct mpcgs_opt_t *options)
{

    struct ms_tab *data;
    int i, err;

    const char *err_name = "estimating theta";
    char *err_str;

    if(!options->gdatfile) {
        err = EINVAL;
        err_str = "no sequence file given";
    }
    data = init_ms_tab(options->gdatfile);
    
    //init tree

    for(i = 0; i < options->nchain; i++) {
        //do burnin

        //do chain

        //do gradient ascent

    }  

    free_ms_tab(data);

errout:
    err_out(err_name, err_str, -err);

}

