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


#include "mpcgs.h"
#include "phylip.h"

void mpcgs_estimate(mpcgs_opt_t *options)
{

    struct ms_tab *data;
    int i;

    if(!filename) {
        //err
    }
    data = init_ms_tab(options->gdatfile);
    
    //init tree

    for(i = 0; i < options->nchain; i++) {
        //do burnin

        //do chain

        //do gradient ascent

    }  


}

