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

#ifndef MPCGS_H
#define MPCGS_H

#include<stdlib.h>


struct mpcgs_opt_t {
    char *gdatfile;
    size_t niter;
    size_t nchain;
    size_t nburn;
    double init_theta;
    long seed;
}; 

void mpcgs_estimate(struct mpcgs_opt_t *options);

#endif /* MPCGS_H */
