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

#include <stdlib.h>

#include "tree.h"

#define EPSILON .0005
#define MAXITER 1000
#define DELTA 2.2204460492503131e-6
#define LNDELTA -13.0178024591766992
#define LN4DELTA -11.6315080980568086
#define MAXJUMP 20.0

struct mpcgs_opt_t
{
    char *gdatfile;
    size_t niter;
    size_t nchain;
    size_t nburn;
    double init_theta;
    long seed;
};

struct chain_param
{
    size_t nburnin;
    size_t nsummaries;
    size_t sum_freq;
};

struct mp_param
{
    size_t nproposals;
    size_t npicks;
    unsigned sampling;
};

struct multi_proposal
{
    struct mp_param *mparam;
    struct gene_tree *proposals;
    float *trans_mtx;
    unsigned int curr_idx;
};

struct chain
{
    float theta;
    struct chain_param *cparam;
    struct multi_proposal mp;
    struct gtree_summary_set *sum_set;
};

void mpcgs_estimate(struct mpcgs_opt_t *options);

#endif /* MPCGS_H */
