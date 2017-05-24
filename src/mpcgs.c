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
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "debug.h"
#include "mpcgs.h"
#include "phylip.h"
#include "tree.h"

#include "SFMT.h"

static float lkhood_gradient_ascent_step(struct gtree_summary_set *summary_set,
                                         float theta)
{

    float thetalow, thetahigh, thetanext;
    float llkhoodlow, llkhood, llkhoodhigh, llkhoodnext;
    float half_diff, half_sum, lndiff, lslope, jump;
    if (!summary_set) {
        // TODO: handle error
    }

    if (theta <= 0.0) {
        // TODO: handle error
    }

    thetalow = theta - (2 * DELTA);
    thetahigh = theta + (2 * DELTA);

#ifdef MPCGS_NOGPU
    llkhoodlow = gtree_summary_set_llkhood(summary_set, thetalow);
    llkhood = gtree_summary_set_llkhood(summary_set, theta);
    llkhoodhigh = gtree_summary_set_llkhood(summary_set, thetahigh);
#else
    llkhoodlow = gtree_summary_set_llkhood_gpu(summary_set, thetalow);
    llkhood = gtree_summary_set_llkhood_gpu(summary_set, theta);
    llkhoodhigh = gtree_summary_set_llkhood_gpu(summary_set, thetahigh);
#endif /* MPCGS_NOGPU */
    half_diff = (llkhoodhigh - llkhoodlow) / 2.0;
    half_sum = (llkhoodhigh + llkhoodlow) / 2.0;

    errno = 0;
    if (half_diff > 0.0) {
        lslope = log(exp(half_diff) - exp(-half_diff)) + half_sum - LN4DELTA;
        jump = exp(lslope);
    } else if (half_diff < 0.0) {
        lslope = log(exp(-half_diff) - exp(half_diff)) + half_sum - LN4DELTA;
        jump = -exp(lslope);
    } else {
        // Curve is flat...gradient ascent is either done or will fail.
        return (theta);
    }

    if (errno == ERANGE) {
        printf("WARNING: an overflow error occurred.\n");
    }

    lndiff = FLT_MAX;
    do {
        jump = fmin(MAXJUMP, jump);
        jump /= 2.0;
        thetanext = theta + jump;
        if (thetanext <= 0.0) {
            continue;
        }
#ifdef MPCGS_NOGPU
        llkhoodnext = gtree_summary_set_llkhood(summary_set, thetanext);
#else
        llkhoodnext = gtree_summary_set_llkhood_gpu(summary_set, thetanext);
#endif /* MPCGS_NOGPU */
        if (llkhoodnext >= llkhood) {
            break;
        }
        half_diff = (llkhood - llkhoodnext) / 2.0;
        half_sum = (llkhood + llkhoodnext) / 2.0;
        lndiff = log(exp(half_diff) - exp(-half_diff)) + half_sum;
    } while (lndiff > LNDELTA);

    return (thetanext);
}

static float lkhood_by_gradient_ascent(struct gtree_summary_set *summary_set,
                                       float drv_theta)
{

    float thn, thnext;
    int iteration;

    if (!summary_set) {
        // TODO: handle error
    }

    if (drv_theta <= 0.0) {
        // TODO: handle error
    }

#ifdef MPCGS_NOGPU
    gtree_summary_set_base_lposteriors(summary_set, drv_theta);
#else
    gtree_summary_set_base_lposteriors_gpu(summary_set, drv_theta);
#endif /* MPCGS_NOGPU */

    thnext = drv_theta;
    iteration = 0;
    do {
        thn = thnext;
        thnext = lkhood_gradient_ascent_step(summary_set, thn);
        iteration++;
    } while ((fabs(thnext - thn) > EPSILON && (iteration < MAXITER)));

    return (thnext);
}

static void do_multi_proposal(struct chain *ch)
{

    struct chain_param *cparam;
    struct mp_param *mparam;
    struct gene_tree *curr_tree, *proposal;
    struct gtree_summary *summary;
    struct gtree_summary_set *sum_set;
    unsigned int tgtidx, pick;
    float llhood_sum;
    int trans_mtx_pos;
    float trans_mtx_val, trans_mtx_sum;
    float *trans_dist;
    float randf;
    int i, j;

    if (!ch) {
        // TODO: handle error
    }

    cparam = ch->cparam;
    mparam = ch->mp->mparam;

    sum_set = ch->sum_set;
    summary = &sum_set->summaries[sum_set->nsummaries];
    curr_tree = &ch->mp->proposals[ch->mp->curr_idx];
    // TODO: refactor
    tgtidx =
      sfmt_genrand_uint32(&(ch->mp->sfmt)) % ((curr_tree->nnodes + curr_tree->ntips) - 1);
    if (tgtidx == curr_tree->root->idx) {
        tgtidx++;
    }
    proposal = ch->mp->proposals;
    llhood_sum = 0;
    for (i = 0; i < mparam->nproposals; i++) {
        if (i != ch->mp->curr_idx) {
            gtree_propose_fixed_target(
              curr_tree, proposal, ch->theta, tgtidx, &(ch->mp->sfmt));
#ifndef MPCGS_NOGPU
            gtree_set_llhood_gpu(proposal);
#else
            gtree_set_llhood(proposal);
#endif /* MPCGS_NOGPU */
            llhood_sum += proposal->llhood;
        }

        proposal++;
    }

    proposal = ch->mp->proposals;
    for (i = 0; i < mparam->nproposals; i++) {
        trans_mtx_sum = 0.0;
        for (j = 0; j < mparam->nproposals; j++) {
            trans_mtx_pos = (i * mparam->nproposals) + j;
            if (i != j) {
                if (ch->mp->proposals[j].llhood >= proposal->llhood) {
                    trans_mtx_val = 1.0;
                } else {
                    trans_mtx_val =
                      exp(ch->mp->proposals[j].llhood - proposal->llhood);
                }
                trans_mtx_val /= (float)(mparam->nproposals - 1);
                trans_mtx_sum += trans_mtx_val;
                ch->mp->trans_mtx[trans_mtx_pos] = trans_mtx_val;
            }
        }
        trans_mtx_pos = (i * mparam->nproposals) + i;
        ch->mp->trans_mtx[trans_mtx_pos] = 1.0 - trans_mtx_sum;
        proposal++;
    }

    pick = ch->mp->curr_idx;
    for (i = 1; i <= (mparam->npicks * cparam->sum_freq); i++) {
        trans_dist = &ch->mp->trans_mtx[pick * mparam->nproposals];
        randf = 1.0 - sfmt_genrand_real2(&(ch->mp->sfmt));
        pick =
          (unsigned)weighted_pick(trans_dist, mparam->nproposals, 1.0, randf);
        if (mparam->sampling && (i % cparam->sum_freq) == 0) {
            gtree_digest(&ch->mp->proposals[pick], summary);
            ch->sum_set->nsummaries++;
            summary++;
            //printf("proposal->llhood = %f\n", ch->mp->proposals[pick].llhood);
        }
    }

    ch->mp->curr_idx = pick;

}

static float run_chain(struct chain *ch, sfmt_t *sfmt)
{

    struct gene_tree *next_tree;
    struct chain_param *cparam;
    struct gtree_summary *summary;
    struct gene_tree *curr_tree;
    int total_steps;
    float accept;
    float estimate;
    unsigned tgtidx;
    int i;

    if (!ch || !sfmt) {
        // TODO: handle error
    }

    cparam = ch->cparam;
    summary = ch->sum_set->summaries;
    total_steps = cparam->nburnin + (cparam->nsummaries * cparam->sum_freq);
    curr_tree = &ch->mp->proposals[ch->mp->curr_idx];
    for (i = 0; i < total_steps; i++) {
        tgtidx = sfmt_genrand_uint32(sfmt) %
                 ((curr_tree->nnodes + curr_tree->ntips) - 1);
        if (tgtidx == curr_tree->root->idx) {
            tgtidx++;
        }
        next_tree = malloc(sizeof(*next_tree));
        gtree_propose_fixed_target(
          curr_tree, next_tree, ch->theta, tgtidx, sfmt);
        gtree_set_llhood(next_tree);

        if (next_tree->llhood > curr_tree->llhood) {
            curr_tree = next_tree;
        } else {
            accept = exp(next_tree->llhood - curr_tree->llhood);
            if (sfmt_genrand_real2(sfmt) < accept) {
                curr_tree = next_tree;
            }
        }

        if ((i % cparam->sum_freq) == 0 && i > cparam->nburnin) {
            gtree_digest(curr_tree, summary);
            // TODO: encapsulate summary set operations:
            summary++;
            ch->sum_set->nsummaries++;
        }
    }

    estimate = lkhood_by_gradient_ascent(ch->sum_set, ch->theta);
    ch->sum_set->nsummaries = 0;
    return (estimate);
}

static float run_chain_with_multi_proposal(struct chain *ch)
{

    struct chain_param *cparam;
    struct gene_tree *curr_tree;
    struct mp_param burnin_param, sampling_param;
    float estimate;
    int to_pick;

    if (!ch) {
        // TODO: handle error
    }

    cparam = ch->cparam;

    burnin_param.nproposals = 100;
    burnin_param.npicks = 100;
    burnin_param.sampling = 0;

    sampling_param.nproposals = 100;
    sampling_param.npicks = 100;
    sampling_param.sampling = 1;

    // TODO: refactor
    to_pick = cparam->nburnin;
    ch->mp->mparam = &burnin_param;
    if (to_pick % burnin_param.npicks) {
        // TODO: warn
    }

    while (to_pick > 0) {
        do_multi_proposal(ch);
        to_pick -= burnin_param.npicks;
    }

    to_pick = cparam->nsummaries;
    ch->mp->mparam = &sampling_param;
    if (to_pick % sampling_param.npicks) {
        // TODO: err
    }

    while (to_pick > 0) {
        do_multi_proposal(ch);
        to_pick -= sampling_param.npicks;
    }

    estimate = lkhood_by_gradient_ascent(ch->sum_set, ch->theta);
    ch->sum_set->nsummaries = 0;
    return (estimate);
}

static unsigned multi_prop_init(struct multi_proposal **mp,
                                struct ms_tab *data,
                                unsigned nproposal,
                                float theta,
                                unsigned seed)
{

    struct gene_tree *init_tree;
    sfmt_t *sfmt;

    if (!mp || !data) {
        // TODO: handle error
    }

    *mp = malloc(sizeof(**mp));
    (*mp)->curr_idx = 0;
    (*mp)->proposals = calloc(nproposal, sizeof(*((*mp)->proposals)));
    (*mp)->trans_mtx =
      malloc(nproposal * nproposal * sizeof(*((*mp)->trans_mtx)));

    sfmt = &((*mp)->sfmt);
    sfmt_init_gen_rand(sfmt, seed);

    init_tree = &(*mp)->proposals[(*mp)->curr_idx];
    gtree_init(theta, data->len, init_tree, sfmt);
    gtree_add_seqs_to_tips(init_tree, data);
    gtree_set_exp(init_tree);
    gtree_print_newick(init_tree);
    gtree_set_llhood(init_tree);
    log_debug("init tree root time: %f\n", init_tree->root->time);
    log_debug("init tree log likelihood: %f\n", init_tree->llhood);

    return (init_tree->nnodes);
}

void mpcgs_estimate(struct mpcgs_opt_t *options)
{

    struct timeval currtime;
    unsigned long seed;
    struct ms_tab *data;
    unsigned num_nodes, nmpproposals;
    struct chain_param small_chain_param, big_chain_param;
    struct chain ch;
    int i, err;

    const char *err_name = "estimating theta";
    char *err_str;

    if (!options->gdatfile) {
        err = EINVAL;
        err_str = "no sequence file given";
    }

    if (options->seed != 0) {
        seed = options->seed;
    } else {
        gettimeofday(&currtime, NULL);
        seed = currtime.tv_usec;
    }
    log_debug("rseed = %li\n", seed);

    data = init_ms_tab(options->gdatfile);

    small_chain_param.nburnin = big_chain_param.nburnin = 1000;
    small_chain_param.nsummaries = 500;
    big_chain_param.nsummaries = 10000;
    small_chain_param.sum_freq = big_chain_param.sum_freq = 20;

    ch.theta = options->init_theta;
    ch.cparam = &small_chain_param;
    nmpproposals = 100;

#ifdef MPCGS_NOGPU
    num_nodes = multi_prop_init(&(ch.mp), data, nmpproposals, ch.theta, seed);
    gtree_summary_set_create(
      &ch.sum_set, big_chain_param.nsummaries, num_nodes);
#else
    num_nodes = multi_prop_init_gpu(&(ch.mp), data, nmpproposals, ch.theta, seed);
    gtree_summary_set_create_gpu(
      &ch.sum_set, big_chain_param.nsummaries, num_nodes);
#endif
    for (i = 0; i < 10; i++) {
        // ch.theta = run_chain(&ch, &sfmt);
        ch.theta = run_chain_with_multi_proposal(&ch);
        printf("Theta estimate after iteration %i: %f\n", (i + 1), ch.theta);
    }

    ch.cparam = &big_chain_param;
    for (i = 0; i < 2; i++) {
        // ch.theta = run_chain(&ch, &sfmt);
        ch.theta = run_chain_with_multi_proposal(&ch);
        printf(
          "Theta estimate after long iteration %i: %f\n", (i + 1), ch.theta);
    }

    // TODO: free memory

    free_ms_tab(data);

    return;

errout:
    err_out(err_name, err_str, -err);
}
