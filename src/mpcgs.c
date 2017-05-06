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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "debug.h"
#include "mpcgs.h"
#include "phylip.h"
#include "tree.h"

#include "SFMT.h"

//float gradient_ascent()

static float run_chain(struct chain *ch,
						struct gtree_summary_set *summary_set,
						sfmt_t *sfmt)
{

	struct gene_tree *next_tree;
	struct chain_param *param;
	struct gtree_summary *summary;
	int total_steps;
	float accept;
	int i;

	param = ch->param;
	summary = summary_set->summaries;
	total_steps = param->nburnin + (param->nsummaries * param->sum_freq);
	for(i = 0; i < total_steps; i++) {

		next_tree = gtree_propose(ch->curr_tree, ch->theta, sfmt);
		gtree_set_llhood(next_tree);

		if(next_tree->llhood > ch->curr_tree->llhood) {
			ch->curr_tree = next_tree;
		} else {
			accept = exp(next_tree->llhood - ch->curr_tree->llhood);
			if(sfmt_genrand_real2(sfmt) < accept) {
				ch->curr_tree = next_tree;
			}
		}

		if((i % param->sum_freq) == 0 && i > param->nburnin) {
			gtree_digest(ch->curr_tree, summary);
			summary++;
		}

	}

}

void mpcgs_estimate(struct mpcgs_opt_t *options)
{

	struct timeval currtime;
	unsigned long seed;
	sfmt_t sfmt;
    struct ms_tab *data;
	struct gene_tree *curr_tree;
	float drv_theta;
	struct gtree_summary_set summary_set;
	struct chain_param small_chain_param, big_chain_param;
	struct chain ch;
    int i, err;

    const char *err_name = "estimating theta";
    char *err_str;

    if(!options->gdatfile) {
        err = EINVAL;
        err_str = "no sequence file given";
    }


    if(options->seed != 0) {
    	seed = options->seed;
    } else {
    	gettimeofday(&currtime, NULL);
    	seed = currtime.tv_usec;
    }
	log_debug("rseed = %li\n", seed);
	sfmt_init_gen_rand(&sfmt, seed);
	
	small_chain_param.nburnin = big_chain_param.nburnin = 1000;
	small_chain_param.nsummaries = 500;
	big_chain_param.nsummaries = 10000;
	small_chain_param.sum_freq = big_chain_param.sum_freq = 20;

	drv_theta = options->init_theta;
    data = init_ms_tab(options->gdatfile);
    curr_tree = gtree_init(drv_theta, data->len, &sfmt);
	gtree_add_seqs_to_tips(curr_tree, data);
	gtree_set_exp(curr_tree);
	gtree_print_newick(curr_tree);
	gtree_set_llhood(curr_tree);
	log_debug("init tree root time: %f\n", curr_tree->root->time);
	log_debug("init tree log likelihood: %f\n", curr_tree->llhood);
	

	gtree_summary_set_create(&summary_set, big_chain_param.nsummaries, curr_tree->nnodes);
	ch.curr_tree = curr_tree;
	ch.theta = drv_theta;
	ch.param = &small_chain_param;


	run_chain(&ch, &summary_set, &sfmt);


    //init tree


        //do burnin

        //do chain

        //do gradient ascent



	//TODO: free tree
	
    free_ms_tab(data);
	
	return;

errout:
    err_out(err_name, err_str, -err);

}

