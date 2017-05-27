/*
        tree.h defines the structures and functions for a phylogenetic tree.

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

#ifndef MPCGS_TREE_H
#define MPCGS_TREE_H

#include <stdlib.h>
#include <string.h>

#include "SFMT.h"
#include "phylip.h"

#define MPCGS_NUM_FREQ_TERM 8
#define NUM_BASE 4 //TODO: get rid of static references to 4 throughout code

enum
{
    FREQ_A = 0,
    FREQ_T,
    FREQ_C,
    FREQ_G,
    FREQ_AR,
    FREQ_GR,
    FREQ_CY,
    FREQ_TY
};

struct gene_node
{
    int order;
    int idx;
    struct gene_node *parent;
    struct gene_node *child1;
    struct gene_node *child2;
    struct gene_node *prev;
    struct gene_node *next;
    struct gene_tree *tree;
    float time;
    int exp_valid;
    float lexpA;
    float lexpB;
    float lexpC;
    struct mol_seq *mseq;
#ifndef MPCGS_NOGPU
    float *tip_llhoods;
#endif /* MPCGS_NOGPU */
};

struct gnode_list
{
    struct gene_node **gnodes;
    int head;
    int tail;
};

struct gene_tree
{
    struct ms_tab *mstab;
    struct gene_node *nodes;
    struct gene_node *tips;
    size_t nnodes;
    size_t ntips;
    struct gene_node *root;
    struct gene_node *last;
    float lfreq[MPCGS_NUM_FREQ_TERM];
    float xrate;
    float yrate;
#ifndef MPCGS_NOGPU
    float *block_scratch;
    //Some pre-computed likelihood kernel parameters:
    size_t block_size;
    size_t num_blocks;
    size_t shared_size;
    struct gene_node **node_list_scratch;
    float *rand_scratch;
#endif /* MPCGS_NOGPU */
    float llhood;
};

struct gtree_summary
{
    float *intervals;
    size_t nintervals;
    float ldrv_posterior;
    float ltmp_lkhood_comp;
};

struct gtree_summary_set
{
    struct gtree_summary *summaries;
    size_t nsummaries;
    size_t szintervals;
#ifndef MPCGS_NOGPU
    float *block_scratch;
#endif /* MPCGS_NOGPU */
};

void gtree_simulate_tree(struct gene_tree *gtree, float theta, sfmt_t *sfmt);
void gtree_init(float theta,
                size_t ntips,
                struct gene_tree *gtree,
                sfmt_t *sfmt);
void gtree_add_seqs_to_tips(struct gene_tree *gtree, struct ms_tab *mstab);
void gtree_set_exp(struct gene_tree *gtree);
void gtree_set_llhood(struct gene_tree *gtree);
void gtree_print_newick(struct gene_tree *gtree);
struct gene_tree *gtree_propose(struct gene_tree *current,
                                float theta,
                                sfmt_t *sfmt);
void gtree_propose_fixed_target(struct gene_tree *current,
                                             struct gene_tree *proposal,
                                             float theta,
                                             unsigned int tgtidx,
                                             sfmt_t *sfmt);
void gtree_digest(struct gene_tree *gtree, struct gtree_summary *digest);
void gtree_summary_set_create(struct gtree_summary_set **sum_set,
                              size_t count,
                              size_t nintervals);
void gtree_summary_set_base_lposteriors(struct gtree_summary_set *sum_set,
                                        float drv_theta);
float gtree_summary_set_llkhood(struct gtree_summary_set *summary_set,
                                float theta);
void gtree_summary_set_print_lkhoods(struct gtree_summary_set *summary_set,
                                     float start,
                                     float stop,
                                     float incr);
size_t weighted_pick(float *dist,
                     size_t num_picks,
                     float sum_dist,
                     float randf);

#ifndef MPCGS_NOGPU
void gtree_add_seqs_to_tips_gpu(struct gene_tree *gtree, struct ms_tab *mstab);
void gtree_nodes_init_gpu(struct gene_tree *gtree, size_t ntips, size_t seq_len);
void gtree_init_gpu(float theta,
                size_t ntips,
                struct gene_tree *gtree,
                sfmt_t *sfmt);
void gtree_summary_set_create_gpu(struct gtree_summary_set **sum_set,
                                  size_t count,
                                  size_t nintervals);
void gtree_summary_set_base_lposteriors_gpu(struct gtree_summary_set *sum_set,
                                            float drv_theta);
float gtree_summary_set_llkhood_gpu(struct gtree_summary_set *summary_set,
                                    float theta);
void gtree_set_llhood_gpu(struct gene_tree *gtree);

#endif /* MPCGS_NOGPU */

#endif /* MPCGS_TREE_H */
