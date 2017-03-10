// genealogy.h
// Philip Davis, 2015
// 
//

#ifndef GENEALOGY_H
#define GENEALOGY_H

#include"sfmt/SFMT.h"
#include<stdlib.h>

#include "sequence.h"

typedef struct g7y {
//leaves will have seq, NULL for a node (branch)
	unsigned isleaf;
	seq *seq;
	struct g7y *parent;
	struct g7y *left;
	struct g7y *right;
//an order is maintained by depth, with tarray being a pointer
//to the node's spot on an array of pointers.
	struct g7y *prev;
	struct g7y *next;
	struct g7y **tarray;
	unsigned order;
        double depth;
//array of size 4 * seqlen, keeping the likelihood of each nucleotide at each
//site.
        double *lk;
//does the likelihood of a given site need updating (could be made more efficient)
        char *lkset;
//the overall likelihood of this node being in this configuration
	double llk;
//is llk set?
        unsigned llkset;
//a set of temporary values for holding a prosposal of a new configuration
	struct g7y *prparent;
	struct g7y *prleft;
	struct g7y *prright;
	struct g7y *prprev;
	struct g7y *prnext;
	unsigned prorder;
	unsigned active;

} g7y;

//create a new node at depth with the given children
g7y *new_node(g7y *lchild, g7y *rchild, double depth);

//create a new leaf holding the given sequence
g7y *new_leaf(seq *s);

//copy the current "production" values into the proposal
void copy_proposal(g7y *g);

//remove the parent from the proposal (the the chold from the former parent)
void prop_orphan(g7y *g);

//delete the node from the proposal order
void prop_order_removal(g7y *g);

//remove the target and parent from the tree proposal
void erase_active_prop(g7y *target, g7y *parent, g7y **child);

//print the genealogy in the newick format
void print_newick(g7y *g);

//find the sum of the distances between each of the leaves in two clusters
double cluster_dsum_g7y(g7y *c1, g7y *c2, double **dmtx);

//find the distance between two clusters
double cluster_dist_g7y(g7y **cl, size_t *cllen, unsigned i, 
					unsigned j, double **dmtx);

//join the closest two clusters in an array of clusters
g7y *upgma_join_g7y(g7y **clusters, size_t ncl, size_t *cllen, double **dmtx);

//build a upgma genealogy of a sequence set.
g7y *sample_upgma_g7y(seqset *sset);

double prob_sub(nuc_t i, nuc_t j, double t, double *nucfreq);

double *upd_like_g7y(g7y *g, unsigned pos, double *nucfreq);

double get_llike_g7y(g7y *root, seqset *sset);

void get_ordered_children(g7y *target, g7y *parent, g7y **child)

double prob_no_coalesc(unsigned x, double t, unsigned z, double theta);

double prob_one_coalesc(unsigned x, double t, unsigned z, double theta);

double prob_two_coalesc(unsigned x, double t, unsigned z, double theta)

double prob_act_coalesc(unsigned x, unsigned y, double t,
                                        unsigned z, double theta);

g7y *sample_posterior(g7y *root, seqset *sset, sfmt_t *sfmt, double theta);

#endif
