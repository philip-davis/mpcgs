/*
 * genealogy.h
 *
 * Philip Davis, 2015-2016
 * ------------------
 *
 * Provides the gene tree data structure genealogy.
 */
 
#ifndef MLGENE_GENEALOGY_H
#define MLGENE_GENEALOGY_H


#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include <curand_kernel.h>

#include "aligned_seq.h"

/* struct: gene_node
 * -----------------
 * A node in a genealogical tree.
 *
 */
typedef struct gene_node {
    int order;
    int idx;
    struct gene_node *parent;
    struct gene_node *child1;
    struct gene_node *child2;
    struct gene_node *prev;
    struct gene_node *next;
    char isseq;
	char *seq;
    float time;
    //char dirty;
	//float *lkspace;
} gene_node;

/* struct: genealogy
 * -----------------
 * A genealogical tree.
 *
 */
typedef struct genealogy {
    gene_node *nodes;
	gene_node *seqs;
    gene_node *root;
	gene_node *lastnode;
	//float *lkspace;
	float *blksum;
	float llike;
} genealogy;

typedef struct neighborhood {
	gene_node *target;
	gene_node *parent;
	gene_node *ancestor;
	gene_node *child[3];
} neighborhood;

__device__ float active_no_coal_prob(unsigned j, float t, 
										unsigned z, float theta);
__device__ float active_one_coal_prob(unsigned j, float t, 
										unsigned z, float theta);
__device__ float active_two_coal_prob(unsigned j, float t,
										unsigned z, float theta);
void build_gene_upgma(aligned_seqs *seq_dat, genealogy *g, size_t num_seq, 											size_t seq_len, float thetascale);
__device__ void dislocate_node(genealogy *g, gene_node *node);
__device__ void duplicate_genealogy(genealogy *src, genealogy *dst, 
											size_t seq_len);
float find_dist(gene_node *node1, gene_node *node2, size_t seq_len);
__device__ void fix_order(genealogy *g, neighborhood *nei);
__device__ void foul_likelihood(gene_node *node);
__device__ float gen_coal_length1(float t, size_t j, size_t z,
								float theta, float randf);
__device__ float gen_coal_length2(float t, size_t j, size_t z,
								float theta, float randf);
__device__ void get_coalescent_intervals(genealogy *proposal, neighborhood *nei, 
									float *active_prob, float *coal_prob, 
									size_t num_intervals, gene_node **coal_int, 
									float *rand_scratch);
__device__ float get_init_coal_rv_prob(float x, float t, size_t z, float theta);
__device__ void get_neighborhood(gene_node *target, neighborhood *nei);
__device__ size_t get_next_active(size_t num_active, float *active_prob,
									float *coal_prob, unsigned child_term,
									float randf);
__device__ size_t get_num_inactive(neighborhood *nei, gene_node *node);
__device__ size_t get_num_intervals(genealogy *g, neighborhood *nei,
										size_t num_seq);
size_t get_tree_size(gene_node *node);
void init_gene_node(gene_node *node, unsigned idx, char isseq);
gene_node **init_trees(aligned_seqs *seq_dat, gene_node *seqs, size_t num_seq);
genealogy *new_genealogy(size_t num_seq, size_t seq_len);
__device__ void order_3_nodes(gene_node **n1, gene_node **n2, gene_node **n3);
__device__ void populate_a_c_probs(genealogy *proposal, neighborhood *nei, 
									float theta, float *aprob, 
									float *cprob, size_t num_intervals);
void print_newick(gene_node *node);
void reduce_genealogy(genealogy *g, float *sample);
__device__ void relocate_nodes(genealogy *proposal, neighborhood *nei, 
								gene_node **coal_int, float theta, 
								float *rand_scratch, unsigned randui);
void upgma_merge(gene_node *node, gene_node **trees, gene_node *next,
						size_t num_trees, size_t num_seq);
__device__ __host__ size_t weighted_pick(float *prob, size_t num, float total,
										float randf);
						
#endif