#include<float.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include <cuda.h>

#include "aligned_seq.h"
#include "genealogy.h"

//extern __constant__ gene_node seqnodes[12];

//return the probability of no coalescent events involving the active
//times over the course of an time
//j - active times
//t - length of time
//z - inactive times
__device__ float active_no_coal_prob(unsigned j, float t, 
										unsigned z, float theta)
{

	float e1, x;
	
	e1 = (j * (j - 1)) + (2 * j * z);
	x = exp(-e1 * (t / theta));
	
	return(x);

}

__device__ float active_no_coal_lprob(unsigned j, float t, 
										unsigned z, float theta)
{

	float e1, x;
	
	e1 = (j * (j - 1)) + (2 * j * z);
	x = -e1 * (t / theta);
	
	return(x);

}

//return the probability of a single coalescent event involving the active
//times over the course of an time
//j - active times
//t - length of time
//z - inactive times
//Note: I am assuming an error in the original paper/coalesce1.5 code. The
// numerator of the factor in the paper and code does not include 2jz.
__device__ float active_one_coal_prob(unsigned j, float t, 
										unsigned z, float theta)
{

  float coeff, e1, e2, x;
  
  coeff = (float)((j * (j - 1)) + (2 * j * z)) / (float)((2 * (j - 1)) + (2 * z));
  e1 = ((j - 1) * (j - 2)) + (2 * (j - 1) * z);
  e2 = (j * (j - 1)) + (2 * j * z);
  x = coeff * (exp(-e1 * (t / theta)) - exp(-e2 * (t / theta)));
  
  return(x);

}

__device__ float active_one_coal_lprob(unsigned j, float t, 
										unsigned z, float theta)
{

  float coeff, e1, e2, x;
  
  coeff = ((float)j / 2.0) * (1.0 + ((float)z / (float)(j - 1.0 + z)));
  e1 = ((j - 1.0) * (j - 2.0)) + (2.0 * (j - 1.0) * z);
  e2 = 1 - expf(-(float)((2.0 * (j - 1.0)) + (2.0 * z)) * (t / theta));
  
  x = logf(coeff) - ((t / theta) * e1) + logf(e2);
  
  return(t==0.0?-FLT_MAX:x);

}


//return the probability of a single coalescent event involving the active
//times over the course of an time
//j - active times
//t - length of time
//z - inactive times
//Note: I am assuming an error in the original paper/coalesce1.5 code. The
// equation in the paper is very different from what's in the original code,
// which in turn has a somewhat different coefficient.
//
// NOTE!!! This actually ignores j, assumingfloat active_two_coal_prob(unsigned j, float t, unsigned z, float theta) that j=3.
__device__ float active_two_coal_prob(unsigned j, float t,
										unsigned z, float theta)
{

	float coeff2, coeff3, e1, e2a, e2b, e3, x;
	
	e1 = (2 * z);
	coeff2 = 1.0 / (float)(2 + z);
	e2a = (2 + (2 * z));
	e2b = (4 + (2 * z));
	coeff3 = 1.0 / (float)(3 + (2 * z));
	e3 = (6 + (4 * z));
	x = (3 + (6 * z)) * exp(-e1 * (t / theta)) * 
			((coeff2 * exp(-e2a * (t / theta)) * (exp(-e2b * (t / theta)) - 1)) -
			 (coeff3 * (exp(-e3 * (t / theta)) - 1)));

	return(x);
	
}

//return the probability of a single coalescent event involving the active
//times over the course of an time
//j - active times
//t - length of time
//z - inactive times
//Note: I am assuming an error in the original paper/coalesce1.5 code. The
// equation in the paper is very different from what's in the original code,
// which in turn has a somewhat different coefficient.
//
// NOTE!!! This actually ignores j, assumingfloat active_two_coal_prob(unsigned j, float t, unsigned z, float theta) that j=3.
__device__ float active_two_coal_lprob(unsigned j, float t,
										unsigned z, float theta)
{

	float coeff1, coeff2, e1a, e1b, e1c, e2a, e2b, x;
	
	coeff1 = 3.0 + (6.0 * z);
	coeff2 = (2 * z) * (t / theta);
	e1a = 1.0 / (float)(2.0 + z);
	e1b = expf(-(float)(2.0 + (2.0 * z)) * (t / theta));
	e1c = expf(-(float)(4.0 + (2.0 * z)) * (t / theta)) - 1.0;
	e2a = 1.0 / (float)(3.0 + (2.0 * z));
	e2b = expf(-(float)(6.0 + (4.0 * z)) * (t / theta)) - 1.0;
	
	x = (logf(coeff1) - coeff2) + logf((e1a * e1b * e1c) - (e2a * e2b));
	
	return(t==0.0?-FLT_MAX:x);
	
}

void build_gene_upgma(aligned_seqs *seq_dat, genealogy *g, size_t num_seq, 											size_t seq_len, float thetascale)
{
	
	size_t num_trees;
	gene_node **trees;
	int i;
	gene_node *node, *next;
	
	if(NULL == g) {
		fprintf(stderr, "Tried to build a UPGMA in a NULL genealogy.\n");
    }

	num_trees = num_seq;
	trees = init_trees(seq_dat, g->seqs, num_seq);
	node = NULL;
	for(i = (num_seq - 2); i >= 0; i--) {
		next = node;
		node = &g->nodes[i];
		upgma_merge(node, trees, next, num_trees, seq_len);
		num_trees--;
	}
	g->root = g->nodes;
	g->lastnode = &g->nodes[num_seq - 2];
	free(trees);
	for(i = 0; i < num_seq - 1; i++) {
		g->nodes[i].time *= thetascale;
	}
	
}

__device__ void dislocate_node(genealogy *g, gene_node *node) 
{
	
	if(NULL != node->prev) {
		(node->prev)->next = node->next;
	}
	if(NULL != node->next) {
		(node->next)->prev = node->prev;
	}
	if(g->lastnode == node) {
		g->lastnode = node->prev;
	}

}


extern __device__ void duplicate_genealogy(genealogy *src, genealogy *dst, 
											size_t num_node)
{
	
	int i;
	gene_node *snode, *dnode;
	
	for(i = 0; i < num_node; i++) {
		snode = &src->nodes[i];
		dnode = &dst->nodes[i];
		*dnode = *snode;
		//This all could be replaced by calculating an offset once, then
		//just incrementing the pointer, but I don't know if I can trust
		//the ptr arithmatic to work that way here.
		if(snode->parent != NULL) {
			dnode->parent = &dst->nodes[(snode->parent - src->nodes)];
		}
		if(snode->child1 != NULL && !snode->child1->isseq) {
			dnode->child1 = &dst->nodes[(snode->child1 - src->nodes)];
		}
		if(snode->child2 != NULL && !snode->child2->isseq) {
			dnode->child2 = &dst->nodes[(snode->child2 - src->nodes)];
		}
		if(snode->prev != NULL) {
			dnode->prev = &dst->nodes[(snode->prev - src->nodes)];
		}
		if(snode->next != NULL) {
			dnode->next = &dst->nodes[(snode->next - src->nodes)];
		}
	}
	dst->root = &dst->nodes[(src->root - src->nodes)];
	dst->lastnode = &dst->nodes[(src->lastnode - src->nodes)];
	
}

float find_dist(gene_node *node1, gene_node *node2, size_t seq_len)
{
    //the flow of this function is to perform a depth first traversal of
    //the tree rooted on node1. At each leaf of node1, in turn perform a
    //depth-first traversal of the tree rooted on node2, summing the
    //distances between the node1 leaf and the node2 leaves. At the root of
    //the recursive function, divide the summation across all the leaves of
    //the node1-rooted tree by the product of the size of the nod1 and node2
    //trees and return that value. 

    float dist;
    unsigned i;
	
	if((NULL == node1) || (NULL == node2)) {
		fprintf(stderr, "Attempted to find distance to a null tree.\n");
		return(0);
    }
	
	if(((NULL == node1->seq) 
			&& ((NULL == node1->child1) || (NULL == node1->child2)))
			|| ((NULL == node2->seq)
			&& ((NULL == node2->child1) || (NULL == node2->child2)))) {
		fprintf(stderr, "Found a node with neither children nor a sample.\n");
		return(0);
    }
	dist = 0;
    if((NULL != node1->seq) && (NULL != node2->seq)) {
	//If both nodes are leaves, just calculate the distance based on
	//the number of different elements of the sample.
		for(i = 0; i < seq_len; i++) {
			if(node1->seq[i] != node2->seq[i]) {
				dist++;
			}
		}
    } else {
	//At least one of node1 and nod2 are not leaves
		if(NULL != node1->seq) {
			//node1 is a leaf, so node2 is not
			dist += find_dist(node1, node2->child1, seq_len);
			dist += find_dist(node1, node2->child2, seq_len);
		} else {
		//node1 is not a leaf
			dist += find_dist(node1->child1, node2, seq_len);
			dist += find_dist(node1->child2, node2, seq_len);
		}
    }
	if((NULL == node1->parent) && (NULL == node2->parent)) {
	//this is not a recursive call, this is the root call, and
	//we are about to finish. Normalize by the number of nodes.
        dist /= (float)(get_tree_size(node1) * get_tree_size(node2));
    }

    return(dist);
}

__device__ void fix_order(genealogy *g, neighborhood *nei)
{
	
	gene_node *first, *last, *node;
	
	first = nei->ancestor;
	if(NULL == first) {
		first = nei->parent;
	}
	last = nei->child[1];
	for(node = first; node != last && NULL != node; node = node->next) {
		if(node->prev == NULL) {
			node->order = 0;
		} else {
			node->order = (node->prev)->order + 1;
		}
		if(NULL == node->next) {
			//A little hidden here
			g->lastnode = node;
		}
	}
	
}

__device__ void foul_likelihood(gene_node *node)
{
	
	while(NULL != node) {
		//node->dirty = 1;
		node = node->parent;
	}
	
}


__device__ float gen_coal_length1(float t, size_t j, size_t z,
								float theta, float randf)
{ 	
	float coeff1, e1, x;
	
	coeff1 = -theta / (float)((2 * (j - 1)) + (2 * z));
	e1 = (float)((2 * (j - 1)) + (2 * z));
	
	x = coeff1 * logf(1.0 - ((1.0 - randf) * 
							(1.0 - expf(-e1 * (t / theta)))));

	return(x);
	
}

//NOTE: ignores j, assumes j==3.
//binary search is a modified version of that found in coalesce v1.50 
// by Felsenstein et al.
__device__ float gen_coal_length2(float t, size_t j, size_t z,
								float theta, float randf)
{
	
	float r;
	float x, xmin, xmax;
	int i;
	
	r = randf;
	xmin = 0;
	xmax = t;
	//Beyond 22, float precision will be exhausted.
	for(i = 0; i < 22; i++) {
		x = (xmin + xmax) / 2.0;
		if(get_init_coal_rv_prob(x, t, z, theta) > r) {
			xmax = x;
		} else {
			xmin = x;
		}
	}
	
	return(x);
	
}

__device__ void get_coalescent_intervals(genealogy *proposal, neighborhood *nei, 
									float *active_prob, float *coal_prob, 
									size_t num_intervals, gene_node **coal_int, 
									float *rand_scratch)
{
	
	gene_node *node, *next;
	size_t num_active, next_active, num_set;
	int i;
	float *apiter, *cpiter;
	unsigned child_term;
	
	coal_int[0] = NULL;
	coal_int[1] = NULL;
	
	node = nei->ancestor;
	if(NULL == node) {
		if((nei->parent)->next == nei->target) {
			next = (nei->target)->next;
		} else {
			next = (nei->parent)->next;
		}
	} else {
		next = node->next;
	}
	num_active = 1;
	num_set = 0;
	apiter = &active_prob[3];
	cpiter = coal_prob;
	for(i = 0; i < num_intervals; i++) {
		child_term = (node == nei->child[0]);
		next_active = get_next_active(num_active, apiter, 
										cpiter, child_term, 
										rand_scratch[i]);
		if(((num_active + 1) == next_active) || 
			(child_term && (num_active == next_active))) {
			num_set++;
			if(1 == num_set) {
				coal_int[0] = node;
			} else {
				coal_int[1] = node;
				break;
			}
		} else if(((num_active + 2) == next_active) ||
					(child_term && ((num_active + 1) == next_active))) {
			coal_int[0] = node;
			coal_int[1] = node;
			num_set += 2;
			break;
		}
		num_active = next_active;
		node = next;
		next = node->next;
		
		apiter += 3;
		cpiter += 6;
	}
	
}									

__device__ float get_init_coal_rv_prob(float x, float t, size_t z, float theta)
{
	
	float coeff1, coeff2;
	float e1a, e1b, e2;
	float num, den;
	
	coeff1 = 1.0 / (float)(2.0 + z);
	coeff2 = 1.0 / (float)(3.0 + (2.0 * z));
	e1a = 2.0 + (2.0 * (float)z);
	e1b = 4.0 + (2.0 * z);
	e2 = 6.0 + (4.0 * z);
	num = (coeff1 * exp(-e1a * (t / theta)) * (exp(-e1b * (x / theta)) - 1)) -
			(coeff2 * (exp(-e2 * (x / theta)) - 1));
	den = (coeff1 * exp(-e1a * (t / theta)) * (exp(-e1b * (t / theta)) - 1)) -
			(coeff2 * (exp(-e2 * (t / theta)) - 1));
	return(num / den);
	
}

__device__ void get_neighborhood(gene_node *target, neighborhood *nei)
{

	nei->target = target;
	nei->parent = target->parent;
	nei->ancestor = nei->parent->parent;
	if(target != (nei->parent)->child1) {
		nei->child[0] = (nei->parent)->child1;
	} else {
		nei->child[0] = (nei->parent)->child2;
    }
    nei->child[1] = target->child1;
    nei->child[2] = target->child2;
    order_3_nodes(&nei->child[0],
				  &nei->child[1],
				  &nei->child[2]);
	
}

__device__ size_t get_next_active(size_t num_active, float *active_prob,
									float *coal_prob, unsigned child_term,
									float randf)
{

	unsigned base;
	int i;
	float lprob[3], prob[3];
	float normalizer;
	size_t next_active;
	unsigned idx;
	float scale;
	
	if(child_term) {
		num_active--;
	}
	base = 3;
	if(1 == num_active) {
		base = 0;
	}
	for(i = 0; i < 3; i++) {
		prob[i] = 0.0;
	}
	scale = -FLT_MAX;
	//---this will cause some execution divergence---
	for(i = base; i < (base + (4 - num_active)); i++) {
		idx = (num_active - 1) + (i - base);
		lprob[idx] = active_prob[idx] + coal_prob[i];
		scale = fmaxf(scale, lprob[idx]);
	}
	normalizer = 0.0;
	for(i = base; i < (base + (4 - num_active)); i++) {
		idx = (num_active - 1) + (i - base);
		prob[idx] = expf(lprob[idx] - scale);
		normalizer += prob[idx];
	}
	//---
	next_active = weighted_pick(prob, 3, normalizer, randf) + 1;
	
	return(next_active);
	
}

__device__ size_t get_num_inactive(neighborhood *nei, gene_node *node)
{
	
	unsigned num_inactive;
	
	num_inactive = node->order + 2;
	if(node->order >= (nei->parent)->order) {
		num_inactive -= 2;
	}
	if(node->order >= (nei->target)->order) {
		num_inactive--;
	}
	if(node->order >= (nei->child[0])->order) {
		num_inactive++;
	}
	
	return(num_inactive);
	
}

__device__ size_t get_num_intervals(genealogy *g, neighborhood *nei,
										size_t num_seq)
{
	
	size_t num_interval;
	
	if(!nei->child[1]->isseq) {
		num_interval = (nei->child[1])->order - 1;
	} else {
		num_interval = num_seq - 2;
	}
	if(NULL != nei->ancestor) {
		num_interval -= (nei->ancestor)->order + 1;
	}
	
	return(num_interval);
	
}

size_t get_tree_size(gene_node *node)
{

    size_t tree_size;

    if(NULL == node) {
		return(0);
    } else if(NULL != node->seq) {
    	return(1);
	}
	
	tree_size = 0;
    if(NULL != node->child1) {
		tree_size += get_tree_size(node->child1);
    }
    if(NULL != node->child2) {
		tree_size += get_tree_size(node->child2);
    }
	
    return(tree_size);   

}

void init_gene_node(gene_node *node, unsigned idx, char isseq)
{

    if(NULL != node) {
		node->idx = idx;
		node->order = -1;
		node->parent = NULL;
		node->child1 = NULL;
		node->child2 = NULL;
		node->prev = NULL;
		node->next = NULL;
		if(!isseq) {
			//node->dirty = 1;
			node->seq = NULL;

		} else {
			//node->dirty = 0;
			node->time = 0;
		}
		node->isseq = isseq;
    }

}

gene_node **init_trees(aligned_seqs *seq_dat, gene_node *seqs, size_t num_seq)
{
	
	gene_node *node, **trees;
	int i;
	
	trees = (gene_node **)malloc(num_seq * sizeof(void *));
	for(i = 0; i < num_seq; i++) {
		node = &seqs[i];
		init_gene_node(node, (num_seq - 1) + i, 1);
		node->seq = seq_dat->seqs[i];
		trees[i] = node;
	}
	
	return(trees);
	
}

__device__ void order_3_nodes(gene_node **n1, gene_node **n2, gene_node **n3)
{

    gene_node *t;

    if((((*n2)->order > (*n3)->order) && (-1 != (*n3)->order)) || 
						(-1 == (*n2)->order)) {
		t = *n2;
		*n2 = *n3;
		*n3 = t;
    }
    if((((*n1)->order > (*n2)->order) && (-1 != (*n2)->order)) || 
						(-1 == (*n1)->order))  {
		t = *n1;
		*n1 = *n2;
		*n2 = t;
    }
    if((((*n2)->order > (*n3)->order) && (-1 != (*n3)->order)) || 
						(-1 == (*n2)->order)) {
        t = *n2;
        *n2 = *n3;
        *n3 = t;
    }

}

__device__ void populate_a_c_probs(genealogy *proposal, neighborhood *nei, 
									float theta, float *aprob, 
									float *cprob, size_t num_intervals)
{
	
	float *apiter, *cpiter;
	gene_node *prev, *node, *next;
	int i;
	size_t num_inactive;
	float time;
	float prod1, prod2, prod3, shift;
	
	apiter = &aprob[3 * num_intervals];
	if(nei->child[0]->isseq) {
		apiter[0] = -FLT_MAX;
		apiter[1] = -FLT_MAX;
		apiter[2] = 0.0;
	} else {
		apiter[0] = -FLT_MAX;
		apiter[1] = 0.0;
		apiter[2] = -FLT_MAX;
	}
	apiter -= 3;
	cpiter = &cprob[6 * (num_intervals - 1)];
	
	prev = (nei->child[1])->prev;
	if(NULL == prev) {
		prev = proposal->lastnode;
	}
	
	for(i = num_intervals - 1; i > 0; i--) {
		node = prev;
		prev = node->prev;
		next = node->next;
		num_inactive = get_num_inactive(nei, node);
		if(next == NULL) {
			time = node->time;
		} else {
			time = node->time - next->time;
		}
		cpiter[0] = active_no_coal_lprob(1, time, num_inactive, theta);
		cpiter[1] = active_one_coal_lprob(2, time, num_inactive, theta);
		cpiter[2] = active_two_coal_lprob(3, time, num_inactive, theta);
		cpiter[3] = active_no_coal_lprob(2, time, num_inactive, theta);
		cpiter[4] = active_one_coal_lprob(3, time, num_inactive, theta);
		cpiter[5] = active_no_coal_lprob(3, time, num_inactive, theta);
		if(node != nei->child[0]) {
			prod1 = apiter[3] + cpiter[0];
			prod2 = apiter[4] + cpiter[1];
			prod3 = apiter[5] + cpiter[2];
			shift = fmaxf(prod1, prod2);
			shift = fmaxf(shift, prod3);
			apiter[0] = expf((apiter[3] + cpiter[0]) - shift);
			apiter[0] += expf((apiter[4] + cpiter[1]) - shift);
			apiter[0] += expf((apiter[5] + cpiter[2]) - shift);
			apiter[0] = logf(apiter[0]) + shift;
			
			prod1 = apiter[4] + cpiter[3];
			prod2 = apiter[5] + cpiter[4];
			shift = fmaxf(prod1, prod2);
			apiter[1] = expf((apiter[4] + cpiter[3]) - shift);
			apiter[1] += expf((apiter[5] + cpiter[4]) - shift);
			apiter[1] = logf(apiter[1]) + shift;
			
			apiter[2] = apiter[5] + cpiter[5];
		} else {
			apiter[0] = -FLT_MAX;
			
			prod1 = apiter[3] + cpiter[0];
			prod2 = apiter[4] + cpiter[1];
			shift = fmaxf(prod1, prod2);
			apiter[1] = expf((apiter[3] + cpiter[0]) - shift);
			apiter[1] += expf((apiter[4] + cpiter[1]) - shift);
			apiter[1] = logf(apiter[1]) + shift;
			
			apiter[2] = apiter[4] + cpiter[3];
		}
		apiter -= 3;
		cpiter -= 6;
	}
	cpiter[0] = 0.0;
	cpiter[1] = 0.0;
	cpiter[2] = 0.0;
	cpiter[3] = -FLT_MAX;
	cpiter[4] = -FLT_MAX;
	cpiter[5] = -FLT_MAX;
	apiter[0] = 0.0;
	apiter[1] = -FLT_MAX;
	apiter[2] = -FLT_MAX;
	
}

void print_newick(gene_node *node) 
{
	
	if(NULL != node->child1) {
		printf("(");
		print_newick(node->child1);
	}
	if(NULL != node->child2) {
		printf(",");
		print_newick(node->child2);
		printf(")");
	}
	printf("%i", node->idx);
	printf(":%f",node->time);
	if(NULL == node->parent) {
		printf(";\n");
	}
	
}

void reduce_genealogy(genealogy *g, float *sample)
{
	
	gene_node *node;
	float *smpiter;
	
	smpiter = sample;
	for(node = g->root; node->next != NULL; node = node->next) {
		*smpiter = node->time - (node->next)->time;
		smpiter++;
	}
	*smpiter = node->time;
	
}

__device__ void relocate_nodes(genealogy *proposal, neighborhood *nei, 
								gene_node **coal_int, float theta, 
								float *rand_scratch, unsigned randui)
{
	
	gene_node *prev[2], *next[2];
	float length[2];
	float interval;
	size_t num_inactive;
	size_t num_active;
	unsigned uncle_idx;
	gene_node **uncle, **seluncle;
	gene_node *swap;
	
	prev[0] = coal_int[0];
	if(coal_int[0] == coal_int[1]) {
		prev[1] = nei->parent;
		next[0] = nei->target;
		if(NULL == coal_int[0]) {
			if((nei->parent)->next == nei->target) {
				next[1] = (nei->target)->next;
			} else {
				next[1] = (nei->parent)->next;
			}
			length[0] = -(theta / 2.0) * log(rand_scratch[0]);
			length[1] = -(theta / 6.0) * log(rand_scratch[1]);
		}
		else {
			next[1] = coal_int[1]->next;
			if(NULL == coal_int[0]->next) {		
				interval = coal_int[0]->time;
			} else {
				interval = coal_int[0]->time - (coal_int[0]->next)->time;
			}
			num_inactive = get_num_inactive(nei, coal_int[0]);
			
			length[1] = gen_coal_length2(interval, 3, num_inactive,
											theta, rand_scratch[1]);
			interval -= length[1];
			length[0] = gen_coal_length1(interval, 2, num_inactive,
											theta, rand_scratch[0]);
		}
	} else {
		prev[1] = coal_int[1];
		if(NULL == coal_int[0]) {
			if((nei->parent)->next == nei->target) {
				next[0] = (nei->target)->next;
			} else {
				next[0] = (nei->parent)->next;
			}
			length[0] = -(theta / 2.0) * log(rand_scratch[0]);
		} else {
			next[0] = coal_int[0]->next;
			if(NULL == coal_int[0]->next) {	
				interval = coal_int[0]->time;
			} else {
				interval = coal_int[0]->time - (coal_int[0]->next)->time;
			}
			num_inactive = get_num_inactive(nei, coal_int[0]);
			
			length[0] = gen_coal_length1(interval, 2, num_inactive,
											theta, rand_scratch[0]);
			
		}
		if(NULL == coal_int[1]->next) {
				interval = coal_int[1]->time;
		} else {
				interval = coal_int[1]->time - (coal_int[1]->next)->time;
		}	
		num_inactive = get_num_inactive(nei, coal_int[1]);
		num_active = 2;
		if((nei->child[1])->time > coal_int[1]->time) {
			num_active++;
		}
		
		length[1] = gen_coal_length1(interval, num_active, num_inactive,
											theta, rand_scratch[1]);
		next[1] = coal_int[1]->next;
		
	}
	(nei->parent)->prev = prev[0];
	(nei->parent)->next = next[0];
	if(NULL != prev[0]) {
		prev[0]->next  = nei->parent;
	}	
	next[0]->prev = nei->parent;
	(nei->target)->prev = prev[1];
	(nei->target)->next = next[1];
	prev[1]->next = nei->target;
	if(NULL != next[1]) {
		next[1]->prev = nei->target;
	}
	if(NULL != next[1]) {
		(nei->target)->time = ((nei->target)->next)->time + length[1];
	} else {
		(nei->target)->time = length[1];
	}
	(nei->parent)->time = ((nei->parent)->next)->time + length[0];
	if(nei->target != (nei->parent)->child1) {
		uncle = &((nei->parent)->child1);
	} else {
		uncle = &((nei->parent)->child2);
	}
//
	if((nei->target)->time > (nei->child[0])->time) {
		uncle_idx = randui % 3;
		seluncle = uncle;
		if(uncle_idx == 1) {
			seluncle = &((nei->target)->child1);
		} else if(uncle_idx == 2) {
			seluncle = &((nei->target)->child2);
		}
	} else if((nei->child[0])->parent == nei->target) {
		if(nei->child[0] == (nei->target)->child1) {
			seluncle = &((nei->target)->child1);
		} else {
			seluncle = &((nei->target)->child2);
		}
	} else {
		seluncle = uncle;
	}
	
	(*uncle)->parent = nei->target;
	(*seluncle)->parent = nei->parent;
	swap = *uncle;
	*uncle = *seluncle;
	*seluncle = swap;
	
}

void upgma_merge(gene_node *node, gene_node **trees, gene_node *next,
						size_t num_trees, size_t seq_len)
{
	
	float dist, min_dist;
	unsigned min_idx1, min_idx2;
	int i, j;
	
	min_dist = DBL_MAX;
	for(i = 0; i < (num_trees - 1); i++) {
		for(j = (i + 1); j < num_trees; j++) {
			dist = find_dist(trees[i], trees[j], seq_len);
			if(dist < min_dist) {
					min_idx1 = i;
					min_idx2 = j;
					min_dist = dist;
			}
		}
	}
	init_gene_node(node, (num_trees - 2), 0);
	node->order = num_trees - 2;
	node->child1 = trees[min_idx1];
	trees[min_idx1]->parent = node;
	node->child2 = trees[min_idx2];
	trees[min_idx2]->parent = node;
	if(NULL != next) {
		node->next = next;
		next->prev = node;
	}
	node->time = min_dist / 2.0;
	trees[min_idx1] = node;
	trees[min_idx2] = trees[(num_trees - 1)];
	
}

__device__ __host__ size_t weighted_pick(float *prob, size_t num, float total, 
										float randf)
{
	
	float threshold;
	int i;
	
	threshold = (1.0 - randf) * total;
	for(i = 0; i < num; i++) {
		if(prob[i] > threshold) {
			return(i);
		}
		threshold -= prob[i];
	}
	
	for(i = 0; i < num; i++) {
		if(prob[i] > 0) {
			return(i);
		}
	}
	
	return(num - 1);
	
}