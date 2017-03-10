#ifndef UPGMA_H
#define UPGMA_H

#include<stdlib.h>

typedef struct btree {

	unsigned isleaf;
	void *leaf;
	unsigned idx;
	struct btree *parent;
	struct btree *left;
	struct btree *right;
	double depth;

} btree;

btree *new_leaf(void *ptr);

size_t leaf_count(btree *b);

void fill_merge_points(btree *b, double *merge);

void fill_node_list(btree *b, btree **nlist);

void print_newick(btree *b, char *(*l_name)(void *p));

btree *new_node(btree *left, btree *right, double depth);

double cluster_dsum(double (*dist_f)(void *p1, void *p2),
                        btree *c1, btree *c2);

double cluster_dist(double (*dist_f)(void *p1, void *p2),
                        btree *c1, btree *c2);

void upgma_join(double (*dist_f)(void *p1, void *p2), btree **cl,
                size_t nclust);

btree *upgma(double (*dist_f)(void *p1, void *p2), size_t nleaf,
                void **lptr);

#endif
