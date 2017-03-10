#include "upgma.h"

#include<float.h>
#include<stdio.h>
#include<stdlib.h>

btree *new_leaf(void *ptr)
{

	btree *leaf;

	leaf = (btree *)malloc(sizeof(btree));
	leaf->isleaf = 1;
	leaf->leaf = ptr;
	leaf->parent = NULL;
	leaf->left = NULL;
	leaf->right = NULL;
	leaf->depth = 0;

	return(leaf);

}

//TODO: error handling
btree *new_node(btree *left, btree *right, double depth)
{

	btree *node;

	node = (btree *)malloc(sizeof(btree));
	node->isleaf = 0;
	node->leaf = NULL;
	node->parent = NULL;
	node->left = left;
	node->right = right;
	node->depth = depth;

	return(node);

}

//TODO: error handling
size_t leaf_count(btree *b)
{

	if(1 == b->isleaf) {
		return(1);
	} else {
		return(leaf_count(b->left) + leaf_count(b->right));
	}

}

void fill_merge_points(btree *b, double *merge)
{

	double depth, swap;
	unsigned i;
	btree *lchild, *hchild;

	if(0 == b->isleaf) {
		depth = b->depth;
		for(i = 0; 0 != depth; i++) {
			if(depth > merge[i]) {
				swap = merge[i];
				merge[i] = depth;
				depth = swap;
			}		
		}
		if(b->left->depth < b->right->depth) {
			lchild = b->left;
			hchild = b->right;
		} else {
			lchild = b->right;
			hchild = b->left;
		}
		fill_merge_points(hchild, merge+1);
		fill_merge_points(lchild, merge+2);
	}
	

}
/*
unsigned fill_node_list(btree *b, btree **nlist)
{

        unsigned idx;

        if(1 == b->isleaf) {
                return(0);
        }

        idx = 0;
        if(NULL != b->parent) {
                nlist[0] = b;
                idx++;
        }
        idx += fill_node_list(b->left, (nlist + idx));
        idx += fill_node_list(b->right, (nlist + idx));

        return(idx);

}
*/

void fill_node_list(btree *b, btree **nlist)
{

	unsigned i;
	btree *lchild, *hchild;
	btree *swap;

	if(0 == b->isleaf) {
		if(b->left->depth < b->right->depth) {
                        lchild = b->left;
                        hchild = b->right;
                } else {
                        lchild = b->right;
                        hchild = b->left;
                }
		for(i = 0; NULL != b; i++) {
			if(NULL == nlist[i]) {
				nlist[i] = b;
				break;
			} else {
				if(nlist[i]->depth < b->depth) {
					swap = nlist[i];
					nlist[i] = b;
					b = swap;
				}
			}
		}
		fill_node_list(hchild, nlist+1);
                fill_node_list(lchild, nlist+2);
	}

}

void print_newick(btree *b, char *(*l_name)(void *p))
{

	double len;
	btree *parent;

	if(1 == b->isleaf) {
		printf("%s", l_name(b->leaf));
	} else {
		printf("(");
		print_newick(b->left, l_name);
		printf(",");
		print_newick(b->right, l_name);
		printf(")");
	}
	if(NULL != b->parent) {
		parent = b->parent;
                len = parent->depth - b->depth;
		printf(":%.1f", len);
	} else {
		printf(";\n");
	}
		
}

//TODO: error handling
double cluster_dsum(double (*dist_f)(void *p1, void *p2),
			btree *c1, btree *c2)
{

	double dsum;

	dsum = 0;
	if(1 == c1->isleaf) {
		if(1 == c2->isleaf) {
			dsum = dist_f(c1->leaf, c2->leaf);
		} else {
			dsum += cluster_dsum(dist_f, c1, c2->left);
			dsum += cluster_dsum(dist_f, c1, c2->right);
		}
	} else {
		dsum += cluster_dsum(dist_f, c1->left, c2);
		dsum += cluster_dsum(dist_f, c1->right, c2);
	}

	return(dsum);

}

//TODO: error handling
double cluster_dist(double (*dist_f)(void *p1, void *p2),
                        btree *c1, btree *c2)
{

	double dsum;
	size_t npair;

	dsum = cluster_dsum(dist_f, c1, c2);
	npair = leaf_count(c1) * leaf_count(c2);
	
	return(dsum / (double)npair);

}

void upgma_join(double (*dist_f)(void *p1, void *p2), btree **cl,
		size_t nclust)
{

	unsigned i, j;
	double dist, lowdist;
	unsigned lowi, lowj;
	btree *node;

	lowdist = DBL_MAX;
	for(j = 1; j < nclust; j++) {
		for(i = 0; i < j; i++) {
			dist = cluster_dist(dist_f, cl[i], cl[j]);
			if(dist < lowdist) {
				lowdist = dist;
				lowi = i;
				lowj = j;
			}
		}
	}
	node = new_node(cl[lowi], cl[lowj], (lowdist / 2.0));
	cl[lowi]->parent = node;
	cl[lowj]->parent = node;
	cl[lowi] = node;
	cl[lowj] = cl[nclust - 1];

}

btree *upgma(double (*dist_f)(void *p1, void *p2), size_t nleaf, 
		void **lptr)
{

	btree **clusters;
	unsigned i;
	size_t nclust;
	btree *tree;

	clusters = (btree **)malloc(nleaf * sizeof(void *));
	for(i = 0; i < nleaf; i++) {
		clusters[i] = new_leaf(lptr[i]);
	}
	nclust = nleaf;
	while(nclust > 1) {
		upgma_join(dist_f, clusters, nclust);
		nclust--;
	}
	tree = clusters[0];
	free(clusters);

	return(tree);

}
