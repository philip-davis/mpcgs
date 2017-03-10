// genealogy.c
// Philip Davis, 2015
//
//

#include<float.h>
#include<math.h>
#include"sfmt/SFMT.h"
#include<stdio.h>
#include<stdlib.h>

#include "genealogy.h"
#include "sequence.h"

g7y *new_node(g7y *lchild, g7y *rchild, double depth)
{

	g7y *node;
	size_t seqlen;

	if((lchild == NULL) || (rchild == NULL)) {
		fprintf(stderr, "New node created with NULL child.\n");
		return(NULL);
	}

	node = (g7y *)malloc(sizeof(g7y));
	node->isleaf = 0;
	node->seq = NULL;
	node->parent = NULL;
	node->left = lchild;
	lchild->parent = node;
	node->right = rchild;
	rchild->parent = node;
	node->prev = NULL;
	node->next = NULL;
	node->depth = depth;
	while(0 == lchild->isleaf) {
		lchild = lchild->left;
	}
	seqlen = lchild->seq->len;
	node->lk = calloc((4 * seqlen), sizeof(double));
	node->lkset = calloc(seqlen, sizeof(char));
	node->llkset = 0;

	return(node);

}

g7y *new_leaf(seq *s)
{

	g7y *leaf;
	unsigned i, j;

	leaf = (g7y *)malloc(sizeof(g7y));
	leaf->isleaf = 1;
	leaf->seq = s;
	leaf->parent = NULL;
	leaf->left = NULL;
	leaf->right = NULL;
	leaf->prev = NULL;
	leaf->next = NULL;
	leaf->depth = 0;
	leaf->lk = malloc((4 * s->len) * sizeof(double));
	leaf->lkset = malloc(s->len * sizeof(char));
	for(i = 0; i < s->len; i++) {
		leaf->lkset[i] = 1;
		for(j = 0; j < 4; j++) {
			if(s->nuc[i] == j) {
				leaf->lk[(i * 4) + j] = 2;
			} else {
				leaf->lk[(i * 4) + j] = 0;
			}
		}
	}
	leaf->llkset = 0;

	return(leaf);

}

void copy_proposal(g7y *g)
{

	g->prparent = g->parent;
	g->prleft = g->left;
	g->prright = g->right;
	g->prprev = g->prev;
	g->prnext = g->next;
	g->prorder = g->order;
	g->active = 0;

}

void prop_orphan(g7y *g)
{

	if(NULL == g) {
		fprintf(stderr, "Tried to orphan NULL.\n");
		exit(EXIT_FAILURE);
	}

	if(NULL != g->parent) {
		if(g->parent->left == g) {
			g->parent->prleft = NULL;
		} else {
			g->parent->prright = NULL;
		}
		g->prparent = NULL;	
	}

}

void prop_order_removal(g7y *g)
{

	if(NULL == g) {
		fprintf(stderr, "Tried to remove NULL from order.\n");
		exit(EXIT_FAILURE);
	}

	if(g->prev != NULL) {
		g->prev->prnext = g->next;
	}
	if(g->next != NULL) {
		g->next->prprev = g->prev;
	}
	while(g->prnext != NULL) {
		g = g->prnext;
		g->prorder--;
	}

}

void erase_active_prop(g7y *target, g7y *parent, g7y **child)
{

	unsigned i;

	for(i = 0; i < 3; i++) {
		prop_orphan(child[i]);
	}
	prop_orphan(target);
	prop_orphan(parent);
	prop_order_removal(target);
	prop_order_removal(parent);
	target->active = 1;
	parent->active = 1;

}

void print_newick(g7y *g)
{

	g7y *parent;
	double len;

	if(NULL != g) {
		if(1 == g->isleaf) {
			printf("%s", g->seq->name);
		} else {
			printf("(");
			print_newick(g->left);
			printf(",");
			print_newick(g->right);
			printf(")");
		}
		if(NULL != g->parent) {
			parent = g->parent;
			len = parent->depth - g->depth;
			printf(":%.1f", len);
		} else {
			printf(";\n");
		}
	}

}

double cluster_dsum_g7y(g7y *c1, g7y *c2, double **dmtx)
{

	double dsum;
	g7y *gswap;

	if((NULL == dmtx) || (NULL == c1) || (NULL == c2)) {
		fprintf(stderr, "Passed NULL comparing clusters.\n");
		exit(EXIT_FAILURE);
	}
	
	dsum = 0;
	if(1 == c1->isleaf) {
		if(1 == c2->isleaf) {
			if(c1->seq->idx > c2->seq->idx) {
				gswap = c1;
				c1 = c2;
				c2 = gswap;
			}		
			dsum = seq_dist(c1->seq, c2->seq);	
		} else {
			dsum += cluster_dsum_g7y(c1, c2->left, dmtx);
			dsum += cluster_dsum_g7y(c1, c2->right, dmtx);
		}
	} else {
		dsum += cluster_dsum_g7y(c1->left, c2, dmtx);
		dsum += cluster_dsum_g7y(c1->right, c2, dmtx);
	}

	return(dsum);

}

double cluster_dist_g7y(g7y **cl, size_t *cllen, unsigned i, 
					unsigned j, double **dmtx)
{

	double dsum;

	if((NULL == cl) || (NULL == cllen) || (NULL == dmtx)) {
		fprintf(stderr, "Received NULL comparing clusters.\n");
		exit(EXIT_FAILURE);
	}

	
	dsum = cluster_dsum_g7y(cl[i], cl[j], dmtx);
	return(dsum / (double)(cllen[i] * cllen[j]));	

}

g7y *upgma_join_g7y(g7y **cl, size_t ncl, size_t *cllen, double **dmtx)
{

	unsigned i, j;
	double dist, lowdist;
	unsigned lowi, lowj;
	g7y *root;

	root = NULL;
	if((NULL != dmtx) && (NULL != cl) && (ncl > 1)) {
		lowdist = DBL_MAX;
		for(j = 1; j < ncl; j++) {
			for(i = 0; i < j; i++) {
				dist = cluster_dist_g7y(cl, cllen, i, j, dmtx);
				if(dist < lowdist) {
					lowdist = dist;
					lowi = i;
					lowj = j;
				}
			}
		}
		root = new_node(cl[lowi], cl[lowj], (lowdist / 2.0));
		cl[lowi] = root;
		cllen[lowi] += cllen[lowj];
		cl[lowj] = cl[ncl - 1];
		cllen[lowj] = cllen[ncl - 1];
	}

	return(root);

}

g7y *sample_upgma_g7y(seqset *sset)
{

	g7y **clusters;
	unsigned i, j;
	size_t ncl;
	size_t *cllen;
	double **dmtx;
	g7y **tarray;
	g7y *prev, *next;
	g7y *root;

	if(NULL == sset) {
		return(NULL);
	}
	ncl = sset->nseq;
	clusters = (g7y **)malloc(ncl * sizeof(void *));
	cllen  = (size_t *)malloc(ncl * sizeof(size_t));
	dmtx = (double **)malloc(ncl * sizeof(void *));
	for(i = 0; i < ncl; i++) {
		clusters[i] = new_leaf(sset->seq[i]);
		cllen[i] = 1;
		dmtx[i] = (double *)malloc(i * sizeof(double));
		for(j = 0; j < i; j++) {
			dmtx[i][j] = seq_dist(sset->seq[i], sset->seq[j]);
		}
	}
	tarray = (g7y **)malloc((ncl - 1) * sizeof(void *));
	next = NULL;
	while(ncl > 1) {
		prev = upgma_join_g7y(clusters, ncl, cllen, dmtx);
		prev->next = next;
		if(NULL != next) {
			next->prev = prev;
			next->prprev = prev;
		}
		next = prev;
		tarray[ncl - 1] = prev;
		prev->tarray = &tarray[ncl - 1];
		prev->order = ncl - 1;
		copy_proposal(prev);
		ncl--;	
	}

	root = clusters[0];
#ifdef DEBUG
	next = root;
	while(NULL != next) {
		printf("%f,", next->depth);
		next = next->next;
	}
	printf("\n");
	printf("%f;%f\n", root->tarray[1]->depth, root->next->depth);
#endif
	free(clusters);
	free(cllen);
	free(dmtx);

	return(root);

}

//Felsenstein 1981, eq7
double prob_sub(nuc_t i, nuc_t j, double t, double *nucfreq)
{

        double p;

        p = (i == j) ? exp(-t * 200) : 0.0;
        p += (1.0 - exp(-t * 200)) * nucfreq[j];

        return(p);

}


double *upd_like_g7y(g7y *g, unsigned pos, double *nucfreq)
{

	unsigned i, j;
	double *lk;
	double *llk, *rlk;
	double llksum, rlksum;
	double llen, rlen;
	double psub;

	if((NULL == g) || (nucfreq == NULL)) {
		fprintf(stderr, "Received a NULL value for P(D|G).\n");
		exit(EXIT_FAILURE);
	}

	lk = &g->lk[4 * pos];
	if(0 == g->lkset[pos]) {
		if(1 == g->isleaf) {
			for(i = 0; i < 4; i++) {
				if(g->seq->nuc[pos] == i) {
					//HACK
					lk[i] = 2;
				} else {
					lk[i] = 0;
				}
			}	
		} else {
			llk = upd_like_g7y(g->left, pos, nucfreq);
			rlk = upd_like_g7y(g->right, pos, nucfreq);
			llen = g->depth - g->left->depth;
			rlen = g->depth - g->right->depth;
			for(i = 0; i < 4; i++) {
				llksum = 0;
				rlksum = 0;
				for(j = 0; j < 4; j++) {
					psub = prob_sub(i, j, llen, nucfreq);
					llksum += psub * llk[j];
					psub = prob_sub(i, j, rlen, nucfreq);
					rlksum += psub * rlk[j];
				}
				lk[i] = llksum * rlksum;			
			}
		}
	}

	return(lk);
	
}

double get_llike_g7y(g7y *root, seqset *sset)
{

	double *nucfreq;
	size_t seqlen;
	unsigned i, j;
	double *lk;
	double lksum;
	double llk;

	if((NULL == root) || (NULL == sset)) {
		fprintf(stderr, "Received a NULL value for P(D|G).\n");
                exit(EXIT_FAILURE);	
	}
	if(1 == root->llkset) {
		return(root->llk);
	}

	nucfreq = sset->nucfreq;
	seqlen = sset->seq[0]->len;
	llk = 0;
	for(i = 0; i < seqlen; i++) {
		lk = upd_like_g7y(root, i, nucfreq);
		lksum = 0;
		for(j = 0; j < 4; j++) {
			lksum += nucfreq[j] * lk[j];
		}
		llk += log(lksum);
		printf("%f\n", lksum);
	}

	//HACK
	llk -= sset->nseq * log(2);
	//END HACK
	root->llkset = 1;
	root->llk = llk;	

	return(llk);

}

//Gross, but there will only ever be three children, so...
void get_ordered_children(g7y *target, g7y *parent, g7y **child)
{

	if((NULL == target) || (NULL == parent) || (NULL == child)) {
		fprintf(stderr, "Received NULL collecting children.\n");
		exit(EXIT_FAILURE);
	}

	if(target == parent->left) {
		child[2] = parent->right;
	} else {
		child[2] = parent->left;
	}
	
	if(target->left->depth <= child[2]->depth) {
		child[1] = target->left;
	} else {
		child[1] = child[2];
		child[2] = target->left;
	}
	if(target->right->depth <= child[1]->depth) {
		child[0] = target->right;
	} else {
		child[0] = child[1];
		if(target->right->depth <= child[2]->depth) {
			child[1] = target->right;
		} else {
			child[1] = child[2];
			child[2] = target->right;
		}
	}

}

//Kuhner 1995, eq 8
double prob_no_coalesc(unsigned x, double t, unsigned z, double theta)
{

	int expcoeff;
	double expmult;
	double expon;

	expcoeff = (x * (x - 1)) + (2 * x * z);
	expmult = -(t / theta);
	expon = (double)expcoeff * expmult;	

	return(exp(expon)); 

}

//Kuhner 1995, eq 9
double prob_one_coalesc(unsigned x, double t, unsigned z, double theta)
{

	int y;
	int coeffnum, coeffden;
	double coeff, mult;
	int expcoeff1, expcoeff2;
	double expmult;
	double expon1, expon2;

	y = x - 1;
	coeffnum = x * y;
	coeffden = (2 * z) + (2 * y);
	coeff = (double)coeffnum / (double)coeffden;
	expcoeff1 = ((2 * z) * y) + (y * (y - 1));
	expcoeff2 = ((2 * z) * x) + (x * y);
	expmult = -(t / theta);
	expon1 = expcoeff1 * expmult;
	expon2 = expcoeff2 * expmult;
	mult = exp(expon1) - exp(expon2);

	return(coeff * mult);

}

//Kuhner 1995, eq 10
double prob_two_coalesc(unsigned x, double t, unsigned z, double theta)
{

	int x1, x2, x3;
	int coeffnum, coeffden;
	double coeff, coeff1, coeff2;
	double expcoeff1, expcoeff2, expcoeff3, expcoeff4;
	double expmult;
	double expon1, expon2, expon3, expon4;
	double mult1, mult2;

	x1 = x - 1;
	x2 = x - 2;
	x3 = x - 3;
	coeffnum = x * (x1 * x1) * x2;
	coeffden = (2 * z) + (2 * x2);
	coeff = (double)coeffnum / (double)coeffden;
	coeff1 = 1.0 / (double)((4 * z) + (4 * x) - 6);
	coeff2 = 1.0 / (double)((2 * z) + (2 * x1));
	expcoeff1 = (2 * z * x2) + (x2 * x3);
	expcoeff2 = (2 * z * x) + (x * x1);
	expcoeff3 = (2 * z * x1) + (x1 * x2);
	expcoeff4 = (2 * z * x) + (x * x1);
	expmult = -(t / theta);
	expon1 = expcoeff1 * expmult;
	expon2 = expcoeff2 * expmult;
	expon3 = expcoeff3 * expmult;
	expon4 = expcoeff4 * expmult;
	mult1 = coeff1 * (exp(expon1) - exp(expon2));
	mult2 = coeff2 * (exp(expon3) - exp(expon4));

	return(coeff * (mult1 - mult2));

}

double prob_act_coalesc(unsigned x, unsigned y, double t, 
					unsigned z, double theta)
{

	switch(x - y) {
		case 0: return(prob_no_coalesc(x, t, z, theta));
		case 1: return(prob_one_coalesc(x, t, z, theta));
		case 2: return(prob_two_coalesc(x, t, z, theta));
		default: return(0);
	}

}

g7y *sample_posterior(g7y *root, seqset *sset, sfmt_t *sfmt, double theta)
{

	g7y *target, *parent, **child;
	unsigned tgtid;
	unsigned n, z, i, j, k;
	unsigned c2int, c3int, baseint;
	double *practlin;
	g7y *giter;
	double *prevint, *practint;
	unsigned nlinpos, linmod;
	double t;

	if((NULL == root) || (NULL == sset) || (NULL == sfmt)) {
		fprintf(stderr, "Received NULL sampling the posterior.\n");
		exit(EXIT_FAILURE);
	}

	tgtid = (sfmt_genrand_uint32(sfmt) % (sset->nseq - 2)) + 1;
	target = root->tarray[tgtid];
	parent = target->parent;
	child = (g7y **)malloc(3 * sizeof(void *));
	get_ordered_children(target, parent, child);
	erase_active_prop(target, parent, child);
	n =  sset->nseq - 2;
	if(1 == child[1]->isleaf) {
		c2int = 0;
		giter = root->tarray[sset->nseq - 2];
		while(1 == giter->active) {
			giter = giter->prev;
		}
	} else {
		c2int = n - child[1]->prorder;
		giter = child[1]->prev;
		while(1 == giter->active) {
			giter = giter->prev;
		}
	}
	if(1 == child[2]->isleaf) {
		c3int = 0;
	} else {
		c3int = n - child[2]->prorder;
	}
	if(parent->parent == NULL) {
		baseint = n;
	} else {
		baseint = n - parent->parent->order;
	}
	practlin = (double *)calloc((3 * (baseint - c2int)), sizeof(void *));
	practint = practlin;
	if(c2int == c3int) {
		practint[2] = 1;
	} else {
		practint[1] = 1;
	}
	for(i = 1; i <= (baseint - c2int); i++) {
		prevint = practint;
		practint += 3;
		if(i < c3int) {
			nlinpos = 2;
			z = n - i + 1;
		} else {
			nlinpos = 3;
			z = n - i;
		}
		if(i == c3int) {
			linmod = 1;
		} else {
			linmod = 0;
		}
		for(j = linmod; j < nlinpos; j++) {
			for(k = 0; k < (3 - linmod); k++) {
				if(giter->prnext == NULL) {
					t = giter->depth;
				} else {
					t = giter->depth - giter->prnext->depth;
				}
				practint[j] += prevint[k] * prob_act_coalesc(((k + 1) + linmod), (j + 1), t, z, theta);		
			}	
		}
		giter = giter->prprev;
	}
#ifdef DEBUG
	printf("Target: %i\n", tgtid);
	practint = practlin;
	for(i = 0; i <= (baseint - c2int); i++) {
		printf("Interval: %2i: %e %e %e\n", i, practlin[0], practlin[1], practlin[2]);
		practlin += 3;
	}
#endif

	return(root);

} 	
