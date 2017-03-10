#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include "sampleset.h"
#include "seqsmp.h"

#include "sfmt/SFMT.h"
#include "upgma.h"

void calc_nuc_dist(sampleset *smpset)
{

	unsigned i, j;
	size_t nnuc;

	for(i = 0; i < 4; i++) {
                smpset->nucdistr[i] = 0;
        }
        nnuc = 0;
	for(i = 0; i < smpset->nsamp; i++) {
                for(j = 0; j < smpset->dat[i]->len; j++) {
			smpset->nucdistr[(unsigned)smpset->dat[i]->seq[j]]++;
                }
		nnuc += smpset->dat[i]->len;
        }
        for(i = 0; i < 4; i++) {
                smpset->nucdistr[i] /= (double)nnuc;
        }

}

//TODO: error handling
void read_sample_set(FILE *fin, sampleset *smpset)
{

	size_t nsamp, seqlen;
	unsigned i;

	fscanf(fin, "%zi", &nsamp);
        fscanf(fin, "%zi", &seqlen);
        while('\n' != getc(fin));
	smpset->nsamp = nsamp;
        smpset->dat = (seqsmp **)malloc(nsamp * sizeof(void *));
        for(i = 0; i < smpset->nsamp; i++) {
                smpset->dat[i] = (seqsmp *)malloc(sizeof(seqsmp));
                read_seq(fin, seqlen, smpset->dat[i]);
        }
	calc_nuc_dist(smpset);

}

void get_sorted_children(btree *target, btree **c1, btree **c2, btree **c3)
{

	if(target == target->parent->left) {
                *c3 = target->parent->right;
        } else {
                *c3 = target->parent->left;
        }

	if(target->left->depth < (*c3)->depth) {
                *c2 = target->left;
        } else {
                *c2 = *c3;
                *c3 = target->left;
        }
        if(target->right->depth < (*c2)->depth) {
                *c1 = target->right;
        } else {
                if(target->right->depth < (*c3)->depth) {
                        *c1 = *c2;
                        *c2 = target->right;
                } else {
                        *c1 = *c2;
                        *c2 = *c3;
                        *c3 = target->right;
                }
        }

}

//Not much of a sample - not stochastic. Change this later.
void sample_upgma(sampleset *smpset)
{

	size_t ns;
#ifdef DEBUG
	unsigned i;
#endif

	smpset->gen = upgma(seq_dist, smpset->nsamp, (void **)smpset->dat);
	ns = smpset->nsamp;
	smpset->nodelist = (btree **)malloc((ns - 1) * sizeof(void *));
	fill_node_list(smpset->gen, smpset->nodelist);
#ifdef DEBUG
	for(i = 0; i < (ns - 1); i++) {
		printf("%.1f\n", smpset->nodelist[i]->depth);
	}
	printf("\n");
#endif

}


double get_no_coalesc_prob(unsigned ncoal, double t, 
				double j, double z, double theta)
{

	unsigned npair;
	double p;


}

//has side-effects
void sample_posterior(sampleset *smpset, sfmt_t *sfmt)
{

	unsigned tidx;
	//TODO: enum or DEFINE
	btree *neigh[5];
	unsigned i;
	unsigned nstlin, nendlin;

	tidx = (sfmt_genrand_uint32(sfmt) % (smpset->nsamp - 2)) + 1;
	neigh[0] = smpset->nodelist[tidx];
	neigh[1] = neigh[0]->parent;
	get_sorted_children(neigh[0], &neigh[2], &neigh[3], &neigh[4]);

#ifdef DEBUG
	printf("%u\n", tidx);
	printf("%f %f %f\n", neigh[2]->depth, neigh[3]->depth, neigh[4]->depth);
#endif

}

//Felsenstein 1981, eq 7
double prob_sub(nuc_t i, nuc_t j, double t, double *ndistr)
{

	double p;

	p = (i == j) ? exp(-t) : 0.0;
	p += (1.0 - exp(-t)) * ndistr[j];

	return(p);

}

//Felsenstein 1981, eq 4
void site_like(double *ndistr, btree *gen, unsigned site, double *lk)
{

	seqsmp *smp;
	double *leftlk, *rightlk;
	double llsum, rlsum;
	double llen, rlen;
	unsigned i, j;
	double psub;

	if(1 == gen->isleaf) {
		smp = (seqsmp*)gen->leaf;
		//THIS BEING 2 IS A HORRIBLE, BAD KLUDGE
		lk[(unsigned)smp->seq[site]] = 2.0;
#ifdef DEBUG
		if(site == 0) {
			printf("%s:%i\n", gen_label(gen->leaf), smp->seq[site]);
		}
#endif
	} else {
		leftlk = (double *)calloc(4, sizeof(double));
		rightlk = (double *)calloc(4, sizeof(double));
		site_like(ndistr, gen->left, site, leftlk);
		site_like(ndistr, gen->right, site, rightlk);
		llen = gen->depth - gen->left->depth;
		rlen = gen->depth - gen->right->depth;
		for(i = 0; i < 4; i++) {
			llsum = 0;
			rlsum = 0;
			for(j = 0; j < 4; j++) {
				psub = prob_sub((nuc_t)i, (nuc_t)j, llen, ndistr);
				llsum += psub * leftlk[j];
				psub = prob_sub((nuc_t)i, (nuc_t)j, rlen, ndistr);
				rlsum += psub * rightlk[j];
			}
			lk[i] = llsum * rlsum;
		}
		free(leftlk);
		free(rightlk);
#ifdef DEBUG
		if(site == 0) {
			for(i = 0; i < 4; i++) {
				printf("%i: %f,", i, lk[i]);
			}
			printf("\n");
		}
#endif
	}	


}

double gen_llike(sampleset *smpset)
{

	size_t seqlen;
	double llk, lksum;
	double *lk;
	unsigned i, j;

	seqlen = smpset->dat[0]->len;
	llk = 0;
	lk = (double *)malloc(4 * sizeof(double));
	
	for(i = 0; i < seqlen; i++) {
		site_like(smpset->nucdistr, smpset->gen, i, lk);	
		lksum = 0;
		for(j = 0; j < 4; j++) {
			lksum += smpset->nucdistr[j] * lk[j];
		}
		llk += log(lksum);
#ifdef DEBUG
		printf("%f\n", lksum);
#endif
	}
	free(lk);

	return(llk);

}
