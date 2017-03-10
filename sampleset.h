#ifndef SAMPLESET_H
#define SAMPLESET_H

#include<stdio.h>
#include<stdlib.h>

#include "seqsmp.h"
#include "sfmt/SFMT.h"
#include "upgma.h"

typedef struct sampleset {
	seqsmp **dat;
	size_t nsamp;
	double nucdistr[4];
	btree *gen;
	double *coal;
	btree **nodelist;
} sampleset;

void calc_nuc_dist(sampleset *smpset);

void read_sample_set(FILE *fin, sampleset *smpset);

void get_sorted_children(btree *target, btree **c1, btree **c2, btree **c3);

void sample_upgma(sampleset *smpset);

void sample_posterior(sampleset *smpset, sfmt_t *sfmt);

double prob_sub(nuc_t i, nuc_t j, double t, double *ndistr);

void site_like(double *ndistr, btree *gen, unsigned site, double *lk);

double gen_llike(sampleset *smpset);

#endif
