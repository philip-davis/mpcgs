// seqsmp.h
// Philip Davis, 2015
/////

#ifndef SEQSMP_H
#define SEQSMP_H

#include<stdio.h>

typedef enum {A, C, G, T} nuc_t;

typedef struct seqsmp {
	char name[10];
	char *seq;
	size_t len;
} seqsmp;

void read_seq(FILE *fin, size_t seqlen, seqsmp *smp);

double seq_dist(void *s1, void *s2);

void dealloc_seq(seqsmp *smp);

char *gen_label(void *s);

#endif
