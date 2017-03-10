// seqsmp.c
// Philip Davis, 2015
//////

#include "seqsmp.h"

#include<stdio.h>
#include<stdlib.h>

//TODO: Support interleaving
//TODO: Better error handling
void read_seq(FILE *fin, size_t seqlen, seqsmp *smp)
{

	unsigned i;

	for(i = 0; i < 10; i++) {
		smp->name[i] = getc(fin);
	}
	smp->seq = (char *)malloc(seqlen * sizeof(char));
	for(i = 0; i < seqlen; i++) {
		switch(getc(fin)) {
			case 'A': smp->seq[i] = A; break;
			case 'T': smp->seq[i] = T; break;
			case 'C': smp->seq[i] = C; break;
			case 'G': smp->seq[i] = G; break;
		}
	}
	smp->len = seqlen;
	//cleanup line
	while('\n' != getc(fin));

}

//TODO: error handling
double seq_dist(void *s1, void *s2)
{

	seqsmp *smp1, *smp2;
	double dist;
	unsigned i;

	smp1 = (seqsmp *)s1;
	smp2 = (seqsmp *)s2;
	dist = 0;
	for(i = 0; i < smp1->len; i++) {
		if(smp1->seq[i] != smp2->seq[i]) {
			dist++;
		}
	}
	
	return(dist);

}

void dealloc_seq(seqsmp *smp)
{

	free(smp->seq);

}

char *gen_label(void *s)
{

	seqsmp *smp;
	char *label;
	unsigned i;

	smp = (seqsmp *)s;
	label = malloc(11 * sizeof(char));
	for(i = 0; i < 10; i++) {
		label[i] = smp->name[i];
	}
	for(i = 10; i > 0; i--) {
		if(' ' != label[i - 1]) {
			break;
		}
	}
	label[i] = '\0';

	return(label);

}
