// sequence.c
// Philip Davis, 2015
//
//

#include "sequence.h"

#include<stdio.h>
#include<stdlib.h>

void err_header_read(FILE *fin)
{

	fprintf(stderr, "Could not read header.\n");
	if(fin != NULL) {
		fclose(fin);
	}
	exit(EXIT_FAILURE);

}

void err_seq_missing(FILE *fin)
{

	fprintf(stderr, "Missing sequence data.\n");
	if(fin != NULL) {
        	fclose(fin);
	}
        exit(EXIT_FAILURE);

}

//TODO: support interleaving/leading garbage
void read_seq(FILE *fin, seq *s)
{

	unsigned i;

	for(i = 0; i < 10; i++) {
		s->name[i] = getc(fin);
		if('\n' == s->name[i] || feof(fin)) {
			err_seq_missing(fin);
		}
	}
	s->name[10] = '\0';
	for(i = 9; i > 0; i--) {
		if(' ' != s->name[i]) {
			break;
		} else {
			s->name[i] = s->name[i + 1]; 
		}
	}
	for(i = 0; i < s->len; i++) {
		switch(getc(fin)) {
			case 'A': s->nuc[i] = A; break;
			case 'T': s->nuc[i] = T; break;
			case 'C': s->nuc[i] = C; break;
			case 'G': s->nuc[i] = G; break;
			default:
				err_seq_missing(fin);
		}
	}
	while(getc(fin) != '\n' && feof(fin) == 0);

}

void read_seq_set(FILE *fin, seqset *sset)
{

	size_t nseq, slen;
	unsigned i, j;
	seq *s;
	nuc_t nuc;

	if(fscanf(fin, "%zi %zi\n", &nseq, &slen) != 2) {
		err_header_read(fin);
	}
	sset->nseq = nseq;
	sset->seq = (seq **)malloc(nseq * sizeof(void *));
	for(i = 0; i < nseq; i++) {
		sset->seq[i] = (seq *)malloc(sizeof(seq));
		s = sset->seq[i];
		s->idx = i;
		s->nuc = (char *)malloc(slen * sizeof(char));
		s->len = slen;
		read_seq(fin, s);
		for(j = 0; j < slen; j++) {
			nuc = s->nuc[j];
			sset->nucfreq[nuc]++;
		}
	}
	for(i = 0; i < 4; i++) {
		sset->nucfreq[i] /= (double)(nseq * slen);
	}

}

double seq_dist(seq *s1, seq *s2)
{

	size_t slen;
	double dist;
	unsigned i;

	if((s1 == NULL) || (s2 == NULL)) {
		fprintf(stderr, "Tried to compare NULL sequence.\n");
                exit(EXIT_FAILURE);
        }
	slen = s1->len < s2->len ? s1->len : s2->len;
	dist = 0;
	for(i = 0; i < slen; i++) {
		if(s1->nuc[i] != s2->nuc[i]) {
			dist += .005;
		}
	}

	return(dist);

}
