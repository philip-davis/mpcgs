// sequence.h
// Philip Davis, 2015
//
//

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include<stdlib.h>
#include<stdio.h>

typedef enum {A, C, G, T} nuc_t;

//A single gene sample
typedef struct seq {
	unsigned idx;
	char name[11];
	char *nuc;
	size_t len;
} seq;

//A collection of gene samples
typedef struct seqset {
	seq **seq;
	size_t nseq;
	double nucfreq[4];
} seqset;

void err_header_read(FILE *fin);

void err_seq_missing(FILE *fin);

//Read a single line (sample) of a .phy file
void read_seq(FILE *fin, seq *s);

//Read a .phy file
void read_seq_set(FILE *fin, seqset *sset);

//Calculate the distance between two different samples as a scaled count
//of the number of different nucleotides.
double seq_dist(seq *s1, seq *s2);

#endif
