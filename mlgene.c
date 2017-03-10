#include "sampleset.h"
#include "seqsmp.h"
#include "upgma.h"

#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#define BAD_USAGE 2
#define FILE_ERROR 3

void print_usage()
{

	printf("usage: mlgene <PHYLIP>\n\n");	

}

int main(int argc, char *argv[])
{

	FILE *fin;
	sampleset *smps;
	sfmt_t sfmt;

	if(argc != 2) {
		fprintf(stderr, "Bad usage.\n");
		print_usage();
		return(BAD_USAGE);
	}

	//TODO: error handling
	fin = fopen(argv[1], "r");
	if(fin == NULL) {
		return(FILE_ERROR);
	}
	smps = (sampleset *)malloc(sizeof(sampleset));
	read_sample_set(fin, smps);
	sample_upgma(smps);
	gen_llike(smps);
	print_newick(smps->gen, gen_label);
	sfmt_init_gen_rand(&sfmt, time(NULL));
	sample_posterior(smps, &sfmt);

	return(0);

}
