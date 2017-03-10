#include<stdio.h>
#include<stdlib.h>
#include"sfmt/SFMT.h"
#include<time.h>

#include "genealogy.h"
#include "sequence.h"

#define BAD_USAGE 2
#define FILE_ERROR 3

void print_usage()
{

        printf("usage: mlgene <PHYLIP>\n\n");

}

int main(int argc, char *argv[])
{

	FILE *fin;
	seqset *sset;
	g7y *root;
	sfmt_t sfmt;

	if(argc != 2) {
                fprintf(stderr, "Bad usage.\n");
                print_usage();
                return(BAD_USAGE);
        }

	fin = fopen(argv[1], "r");
	if(fin == NULL) {
                return(FILE_ERROR);
        }
	sset = (seqset *)malloc(sizeof(seqset));
	read_seq_set(fin, sset);
	root = sample_upgma_g7y(sset);
	print_newick(root);
	get_llike_g7y(root, sset);
	//sfmt_init_gen_rand(&sfmt, time(NULL));
	sfmt_init_gen_rand(&sfmt, 0);
	sample_posterior(root, sset, &sfmt, 1.0); 
	
	free(sset);

	//TODO cleanup
	
	return(0);

}
