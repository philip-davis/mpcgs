/*
 * phylip.c
 *
 * Philip Davis, 2015
 *
 */

#include "phylip.h"

void free_phylip(aligned_seqs *as)
{

    if(NULL != as) {
	if(NULL != as->seq_dat) {
	    free(as->seq_dat);
	}
	if(NULL != as->seqs) {
	    free(as->seqs);
	}
	if(NULL != as->name_dat) {
	    free(as->name_dat);
	}
	if(NULL != as->names) {
	   free(as->names);
	}
	free(as);
    }

}

aligned_seqs *read_phylip(FILE *phylip_file)
{

    size_t nseq, slen, sdatlen, snamelen;
    aligned_seqs *as;
    unsigned i;

    if((NULL == phylip_file) || feof(phylip_file)) {
	fprintf(stderr, "Failed to read from the alignment file.\n");
	return(NULL);
    }
  
    if(fscanf(phylip_file, "%zi %zi\n", &nseq, &slen) != 2) {
	fprintf(stderr, "Could not read header from the alignment file.\n");
    	return(NULL);
    }
    as = (aligned_seqs *)malloc(sizeof(aligned_seqs));
    as->num_seq = nseq;
    as->seq_len = slen;
    sdatlen = as->num_seq * as->seq_len;
    as->seq_dat = (char *)malloc(sdatlen * sizeof(char));
    as->seqs = (char **)malloc(as->num_seq * sizeof(void *));
    snamelen = as->num_seq * (SEQ_NAME_LEN + 1);
    as->name_dat = (char *)malloc(snamelen * sizeof(char));
    as->names = (char **)malloc(as->num_seq * sizeof(void *));
    for(i = 0; i < as->num_seq; i++) {
        as->seqs[i] = &as->seq_dat[(i * as->seq_len)];
        as->names[i] = &as->name_dat[(i * (SEQ_NAME_LEN + 1))];
	if(read_seq_name(phylip_file, as->names[i]) != FILE_READ_OK) {
	    free_phylip(as);
	    return(NULL);
	}
	if(read_seq(phylip_file, as->seqs[i], as->seq_len) != FILE_READ_OK) {
	    free_phylip(as);
	    return(NULL);
	}
        if((getc(phylip_file) != '\n') && !feof(phylip_file)) {
	    fprintf(stderr, "Sequence longer than expected.\n");
	    free_phylip(as);
            return(NULL);
        }
    }

    return(as);

}

int read_seq(FILE *phylip_file, char *seq, const size_t slen)
{

    unsigned i;
    int inp;

    for(i = 0; i < slen; i++) {
	inp = getc(phylip_file);
	if(EOF == inp) {
		fprintf(stderr, "Unexpected end of alignment file reached.\n");
        return(FILE_READ_ERROR);
    } else {
	    switch(inp) {
	    	case 'A': seq[i] = (char)A; break;
			case 'T': seq[i] = (char)T; break;
			case 'G': seq[i] = (char)G; break;
			case 'C': seq[i] = (char)C; break;
			default:
				fprintf(stderr, "Unexpected character in sequence\n");
				return(FILE_READ_ERROR);
			}
    	}
    }

    return(FILE_READ_OK);

}

int read_seq_name(FILE *phylip_file, char *seq_name)
{

    unsigned i;
    int inp;

    for(i = 0; i < SEQ_NAME_LEN; i++) {
	inp = getc(phylip_file);
	if(EOF == inp) {
	    fprintf(stderr, "Unexpected end of alignment file reached.\n");
	    return(FILE_READ_ERROR);
	} else {
	    seq_name[i] = (char)inp;
	}
    }      
    seq_name[SEQ_NAME_LEN] = '\0';
    
    return(FILE_READ_OK);     

}
