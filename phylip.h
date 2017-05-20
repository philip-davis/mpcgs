/*
 * phylip.h
 * 
 * Philip Davis, 2015
 * ---
 * Provide some utility functions for handling PHYLIP-formatted input files.
 * http://evolution.genetics.washington.edu/phylip/doc/sequence.html
 * This is not a robust implementation, and only handles a subset of PHYLIP
 * files, as documented below.
 *
 */

#ifndef MLGENE_PHYLIP_H

#include<stdio.h>
#include<stdlib.h>

#include"aligned_seq.h"

#define SEQ_NAME_LEN 10
#define FILE_READ_ERROR -1
#define FILE_READ_OK 1

/*
 * free_phylip
 * -----------
 * Deallocate memory for an alligned_seqs struct pointer.
 *
 * as: the aligned_seqs to be freed
 */
void free_phylip(aligned_seqs *as);

/*
 * read_phylip
 * -----------
 * Read a PHYLIP-formatted file into a new aligned-seqs struct.
 * Currently, this will choke on the interleaved format, and expects all
 * sequences to be upper-case nucleotides (A/T/G/C)
 *
 * phylip_file: a file stream pointing to the beginning of a valid PHYLIP-
 *              formatted file.
 *
 * If succesful, return a pointer to a newly allocated aligned_seqs struct 
 * containing the data found in the PHYLIP file. Otherwise, return NULL
 *
 * Prints error messages to stderr.
 */
aligned_seqs *read_phylip(FILE *phylip_file);


/*
 * read_seq
 * --------
 * Read a sequence into an array
 *
 * phylip_file: a files stream pointing to the beggining of a sequence
 * seq:         the array into which to read the sequence
 * slen:        the expected length of the sequence
 *
 * Return FILE_READ_OK is successful, FILE_READ_ERROR otherwise.
 *
 * Prints error messages to stderr.
 */
int read_seq(FILE *phylip_file, char *seq, const size_t slen);

/*
 * read_seq_name
 * -------------
 * Read a sequence name into a string.
 * 
 * phylip_file: a files stream pointing to the beggining of a sequence name
 * seq_name:    the array into which to read the name
 *
 * Return FILE_READ_OK is successful, FILE_READ_ERROR otherwise.
 *
 * Prints error messages to stderr.
 */
int read_seq_name(FILE *phylip_file, char *seq_name);

#endif
