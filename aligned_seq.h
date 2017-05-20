/*
 * aligned_seq.h
 *
 * Philip Davis, 2015
 * ------------------
 *
 * Provides the aligned sequence data structure aliged_seqs.
 */


#ifndef MLGENE_ALIGN_SEQ_H
#define MLGENE_ALIGN_SEQ_H

typedef enum { A, T, G, C } element;

/* struct: aligned_seqs
 * --------------------
 * A set of aligned sequences.
 *
 * num_seqs: the number of sequences in the set
 * seq_len:  the length, in chars, of each sequence
 * seq_dat:  the sequence data. This will be an array of chars of size
 *           num_seqs * seq_len
 * seqs:     an array of pointers to the seq_dat. seqs[i] will point to
 *           seq_dat[i * seq_len]
 * name_dat: the sequence name data. This will be an array of chars of size
 *           num_seqs * (SEQ_NAME_LEN + 1) the extra char is for the string
 *           terminator
 * names:    an array of pointers to the name_dat. names[i] will point to
 *           name_dat[i * (SEQ_NAME_LEN + 1)]
 */
typedef struct aligned_seqs {
    size_t num_seq;
    size_t seq_len;
    char *seq_dat;
    char **seqs;
    char *name_dat;
    char **names;
} aligned_seqs;

#endif
