#ifndef MAX_LIKE_H
#define MAX_LIKE_H

#include<stdlib.h>

#include <curand_kernel.h>

#include "genealogy.h"
#include "sfmt/SFMT.h"

#define BLKSZ 256
#define WRPSZ 32

#define CHAINLEN 100000
#define BURNIN 1000
#define EPSILON .0005
#define MAXITER 1000
#define DELTA 2.2204460492503131e-6
#define MAXJUMP 10.0

typedef struct chain {
	float **samples;
	size_t ncoal;
	size_t nsamp;
	float gentheta;
} chain;

typedef struct multiprop {
	size_t slots;
	size_t num_seq;
	size_t seq_len;
	unsigned geni;
	float *transmtx;
	genealogy *proposals;
	gene_node *nodespace;
	//float *lkspace;
	float *act_scratch;
	float *coal_scratch;
	float *rand_scratch;
	float *blksum_scratch;
	//float ndist[4];
} multiprop;

void do_burnin(multiprop *mp, size_t count, float theta,
						curandStateMtgp32 *mtgp, sfmt_t *sfmt);
void do_chain(multiprop *mp, chain *ch, curandStateMtgp32 *mtgp,
												sfmt_t *sfmt);
void do_multi_prop(multiprop *mp, float inittheta, curandStateMtgp32 *mtgp, 
						sfmt_t *sfmt, float **samples, size_t interval,
						size_t count, unsigned sampling);
void do_sampling(multiprop *mp, chain *ch, curandStateMtgp32 *mtgp,
												sfmt_t *sfmt);
void free_chain(chain *ch);
void free_multiprop(multiprop *mp);
__device__ void get_branch_likelihood(float *lk, gene_node *node, 								multiprop *mp, unsigned seq, 
								unsigned num_seq, unsigned seq_len);
__global__ void get_genealogy_posterior(chain *ch, float theta, 
						float inittheta, float *gllikes, float *blkshift);
float get_theta_likelihood(chain *ch, float theta, float *blockshift, 
									float *gllikes, float *shift);
float gradient_ascent_step(chain *ch, float x0);
chain *init_chain(size_t nsamp, size_t ncoal, float gentheta);
multiprop *init_multiprop(aligned_seqs *seq_dat, float inittheta, size_t slots);
__global__ void lamarc_set_propose(multiprop *mp, float theta, genealogy *gen, 
									unsigned tgt, curandStateMtgp32 *mtgp);
void print_curve_data(chain *ch, float start, float stop, float incr);
__global__ void reduce_genealogy_posterior(size_t num_samp, float *gllikes, 
															float shift);
void seed_upgma(multiprop *mp, aligned_seqs *seq_dat, float thetascale);
__global__ void set_genealogy_llikelihood(genealogy *g, multiprop *mp, 
											unsigned num_seq, unsigned seq_len);

#endif