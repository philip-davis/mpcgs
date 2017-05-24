
//TODO: make this a configuration option, and conditionally compile/link

#ifndef MPCGS_NOGPU

#include <cuda.h>
#include <curand_kernel.h>
#include <curand_mtgp32_host.h>

extern "C" {

#include "mpcgs.h"
#include "SFMT.h"
#include "tree.h"
#include "debug.h"
#include "mpcgs-gpu.h"
#include "tree.cuh"

#include<stdlib.h>
#include<sys/time.h>

static void *init_prng_gpu(unsigned seed)
{

	struct timeval curr_time;
	curandStateMtgp32 *devMTGPStates;
	mtgp32_kernel_params *devKernelParams;

	if(seed == 0) {
		gettimeofday(&curr_time, NULL);
		seed = curr_time.tv_usec;
	}

	//Some code adapted from NVidia's cuRAND guide.
	cudaMalloc((void **)&devMTGPStates, sizeof(curandStateMtgp32));
	cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params));
	curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams);
	curandMakeMTGP32KernelState(devMTGPStates,
	                mtgp32dc_params_fast_11213, devKernelParams, 1, seed);

	return(devMTGPStates);

}

__global__ void propose_multiple(struct multi_proposal *mp, float theta, struct gene_tree *curr_tree, unsigned tgt_idx, unsigned num_seq)
{

	unsigned prop_idx;
	struct gene_tree *proposal;
	extern __shared__ float prop_llhoods[];

	prop_idx = threadIdx.x;

	if(prop_idx != mp->curr_idx) {
		proposal = &(mp->proposals[prop_idx]);
	}

}

static void do_multi_proposal_gpu(struct chain *ch)
{

	struct chain_param *cparam;
	struct mp_param *mparam;
	struct gene_tree *curr_tree;
	struct gtree_summary *summary;
	struct gtree_summary_set *sum_set;
	unsigned tgt_idx;
	size_t num_blocks, block_size, shared_size;
	int i;

	if (!ch) {
	        // TODO: handle error
	}

	cparam = ch->cparam;
	mparam = ch->mp->mparam;

	curr_tree = &ch->mp->proposals[ch->mp->curr_idx];
	tgt_idx =
	      sfmt_genrand_uint32(&(ch->mp->sfmt)) % ((curr_tree->nnodes + curr_tree->ntips) - 1);
	if (tgt_idx == curr_tree->root->idx) {
		tgt_idx++;
	}

	num_blocks = 1;
	block_size = mparam->nproposals;
	if(block_size > NUM_MTGP_STATES) {
		//TODO: handle error
	}
	shared_size = block_size * sizeof(float);

	propose_multiple<<<num_blocks, block_size, shared_size>>> (ch->mp, ch->theta, curr_tree, tgt_idx, curr_tree->ntips);


}

unsigned multi_prop_init_gpu(struct multi_proposal **mp,
                             struct ms_tab *data,
                             unsigned nproposal,
							 float theta,
							 unsigned seed)
{

	struct gene_tree *proposal, *init_tree;
	sfmt_t *sfmt;
	int i;

	if (!mp || !data) {
		// TODO: handle error
	}

	*mp = (multi_proposal *)malloc(sizeof(**mp));
	cudaMallocManaged(mp, sizeof(*mp));

	cudaMallocManaged(
			&(*mp)->proposals, nproposal * sizeof(*((*mp)->proposals)));
	proposal = (*mp)->proposals;
	for(i = 0; i < nproposal; i++) {
		gtree_nodes_init_gpu(proposal, data->len, data->seq_len);
		proposal++;
	}

	cudaMallocManaged(
			&(*mp)->trans_mtx,
			nproposal * nproposal * sizeof(*((*mp)->trans_mtx)));


	sfmt = &((*mp)->sfmt);
	sfmt_init_gen_rand(sfmt, seed);
	(*mp)->gtp_states = init_prng_gpu(seed);

	(*mp)->curr_idx = 0;
	init_tree = (*mp)->proposals;
	gtree_simulate_tree(init_tree, theta, sfmt);
	gtree_add_seqs_to_tips_gpu(init_tree, data);
	gtree_set_exp(init_tree);
	gtree_print_newick(init_tree);
	gtree_set_llhood_gpu(init_tree);
	log_debug("init tree root time: %f\n", init_tree->root->time);
	log_debug("init tree log likelihood: %f\n", init_tree->llhood);

	return (init_tree->nnodes);

}

}

#endif
