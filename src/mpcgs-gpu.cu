
//TODO: make this a configuration option, and conditionally compile/link

#ifndef MPCGS_NOGPU

#include <cuda.h>

#include "tree.cuh"

#include <curand_mtgp32_host.h>

extern "C" {

#include "mpcgs.h"
#include "SFMT.h"
#include "tree.h"
#include "debug.h"
#include "mpcgs-gpu.h"

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

__global__ void propose_multiple(struct multi_proposal *mp, float theta, struct gene_tree *curr_tree, unsigned tgt_idx, unsigned num_prop, unsigned num_seq, unsigned seq_len)
{

	unsigned prop_idx;
	struct gene_tree *proposal;
	extern __shared__ float prop_llhoods[];
	unsigned num_blocks, block_size, shared_size;
	cudaStream_t cstream;
	float trans_mtx_val, trans_mtx_sum;
	float *trans_dist;
	float llhood, other_llhood;
	int i;

	prop_idx = threadIdx.x;

	proposal = &(mp->proposals[prop_idx]);
	// pre-prepare the most random floats the proposal mechanism could need. All threads must
	// run this simultaneously, so do this even for the generator thread.
	fill_rand_array(proposal->rand_scratch, (2 * num_seq), (curandStateMtgp32 *)mp->gtp_states);

	if(prop_idx != mp->curr_idx) {
		gtree_propose_fixed_target_gpu(curr_tree, proposal, theta, tgt_idx);
	}

	__syncthreads();

	block_size = curr_tree->block_size;
	num_blocks = curr_tree->num_blocks;
	shared_size = curr_tree->shared_size;
	cudaStreamCreateWithFlags(&cstream, cudaStreamNonBlocking);
	gtree_compute_llhood<<<num_blocks, block_size, shared_size, cstream>>>
				(proposal, num_seq, seq_len);

	cudaDeviceSynchronize();

	llhood = 0;
	for(i = 0; i < num_blocks; i++) {
		llhood += proposal->block_scratch[i];
	}
	prop_llhoods[prop_idx] = llhood;
	__syncthreads();

	trans_dist = &mp->trans_mtx[prop_idx * num_prop];
	trans_mtx_sum = 0.0;
	for(i = 0; i < num_prop; i++) {
		if(i != prop_idx) {
			other_llhood = prop_llhoods[i];
			if(other_llhood >= llhood) {
				trans_mtx_val = 1.0;
			} else {
				trans_mtx_val = expf(other_llhood - llhood);
			}
			trans_mtx_val /= (float)(num_prop - 1);
			trans_mtx_sum += trans_mtx_val;
			trans_dist[i] = trans_mtx_val;
		}
	}
	trans_dist[prop_idx] = 1.0 - trans_mtx_sum;

}

void do_multi_proposal_gpu(struct chain *ch)
{

	struct chain_param *cparam;
	struct mp_param *mparam;
	struct gene_tree *curr_tree;
	struct gtree_summary *summary;
	unsigned tgt_idx, pick;
	size_t num_blocks, block_size, shared_size;
	float *trans_dist;
	float randf;
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

	propose_multiple<<<num_blocks, block_size, shared_size>>> (ch->mp, ch->theta, curr_tree, tgt_idx, mparam->nproposals, curr_tree->ntips, curr_tree->mstab->seq_len);
	cudaDeviceSynchronize();

	pick = ch->mp->curr_idx;
	for (i = 1; i <= (mparam->npicks * cparam->sum_freq); i++) {
		trans_dist = &ch->mp->trans_mtx[pick * mparam->nproposals];
		randf = 1.0 - sfmt_genrand_real2(&(ch->mp->sfmt));
		pick =
				(unsigned)weighted_pick(trans_dist, mparam->nproposals, 1.0, randf);
		if (mparam->sampling && (i % cparam->sum_freq) == 0) {
			gtree_digest(&ch->mp->proposals[pick], summary);
			ch->sum_set->nsummaries++;
			summary++;
			//printf("proposal->llhood = %f\n", ch->mp->proposals[pick].llhood);
		}
	}

	ch->mp->curr_idx = pick;


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

	//*mp = (multi_proposal *)malloc(sizeof(**mp));
	cudaMallocManaged(mp, sizeof(**mp));

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
