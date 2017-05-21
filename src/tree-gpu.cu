
extern "C" {

#ifndef MPCGS_NOGPU

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <string.h>

#include "mpcgs-gpu.h"
#include "tree.h"



__device__ static float sum_calc_lposterior_gpu(struct gtree_summary *sum,
		float theta)
{

	int i;
	int lineages;
	float exp1, coeff;

	exp1 = 0;
	for (i = 0; i < sum->nintervals; i++) {
		lineages = i + 2;
		exp1 += -(float)(lineages * (lineages - 1)) * sum->intervals[i];
	}
	coeff = (float)(lineages - 1) * logf(2.0 / theta);

	return (coeff + (exp1 / theta));

}

__global__ static void set_base_lposteriors(struct gtree_summary_set *sum_set,
        unsigned int nsummaries, float theta)
{

	unsigned int sum_idx;
	struct gtree_summary *sum;

	sum_idx = (blockDim.x * blockIdx.x) + threadIdx.x;

	if(sum_idx < nsummaries) {
		sum = &sum_set->summaries[sum_idx];
		sum->ldrv_posterior = sum_calc_lposterior_gpu(sum, theta);
	}

}

__global__ static void set_llhood_comps(struct gtree_summary_set *sum_set,
		unsigned int nsummaries, float theta)
{

	unsigned sum_idx, delta;
	extern __shared__ float warp_normal[];
	struct gtree_summary *sum;
	float lcomp, other_lcomp;
	int i;

	sum_idx = (blockDim.x * blockIdx.x) + threadIdx.x;

	if(sum_idx < nsummaries) {
		sum = &sum_set->summaries[sum_idx];
		lcomp = sum_calc_lposterior_gpu(sum, theta) - sum->ldrv_posterior;
		sum->ltmp_lkhood_comp = lcomp;
	} else {
		lcomp = -FLT_MAX;
	}
	__syncthreads();
	for(delta = (WARPSZ / 2); delta >= 1; delta /= 2) {
		other_lcomp = __shfl_down(lcomp, delta);
		lcomp = fmaxf(lcomp, other_lcomp);
	}
	if((sum_idx % WARPSZ) == 0) {
		warp_normal[threadIdx.x / WARPSZ] = lcomp;
	}
	__syncthreads();
	if(threadIdx.x == 0) {
		for(i = 0; i < BLKSZ / WARPSZ; i++) {
			lcomp = fmaxf(lcomp, warp_normal[i]);
		}
		sum_set->block_scratch[blockIdx.x] = lcomp;
	}

}

__global__ static void sum_comp_fracs(struct gtree_summary_set *sum_set,
		unsigned int nsummaries, float normal)
{
	unsigned sum_idx, delta;
	extern __shared__ float warp_sum[];
	struct gtree_summary *sum;
	float comp, other_comp, sum_comp;
	int i;

	sum_idx = (blockDim.x * blockIdx.x) + threadIdx.x;

	if(sum_idx < nsummaries) {
		sum = &sum_set->summaries[sum_idx];
		comp = expf(sum->ltmp_lkhood_comp - normal);
	} else {
		comp = 0;
	}
	__syncthreads();

	for(delta = (WARPSZ/2); delta >= 1; delta /= 2) {
		other_comp = __shfl_down(comp, delta);
		comp += other_comp;
	}
	if((sum_idx % WARPSZ) == 0) {
		warp_sum[threadIdx.x / WARPSZ] = comp;
	}
	__syncthreads();

	if(threadIdx.x == 0) {
		sum_comp = 0;
		for(i = 0; i < BLKSZ / WARPSZ; i++) {
			sum_comp += warp_sum[i];
		}
		sum_set->block_scratch[blockIdx.x] = sum_comp;
	}

}

float gtree_summary_set_llkhood_gpu(struct gtree_summary_set *sum_set,
                                float theta)
{

	size_t num_block, shared_size;
	float normal, lkhood;
	unsigned i;

	num_block = (size_t)ceil((float)sum_set->nsummaries/(float)BLKSZ);
	shared_size = (size_t)ceil(BLKSZ/WARPSZ) * sizeof(float);

	set_llhood_comps<<<num_block, BLKSZ, shared_size>>> (sum_set, sum_set->nsummaries, theta);
	cudaDeviceSynchronize();

	normal = -FLT_MAX;
	for(i = 0; i < num_block; i++) {
		normal = fmaxf(normal, sum_set->block_scratch[i]);
	}

	sum_comp_fracs <<<num_block, BLKSZ, shared_size>>> (sum_set, sum_set->nsummaries, normal);
	cudaDeviceSynchronize();

	lkhood = 0;
	for(i = 0; i < num_block; i++) {
		lkhood += sum_set->block_scratch[i];
	}

	return(logf(lkhood) + normal);

}


void gtree_summary_set_base_lposteriors_gpu(struct gtree_summary_set *sum_set,
                                        float drv_theta)
{

	size_t num_block;

    if (!sum_set) {
        // TODO: handle error
    }

    if (drv_theta <= 0.0) {
        // TODO: handle error
    }

    num_block = (size_t)ceil((float)sum_set->nsummaries/(float)BLKSZ);

    set_base_lposteriors<<<num_block, BLKSZ>>> (sum_set, sum_set->nsummaries, drv_theta);
    cudaDeviceSynchronize();

}

void gtree_summary_set_create_gpu(struct gtree_summary_set **sum_set,
                              size_t count,
                              size_t nintervals)
{

    int i;
    struct gtree_summary *summary;
    size_t num_block;

    if (!sum_set) {
        // TODO: handle error
    }

    cudaMallocManaged(sum_set, sizeof(**sum_set));

    (*sum_set)->nsummaries = 0;
    (*sum_set)->szintervals = count;
    cudaMallocManaged(&(*sum_set)->summaries, count * sizeof(*((*sum_set)->summaries)));
    if (!(*sum_set)->summaries) {
        // TODO: handle error
    }

    summary = (*sum_set)->summaries;
    for (i = 0; i < count; i++) {
        summary->nintervals = nintervals;
        cudaMallocManaged(&summary->intervals, nintervals * sizeof(*(summary->intervals)));
        if (!summary->intervals) {
            // TODO: handle error
        }
        summary++;
    }

    num_block = (size_t)ceil((float)count/(float)BLKSZ);
    cudaMallocManaged(&((*sum_set)->block_scratch), num_block * sizeof(float));

}

void gtree_set_llhood_gpu(struct gene_tree *gtree)
{



}


void gtree_init_gpu(float theta,
                size_t ntips,
                struct gene_tree *gtree,
                sfmt_t *sfmt)
{


    memset(gtree, sizeof(*gtree), 0);

    gtree_nodes_init(gtree, ntips);
    gtree_simulate_tree(gtree, theta, sfmt);



}


#endif /* MPCGS_NOGPU */

}
