
extern "C" {

#ifndef MPCGS_NOGPU

#include <cuda.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "phylip.h"
#include "mpcgs-gpu.h"
#include "tree.h"

__constant__ float lfreq[MPCGS_NUM_FREQ_TERM];

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
                                            unsigned int nsummaries,
                                            float theta)
{

    unsigned int sum_idx;
    struct gtree_summary *sum;

    sum_idx = (blockDim.x * blockIdx.x) + threadIdx.x;

    if (sum_idx < nsummaries) {
        sum = &sum_set->summaries[sum_idx];
        sum->ldrv_posterior = sum_calc_lposterior_gpu(sum, theta);
    }
}

//TODO: rename to distinguish theta likelihood from tree likelihood
__global__ static void set_llhood_comps(struct gtree_summary_set *sum_set,
                                        unsigned int nsummaries,
                                        float theta)
{

    unsigned sum_idx, delta;
    extern __shared__ float warp_normal[];
    struct gtree_summary *sum;
    float lcomp, other_lcomp;
    int i;

    sum_idx = (blockDim.x * blockIdx.x) + threadIdx.x;

    if (sum_idx < nsummaries) {
        sum = &sum_set->summaries[sum_idx];
        lcomp = sum_calc_lposterior_gpu(sum, theta) - sum->ldrv_posterior;
        sum->ltmp_lkhood_comp = lcomp;
    } else {
        lcomp = -FLT_MAX;
    }
    __syncthreads();
    for (delta = (WARPSZ / 2); delta >= 1; delta /= 2) {
        other_lcomp = __shfl_down(lcomp, delta);
        lcomp = fmaxf(lcomp, other_lcomp);
    }
    if ((sum_idx % WARPSZ) == 0) {
        warp_normal[threadIdx.x / WARPSZ] = lcomp;
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        for (i = 0; i < DEF_BLKSZ / WARPSZ; i++) {
            lcomp = fmaxf(lcomp, warp_normal[i]);
        }
        sum_set->block_scratch[blockIdx.x] = lcomp;
    }
}

__global__ static void sum_comp_fracs(struct gtree_summary_set *sum_set,
                                      unsigned int nsummaries,
                                      float normal)
{
    unsigned sum_idx, delta;
    extern __shared__ float warp_sum[];
    struct gtree_summary *sum;
    float comp, other_comp, sum_comp;
    int i;

    sum_idx = (blockDim.x * blockIdx.x) + threadIdx.x;

    if (sum_idx < nsummaries) {
        sum = &sum_set->summaries[sum_idx];
        comp = expf(sum->ltmp_lkhood_comp - normal);
    } else {
        comp = 0;
    }
    __syncthreads();

    for (delta = (WARPSZ / 2); delta >= 1; delta /= 2) {
        other_comp = __shfl_down(comp, delta);
        comp += other_comp;
    }
    if ((sum_idx % WARPSZ) == 0) {
        warp_sum[threadIdx.x / WARPSZ] = comp;
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        sum_comp = 0;
        for (i = 0; i < DEF_BLKSZ / WARPSZ; i++) {
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

    num_block = (size_t)ceil((float)sum_set->nsummaries / (float)DEF_BLKSZ);
    shared_size = (size_t)ceil(DEF_BLKSZ / WARPSZ) * sizeof(float);

    set_llhood_comps<<<num_block, DEF_BLKSZ, shared_size>>>(
      sum_set, sum_set->nsummaries, theta);
    cudaDeviceSynchronize();

    normal = -FLT_MAX;
    for (i = 0; i < num_block; i++) {
        normal = fmaxf(normal, sum_set->block_scratch[i]);
    }

    sum_comp_fracs<<<num_block, DEF_BLKSZ, shared_size>>>(
      sum_set, sum_set->nsummaries, normal);
    cudaDeviceSynchronize();

    lkhood = 0;
    for (i = 0; i < num_block; i++) {
        lkhood += sum_set->block_scratch[i];
    }

    return (logf(lkhood) + normal);
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

    num_block = (size_t)ceil((float)sum_set->nsummaries / (float)DEF_BLKSZ);

    set_base_lposteriors<<<num_block, DEF_BLKSZ>>>(
      sum_set, sum_set->nsummaries, drv_theta);
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
    cudaMallocManaged(&(*sum_set)->summaries,
                      count * sizeof(*((*sum_set)->summaries)));
    if (!(*sum_set)->summaries) {
        // TODO: handle error
    }

    summary = (*sum_set)->summaries;
    for (i = 0; i < count; i++) {
        summary->nintervals = nintervals;
        cudaMallocManaged(&summary->intervals,
                          nintervals * sizeof(*(summary->intervals)));
        if (!summary->intervals) {
            // TODO: handle error
        }
        summary++;
    }

    num_block = (size_t)ceil((float)count / (float)DEF_BLKSZ);
    cudaMallocManaged(&((*sum_set)->block_scratch), num_block * sizeof(float));
    cudaDeviceSynchronize();

}

__device__ static void gnode_get_llhood_lcomps_gpu(struct gene_node *gnode,
                                    float *cllike,
                                    float *lcomps)
{

    float sumAllfact, lsumAll, sumPur, lsumPur, sumPyr, lsumPyr;
    float normal, comp;
    int i;

    /************************************************************
     * The following code adapted from LAMARC, (c) 2002
     * Peter Beerli, Mary Kuhner, Jon  Yamato and Joseph Felsenstein
     * TODO: license ref?
     ***********************************************************/
    normal = -FLT_MAX;
    for (i = 0; i < NUM_BASE; i++) {
        normal = fmaxf(normal, (lfreq[i] + cllike[i]));
    }
    sumAllfact = 0;
    for (i = 0; i < NUM_BASE; i++) {
        if (cllike[i] > -FLT_MAX) {
            sumAllfact += exp((lfreq[i] + cllike[i]) - normal);
        }
    }
    lsumAll = gnode->lexpA + log(sumAllfact) + normal;

    normal =
      fmaxf((lfreq[FREQ_AR] + cllike[DNA_A]), (lfreq[FREQ_GR] + cllike[DNA_G]));
    sumPur = 0;
    if (cllike[DNA_A] > -FLT_MAX) {
        sumPur += exp((lfreq[FREQ_AR] + cllike[DNA_A]) - normal);
    }
    if (cllike[DNA_G] > -FLT_MAX) {
        sumPur += exp((lfreq[FREQ_GR] + cllike[DNA_G]) - normal);
    }

    if (sumPur > 0) {
        lsumPur = log(sumPur) + normal;
    } else {
        lsumPur = -FLT_MAX;
    }

    normal =
      fmaxf((lfreq[FREQ_CY] + cllike[DNA_C]), (lfreq[FREQ_TY] + cllike[DNA_T]));
    sumPyr = 0;
    if (cllike[DNA_C] > -FLT_MAX) {
        sumPyr += exp((lfreq[FREQ_CY] + cllike[DNA_C]) - normal);
    }
    if (cllike[DNA_T] > -FLT_MAX) {
        sumPyr += exp((lfreq[FREQ_TY] + cllike[DNA_T]) - normal);
    }

    if (sumPyr > 0) {
        lsumPyr = log(sumPyr) + normal;
    } else {
        lsumPyr = -FLT_MAX;
    }

    for (i = 0; i < NUM_BASE; i++) {
        // TODO: place these components into an array rather than recalc?

        normal = fmaxf(lsumAll, (gnode->lexpB + cllike[i]));
        normal = fmaxf(
          normal,
          (gnode->lexpC + ((i == DNA_A || i == DNA_G) ? lsumPur : lsumPyr)));

        comp = exp(lsumAll - normal);
        comp += exp((gnode->lexpB + cllike[i]) - normal);
        comp += exp(
          (gnode->lexpC + ((i == DNA_A || i == DNA_G) ? lsumPur : lsumPyr)) -
          normal);
        lcomps[i] = log(comp) + normal;
    }
    /***********************************************************/
}

__global__ void gtree_compute_llhood(struct gene_tree *gtree, unsigned nseq, unsigned seq_len)
{

	extern __shared__ float s[];
	float *warp_sum = s;
	float *llhood_scratch = &warp_sum[(int)ceil((float)blockDim.x/(float)WARPSZ)];
	struct gene_node *node, *child1, *child2;
	float *llhoods, *tip_llhoods, *ch1_llhoods, *ch2_llhoods;
	float ch1_llhood_comp[NUM_BASE], ch2_llhood_comp[NUM_BASE];
	float normal, pos_lhood, llhood, other_llhood, sum_llhood;
	unsigned delta;
	int i, j, node_idx, ch1_idx, ch2_idx, seq_idx, scratch_idx, scratch_idx_ch1, scratch_idx_ch2;

	//Each thread is a base-pair position.
	seq_idx = (blockDim.x * blockIdx.x) + threadIdx.x;
	if(seq_idx < seq_len) {
		//First copy in pre-computed likelihood data from the tips (this will never change over the life of the program.
		for(i = 0; i < nseq; i++) {
			node = &gtree->tips[i];
			node_idx = node->idx;

			scratch_idx = NUM_BASE * ((node_idx * blockDim.x) + threadIdx.x); //try to combine reads as much as possible.

			llhoods = &llhood_scratch[scratch_idx];
			tip_llhoods = &node->tip_llhoods[NUM_BASE * seq_idx];
			for(j = 0 ; j < NUM_BASE; j++) {
				llhoods[j] = tip_llhoods[j];
			}
		}

		//Next walk backwards from the last coalescent event
		for(node = gtree->last; node; node = node->prev) {
			//TODO: see if any of these can be precomputed to save memory accesses
			child1 = node->child1;
			child2 = node->child2;
			node_idx = node->idx;
			ch1_idx = child1->idx;
			ch2_idx = child2->idx;

			scratch_idx = NUM_BASE * ((node_idx * blockDim.x) + threadIdx.x);
			scratch_idx_ch1 = NUM_BASE * ((ch1_idx * blockDim.x) + threadIdx.x);
			scratch_idx_ch2 = NUM_BASE * ((ch2_idx * blockDim.x) + threadIdx.x);

			llhoods = &llhood_scratch[scratch_idx];
			ch1_llhoods = &llhood_scratch[scratch_idx_ch1];
			ch2_llhoods = &llhood_scratch[scratch_idx_ch2];
			gnode_get_llhood_lcomps_gpu(child1, ch1_llhoods, ch1_llhood_comp);
			gnode_get_llhood_lcomps_gpu(child2, ch2_llhoods, ch2_llhood_comp);
			for(i = 0; i < NUM_BASE; i++) {
				llhoods[i] = ch1_llhood_comp[i] + ch2_llhood_comp[i];
			}
		}

		//Now calculate the root likelihood at this base position.
		//after the loop, llhoods should point at the root likelihood components.
		normal = -FLT_MAX;
		for(i = 0; i < NUM_BASE; i++) {
			normal = fmaxf(normal, llhoods[i]);
		}
		pos_lhood = 0;
		for(i = 0; i < NUM_BASE; i++) {
			pos_lhood += expf((llhoods[i] + lfreq[i]) - normal);
		}
		llhood = logf(pos_lhood) + normal;
	} else {
		llhood = 0;
	}
	__syncthreads();

	//Now sum the llhoods across the block
	for (delta = (WARPSZ / 2); delta >= 1; delta /= 2) {
		other_llhood = __shfl_down(llhood, delta);
		llhood += other_llhood;
	}
	if ((seq_idx % WARPSZ) == 0) {
		warp_sum[threadIdx.x / WARPSZ] = llhood;
	}
	__syncthreads();

	if (threadIdx.x == 0) {
		sum_llhood = 0;
		for (i = 0; i < (int)ceil((float)blockDim.x/(float)WARPSZ); i++) {
			sum_llhood += warp_sum[i];
		}
		gtree->block_scratch[blockIdx.x] = sum_llhood;
	}

}

void gtree_set_llhood_gpu(struct gene_tree *gtree)
{

	size_t num_blocks, block_size, shared_size;
	int i = 0;
	float llhood;

	if(!gtree) {
		//TODO: handle error
	}

	block_size = gtree->block_size;
	num_blocks = gtree->num_blocks;
	shared_size = gtree->shared_size;

	gtree_compute_llhood<<<num_blocks, block_size, shared_size>>>
			(gtree, gtree->ntips, gtree->mstab->seq_len);
	cudaDeviceSynchronize();

	llhood = 0;
	for(i = 0; i < num_blocks; i++) {
		llhood += gtree->block_scratch[i];
	}

	gtree->llhood = llhood;

}

static void gnode_add_seq(struct gene_node *tip, struct mol_seq *mseq)
{

	size_t tip_llhood_sz;
	float *llhood;
	int i, j;

	if(!tip || !mseq) {
		//TODO: handle error
	}

	if(tip->child1 || tip->child2) {
		//TODO: warn
	}

	tip_llhood_sz = mseq->len * NUM_BASE * sizeof(*tip->tip_llhoods);
	cudaMallocManaged(&tip->tip_llhoods, tip_llhood_sz);

	llhood = tip->tip_llhoods;
	for(i = 0; i < mseq->len; i++) {
		for(j = 0; j < NUM_BASE; j++) {
			if(j == mseq->seq[i]) {
				llhood[j] = 0;
			} else {
				llhood[j] = -FLT_MAX;
			}
		}
		llhood += NUM_BASE;
	}

	tip->mseq = mseq;

	cudaDeviceSynchronize();

}

void gtree_add_seqs_to_tips_gpu(struct gene_tree *gtree, struct ms_tab *mstab)
{

    int i;
    unsigned int mol_counts[PHY_NUM_MOL_T] = { 0 };
    float freqa, freqg, freqc, freqt; // for readability
    float freqar, freqgr, freqcy, freqty;
    float pur, pyr, ag, ct, m, n, fracchange;
    unsigned int nmol;
    //cudaError_t cucode;

    if (!gtree || !mstab) {
        // TODO: handle error
    }

    if (gtree->ntips != mstab->len) {
        // TODO: handle error
    }

    gtree->mstab = mstab;

    for (i = 0; i < mstab->len; i++) {
        gtree->tips[i].mseq = &mstab->mseq[i];
        gnode_add_seq(&(gtree->tips[i]), &(mstab->mseq[i]));
    }
    nmol = get_mol_counts(mstab, mol_counts);
    if (!nmol) {
        // TODO: handle error
    }

    freqa = (float)mol_counts[DNA_A] / (float)nmol;
    freqt = (float)mol_counts[DNA_T] / (float)nmol;
    freqc = (float)mol_counts[DNA_C] / (float)nmol;
    freqg = (float)mol_counts[DNA_G] / (float)nmol;

    gtree->lfreq[FREQ_A] = logf(freqa);
    gtree->lfreq[FREQ_T] = logf(freqt);
    gtree->lfreq[FREQ_C] = logf(freqc);
    gtree->lfreq[FREQ_G] = logf(freqg);

    /************************************************************
     * The following code adapted from LAMARC, (c) 2002
     * Peter Beerli, Mary Kuhner, Jon  Yamato and Joseph Felsenstein
     * TODO: license ref?
     ***********************************************************/
    pur = freqa + freqg;
    pyr = freqc + freqt;
    if (!pur || !pyr) {
        // TOO: handle error
    }
    freqar = freqa / pur;
    freqgr = freqg / pur;
    freqcy = freqc / pyr;
    freqty = freqt / pyr;
    gtree->lfreq[FREQ_AR] = log(freqar);
    gtree->lfreq[FREQ_GR] = log(freqgr);
    gtree->lfreq[FREQ_CY] = log(freqcy);
    gtree->lfreq[FREQ_TY] = log(freqty);
    ag = freqa * freqg;
    ct = freqc * freqt;
    m = (2.0 * pur * pyr) - (ag + ct);
    n = (ag / pur) + (ct / pyr);
    gtree->yrate = m / (m + n);
    gtree->xrate = 1.0 - gtree->yrate;
    fracchange = gtree->yrate * (2.0 * freqa * freqgr + 2.0 * freqc * freqty) +
                 gtree->xrate * (1.0 - freqa * freqa - freqc * freqc -
                                 freqg * freqg - freqt * freqt);
    gtree->xrate /= -(fracchange);
    gtree->yrate /= -(fracchange);

    /***********************************************************/

    //TODO: error checking (for this and other CUDA calls.)
    cudaMemcpyToSymbol(lfreq, gtree->lfreq, MPCGS_NUM_FREQ_TERM * sizeof(float));

}

void gtree_nodes_init_gpu(struct gene_tree *gtree, size_t ntips, size_t seq_len)
{

    struct gene_node *nodes;
    struct gene_node *prev = NULL;
    size_t num_nodes, nodesSz;
    int i;

    gtree->ntips = ntips;
    gtree->nnodes = ntips - 1;

    num_nodes = gtree->nnodes + gtree->ntips;
    nodesSz = num_nodes * sizeof(*nodes);
    cudaMallocManaged(&nodes, nodesSz);
    memset(nodes, 0, nodesSz);

    //Precompute some likelihood kernel sizing parameters
    gtree->block_size = 2 * DEF_BLKSZ;
    do {
    	gtree->block_size /= 2;
    	gtree->shared_size = (gtree->block_size * num_nodes * NUM_BASE * sizeof(float)) +
    			ceil(((float)gtree->block_size / (float)WARPSZ) * sizeof(float));
    }while(gtree->shared_size > MAX_SHARED_SZ);

    gtree->num_blocks = ceil((float)seq_len/(float)gtree->block_size);

    cudaMallocManaged(&gtree->block_scratch, gtree->num_blocks * sizeof(float));
    cudaDeviceSynchronize();


    for (i = 0; i < gtree->nnodes; i++) {
    	nodes[i].order = nodes[i].idx = i;
    	nodes[i].prev = prev;
    	nodes[i].tree = gtree;
    	if (prev) {
            prev->next = &nodes[i];
        }
        gtree->last = prev = &nodes[i];
    }

    for (i = gtree->nnodes; i < num_nodes; i++) {
        nodes[i].order = -1;
        nodes[i].idx = i;
        nodes[i].tree = gtree;
    }

    gtree->nodes = nodes;
    gtree->tips = &nodes[gtree->nnodes];
    gtree->root = &nodes[0];

}

#endif /* MPCGS_NOGPU */
}
