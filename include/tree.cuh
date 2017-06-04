
#ifndef MPCGS_TREE_CUH
#define MPCGS_TREE_CUH

#include <curand_kernel.h>

extern "C" {

#include "tree.h"

__device__ void gtree_propose_fixed_target_gpu(struct gene_tree *current,
                                             struct gene_tree *proposal,
                                             float theta,
                                             unsigned int tgtidx);

__device__ void fill_rand_array(float *rand_array, unsigned count, curandStateMtgp32 *mtgp);
__global__ void gtree_compute_llhood(struct gene_tree *gtree, unsigned nseq, unsigned seq_len);

}

#endif /* MPCGS_TREE_CUH */
