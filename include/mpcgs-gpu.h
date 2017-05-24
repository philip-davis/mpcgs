/*
 * mpcgs-gpu.h
 *
 */

#ifndef MPCGS_GPU_H_
#define MPCGS_GPU_H_

#define DEF_BLKSZ 256 // default blocksize for gradient ascent and multi_prop
                  // kernels
#define WARPSZ 32
#define MAX_SHARED_SZ 49152
#define NUM_MTGP_STATES 256

#ifdef __CDT_PARSER__
#define __global__
#define __device__
#define __shared__
#endif /* __CDT_PARSER__ */

#endif /* MPCGS_GPU_H_ */
