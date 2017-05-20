/*
 * mpcgs-gpu.h
 *
 */

#ifndef MPCGS_GPU_H_
#define MPCGS_GPU_H_

#define BLKSZ 256
#define WARPSZ 32

#ifdef __CDT_PARSER__
#define __global__
#define __device__
#define __shared__
#endif /* __CDT_PARSER__ */

#endif /* MPCGS_GPU_H_ */
