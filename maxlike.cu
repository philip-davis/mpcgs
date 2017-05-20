#include<float.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>

#include<cuda.h>

#include "aligned_seq.h"
#include "genealogy.h"
#include "maxlike.h"

__constant__ float ndist[4];
__constant__ unsigned long seqdat[2000];
//__constant__ gene_node seqnodes[12];

void do_burnin(multiprop *mp, size_t count, float theta,
						curandStateMtgp32 *mtgp, sfmt_t *sfmt)
{
	
	size_t nsample, samplestep;
	
	nsample = count;
	while(nsample > 0) {
		if(nsample >= mp->slots) {
			samplestep = mp->slots;
		} else {
			samplestep = nsample;
		}
		do_multi_prop(mp, theta, mtgp, sfmt, NULL, 1, samplestep, 0);
		nsample -= samplestep;
	}
	
}

void do_chain(multiprop *mp, chain *ch, curandStateMtgp32 *mtgp,
												sfmt_t *sfmt)
{
	
	float xn, xnext;
	unsigned niter;
	struct timeval start, stop;
	
	
	do_burnin(mp, BURNIN, ch->gentheta, mtgp, sfmt);
	do_sampling(mp, ch, mtgp, sfmt);
#ifdef DEBUG
	//print_curve_data(ch, 20.0, 60.0, .01);
#endif
	xnext = ch->gentheta;
	niter = 0;
	do {
		xn = xnext;
#ifdef DEBUGGA
		printf("xnext: %f\n", xnext);
#endif
		gettimeofday(&start, NULL);
		xnext = gradient_ascent_step(ch, xn);
		gettimeofday(&stop, NULL);
		//printf("Gradient ascent time: %u us.\n", stop.tv_usec - start.tv_usec);
		niter++;
	}while((fabs(xnext - xn) > EPSILON) && (niter < MAXITER));
	ch->gentheta = xnext;
#ifdef DEBUGGA
	printf("Chain maxlike: %f\n", xnext);
#endif
	
}

void do_multi_prop(multiprop *mp, float inittheta, curandStateMtgp32 *mtgp, 
						sfmt_t *sfmt, float **samples, size_t interval,
						size_t count, unsigned sampling)
{
	
	genealogy *generator;
	unsigned target;
	struct timeval start, stop;
	size_t Ns;
	size_t num_slots;
	int i;
	float *idist;
	unsigned pick;
	float randf;
	float **smpptriter;
	
	generator = &(mp->proposals[mp->geni]);
#ifdef THROWWRENCH
	target = THROWWRENCH;
#else
	target = sfmt_genrand_uint32(sfmt) % (mp->num_seq - 2);
#endif
#ifdef DEBUG
	print_newick(generator->root);
	printf("Target idx: %u\n", target);
#endif
	num_slots = mp->slots;
	Ns = num_slots * sizeof(float);
	gettimeofday(&start, NULL);
	lamarc_set_propose<<<1, num_slots, Ns>>>(mp, inittheta, generator, target, mtgp);
	cudaDeviceSynchronize();
	gettimeofday(&stop, NULL);
	//printf("Multiprop Time: %i us\n", stop.tv_usec - start.tv_usec); 
	pick = mp->geni;
	//printf("Init reject: %e\n", mp->transmtx[pick * (num_slots + 1)]);
	smpptriter = samples;
	
#ifdef DEBUG
	float likesum = 0;
	for(i = 0; i < num_slots; i++) {
		if(i == mp->geni) {
			printf("mp->geni (%u) llike = %f\n", mp->geni, mp->proposals[i].llike);
		} else {
			likesum += mp->proposals[i].llike;
			printf("Proposal %i like: %f\n", i, mp->proposals[i].llike);
		}
	}
	printf("Average proposal llike = %f\n", (likesum / (double)(mp->slots - 1)));
#endif
	
	for(i = 1; i <= (count * interval); i++) {
		idist = &(mp->transmtx[num_slots * pick]);
		randf = 1.0 - sfmt_genrand_real2(sfmt);
		pick = (unsigned)weighted_pick(idist, num_slots, 1.0, randf);
		if(((i % interval) == 0) && sampling) {
			reduce_genealogy(&mp->proposals[pick], *smpptriter);
			smpptriter++;
		}
#ifdef DEBUG
		printf("%f, %u\n", randf, pick);
#endif
	}
	mp->geni = pick;
	gettimeofday(&stop, NULL);
	
}

void do_sampling(multiprop *mp, chain *ch, curandStateMtgp32 *mtgp,
												sfmt_t *sfmt)
{
	
	size_t nsample, samplestep;
	float **smpptriter;
	
	nsample = ch->nsamp;
	smpptriter = ch->samples;
	while(nsample > 0) {
		if(nsample >= mp->slots) {
			samplestep = mp->slots;
		} else {
			samplestep = nsample;
		}
		do_multi_prop(mp, ch->gentheta, mtgp, sfmt, smpptriter, 1, samplestep, 1);
		nsample -= samplestep;
		smpptriter += samplestep;
	}
	
}

void free_chain(chain *ch)
{
	
	cudaDeviceSynchronize();
	cudaFree(ch->samples[0]);
	cudaFree(ch->samples);
	cudaFree(ch);
	
}

void free_multiprop(multiprop *mp)
{
	
	cudaDeviceSynchronize();
	cudaFree(mp->coal_scratch);
	cudaFree(mp->act_scratch);
	cudaFree(mp->rand_scratch);
	cudaFree(mp->blksum_scratch);
	//cudaFree(mp->lkspace);
	cudaFree(mp->nodespace);
	cudaFree(mp->proposals);
	cudaFree(mp->transmtx);
	cudaFree(mp);
	
}

__device__ void get_branch_likelihood(float *lk, gene_node *node, 								multiprop *mp, unsigned seq, 
								unsigned num_seq, unsigned seq_len)
{
	
	int i;
	gene_node *c1, *c2;
	float c1lk[4], c2lk[4];
	float time, c1time, c2time;
	float c1fact, c2fact;
	float c1base, c2base;
	unsigned long seqlk;
	unsigned seqdatidx;
	
	if(node->isseq) {
		for(i = 0; i < 4; i++) {
			lk[i] = 0.0;
		}
		seqdatidx = (seq_len * (node->idx - (num_seq  - 1))) + seq;
		seqlk = seqdat[seqdatidx / 32];
		seqlk >>= 2 * (seqdatidx % 32);
		seqlk &= 0x3;
		lk[seqlk] = 1.0;
	} else {
		c1 = node->child1;
		c2 = node->child2;
		get_branch_likelihood(c1lk, c1, mp, seq, num_seq, seq_len);
		get_branch_likelihood(c2lk, c2, mp, seq, num_seq, seq_len);
		time = node->time;
		c1time = time - c1->time;
		c2time = time - c2->time;
		c1fact = 1.0 - expf(-c1time);
		c2fact = 1.0 - expf(-c2time);
		c1base = c1fact * 
							((ndist[A] * c1lk[A]) + 
							(ndist[C] * c1lk[C]) + 
							(ndist[G] * c1lk[G]) + 
							(ndist[T] * c1lk[T]));
		c2base = c2fact * 
							((ndist[A] * c2lk[A]) + 
							(ndist[C] * c2lk[C]) + 
							(ndist[G] * c2lk[G]) + 
							(ndist[T] * c2lk[T]));
		for(i = 0; i < 4; i++) {
			lk[i] = (c1base + ((1.0 - c1fact) * c1lk[i])) *
						(c2base + ((1.0 - c2fact) * c2lk[i]));
		}
	}
	
}

__global__ void get_genealogy_posterior(chain *ch, float theta, 
						float inittheta, float *gllikes, float *blkshift)
{
	
	unsigned smpidx;
	extern __shared__ float mpshift[];
	float *sample;
	size_t num_coal;
	int i;
	unsigned lineages;
	float coeff0, coeff1, exp1;
	float lpost, rlpost;
	unsigned delta;
	
	smpidx = (blockDim.x * blockIdx.x) + threadIdx.x;
	num_coal = ch->ncoal;
	if(smpidx < ch->nsamp) {
		exp1 = 0.0;
		sample = ch->samples[smpidx];
		for(i = 0; i < num_coal; i++) {
			lineages = i + 2;
			exp1 += -(float)(lineages * (lineages - 1)) * sample[i];
		}
		coeff0 = powf((2.0 / inittheta), (float)(lineages - 1));
		coeff1 = powf((2.0 / theta), (float)(lineages - 1));
		//lpost = logf(coeff1) + (exp1 / theta);
		lpost = (logf(coeff1) - logf(coeff0)) + (exp1 / theta) - 
					(exp1 / inittheta);
		gllikes[smpidx] = lpost;
	} else {
		lpost = -FLT_MAX;
	}
	__syncthreads();
	for(delta = (WRPSZ / 2); delta >= 1; delta /= 2) {
		rlpost = __shfl_down(lpost, delta);
		lpost = fmaxf(lpost, rlpost);
	}
	if((smpidx % WRPSZ) == 0) {
		mpshift[threadIdx.x / WRPSZ] = lpost;
	}
	__syncthreads();
	if(threadIdx.x == 0) {
		for(i = 0; i < BLKSZ / WRPSZ; i++) {
			lpost = fmaxf(lpost, mpshift[i]);
		}
		blkshift[blockIdx.x] = lpost;
	}
	
}

float get_theta_likelihood(chain *ch, float theta, float *blockshift, 
									float *gllikes, float *shift)
{
	
	float gshift;
	size_t num_samp, numblk, Ns;
	float likelihood;
	int i;
	
	num_samp = ch->nsamp;
	numblk = (size_t)ceil((float)num_samp/(float)BLKSZ);
	Ns = ceil((float)BLKSZ / (float)WRPSZ) * sizeof(float);
	get_genealogy_posterior<<<numblk, BLKSZ, Ns>>>
			(ch, theta, ch->gentheta, gllikes, blockshift);
	cudaDeviceSynchronize();
	gshift = *shift;
	for(i = 0; i < numblk; i++) {
		gshift = max(gshift, blockshift[i]);
	}
	reduce_genealogy_posterior<<<numblk, BLKSZ, Ns>>>
			(ch->nsamp, gllikes, gshift);
	cudaDeviceSynchronize();
	likelihood = 0;
	for(i = 0; i < numblk; i++) {
		likelihood += gllikes[i * BLKSZ];
	}
	likelihood /= (float)num_samp;
	
	*shift = gshift;
	return(likelihood);
	
}

float gradient_ascent_step(chain *ch, float x0)
{
	
	float *blockshift, *gllikes;
	size_t num_samp, numblk;
	float lkxpos[3], shiftset[3];
	float shift;
	int i;
	float theta;
	float diff, xnext, lkxnext;
	
	num_samp = ch->nsamp;
	numblk = (size_t)ceil((float)num_samp/(float)BLKSZ);
	cudaMallocManaged(&blockshift, sizeof(float) * numblk);
	cudaMallocManaged(&gllikes, sizeof(float) * num_samp);
	cudaDeviceSynchronize();
	
	shift = -FLT_MAX;
	for(i = 0; i < 3; i++) {
		theta = x0 + ((float)(2 * (i - 1)) * DELTA);
		lkxpos[i] = get_theta_likelihood(ch, theta, blockshift, 
												gllikes, &shift);
		shiftset[i] = shift;
	}
	for(i = 0; i < 3; i++) {
		lkxpos[i] /= exp(shift - shiftset[i]);
		shiftset[i] = shift;
	}
	diff = (lkxpos[2] - lkxpos[0]) / (4.0 * DELTA);
	do {
#ifdef DEBUGGA
		printf(" diff: %e\n", diff);
#endif
		diff = min(20.0, diff);
		diff /= 2.0;
		xnext = x0 + diff;
		if(xnext <= 0) {
			continue;
		}			
		lkxnext = get_theta_likelihood(ch, xnext, blockshift, gllikes, &shift);
		lkxpos[1] /= exp(shift - shiftset[1]);
		shiftset[1] = shift;
	}while((lkxpos[1] - lkxnext) > DELTA);
	
	cudaFree(blockshift);
	cudaFree(gllikes);
	
	return(xnext);
	
}

chain *init_chain(size_t nsamp, size_t ncoal, float gentheta)
{
	
	chain *ch;
	float *sampledat;
	int i;
	
	cudaMallocManaged(&ch, sizeof(chain));
	cudaDeviceSynchronize();
	cudaMallocManaged(&ch->samples, sizeof(void *) * nsamp);
	ch->ncoal = ncoal;
	ch->nsamp = nsamp;
	ch->gentheta = gentheta;
	cudaMallocManaged(&sampledat, sizeof(float) * ncoal * nsamp);
	cudaDeviceSynchronize();
	for(i = 0; i < nsamp; i++) {
		ch->samples[i] = &sampledat[i * ncoal];
	}
	
	return(ch);
	
}

multiprop *init_multiprop(aligned_seqs *seq_dat, float inittheta, size_t slots)
{
	
	multiprop *mp;
	size_t num_seq, seq_len, nblk;
	size_t nspsz, actscsz, coalscsz, rscsz, bsumsz, slkspsz;
	unsigned long *seqlkspace;
	int i, j;
	gene_node *sharedseq;
	float myndist[4];
	unsigned long lkpos;
	cudaError_t cucode;
	size_t limit;
	
	num_seq = seq_dat->num_seq;
	seq_len = seq_dat->seq_len;
	cucode = cudaDeviceGetLimit(&limit, cudaLimitStackSize);
	if(cucode != cudaSuccess) {
		fprintf(stderr, "An error occurred checking the stack size: %s\n",  			cudaGetErrorString(cucode)); 
	}
	printf("cudaLimitStackSize: %u\n", (unsigned)limit);
	limit = 2048;
	cucode = cudaDeviceSetLimit(cudaLimitStackSize, limit);
	if(cucode != cudaSuccess) {
		fprintf(stderr, "An error occurred checking the stack size: %s\n",  			cudaGetErrorString(cucode)); 
	}
    cucode = cudaDeviceGetLimit(&limit, cudaLimitStackSize);
	if(cucode != cudaSuccess) {
		fprintf(stderr, "An error occurred checking the stack size: %s\n",  			cudaGetErrorString(cucode)); 
	}
	nspsz = ((num_seq - 1) * slots) + num_seq;
	actscsz = (num_seq - 1) * slots * 3;
	coalscsz = (num_seq - 2) * slots * 6;
	rscsz = (num_seq - 2) * slots;
	nblk = ceil((float)seq_len / (float)BLKSZ);
	bsumsz = nblk * slots;
	slkspsz = (num_seq * seq_len) / 4;
	if((slkspsz % sizeof(long)) != 0) {
		slkspsz += sizeof(long) - (slkspsz % sizeof(unsigned long));
	}
	
	if(slkspsz > (300 * sizeof(unsigned long))) {
		fprintf(stderr, "Too much sequence data - will crash due to kludge.\n");
	}
	
	cudaMallocManaged(&mp, sizeof(multiprop));
	cudaDeviceSynchronize();
	mp->slots = slots;
	mp->num_seq = num_seq;
	mp->seq_len = seq_len;
	mp->geni = 0;
	cudaMallocManaged(&(mp->transmtx), sizeof(float) * slots * slots);
	cudaMallocManaged(&(mp->proposals), sizeof(genealogy) * slots);
	cudaMallocManaged(&(mp->nodespace), sizeof(gene_node) * nspsz);
	cudaMallocManaged(&(mp->act_scratch), sizeof(float) * actscsz);
	cudaMallocManaged(&(mp->coal_scratch), sizeof(float) * coalscsz);
	cudaMallocManaged(&(mp->rand_scratch), sizeof(float) * rscsz);
	cudaMallocManaged(&(mp->blksum_scratch), sizeof(float) * bsumsz);
	cudaDeviceSynchronize();
	
	sharedseq = &(mp->nodespace)[slots * (num_seq - 1)];
	for(i = 0; i < slots; i++) {
		(mp->proposals[i]).nodes = &(mp->nodespace)[i * (num_seq - 1)];
		(mp->proposals[i]).seqs = sharedseq;
		(mp->proposals[i]).blksum = &mp->blksum_scratch[i * nblk];
	}
	
	for(i = 0; i < 4; i++) {
		myndist[i] = 0;
	}
	seqlkspace = (unsigned long *)calloc(1, slkspsz);
	for(i = 0; i < num_seq; i++) {
		for(j = 0; j < seq_len; j++) {
			lkpos = seq_dat->seqs[i][j];
			lkpos <<= 2 * (((i * seq_len) + j) % (sizeof(long) * 4));
			seqlkspace[((i * seq_len) + j) / (sizeof(long) * 4)] ^= lkpos;
			myndist[seq_dat->seqs[i][j]]++;
		}
	}
	for(i = 0; i < 4; i++) {
		myndist[i] /= (float)(num_seq * seq_len);
	}
	cucode = cudaMemcpyToSymbol(ndist, myndist, 4 * sizeof(float));
	if(cucode != cudaSuccess) {
		fprintf(stderr, "An error occurred copying to the const space: %s\n",  cudaGetErrorString(cucode)); 
	}
	cudaMemcpyToSymbol(seqdat, seqlkspace, slkspsz);
	cudaDeviceSynchronize();
	free(seqlkspace);
	
	return(mp);
	
}

/*
multiprop *init_multiprop(aligned_seqs *seq_dat, float inittheta, size_t slots)
{
	
	multiprop *mp;
	size_t num_seq, seq_len, nblk;
	size_t lksz, nspsz, actscsz, coalscsz, rscsz, lkidx, bsumsz;
	int i, j;
	float *seqlkspace;
	gene_node *sharedseq;
	
	num_seq = seq_dat->num_seq;
	seq_len = seq_dat->seq_len;
	nspsz = ((num_seq - 1) * slots) + num_seq;
	lksz = nspsz * seq_len * 4;
	actscsz = (num_seq - 1) * slots * 3;
	coalscsz = (num_seq - 2) * slots * 6;
	rscsz = (num_seq - 2) * slots;
	nblk = ceil((float)seq_len / (float)BLKSZ);
	bsumsz = nblk * slots;
	
	cudaMallocManaged(&mp, sizeof(multiprop));
	cudaDeviceSynchronize();
	mp->slots = slots;
	mp->geni = 0;
	mp->num_seq = num_seq;
	mp->seq_len = seq_len;
	
	cudaMallocManaged(&(mp->transmtx), sizeof(float) * slots * slots);
	cudaMallocManaged(&(mp->proposals), sizeof(genealogy) * slots);
	cudaMallocManaged(&(mp->nodespace), sizeof(gene_node) * nspsz);
	cudaMallocManaged(&(mp->lkspace), sizeof(float) * lksz);
	cudaMallocManaged(&(mp->act_scratch), sizeof(float) * actscsz);
	cudaMallocManaged(&(mp->coal_scratch), sizeof(float) * coalscsz);
	cudaMallocManaged(&(mp->rand_scratch), sizeof(float) * rscsz);
	cudaMallocManaged(&(mp->blksum_scratch), sizeof(float) * bsumsz);
	seqlkspace = &(mp->lkspace[(num_seq - 1) * slots * seq_len * 4]);
	cudaDeviceSynchronize();
	cudaMemset(seqlkspace, 0x0, (sizeof(float) * 4 * seq_len * num_seq));
	cudaDeviceSynchronize();
	
	sharedseq = &(mp->nodespace)[slots * (num_seq - 1)];
	for(i = 0; i < slots; i++) {
		(mp->proposals[i]).nodes = &(mp->nodespace)[i * (num_seq - 1)];
		(mp->proposals[i]).seqs = sharedseq;
		lkidx = ((num_seq - 1) * seq_len * 4) * i;
		(mp->proposals[i]).lkspace = &mp->lkspace[lkidx];
		(mp->proposals[i]).blksum = &mp->blksum_scratch[i * nblk];
	}
	for(i = 0; i < 4; i++) {
		mp->ndist[i] = 0;
	}
	for(i = 0; i < seq_len; i++) {
		for(j = 0; j < num_seq; j++) {
			seqlkspace[(4 * ((num_seq * i) + j)) + (seq_dat->seqs[j])[i]] = 1.0;
			mp->ndist[(seq_dat->seqs[j])[i]]++;
		}
	}
	for(i = 0; i < 4; i++) {
		mp->ndist[i] /= (num_seq * seq_len);
	}
	cudaMemcpyToSymbol(ndist, mp->ndist, 4 * sizeof(float));
	//cudaMemcpyToSymbol(seqdat, seqlkspace, mp->num_seq * mp->seq_len * sizeof
	for(i = 0; i < num_seq; i++) {
		sharedseq[i].lkspace = seqlkspace;
	}
	
	return(mp);
	
}
*/

__global__ void lamarc_set_propose(multiprop *mp, float theta, genealogy *gen, 
									unsigned tgt, curandStateMtgp32 *mtgp)
{
	
	unsigned propid;
	genealogy *proposal;
	gene_node *target;
	neighborhood nei;
	float *aprob, *cprob, *rand_scratch;
	unsigned randui;
	size_t num_seq, num_intervals, num_slots;
	gene_node *coal_int[2];
	int i;
	size_t nblk, Ns;
	cudaStream_t cstream;
	float tmtxentry;
	extern __shared__ float proplike[];
	float distsum;
	
	propid = threadIdx.x + blockIdx.x * BLKSZ;
	num_seq = mp->num_seq;
	num_slots = mp->slots;
	
	if(propid != mp->geni) {
		proposal = &(mp->proposals[propid]);
		num_seq = mp->num_seq;
		
		duplicate_genealogy(gen, proposal, num_seq - 1);
		//------------
		//Nothing in this block couldn't be done outside of the kernel.
		//However, if the proposal mechanism changes, that may no longer
		//be true.
		
		target = &proposal->nodes[tgt];
		if(target >= proposal->root) {
			target++;
		}
		
		get_neighborhood(target, &nei);
		
		dislocate_node(proposal, nei.target);
		dislocate_node(proposal, nei.parent);
		aprob = &mp->act_scratch[3 * (num_seq - 1) * propid];
		cprob = &mp->coal_scratch[6 * (num_seq - 2) * propid];
		
		num_intervals = get_num_intervals(proposal, &nei, num_seq);
		populate_a_c_probs(proposal, &nei, theta, aprob, cprob, num_intervals);
		
		//------------
	} else {
		proposal = gen;
	}
	
	__syncthreads();
	rand_scratch = &mp->rand_scratch[((num_seq - 2) * propid)];
	for(i = 0; i < (num_seq - 2); i++) {
		rand_scratch[i] = curand_uniform(mtgp);
	}
	
	if(propid != mp->geni) {
		get_coalescent_intervals(proposal, &nei, aprob, cprob, num_intervals,
												coal_int, rand_scratch);
	}
	
	__syncthreads();
	rand_scratch[0] = curand_uniform(mtgp);
	rand_scratch[1] = curand_uniform(mtgp);
	randui = curand(mtgp);
	if(propid != mp->geni) {	
#ifdef DEBUG
		//printf("propid: %u, randui: %u\n", propid, randui);
#endif
		relocate_nodes(proposal, &nei, coal_int, theta, 
							rand_scratch, randui);
		fix_order(proposal, &nei);
		//foul_likelihood(nei.target);
		__syncthreads();
		
		nblk = ceil((float)(mp->seq_len) / (float)BLKSZ);
		Ns = (BLKSZ / WRPSZ) * sizeof(float);
		cudaStreamCreateWithFlags(&cstream, cudaStreamNonBlocking);
#ifdef SGL2
		set_genealogy_llikelihood<<< nblk, BLKSZ, Ns, cstream >>>(proposal, mp);
#else
		set_genealogy_llikelihood<<< nblk, BLKSZ, Ns, cstream >>>
						(proposal, mp, mp->num_seq, mp->seq_len);
#endif
		cudaDeviceSynchronize();
		proposal->llike = 0;
		for(i = 0; i < nblk; i++) {
			proposal->llike += proposal->blksum[i];
		}	
	}
	
	proplike[propid] = proposal->llike - gen->llike;
	__syncthreads();
	distsum = 0;
	for(i = 0; i < num_slots; i++) {
		if(propid != i) {
			if(proplike[i] >= proplike[propid]) {
				tmtxentry = 1.0;
			} else {
				tmtxentry = expf(proplike[i] - proplike[propid]);
			}
			tmtxentry /= (num_slots - 1);
			distsum += tmtxentry;
			mp->transmtx[(num_slots * propid) + i] = tmtxentry;
		}
	}
	mp->transmtx[(num_slots * propid) + propid] = 1.0 - distsum;
	
}


/*
__global__ void lamarc_set_propose(multiprop *mp, float theta, genealogy *gen, 
									unsigned tgt, curandStateMtgp32 *mtgp)
{
	
	unsigned propid;
	genealogy proposal;
	genealogy *propslot;
	gene_node *target;
	neighborhood nei;
	float *aprob, *cprob, *rand_scratch;
	unsigned randui;
	size_t num_seq, num_intervals, num_slots;
	gene_node *coal_int[2];
	int i;
	size_t nblk, Ns;
	cudaStream_t cstream;
	float tmtxentry;
	extern __shared__ float proplike[];
	float distsum;
	float llike;
	
	propid = threadIdx.x + blockIdx.x * BLKSZ;
	num_seq = mp->num_seq;
	num_slots = mp->slots;
	
	if(propid != mp->geni) {
		proposal = mp->proposals[propid];
		num_seq = mp->num_seq;
		
		duplicate_genealogy(gen, &proposal, num_seq - 1);
		//------------
		//Nothing in this block couldn't be done outside of the kernel.
		//However, if the proposal mechanism changes, that may no longer
		//be true.
		
		target = &proposal.nodes[tgt];
		if(target >= proposal.root) {
			target++;
		}
		
		get_neighborhood(target, &nei);
		
		dislocate_node(&proposal, nei.target);
		dislocate_node(&proposal, nei.parent);
		aprob = &mp->act_scratch[3 * (num_seq - 1) * propid];
		cprob = &mp->coal_scratch[6 * (num_seq - 2) * propid];
		
		num_intervals = get_num_intervals(&proposal, &nei, num_seq);
		populate_a_c_probs(&proposal, &nei, theta, aprob, cprob, num_intervals);
		
		//------------
	} else {
		proposal = *gen;
	}
	
	__syncthreads();
	rand_scratch = &mp->rand_scratch[((num_seq - 2) * propid)];
	for(i = 0; i < (num_seq - 2); i++) {
		rand_scratch[i] = curand_uniform(mtgp);
	}
	
	if(propid != mp->geni) {
		get_coalescent_intervals(&proposal, &nei, aprob, cprob, num_intervals,
												coal_int, rand_scratch);
	}
	
	__syncthreads();
	rand_scratch[0] = curand_uniform(mtgp);
	rand_scratch[1] = curand_uniform(mtgp);
	randui = curand(mtgp);
	if(propid != mp->geni) {	
		relocate_nodes(&proposal, &nei, coal_int, theta, 
							rand_scratch, randui);
		fix_order(&proposal, &nei);
		foul_likelihood(nei.target);
		__syncthreads();
		nblk = ceil((float)(mp->seq_len) / (float)BLKSZ);
		Ns = (BLKSZ / WRPSZ) * sizeof(float);
		cudaStreamCreateWithFlags(&cstream, cudaStreamNonBlocking);
		propslot = &mp->proposals[propid];
		*propslot = proposal;
		
		set_genealogy_llikelihood<<< nblk, BLKSZ, Ns, cstream >>>(propslot, mp);
		llike = 0;
		for(i = 0; i < nblk; i++) {
			llike += propslot->blksum[i];
		)
		
	}
	
	__syncthreads();
	proplike[propid] = llike - gen->llike;
	distsum = 0;
	for(i = 0; i < num_slots; i++) {
		if(propid != i) {
			if(proplike[i] >= proplike[propid]) {
				tmtxentry = 1.0;
			} else {
				tmtxentry = expf(proplike[i] - proplike[propid]);
			}
			tmtxentry /= num_slots;
			distsum += tmtxentry;
	//		mp->transmtx[(num_slots * propid) + i];
		}
	}
	
}
*/

void print_curve_data(chain *ch, float start, float stop, float incr)
{
	
	float *blockshift, *gllikes;
	size_t num_samp, numblk;
	size_t num_data;
	float *ldata;
	float shift;
	float *shiftset;
	int i;
	float theta;
	
	num_samp = ch->nsamp;
	numblk = (size_t)ceil((float)num_samp/(float)BLKSZ);
	cudaMallocManaged(&blockshift, sizeof(float) * numblk);
	cudaMallocManaged(&gllikes, sizeof(float) * num_samp);
	cudaDeviceSynchronize();
	num_data = (size_t)((stop - start) / incr) + 1;
	ldata = (float *)malloc(sizeof(float) * num_data);
	shiftset = (float *)malloc(sizeof(float) * num_data);
	shift = -FLT_MAX;
	for(i = 0; i < num_data; i++) {
		theta = start + (incr * (float)i);
		ldata[i] = get_theta_likelihood(ch, theta, blockshift, gllikes, &shift);
		shiftset[i] = shift;
	}
	for(i = 0; i < num_data; i++) {
		ldata[i] /= exp(shift - shiftset[i]);
		//printf("%f,%f\n", theta, ldata[i]);
		printf("%f\n", ldata[i]);
	}
	
	cudaFree(blockshift);
	cudaFree(gllikes);
	free(ldata);
	free(shiftset);
	
}

__global__ void reduce_genealogy_posterior(size_t num_samp, float *gllikes, float shift)
{
	
	unsigned smpidx;
	float posterior, blocksum;
	extern __shared__ float mppost[];
	unsigned delta;
	int i;
	
	smpidx = (blockDim.x * blockIdx.x) + threadIdx.x;
	
	if(smpidx < num_samp) {
		posterior = expf(gllikes[smpidx] - shift);
	} else {
		posterior = 0;
	}
	
	for(delta = (WRPSZ / 2); delta >= 1; delta /= 2) {
		posterior += __shfl_down(posterior, delta);
	}
	if((smpidx % WRPSZ) == 0) {
		mppost[threadIdx.x / WRPSZ] = posterior;
	}
	blocksum = 0;
	if(threadIdx.x == 0) {
		for(i = 0; i < (BLKSZ / WRPSZ); i++) {
			blocksum += mppost[i];
		}
		gllikes[smpidx] = blocksum;
	}
	
}

void seed_upgma(multiprop *mp, aligned_seqs *seq_dat, float thetascale)
{
	
	genealogy *generator;
	size_t nblk, Ns;
	struct timeval start, stop;
	int i;
	
	generator = &mp->proposals[mp->geni];
	build_gene_upgma(seq_dat, generator, mp->num_seq, mp->seq_len, thetascale);
	nblk = ceil((float)(mp->seq_len) / (float)BLKSZ);
	Ns = (BLKSZ / WRPSZ) * sizeof(float);
	print_newick(generator->root);
	gettimeofday(&start, NULL);
#ifdef SBL2
	set_genealogy_llikelihood<<< nblk, BLKSZ, Ns >>>(generator, mp);
#else
	set_genealogy_llikelihood<<< nblk, BLKSZ, Ns >>>
							(generator, mp, mp->num_seq, mp->seq_len);
#endif
	cudaDeviceSynchronize();
	gettimeofday(&stop, NULL);
	printf("Likelihood Time: %i us\n", stop.tv_usec - start.tv_usec); 
	generator->llike = 0;
	for(i = 0; i < nblk; i++) {
		generator->llike += generator->blksum[i];
	}
	
	printf("Upgma likelihood: %f\n", generator->llike);
	
}


#ifdef SGL2
__global__ void set_genealogy_llikelihood(genealogy *g, multiprop *mp)
{
	
	unsigned seq;
	gene_node *node, *c1, *c2;
	float *nlksp, *c1lksp, *c2lksp;
	float *nlk, *c1lk, *c2lk;
	size_t num_seq, num_nodes;
	float c1time, c2time;
	float c1fact, c2fact;
	float c1base, c2base;
	int i;
	float llike;
	unsigned delta;
	unsigned seqidx;
	extern __shared__ float blkllike[];
	
	seq = (blockDim.x * blockIdx.x) + threadIdx.x;
	
	if(seq < mp->seq_len) {
		num_seq = mp->num_seq;
		num_nodes = mp->num_seq - 1;
	
		for(node = g->lastnode; node != NULL; node = node->prev) {
			if(node->dirty) {
				nlksp = g->lkspace;
				if(seq == 0) {
					node->lkspace = nlksp;
				}
				nlk = &nlksp[((num_nodes * seq) + node->idx) * 4]; 
				c1 = node->child1;
				c2 = node->child2;
				c1lksp = c1->lkspace;
				c2lksp = c2->lkspace;
				if(c1->isseq) {
					seqidx = c1->idx - (num_seq - 1);
					c1lk = &c1lksp[((num_seq * seq) + seqidx) * 4];
				} else {
					c1lk = &c1lksp[((num_nodes * seq) + c1->idx) * 4];
				}
				if(c2->isseq) {
					seqidx = c2->idx - (num_seq - 1);
					c2lk = &c2lksp[((num_seq * seq) + seqidx) * 4];
				} else {
					c2lk = &c2lksp[((num_nodes * seq) + c2->idx) * 4];
				}
				c1time = node->time - c1->time;
				c2time = node->time - c2->time;
				c1fact = 1.0 - expf(-c1time);
				c2fact = 1.0 - expf(-c2time);
				c1base = c1fact * 
							((ndist[A] * c1lk[A]) + 
							(ndist[C] * c1lk[C]) + 
							(ndist[G] * c1lk[G]) + 
							(ndist[T] * c1lk[T]));
				c2base = c2fact * 
							((ndist[A] * c2lk[A]) + 
							(ndist[C] * c2lk[C]) + 
							(ndist[G] * c2lk[G]) + 
							(ndist[T] * c2lk[T]));
				for(i = 0; i < 4; i++) {
					nlk[i] = (c1base + ((1.0 - c1fact) * c1lk[i])) *
								(c2base + ((1.0 - c2fact) * c2lk[i]));
				}
				node->lkspace = nlksp;
				if(seq == 0) {
					node->dirty = 0;
				}
			}
		}
		llike = 0;
		for(i = 0; i < 4; i++) {
			llike += nlk[i] * ndist[i];
		}
		llike = logf(llike);
	} else {
		llike = 0.0;
	}
	for(delta = (WRPSZ / 2); delta >= 1; delta /= 2) {
		llike += __shfl_down(llike, delta);
	}
	if((seq % WRPSZ) == 0) {
		blkllike[seq / WRPSZ] = llike;
	}
	if(threadIdx.x == 0) {
		for(i = 1; i < BLKSZ / WRPSZ; i++) {
			llike += blkllike[i];
		}
		g->blksum[blockIdx.x] = llike;
	}
	
}
#else
__global__ void set_genealogy_llikelihood(genealogy *g, multiprop *mp, 
											unsigned num_seq, unsigned seq_len)
{
	
	unsigned seq;
	float lk[4];
	float llike;
	int i;
	unsigned delta;
	extern __shared__ float blkllike[];
	
	seq = (blockDim.x * blockIdx.x) + threadIdx.x;
	
	if(seq < seq_len) {
		get_branch_likelihood(lk, g->root, mp, seq, num_seq, seq_len);
		llike = 0;
		for(i = 0; i < 4; i++) {
			llike += lk[i] * ndist[i];
		}
		llike = logf(llike);
	} else {
		llike = 0.0;
	}
	for(delta = (WRPSZ / 2); delta >= 1; delta /= 2) {
		llike += __shfl_down(llike, delta);
	}
	if((seq % WRPSZ) == 0) {
		blkllike[seq / WRPSZ] = llike;
	}
	__syncthreads();
	if(threadIdx.x == 0) {
		for(i = 1; i < (BLKSZ / WRPSZ); i++) {
			llike += blkllike[i];
		}
		g->blksum[blockIdx.x] = llike;
	}
	
}
	
#endif