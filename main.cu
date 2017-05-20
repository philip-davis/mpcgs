#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<sys/time.h>

#include "aligned_seq.h"
#include "genealogy.h"
#include "maxlike.h"
#include "phylip.h"
#include "sfmt/SFMT.h"

#include <curand_kernel.h>
#include <curand_mtgp32_host.h>

#define BAD_USAGE 1
#define FILE_ERROR 2

#define NUMCHAINS 1
#define CALDSLOTS 100

unsigned long time_diff(struct timeval *start, struct timeval *stop)
{
	
	unsigned long usec;
	
	usec = 1000000 * (stop->tv_sec - start->tv_sec);
	usec += stop->tv_usec;
	usec -= start->tv_usec;
	
	printf("start: %zi, %zi\n", start->tv_sec, start->tv_usec);
	printf("stop: %zi, %zi\n", stop->tv_sec, stop->tv_usec);
	
	return(usec);
	
}

float estimate_theta(aligned_seqs *seq_dat, float inittheta)
{
	
	chain *ch;
	multiprop *mp;
	curandStateMtgp32 *devMTGPState;
    mtgp32_kernel_params *devKernelParams;
	sfmt_t sfmt;
	struct timeval seed;
	unsigned rseed;
	float theta;
	int i;
	
	gettimeofday(&seed, NULL);
	//Some code copied from cuRAND guide.
	
	cudaMalloc((void **)&devMTGPState, sizeof(curandStateMtgp32));
	cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params));
	curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams);
#ifndef RSEED
	rseed = seed.tv_usec;
#else
	rseed = RSEED;
#endif
	curandMakeMTGP32KernelState(devMTGPState, 
                mtgp32dc_params_fast_11213, devKernelParams, 1, rseed);
	sfmt_init_gen_rand(&sfmt, rseed);
	printf("Random seed: %u\n", rseed);
	mp = init_multiprop(seq_dat, inittheta, CALDSLOTS);
	seed_upgma(mp, seq_dat, inittheta);
	ch = init_chain(CHAINLEN,(mp->num_seq - 1), inittheta);
	for(i = 0; i < NUMCHAINS; i++) {
		do_chain(mp, ch, devMTGPState, &sfmt);
		theta = ch->gentheta;
	}
	free_multiprop(mp);
	free_chain(ch);
	
	return(theta);
	
}

void print_usage()
{

    printf("Usage: Placeholder.\n");    

}

int main(int argc, const char *argv[])
{

    FILE *dat_file;
    aligned_seqs *seq_dat;
	float inittheta, theta;
	struct timeval start, stop;
	unsigned long time;
	
	if(argc != 3) {
		print_usage();
		exit(BAD_USAGE);
    }
	
	dat_file = fopen(argv[1], "r");
    if(NULL == dat_file) {
		fprintf(stderr, "Couldn't open %s.\n", argv[1]);
		exit(FILE_ERROR);
    }
    seq_dat = read_phylip(dat_file);
    if(NULL == seq_dat) {
		fprintf(stderr, "Couldn't understand %s.\n", argv[1]);
		exit(FILE_ERROR);
    }
    fclose(dat_file);
	
	gettimeofday(&start, NULL);
	inittheta = atof(argv[2]);
	
	theta = estimate_theta(seq_dat, inittheta);
	gettimeofday(&stop, NULL);
	free_phylip(seq_dat);
	
	printf("Total time for estimation: %u usec.\n", time_diff(&start, &stop));
	printf("Estimated theta: %f\n", theta);
	
	return(0);
	
}