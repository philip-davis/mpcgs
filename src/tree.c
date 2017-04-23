/*
    Copyright (C) 2017  Philip Davis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include<errno.h>
#include<float.h>
#include<math.h>
#include<sys/time.h>

#include "phylip.h"
#include "tree.h"
#include "debug.h"

#include "SFMT.h"

float get_next_coal_time(unsigned int active, unsigned int inactive, 
						unsigned int act_coal, float theta, sfmt_t *sfmt)
{
	
	float r, time, denom;
	
	debug_var_decl("allocating tree nodes");
	
	param_chk(sfmt && act_coal <= 2 && act_coal > 0, errout);
	
	if(active < act_coal || inactive < 2 - act_coal) {
		return(FLT_MAX);
	}
	
	r = sfmt_genrand_real3(sfmt); // Uniform on interval (0,1)
	if(2 == act_coal) {
		denom = ((float)active * (float)(active - 1)) / 2.0; 
	} else if(1 == act_coal) {
		denom = (float)active * (float)inactive;
	} else {
		//show be caught by param_chk
		return(FLT_MAX);
	}
	
	errno = 0;
	time = -(log(r) * theta) / (2.0 * denom);
	if(errno || time == 0.0) {
		err_warn("Random number generator returned an invalid value.\n");
	}
	return(time);
	
errout:
	debug_err_out();
	
}

static void gtree_nodes_init(struct gene_tree *gtree)
{
	
	struct gene_node *nodes;
	struct gene_node *prev = NULL;
	int i;
	
	debug_var_decl("allocating tree nodes");
	
	param_chk(gtree, errout);
	
	gtree->nnodes = gtree->ntips - 1;
	
	nodes = calloc(gtree->nnodes, sizeof(*nodes)); 
	alloc_chk(nodes, errout);
	
	for(i = 0; i < gtree->nnodes; i++) {
		nodes[i].order = nodes[i].idx = i;
		nodes[i].prev = prev;
		nodes[i].tree = gtree;
		if(prev) {
			prev->next = &nodes[i];
		}
		gtree->last = prev = &nodes[i];
		
	}
	
	gtree->nodes = nodes;
	gtree->root = &nodes[0];
	
	return;
	
errout:
	debug_err_out();
	
}

static void gtree_tips_init(struct gene_tree *gtree, size_t ntips)
{
	
	struct gene_node *tips;
	int i;
	
	debug_var_decl("allocating tree tips");
	
	param_chk(gtree && ntips > 0, errout);
	
	tips = calloc(ntips, sizeof(*tips));
	alloc_chk(tips, errout);
	
	for(i = 0; i < ntips; i++) {
		tips[i].order = -1;
		tips[i].idx = (ntips - 1) + i; // tips come after (ntips - 1) nodes
		tips[i].tree = gtree;
	}
	
	gtree->tips = tips;
	gtree->ntips = ntips;
	
	return;
	
errout:
	debug_err_out();
	
}

static void add_child(struct gene_node *parent, struct gene_node *child) {
	
	if(!parent || !child) {
		//TODO: handle error
	}
	
	if(!parent->child1) {
		parent->child1 = child;
	} else if(!parent->child2) {
		parent->child2 = child;
	} else {
		//TODO: handle error
	}
	
	child->parent = parent;
	child->exp_valid = 0;
	
}

static void gtree_simulate_tree(struct gene_tree *gtree, float theta)
{
	
	sfmt_t sfmt;
	struct timeval seed;
	struct gene_node **branches;
	struct gene_node *node, **child1, **child2;
	int nbranch;
	float current = 0, coal_time;
	int i;
	
	debug_var_decl("simulating gene tree");
	
	param_chk(gtree && theta > 0.0, errout);

	gettimeofday(&seed, NULL);
	log_debug("rseed = %li\n", seed.tv_usec);
	sfmt_init_gen_rand(&sfmt, seed.tv_usec);
	
	nbranch = gtree->ntips;
	branches = calloc(gtree->ntips, sizeof(*branches));
	alloc_chk(branches, errout);
	for(i = 0; i < nbranch; i++) {
		branches[i] = &gtree->tips[i];
	}
	//TODO: refactor?
	for(node = gtree->last; node; node = node->prev) {
		coal_time = get_next_coal_time(nbranch, 0, 2, theta, &sfmt);
		if(coal_time == FLT_MAX) {
			//TODO: handle error
		}
		current += coal_time;
		node->time = current;
		child1 = &branches[sfmt_genrand_uint32(&sfmt) % nbranch];
		child2 = &branches[sfmt_genrand_uint32(&sfmt) % (nbranch-1)];
		if(child1 <= child2) {
			child2++;
		}
		add_child(node, *child1);
		add_child(node, *child2);
		if(child1 < child2) {
			*child1 = node;
			*child2 = branches[nbranch-1];
		} else {
			*child1 = branches[nbranch-1];
			*child2 = node;
		}
		nbranch--;
	}
	
	//TODO: handle error?
	free(branches);
	
	return;
	
errout:
	debug_err_out();
	
}

struct gene_tree *gtree_init(float theta, size_t ntips)
{
	
	struct gene_tree *gtree;
	
	debug_var_decl("creating new gene tree");
	
	param_chk(0 != ntips && theta > 0.0, errout);
	
	gtree = calloc(1, sizeof(*gtree));
	alloc_chk(gtree, errout);
	
	gtree_tips_init(gtree, ntips);
	gtree_nodes_init(gtree);
	gtree_simulate_tree(gtree, theta);

	return(gtree);
	
errout:
	debug_err_out();
	
}

void gtree_add_seqs_to_tips (struct gene_tree *gtree, struct ms_tab *mstab)
{
	
	int i;
	unsigned int mol_counts[PHY_NUM_MOL_T] = {0};
	float freqa, freqg, freqc, freqt; //for readability
	float freqar, freqgr, freqcy, freqty;
	float pur, pyr, ag, ct, m, n, fracchange;
	unsigned int nmol;
	
	if(!gtree || !mstab) {
		//TODO: handle error
	}
	
	if(gtree->ntips != mstab->len) {
		//TODO: handle error
	}
	
	gtree->mstab = mstab;
	
	for(i = 0; i < mstab->len; i++) {
		gtree->tips[i].mseq = &mstab->mseq[i];
	}
	nmol = get_mol_counts(mstab, mol_counts);
	if(!nmol) {
		//TODO: handle error
	}
	
	freqa = (float)mol_counts[DNA_A]/(float)nmol;
	freqt = (float)mol_counts[DNA_T]/(float)nmol;
	freqc = (float)mol_counts[DNA_C]/(float)nmol;
	freqg = (float)mol_counts[DNA_G]/(float)nmol;
	
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
	if(!pur || !pyr) {
		//TOO: handle error
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
	fracchange = gtree->yrate * (2.0 * freqa * freqgr + 2.0 * freqc * freqty) + gtree->xrate * (1.0 - freqa * freqa - freqc * freqc - freqg * freqg - freqt * freqt);
	gtree->xrate /= -(fracchange);
	gtree->yrate /= -(fracchange);
	
	/***********************************************************/
	
}

static void gnode_set_exp(struct gene_node *gnode, float xrate, float yrate)
{
	
	float length, n, n1, n2;
	struct gene_node *parent;
	
	if(!gnode) {
		//TODO: handle error
	}
	if(gnode->parent && !gnode->exp_valid) {
		parent = gnode->parent;
		length = parent->time - gnode->time;
		/********************************************************
		 * The following code adapted from LAMARC, (c) 2002 
		 * Peter Beerli, Mary Kuhner, Jon  Yamato and Joseph Felsenstein
		 * TODO: license ref?
		********************************************************/
		errno = 0;
		n = exp(length * xrate);
		n1 = (length * xrate) / 2.0;
		n2 = (length * yrate) / 2.0;
		if(errno) {
			err_debug("Value overflow setting n\n");
		}
		if(n == 1.0) {
			err_debug("n is 1.0, this will cause problems later...\n");
		}
		gnode->expA = 1.0 - n;
		gnode->lexpA = n1 + log(exp(-n1) - exp(n1));
		gnode->expB = n * exp(length * yrate);
		gnode->lexpB = (2.0 * n1) + (length * yrate);
		gnode->expC = n - gnode->expB;
		gnode->lexpC = (2.0 * n1) + n2 + log(exp(-n2) - exp(n2));
		/*******************************************************/
		gnode->exp_valid = 1;
	}
	if(gnode->child1) {
		gnode_set_exp(gnode->child1, xrate, yrate);
	}
	if(gnode->child2) {
		gnode_set_exp(gnode->child2, xrate, yrate);
	}
	
}

void gtree_set_exp(struct gene_tree *gtree)
{
	
	if(!gtree) {
		//TODO: handle error
	}
	
	gnode_set_exp(gtree->root, gtree->xrate, gtree->yrate);
	
}

static void gnode_get_llhood_lcomps(struct gene_node *gnode, float *cllike, 
							float *lfreq, float *lcomps)
{
	
	float sumAllfact, lsumAll, sumPur, lsumPur, sumPyr, lsumPyr;
	float normal, comp;
	int i;
	
	if(!gnode || !cllike || !lfreq || !lcomps) {
		//TODO: handle error
	}
	
	//log_debug("Setting components for node %i:\n", gnode->idx);
	
	/************************************************************
	 * The following code adapted from LAMARC, (c) 2002 
	 * Peter Beerli, Mary Kuhner, Jon  Yamato and Joseph Felsenstein
	 * TODO: license ref?
	 ***********************************************************/
	normal = -FLT_MAX;
	for(i = 0; i < NUM_BASE; i++) {
		normal = fmaxf(normal, (lfreq[i] + cllike[i]));
	}
	sumAllfact = 0;
	errno = 0;
	for(i = 0; i < NUM_BASE; i++) {
		if(cllike[i] > -FLT_MAX) {
			sumAllfact += exp((lfreq[i] + cllike[i]) - normal);
		}
	}
	if(sumAllfact <= 0.0 || errno) {
		err_debug("Intermediate sumAllfact is invalid.\n");
	}
	lsumAll = gnode->lexpA + log(sumAllfact) + normal;
	
	normal = fmaxf((lfreq[FREQ_AR] + cllike[DNA_A]), 
				   (lfreq[FREQ_GR] + cllike[DNA_G]));
	errno = 0;
	sumPur = 0;
	if(cllike[DNA_A] > -FLT_MAX) {
		sumPur += exp((lfreq[FREQ_AR] + cllike[DNA_A]) - normal);
	}
	if(cllike[DNA_G] > -FLT_MAX) {
		sumPur += exp((lfreq[FREQ_GR] + cllike[DNA_G]) - normal);
	}
	if(sumPur < 0 || errno) {
		err_debug("Intermediate sumPur is invalid.\n");
	}
	
	if(sumPur > 0) {
		lsumPur = log(sumPur) + normal;
	} else {
		lsumPur = -FLT_MAX;
	}
	
	normal = fmaxf((lfreq[FREQ_CY] + cllike[DNA_C]), 
				   (lfreq[FREQ_TY] + cllike[DNA_T]));
	errno = 0;
	sumPyr = 0;
	if(cllike[DNA_C] > -FLT_MAX) {
		sumPyr += exp((lfreq[FREQ_CY] + cllike[DNA_C]) - normal);
	}
	if(cllike[DNA_T] > -FLT_MAX) {
		sumPyr += exp((lfreq[FREQ_TY] + cllike[DNA_T]) - normal);
	}
	
	if(sumPyr < 0.0 || errno) {
		err_debug("Intermediate sumPyr is invalid.\n");
	}
	
	if(sumPyr > 0) {
		lsumPyr = log(sumPyr) + normal;
	} else {
		lsumPyr = -FLT_MAX;
	}
	
	//log_debug("sumAllfact = %f, lsumAll = %f, sumPur = %f, lsumPur = %f, sumPyr = %f, lsumPyr = %f\n", sumAllfact, lsumAll, sumPur, lsumPur, sumPyr, lsumPyr);
	
	for(i = 0; i < NUM_BASE; i++) {
		//TODO: place these components into an array rather than recalc?
	
		normal = fmaxf(lsumAll, (gnode->lexpB + cllike[i]));
		normal = fmaxf(normal, (gnode->lexpC + ((i==DNA_A||i==DNA_G)?lsumPur:lsumPyr)));

		comp = exp(lsumAll - normal);
		comp += exp((gnode->lexpB + cllike[i]) - normal);
		comp += exp((gnode->lexpC + ((i==DNA_A||i==DNA_G)?lsumPur:lsumPyr)) - 			normal);
		lcomps[i] = log(comp) + normal;
		//log_debug("comp[%i] = %f, lcomps[%i] = %f\n", i, comp, i, lcomps[i]);
	}
	/***********************************************************/
	
}

static void gnode_get_llhoods(struct gene_node *gnode, float *llhoods, 
								int pos, float *lfreq)
{
	
	float c1llike[NUM_BASE];
	float c2llike[NUM_BASE];
	float lcomp1[NUM_BASE];
	float lcomp2[NUM_BASE];
	struct gene_node *child1, *child2;
	struct mol_seq *mseq;
	int i;
	
	if(!gnode || !llhoods || !lfreq) {
		//TODO: handle error
	}
	
	if(gnode->time == 0) {
		mseq = gnode->mseq;
		if(!mseq || (pos < 0 || pos > mseq->len)) {
			//TODO: handle error
		}
		for(i = 0; i < NUM_BASE; i++) {
			llhoods[i] = -FLT_MAX;
		}
		llhoods[mseq->seq[pos]] = 0;
		return;
	}
	
	child1 = gnode->child1;
	child2 = gnode->child2;
	gnode_get_llhoods(child1, c1llike, pos, lfreq);
	gnode_get_llhood_lcomps(child1, c1llike, lfreq, lcomp1);
	gnode_get_llhoods(child2, c2llike, pos, lfreq);
	gnode_get_llhood_lcomps(child2, c2llike, lfreq, lcomp2);
	
	for(i = 0; i < NUM_BASE; i++) {
		llhoods[i] = lcomp1[i] + lcomp2[i];
	}
	
}


void gtree_set_llhood(struct gene_tree *gtree)
{
	
	float pos_llhood_vals[4];
	float lhood, llhood, normal;
	struct ms_tab *mstab;
	int i, j;
	
	if(!gtree || !gtree->root) {
		//TODO: handle error
	}
	
	mstab = gtree->mstab;
	llhood = 0;
	
	for(i = 0; i < mstab->seq_len; i++) {
		lhood = 0;
		normal = -FLT_MAX;
		gnode_get_llhoods(gtree->root, pos_llhood_vals, i, gtree->lfreq);
		//Jumping through this hoop to prevent under-run.
		for(j = 0; j < 4; j++) {
			normal = fmaxf(normal, (pos_llhood_vals[j] + gtree->lfreq[j]));
		}
		for(j = 0; j < 4; j++) {
			lhood += exp((pos_llhood_vals[j] + gtree->lfreq[j]) - normal);
		}
		printf("%f,", log(lhood) + normal);
		llhood += log(lhood) + normal;
	}
	printf("\n");
	gtree->llhood = llhood;
	
}

static void gnode_print_newick(struct gene_node *gnode)
{
	
	if(!gnode) {
		//TODO: handle error
	}
	
	if(gnode->child1) {
		printf("(");
        gnode_print_newick(gnode->child1);
	}
	if(NULL != gnode->child2) {
        printf(",");
        gnode_print_newick(gnode->child2);
        printf(")");
    }
	printf("%i:%f", gnode->idx, gnode->time);
    if(!gnode->parent) {
        printf(";\n");
    }
}

void gtree_print_newick(struct gene_tree *gtree)
{
	
	if(!gtree) {
		//TODO: handle error
	}
	
	gnode_print_newick(gtree->root);
	
}