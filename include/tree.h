/*
	tree.h defines the structures and functions for a phylogenetic tree.
	
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

#ifndef MPCGS_TREE_H
#define MPCGS_TREE_H

#include<stdlib.h>
#include<string.h>

#include "phylip.h"

enum {
		FREQ_A = 0,
		FREQ_T,
		FREQ_C,
		FREQ_G,
		FREQ_AR,
		FREQ_GR,
		FREQ_CY,
		FREQ_TY
};

struct gene_node {
	int order;
	int idx;
	struct gene_node *parent;
	struct gene_node *child1;
	struct gene_node *child2;
	struct gene_node *prev;
	struct gene_node *next;
	struct mol_seq *seq;
	float time;
	float expA, expB, expC;
};

struct gene_tree {
	struct gene_node *nodes;
	struct gene_node *tips;
	size_t nnodes;
	size_t ntips;
	struct gene_node *root;
	struct gene_node *lastnode;
	float freq[8];
};

struct gene_tree *create_new_tree(float theta, size_t ntips);

#endif /* MPCGS_TREE_H */