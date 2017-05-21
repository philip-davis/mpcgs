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

#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

#include "debug.h"
#include "phylip.h"
#include "tree.h"

#include "SFMT.h"

static void gnode_list_create(size_t list_size, struct gnode_list *list)
{

    if (!list) {
        // TODO: handle error
    }

    list->gnodes = calloc(list_size, sizeof(*list->gnodes));
    if (!list->gnodes) {
        // TODO: handle error
    }
    list->head = 0;
    list->tail = 0;
}

static void gnode_list_destroy(struct gnode_list *list)
{

    free(list->gnodes);
}

static unsigned int gnode_list_get_size(struct gnode_list *list)
{

    if (!list) {
        // TODO: warn
        return (0);
    }

    return (list->head - list->tail);
}

static void gnode_list_enqueue(struct gnode_list *list, struct gene_node *node)
{

    if (!list) {
        // TODO: handle error
    }

    list->gnodes[list->head++] = node;
}

static void gnode_list_collate_head(struct gnode_list *list)
{

    struct gene_node *collate_target;
    struct gene_node **target_pos;
    float parent_time;

    if (!list) {
        // TODO: handle error
    }

    if (gnode_list_get_size(list) <= 1) {
        return;
    }

    collate_target = list->gnodes[list->head - 1];
    if (collate_target->parent) {
        parent_time = collate_target->parent->time;
    } else {
        parent_time = FLT_MAX;
    }

    target_pos = &list->gnodes[list->head - 1];
    while (target_pos > &list->gnodes[list->tail]) {
        *target_pos = *(target_pos - 1);
        if (parent_time > (*target_pos)->parent->time) {
            break;
        }
        target_pos--;
    }

    *target_pos = collate_target;
}

static struct gene_node *gnode_list_dequeue(struct gnode_list *list)
{

    struct gene_node *node;

    if (!list) {
        // TODO: handle error
    }

    if (list->head == list->tail) {
        return (NULL);
    }

    node = list->gnodes[list->tail++];

    return (node);
}

static struct gene_node *gnode_list_get_tail(struct gnode_list *list)
{

    if (!list) {
        // TODO: handle error
    }

    return (list->gnodes[list->tail]);
}

static int gnode_list_empty(struct gnode_list *list)
{

    if (!list) {
        // TODO: handle error
    }

    return (list->head == list->tail);
}

static struct gene_node *gnode_list_get_random(struct gnode_list *list,
                                               sfmt_t *sfmt)
{

    if (!list || !sfmt) {
        // TODO: handle error
    }

    if (gnode_list_get_size(list) == 0) {
        return (NULL);
    }

    return (list->gnodes[list->tail + (sfmt_genrand_uint32(sfmt) %
                                       gnode_list_get_size(list))]);
}

float get_next_coal_time(unsigned int active,
                         unsigned int inactive,
                         unsigned int act_coal,
                         float theta,
                         sfmt_t *sfmt)
{

	float r, time, denom;

    debug_var_decl("allocating tree nodes");

    param_chk(sfmt && act_coal <= 2 && act_coal > 0, errout);

    if (active < act_coal || inactive < 2 - act_coal) {
        return (FLT_MAX);
    }

    do {
    	r = sfmt_genrand_real3(sfmt); // Uniform on interval (0,1)
    }while(r >= 1.0 || r <= 0.0);

    if (2 == act_coal) {
        denom = ((float)active * (float)(active - 1)) / 2.0;
    } else if (1 == act_coal) {
        denom = (float)active * (float)inactive;
    } else {
        // show be caught by param_chk
        return (FLT_MAX);
    }

    errno = 0;
    time = -(log(r) * theta) / (2.0 * denom);
    if (errno || time == 0.0) {
        err_warn("Random number generator returned an invalid value.\n");
    }
    return (time);

errout:
    debug_err_out();
}

void gtree_nodes_init(struct gene_tree *gtree, size_t ntips)
{

    struct gene_node *nodes;
    struct gene_node *prev = NULL;
    int i;

    debug_var_decl("allocating tree nodes");

    param_chk(gtree, errout);

    gtree->ntips = ntips;
    gtree->nnodes = ntips - 1;

    nodes = calloc(gtree->nnodes + gtree->ntips, sizeof(*nodes));
    alloc_chk(nodes, errout);

    for (i = 0; i < gtree->nnodes; i++) {
        nodes[i].order = nodes[i].idx = i;
        nodes[i].prev = prev;
        nodes[i].tree = gtree;
        if (prev) {
            prev->next = &nodes[i];
        }
        gtree->last = prev = &nodes[i];
    }

    for (i = gtree->nnodes; i < gtree->nnodes + gtree->ntips; i++) {
        nodes[i].order = -1;
        nodes[i].idx = i;
        nodes[i].tree = gtree;
    }

    gtree->nodes = nodes;
    gtree->tips = &nodes[gtree->nnodes];
    gtree->root = &nodes[0];

    return;

errout:
    debug_err_out();
}

static void gnode_add_child(struct gene_node *parent, struct gene_node *child)
{

    if (!parent || !child) {
        // TODO: handle error
    }

    if (!parent->child1) {
        parent->child1 = child;
    } else if (!parent->child2) {
        parent->child2 = child;
    } else {
        // TODO: handle error
    }

    child->parent = parent;
    child->exp_valid = 0;
}

void gtree_simulate_tree(struct gene_tree *gtree,
                                float theta,
                                sfmt_t *sfmt)
{

    struct gene_node **branches;
    struct gene_node *node, **child1, **child2;
    int nbranch;
    float current = 0, coal_time;
    int i;

    debug_var_decl("simulating gene tree");

    param_chk(gtree && theta > 0.0, errout);

    nbranch = gtree->ntips;
    branches = calloc(gtree->ntips, sizeof(*branches));
    alloc_chk(branches, errout);
    for (i = 0; i < nbranch; i++) {
        branches[i] = &gtree->tips[i];
    }
    // TODO: refactor?
    for (node = gtree->last; node; node = node->prev) {
        coal_time = get_next_coal_time(nbranch, 0, 2, theta, sfmt);
        if (coal_time == FLT_MAX) {
            // TODO: handle error
        }
        current += coal_time;
        node->time = current;
        child1 = &branches[sfmt_genrand_uint32(sfmt) % nbranch];
        child2 = &branches[sfmt_genrand_uint32(sfmt) % (nbranch - 1)];
        if (child1 <= child2) {
            child2++;
        }
        gnode_add_child(node, *child1);
        gnode_add_child(node, *child2);
        if (child1 < child2) {
            *child1 = node;
            *child2 = branches[nbranch - 1];
        } else {
            *child1 = branches[nbranch - 1];
            *child2 = node;
        }
        nbranch--;
    }

    // TODO: handle error?
    free(branches);

    return;

errout:
    debug_err_out();
}

void gtree_init(float theta,
                size_t ntips,
                struct gene_tree *gtree,
                sfmt_t *sfmt)
{

    debug_var_decl("creating new gene tree");

    param_chk(0 != ntips && theta > 0.0 && gtree, errout);

    memset(gtree, sizeof(*gtree), 0);

    gtree_nodes_init(gtree, ntips);
    gtree_simulate_tree(gtree, theta, sfmt);

    return;

errout:
    debug_err_out();
}

void gtree_add_seqs_to_tips(struct gene_tree *gtree, struct ms_tab *mstab)
{

    int i;
    unsigned int mol_counts[PHY_NUM_MOL_T] = { 0 };
    float freqa, freqg, freqc, freqt; // for readability
    float freqar, freqgr, freqcy, freqty;
    float pur, pyr, ag, ct, m, n, fracchange;
    unsigned int nmol;

    if (!gtree || !mstab) {
        // TODO: handle error
    }

    if (gtree->ntips != mstab->len) {
        // TODO: handle error
    }

    gtree->mstab = mstab;

    for (i = 0; i < mstab->len; i++) {
        gtree->tips[i].mseq = &mstab->mseq[i];
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
}

static void gnode_set_exp(struct gene_node *gnode, float xrate, float yrate)
{

    float length, n1, n2;
    struct gene_node *parent;

    if (!gnode) {
        // TODO: handle error
    }
    if (gnode->parent && !gnode->exp_valid) {
        parent = gnode->parent;
        length = parent->time - gnode->time;
        /********************************************************
         * The following code adapted from LAMARC, (c) 2002
         * Peter Beerli, Mary Kuhner, Jon  Yamato and Joseph Felsenstein
         * TODO: license ref?
        ********************************************************/
        n1 = (length * xrate) / 2.0;
        n2 = (length * yrate) / 2.0;
        gnode->lexpA = n1 + log(exp(-n1) - exp(n1));
        gnode->lexpB = (2.0 * n1) + (length * yrate);
        gnode->lexpC = (2.0 * n1) + n2 + log(exp(-n2) - exp(n2));
        /*******************************************************/
        gnode->exp_valid = 1;
    }
    if (gnode->child1) {
        gnode_set_exp(gnode->child1, xrate, yrate);
    }
    if (gnode->child2) {
        gnode_set_exp(gnode->child2, xrate, yrate);
    }
}

void gtree_set_exp(struct gene_tree *gtree)
{

    if (!gtree) {
        // TODO: handle error
    }

    gnode_set_exp(gtree->root, gtree->xrate, gtree->yrate);
}

static void gnode_get_llhood_lcomps(struct gene_node *gnode,
                                    float *cllike,
                                    float *lfreq,
                                    float *lcomps)
{

    float sumAllfact, lsumAll, sumPur, lsumPur, sumPyr, lsumPyr;
    float normal, comp;
    int i;

    if (!gnode || !cllike || !lfreq || !lcomps) {
        // TODO: handle error
    }

    // log_debug("Setting components for node %i:\n", gnode->idx);

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
    errno = 0;
    for (i = 0; i < NUM_BASE; i++) {
        if (cllike[i] > -FLT_MAX) {
            sumAllfact += exp((lfreq[i] + cllike[i]) - normal);
        }
        if (errno == ERANGE) {
            printf("WARNING: an overflow error occurred.\n");
        }
    }
    if (sumAllfact <= 0.0 || errno) {
        err_debug("Intermediate sumAllfact is invalid.\n");
    }
    lsumAll = gnode->lexpA + log(sumAllfact) + normal;

    normal =
      fmaxf((lfreq[FREQ_AR] + cllike[DNA_A]), (lfreq[FREQ_GR] + cllike[DNA_G]));
    errno = 0;
    sumPur = 0;
    if (cllike[DNA_A] > -FLT_MAX) {
        sumPur += exp((lfreq[FREQ_AR] + cllike[DNA_A]) - normal);
        if (errno == ERANGE) {
            printf("WARNING: an overflow error occurred.\n");
        }
    }
    if (cllike[DNA_G] > -FLT_MAX) {
        sumPur += exp((lfreq[FREQ_GR] + cllike[DNA_G]) - normal);
        if (errno == ERANGE) {
            printf("WARNING: an overflow error occurred.\n");
        }
    }
    if (sumPur < 0 || errno) {
        err_debug("Intermediate sumPur is invalid.\n");
    }

    if (sumPur > 0) {
        lsumPur = log(sumPur) + normal;
    } else {
        lsumPur = -FLT_MAX;
    }

    normal =
      fmaxf((lfreq[FREQ_CY] + cllike[DNA_C]), (lfreq[FREQ_TY] + cllike[DNA_T]));
    errno = 0;
    sumPyr = 0;
    if (cllike[DNA_C] > -FLT_MAX) {
        sumPyr += exp((lfreq[FREQ_CY] + cllike[DNA_C]) - normal);
        if (errno == ERANGE) {
            printf("WARNING: an overflow error occurred.\n");
        }
    }
    if (cllike[DNA_T] > -FLT_MAX) {
        sumPyr += exp((lfreq[FREQ_TY] + cllike[DNA_T]) - normal);
        if (errno == ERANGE) {
            printf("WARNING: an overflow error occurred.\n");
        }
    }

    if (sumPyr < 0.0 || errno) {
        err_debug("Intermediate sumPyr is invalid.\n");
    }

    if (sumPyr > 0) {
        lsumPyr = log(sumPyr) + normal;
    } else {
        lsumPyr = -FLT_MAX;
    }

    // log_debug("sumAllfact = %f, lsumAll = %f, sumPur = %f, lsumPur = %f,
    // sumPyr = %f, lsumPyr = %f\n", sumAllfact, lsumAll, sumPur, lsumPur,
    // sumPyr, lsumPyr);

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
        // log_debug("comp[%i] = %f, lcomps[%i] = %f\n", i, comp, i, lcomps[i]);
    }
    /***********************************************************/
}

static void gnode_get_llhoods(struct gene_node *gnode,
                              float *llhoods,
                              int pos,
                              float *lfreq)
{

    float c1llike[NUM_BASE];
    float c2llike[NUM_BASE];
    float lcomp1[NUM_BASE];
    float lcomp2[NUM_BASE];
    struct gene_node *child1, *child2;
    struct mol_seq *mseq;
    int i;

    if (!gnode || !llhoods || !lfreq) {
        // TODO: handle error
    }

    if (gnode->time == 0) {
        mseq = gnode->mseq;
        if (!mseq || (pos < 0 || pos > mseq->len)) {
            // TODO: handle error
        }
        for (i = 0; i < NUM_BASE; i++) {
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

    for (i = 0; i < NUM_BASE; i++) {
        llhoods[i] = lcomp1[i] + lcomp2[i];
    }
}

static void gnode_connect(struct gene_node *child, struct gene_node *parent)
{

    if (parent) {
        if (!parent->child1) {
            parent->child1 = child;
        } else if (!parent->child2) {
            parent->child2 = child;
        } else {
            printf("Connecting to a full parent.\n");
        }
    }

    if (child) {
        if (child->parent) {
            printf("Child already has a parent.\n");
        } else {
            child->parent = parent;
        }
    }
}

static void gnode_disconnect(struct gene_node *gnode)
{

    struct gene_node *parent;

    if (!gnode) {
        // TODO: handle error
    }

    parent = gnode->parent;

    if (parent) {
        if (parent->child1 == gnode) {
            parent->child1 = parent->child2;
            parent->child2 = NULL;
        } else {
            parent->child2 = NULL;
        }
    }

    gnode->parent = NULL;
}

static void gnode_extract(struct gene_node *gnode)
{

    if (!gnode) {
        // TODO handle error
    }

    if (gnode->tree->root == gnode) {
        gnode->tree->root = gnode->next;
    }

    if (gnode->tree->last == gnode) {
        gnode->tree->last = gnode->prev;
    }

    if (gnode->prev) {
        gnode->prev->next = gnode->next;
    }
    if (gnode->next) {
        gnode->next->prev = gnode->prev;
    }

    gnode->prev = NULL;
    gnode->next = NULL;
}

static void gnode_insert_after(struct gene_node *gnode, struct gene_node *prev)
{

    if (!gnode) {
        // TODO: handle error
    }

    gnode->prev = prev;
    if (prev) {
        gnode->next = prev->next;
        prev->next = gnode;
    } else {
        gnode->next = gnode->tree->root;
        gnode->tree->root = gnode;
    }

    if (gnode->next) {
        gnode->next->prev = gnode;
    } else {
        gnode->tree->last = gnode;
    }
}

void gtree_set_llhood(struct gene_tree *gtree)
{

    float pos_llhood_vals[4];
    float lhood, llhood, normal;
    struct ms_tab *mstab;
    int i, j;

    if (!gtree || !gtree->root) {
        // TODO: handle error
    }

    mstab = gtree->mstab;
    llhood = 0;

    for (i = 0; i < mstab->seq_len; i++) {
        lhood = 0;
        normal = -FLT_MAX;
        gnode_get_llhoods(gtree->root, pos_llhood_vals, i, gtree->lfreq);
        // Jumping through this hoop to prevent under-run.
        for (j = 0; j < 4; j++) {
            normal = fmaxf(normal, (pos_llhood_vals[j] + gtree->lfreq[j]));
        }
        for (j = 0; j < 4; j++) {
            lhood += exp((pos_llhood_vals[j] + gtree->lfreq[j]) - normal);
        }
        llhood += log(lhood) + normal;
    }
    gtree->llhood = llhood;
}

static void gtree_copy(struct gene_tree *gtree, struct gene_tree *newtree)
{

    struct gene_node *newnodes;
    struct gene_node *gnode;
    size_t nodesSz;
    int i;

    if (!gtree || !newtree) {
        // TODO: handle error
    }

    nodesSz = sizeof(*gtree->nodes) * (gtree->nnodes + gtree->ntips);

    if(newtree->nodes) {
    	newnodes = newtree->nodes;
    } else {
    	newnodes = malloc(nodesSz);
    }
    memcpy(newnodes, gtree->nodes, nodesSz);

    *newtree = *gtree;

    for (i = 0; i < newtree->nnodes + newtree->ntips; i++) {
        // TODO: find a better way to do this
        gnode = &newnodes[i];
        if (gnode->parent) {
            gnode->parent = &newnodes[gnode->parent - gtree->nodes];
        }
        if (gnode->child1) {
            gnode->child1 = &newnodes[gnode->child1 - gtree->nodes];
        }
        if (gnode->child2) {
            gnode->child2 = &newnodes[gnode->child2 - gtree->nodes];
        }
        if (gnode->prev) {
            gnode->prev = &newnodes[gnode->prev - gtree->nodes];
        }
        if (gnode->next) {
            gnode->next = &newnodes[gnode->next - gtree->nodes];
        }
        gnode->tree = newtree;
        if (gnode->order == 0) {
            newtree->root = gnode;
        }
        if (gnode->order == (newtree->nnodes - 1)) {
            newtree->last = gnode;
        }
    }

    newtree->nodes = newnodes;
    newtree->tips = &newnodes[gtree->nnodes];
}

struct gene_tree *gtree_create_copy(struct gene_tree *gtree)
{

    struct gene_tree *newtree;
    struct gene_node *newnodes;
    struct gene_node *gnode;
    size_t nodesSz;
    int i;

    if (!gtree || !newtree) {
        // TODO: handle error
    }

    newtree = malloc(sizeof(*newtree));
    *newtree = *gtree;

    nodesSz = sizeof(*newtree->nodes) * (newtree->nnodes + newtree->ntips);
    newnodes = malloc(nodesSz);
    memcpy(newnodes, gtree->nodes, nodesSz);

    for (i = 0; i < newtree->nnodes + newtree->ntips; i++) {
        // TODO: find a better way to do this
        gnode = &newnodes[i];
        if (gnode->parent) {
            gnode->parent = &newnodes[gnode->parent - gtree->nodes];
        }
        if (gnode->child1) {
            gnode->child1 = &newnodes[gnode->child1 - gtree->nodes];
        }
        if (gnode->child2) {
            gnode->child2 = &newnodes[gnode->child2 - gtree->nodes];
        }
        if (gnode->prev) {
            gnode->prev = &newnodes[gnode->prev - gtree->nodes];
        }
        if (gnode->next) {
            gnode->next = &newnodes[gnode->next - gtree->nodes];
        }
        gnode->tree = newtree;
        if (gnode->order == 0) {
            newtree->root = gnode;
        }
        if (gnode->order == (newtree->nnodes - 1)) {
            newtree->last = gnode;
        }
    }

    newtree->nodes = newnodes;
    newtree->tips = &newnodes[gtree->nnodes];

    return(newtree);

}

static void gtree_fixup_order(struct gene_tree *gtree, struct gene_node *stopat)
{

    struct gene_node *node;

    if (!gtree) {
        // TODO: handle error
    }

    gtree->root->order = 0;
    for (node = gtree->root; node != gtree->last; node = node->next) {
        if (node == stopat) {
            break;
        }
        node->next->order = node->order + 1;
    }
}

struct gene_tree *gtree_propose(struct gene_tree *current,
                                float theta,
                                sfmt_t *sfmt)
{

    struct gene_tree *proposal;
    struct gene_node *target, *parent, *gparent, *newgparent, *node, *tail;
    struct gene_node *child1, *child2, *oldsibling, *sibling, *newnode;
    struct gene_node *ival_end;
    struct gnode_list ival_list;
    float currT, nextT, eventT;

    if (!current || !sfmt) {
        // TODO: handle error
    }

    gnode_list_create(current->nnodes + current->ntips, &ival_list);

    proposal = gtree_create_copy(current);

    target = &proposal->nodes[sfmt_genrand_uint32(sfmt) %
                              ((current->nnodes + current->ntips) - 1)];
    if (target >= proposal->root) {
        target++;
    }
    parent = target->parent;

    gnode_disconnect(target);
    oldsibling = parent->child1;

    if (target->time == 0) {
        // target is a tip
        node = proposal->last;
    } else {
        node = target->prev;
    }
    while (node) {
        child1 = node->child1;
        child2 = node->child2;
        if (child1 && child1->time <= target->time) {
            gnode_list_enqueue(&ival_list, child1);
        }
        if (child2 && child2->time <= target->time) {
            gnode_list_enqueue(&ival_list, child2);
        }
        node = node->prev;
    }

    /********************************************************
     * The following code adapted from LAMARC, (c) 2002
     * Peter Beerli, Mary Kuhner, Jon  Yamato and Joseph Felsenstein
     * TODO: license ref?
     ********************************************************/
    currT = target->time;

    while (1) {
        if (gnode_list_empty(&ival_list)) {
            // TODO: handle error
        }

        ival_end = (gnode_list_get_tail(&ival_list))->parent;
        if (ival_end) {
            nextT = ival_end->time;
            eventT = get_next_coal_time(
              1, gnode_list_get_size(&ival_list), 1, theta, sfmt);
        } else {
            nextT = FLT_MAX;
            eventT = get_next_coal_time(2, 0, 2, theta, sfmt);
        }
        if ((currT + eventT) < nextT) {
            sibling = gnode_list_get_random(&ival_list, sfmt);
            if (sibling == parent) {
                // Parent is a stick at this point, so it can't be a sibling.
                sibling = parent->child1;
            }

            newnode = parent; // for clarity
            gparent = parent->parent;
            newgparent = sibling->parent;

            if (parent != ival_end) {
                gnode_extract(parent);
                gnode_insert_after(newnode, ival_end);
            }
            if (parent != sibling->parent) {
                gnode_disconnect(oldsibling);
                gnode_disconnect(parent);
                gnode_disconnect(sibling);
                gnode_connect(oldsibling, gparent);
                gnode_connect(newnode, newgparent);
                gnode_connect(sibling, newnode);
            }

            gnode_connect(target, newnode);
            newnode->time = currT + eventT;
            newnode->exp_valid = 0;
            target->exp_valid = 0;
            sibling->exp_valid = 0;

            gtree_fixup_order(proposal, target);
            gnode_set_exp(parent, proposal->xrate, proposal->yrate);

            break;

        } else {
            node = gnode_list_dequeue(&ival_list);
            if (!gnode_list_empty(&ival_list)) {
                tail = gnode_list_get_tail(&ival_list);
                if (tail->parent == node->parent) {
                    // TODO: some explanation here
                    gnode_list_dequeue(&ival_list);
                }
            }
            // parent is guaranteed to exist since we would have merged the root
            // otherwise.
            gnode_list_enqueue(&ival_list, node->parent);
            gnode_list_collate_head(&ival_list);
        }

        currT = nextT;
    }

    /***********************************************************/

    gnode_list_destroy(&ival_list);

    return (proposal);
}

struct gene_tree *gtree_propose_fixed_target(struct gene_tree *current,
                                             struct gene_tree *proposal,
                                             float theta,
                                             unsigned int tgtidx,
                                             sfmt_t *sfmt)
{

    struct gene_node *target, *parent, *gparent, *newgparent, *node, *tail;
    struct gene_node *child1, *child2, *oldsibling, *sibling, *newnode;
    struct gene_node *ival_end;
    struct gnode_list ival_list;
    float currT, nextT, eventT;

    if (!current || !sfmt) {
        // TODO: handle error
    }

    if (tgtidx >= (current->nnodes + current->ntips) ||
        tgtidx == proposal->root->idx) {
        // TODO: handle error
    }

    gnode_list_create(current->nnodes + current->ntips, &ival_list);

    gtree_copy(current, proposal);

    target = &proposal->nodes[tgtidx];
    if (target >= proposal->root) {
        target++;
    }
    parent = target->parent;

    gnode_disconnect(target);
    oldsibling = parent->child1;

    if (target->time == 0) {
        // target is a tip
        node = proposal->last;
    } else {
        node = target->prev;
    }
    while (node) {
        child1 = node->child1;
        child2 = node->child2;
        if (child1 && child1->time <= target->time) {
            gnode_list_enqueue(&ival_list, child1);
        }
        if (child2 && child2->time <= target->time) {
            gnode_list_enqueue(&ival_list, child2);
        }
        node = node->prev;
    }

    /********************************************************
     * The following code adapted from LAMARC, (c) 2002
     * Peter Beerli, Mary Kuhner, Jon  Yamato and Joseph Felsenstein
     * TODO: license ref?
     ********************************************************/
    currT = target->time;

    while (1) {
        if (gnode_list_empty(&ival_list)) {
            // TODO: handle error
        }

        ival_end = (gnode_list_get_tail(&ival_list))->parent;
        if (ival_end) {
            nextT = ival_end->time;
            eventT = get_next_coal_time(
              1, gnode_list_get_size(&ival_list), 1, theta, sfmt);
        } else {
            nextT = FLT_MAX;
            eventT = get_next_coal_time(2, 0, 2, theta, sfmt);
        }
        if ((currT + eventT) < nextT) {
            sibling = gnode_list_get_random(&ival_list, sfmt);
            if (sibling == parent) {
                // Parent is a stick at this point, so it can't be a sibling.
                sibling = parent->child1;
            }

            newnode = parent; // for clarity
            gparent = parent->parent;
            newgparent = sibling->parent;

            if (parent != ival_end) {
                gnode_extract(parent);
                gnode_insert_after(newnode, ival_end);
            }
            if (parent != sibling->parent) {
                gnode_disconnect(oldsibling);
                gnode_disconnect(parent);
                gnode_disconnect(sibling);
                gnode_connect(oldsibling, gparent);
                gnode_connect(newnode, newgparent);
                gnode_connect(sibling, newnode);
            }

            gnode_connect(target, newnode);
            newnode->time = currT + eventT;
            newnode->exp_valid = 0;
            target->exp_valid = 0;
            sibling->exp_valid = 0;

            gtree_fixup_order(proposal, target);
            gnode_set_exp(parent, proposal->xrate, proposal->yrate);

            break;

        } else {
            node = gnode_list_dequeue(&ival_list);
            if (!gnode_list_empty(&ival_list)) {
                tail = gnode_list_get_tail(&ival_list);
                if (tail->parent == node->parent) {
                    // TODO: some explanation here
                    gnode_list_dequeue(&ival_list);
                }
            }
            // parent is guaranteed to exist since we would have merged the root
            // otherwise.
            gnode_list_enqueue(&ival_list, node->parent);
            gnode_list_collate_head(&ival_list);
        }

        currT = nextT;
    }

    /***********************************************************/

    gnode_list_destroy(&ival_list);
}

void gtree_digest(struct gene_tree *gtree, struct gtree_summary *digest)
{

    struct gene_node *node;
    float *ival;

    if (!gtree) {
        // TODO: handle error
    }

    if (digest->nintervals != gtree->nnodes) {
        // TODO: handle error
    }

    ival = digest->intervals;
    for (node = gtree->root; node; node = node->next) {
        *ival = node->next ? (node->time - node->next->time) : node->time;
        ival++;
    }
}

void gtree_summary_set_create(struct gtree_summary_set **sum_set,
                              size_t count,
                              size_t nintervals)
{

    int i;
    struct gtree_summary *summary;

    if (!sum_set) {
        // TODO: handle error
    }

    *sum_set = malloc(sizeof(**sum_set));

    (*sum_set)->nsummaries = 0;
    (*sum_set)->szintervals = count;
    (*sum_set)->summaries = malloc(count * sizeof(*((*sum_set)->summaries)));
    if (!(*sum_set)->summaries) {
        // TODO: handle error
    }

    summary = (*sum_set)->summaries;
    for (i = 0; i < count; i++) {
        summary->nintervals = nintervals;
        summary->intervals = malloc(nintervals * sizeof(*(summary->intervals)));
        if (!summary->intervals) {
            // TODO: handle error
        }
        summary++;
    }
}

static float gtree_summary_calc_lposterior(struct gtree_summary *sum, float theta)
{

    int i;
    int lineages;
    float exp1, coeff;

    if (!sum) {
        // TODO: handle error
    }

    if (theta <= 0.0) {
        // TODO: handle error
    }

    exp1 = 0;
    for (i = 0; i < sum->nintervals; i++) {
        lineages = i + 2;
        exp1 += -(float)(lineages * (lineages - 1)) * sum->intervals[i];
    }
    coeff = (float)(lineages - 1) * logf(2.0 / theta);

    return (coeff + (exp1 / theta));
}

void gtree_summary_set_base_lposteriors(struct gtree_summary_set *sum_set,
                                        float drv_theta)
{

    struct gtree_summary *summary;
    int i;

    if (!sum_set) {
        // TODO: handle error
    }

    if (drv_theta <= 0.0) {
        // TODO: handle error
    }

    summary = sum_set->summaries;
    for (i = 0; i < sum_set->nsummaries; i++) {
        summary->ldrv_posterior = gtree_summary_calc_lposterior(summary, drv_theta);
        summary++;
    }
}

float gtree_summary_set_llkhood(struct gtree_summary_set *summary_set,
                                float theta)
{

    struct gtree_summary *summary;
    float lkhood = 0;
    float normal;
    int i;

    if (!summary_set) {
        // TODO: handle error
    }

    if (theta <= 0.0) {
        // TODO: handle error
    }

    normal = -FLT_MAX;

    summary = summary_set->summaries;

    // calculate posterior for each summary tree
    for (i = 0; i < summary_set->nsummaries; i++) {
        summary->ltmp_lkhood_comp =
          gtree_summary_calc_lposterior(summary, theta) - summary->ldrv_posterior;
        if (summary->ltmp_lkhood_comp > normal) {
            normal = summary->ltmp_lkhood_comp;
        }
        summary++;
    }

    summary = summary_set->summaries;
    for (i = 0; i < summary_set->nsummaries; i++) {
        lkhood += exp(summary->ltmp_lkhood_comp - normal);
        summary++;
    }

    return (log(lkhood) + normal);
}

void gtree_summary_set_print_lkhoods(struct gtree_summary_set *summary_set,
                                     float start,
                                     float stop,
                                     float incr)
{

    float theta, llkhood;

    if (!summary_set) {
        // TODO: handle error
    }

    if (start <= 0.0) {
        // TODO: handle error
    }

    for (theta = start; theta <= stop; theta += incr) {
        llkhood = gtree_summary_set_llkhood(summary_set, theta);
        printf("%f,%f\n", theta, llkhood);
    }
}

static void gnode_print_newick(struct gene_node *gnode)
{

    if (!gnode) {
        // TODO: handle error
    }

    if (gnode->child1) {
        printf("(");
        gnode_print_newick(gnode->child1);
    }
    if (NULL != gnode->child2) {
        printf(",");
        gnode_print_newick(gnode->child2);
        printf(")");
    }
    printf("%i:%f", gnode->idx, gnode->time);
    if (!gnode->parent) {
        printf(";\n");
    }
}

void gtree_print_newick(struct gene_tree *gtree)
{

    if (!gtree) {
        // TODO: handle error
    }

    gnode_print_newick(gtree->root);
}

size_t weighted_pick(float *dist, size_t size_dist, float sum_dist, float randf)
{

    float threshold;
    int i;

    threshold = (1.0 - randf) * sum_dist;
    for (i = 0; i < size_dist; i++) {
        if (dist[i] > threshold) {
            return (i);
        }
        threshold -= dist[i];
    }

    // Rounding error (probably). Pick the first non-zero entry
    for (i = 0; i < size_dist; i++) {
        if (dist[i] > 0.0) {
            return (i);
        }
    }

    // all zero. Return last
    return (size_dist - 1);
}
