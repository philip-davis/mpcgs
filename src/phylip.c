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
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "phylip.h"
#include "debug.h"

#define MPCGS_EBADFORMAT 201

//TODO expand supported molecules
static enum mol_t char_to_mol(char c, enum mol_cat mcat)
{

    if(mcat == MOL_DNA) {
        switch(c) {
            case 'A': return DNA_A;
            case 'C': return DNA_C;
            case 'G': return DNA_G;
            case 'T': return DNA_T;
            default:  return MOL_INVALID;
        }
    } else {
        return MOL_INVALID;
    }
}

//TODO refactor

struct ms_tab *init_ms_tab(const char *filename)
{

    FILE *gendatf;
    struct ms_tab *mstab;
    size_t nseq, seq_len, nread;
    struct mol_seq *mseq;
    char scanstr[100];
    char *tmp_buf = NULL;
    int i, err, count;
 
    const char *err_name = "reading sequence data";
    char *err_str;    

    if(!filename) {
        err = EINVAL;
        err_str = "no sequence data file given";
        goto errout;
    } 

    gendatf = fopen(filename, "r");
    if(!gendatf) {
        err = errno;
        err_str = strerror(err);
        goto errout;
    }

    log_debug("Reading sequence data from %s.\n", filename);

    mstab = malloc(sizeof mstab);
    alloc_chk(mstab, errclose, err_str);

    count = fscanf(gendatf, "%zi %zi ", &nseq, &seq_len);
    if(count != PHY_N_HDR_ATTR) {
        err = -MPCGS_EBADFORMAT;
        err_str = "bad file format";
        goto errfree;
    }
    log_debug("nseq = %zi, seqlen = %zi\n", nseq, seq_len);
    
    mstab->len = nseq;
    sprintf(scanstr, "%%10c%%%zic%%n ", seq_len);
    mstab->mseq = malloc(sizeof(*(mstab->mseq)) * mstab->len);
    alloc_chk(mstab->mseq, errfree, err_str);

    tmp_buf = malloc(sizeof(*tmp_buf) * seq_len);
    alloc_chk(tmp_buf, errfree, err_str);

    for(mseq = mstab->mseq; mseq < (mstab->mseq + mstab->len); mseq++) {
        mseq->len = seq_len;
        count = fscanf(gendatf, scanstr, mseq->name, tmp_buf, &nread);
        if(count != PHY_N_SEQ_FIELD || nread != (mseq->len + PHY_NAME_LEN)) {
            err = -MPCGS_EBADFORMAT;
            err_str = "bad file format";
            goto errfree;
        } 
        mseq->name[PHY_NAME_LEN] = '\0';

        mseq->seq = malloc(sizeof(*(mseq->seq)) * mseq->len);
        alloc_chk(mseq->seq, errfree, err_str);
        for(i = 0; i < mseq->len; i++) {
            mseq->seq[i] = char_to_mol(tmp_buf[i], MOL_DNA);
            if(mseq->seq[i] == MOL_INVALID) {
                err = -MPCGS_EBADFORMAT;
                err_str = "unrecognized molecule";
                goto errfree;
            }
        }
    }
    free(tmp_buf);
    
    fclose(gendatf);
    return(mstab);

errfree:
    free_ms_tab(mstab);
    if(tmp_buf) {
        free(tmp_buf);
    }
errclose:
    fclose(gendatf);
errout:
    err_out(err_name, err_str, -err);

}

void free_ms_tab(struct ms_tab *mstab)
{

    int i;

    if(!mstab) {
        err_warn("Tried to free a null ms_tab.\n");
        return;
    }

    if(mstab->mseq) {
        for(i = 0; i < mstab->len; i++) {
            if(mstab->mseq[i].seq) {
                free(mstab->mseq[i].seq);
            } else {
                err_warn("Freeing an ms_tab that is missing some sequence data.\n");
            }
        }
        free(mstab->mseq);
    } else {
        err_warn("Freeing an ms_tab that is missing all sequence data.\n");
    }
    free(mstab);

}

unsigned int get_mol_counts(struct ms_tab *mstab, unsigned int *counts)
{
	
	int i, j;
	struct mol_seq *mseq;
	unsigned int total = 0;
	
	if(!mstab) {
		err_warn("Trying to find molecule counts for a null ms_tab.\n");
		return(0);
	}
	
	if(!counts) {
		err_warn("Trying to write counts into a null buffer.\n");
		return(0);
	}
	
	memset(counts, 0, PHY_NUM_MOL_T);
	
	for(i = 0; i < mstab->len; i++) {
		mseq = &mstab->mseq[i];
		for(j = 0; j < mseq->len; j++) {
				counts[mseq->seq[j]]++;
				total++;
		}
	}
	
	return(total);
	
}
