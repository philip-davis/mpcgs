/*
    Copyright (C) {2017}  {Philip Davis}

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

struct ms_tab *init_ms_tab(const char *filename)
{

    FILE *gendatf;
    struct ms_tab *mstab;
    size_t nseq, seq_len, nread;
    struct mol_seq *mseq;
    char scanstr[100];
    char *tmp_buf;
    int err;
 
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

    if(PHY_N_HDR_ATTR != fscanf(gendatf, "%zi %zi ", &nseq, &seq_len)) {
        err = -MPCGS_EBADFORMAT;
        err_str = "bad file format";
        goto errfree;
    }
    log_debug("nseq = %zi, seqlen = %zi\n", nseq, seq_len);
    
    mstab->len = nseq;
    sprintf(scanstr, "%%10c%%%zic%%n ", seq_len);
    mstab->mseq = malloc(sizeof(*(mstab->mseq)) * mstab->len);
    alloc_chk(mstab->mseq, errfree, err_str);
    for(mseq = mstab->mseq; mseq < mstab->mseq + mstab->len; mseq++) {
        mseq->len = seq_len;
        tmp_buf = malloc(sizeof(*(mseq->seq) * mseq->len));
        alloc_chk(tmp_buf, errfree, err_str);
        if(fscanf(gendatf, scanstr, mseq->name, tmp_buf, &nread) != PHY_N_SEQ_FIELD || nread != (mseq->len + PHY_NAME_LEN)) {
            err = -MPCGS_EBADFORMAT;
            err_str = "bad file format";
            goto errfree;
        } 
        mseq->name[PHY_NAME_LEN] = '\0';
    }

    return(mstab);

errfree:
    free_ms_tab(mstab);
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
