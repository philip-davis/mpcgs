/*
    phylip.h allows parsing of a .phy file. See:
    http://evolution.genetics.washington.edu/phylip/doc/sequence.html

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

#ifndef MPCGS_PHYLIP_H
#define MPCGS_PHYLIP_H

#include<stdlib.h>

//TODO add support for interleaved .phy files.

#define PHY_NAME_LEN 10
#define PHY_N_HDR_ATTR 2
#define PHY_N_SEQ_FIELD 2

enum mol_t { 
    DNA_A, 
    DNA_T, 
    DNA_C, 
    DNA_G,
    DNA_U,
    DNA_PYR,
    DNA_PUR,
    DNA_WEAK,
    DNA_STRONG,
    DNA_KETO,
    DNA_AMINO,
    DNA_NOT_A,
    DNA_NOT_C,
    DNA_NOT_G,
    DNA_NOT_T,
    DNA_UNKNOWN,
    DNA_DEL,
    PRO_ALA,
    PRO_ASX,
    PRO_CYS,
    PRO_ASP,
    PRO_GLU,
    PRO_PHE,
    PRO_GLY,
    PRO_HIS,
    PRO_ILEU,
    PRO_LYS,
    PRO_LEU,
    PRO_MET,
    PRO_ASN,
    PRO_PRO,
    PRO_GLN,
    PRO_ARG,
    PRO_SER,
    PRO_THR,
    PRO_VAL,
    PRO_TRP,
    PRO_UNKNOWN,
    PRO_TYR,
    PRO_GLX,
    PRO_NONSENSE,
    PRO_UNKNOWN_OR_DEL,
    PRO_DEL,
    MOL_INVALID = 255
};

enum mol_cat {
    MOL_DNA,
    MOL_PRO
};

struct mol_seq {
    size_t len;
    char *seq;
    char name[PHY_NAME_LEN + 1];
};

struct ms_tab {
    size_t len;
    struct mol_seq *mseq;
};

struct ms_tab *init_ms_tab(const char *filename);
void free_mol_seq(struct mol_seq *mseq);
void free_ms_tab(struct ms_tab *mstab);

#endif /* MPCGS_PHYLIP_H */
