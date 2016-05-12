/*
 *  mfl_searchrec.c
 *
 *  THE MORPHY FUNCTION LIBRARY
 *  A library for phylogenetic analysis with emphasis on parsimony and
 *  morphology (but someday other methods)
 *
 *  Copyright (C) 2016  by Martin D. Brazeau, Thomas Guillerme,
 *  and Chris Desjardins
 *
 *  Created by Martin Brazeau and Chris Desjardins on 1/26/12.
 *  Some data structs, routines and ideas derived from:
 *      - The PHYLIP package by Joe Felsenstein
 *          <http://evolution.genetics.washington.edu/phylip.html>
 *      - MrBayes by John Huelsenbeck and Fredrik Ronquist
 *          <http://mrbayes.sourceforge.net/>
 *
 *  Any bugs, errors, inefficiences and general amateurish handling are our own
 *  and most likely the responsibility of MDB. We make no guarantees.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "morphy.h"

/* temporary place for prototypes */
void mfl_copy_row_from_partition_into_nodedata(mfl_charstate* target, mfl_datapartition_t* datapart, int row);
void mfl_copy_from_all_partitions_into_node_data(mfl_node_t* n, mfl_partition_set_t* partset);
void mfl_apply_characters_to_tips(mfl_tree_t* t, mfl_handle_s* handle, mfl_partition_set_t* parts);
/**/

void mfl_initialise_searchrec(mfl_searchrec_t* searchrec)
{
    searchrec->sr_best_length = MORPHY_UINTMAX;
}

void mfl_copy_row_from_partition_into_nodedata(mfl_charstate* target, mfl_datapartition_t* datapart, int row)
{
    int i = 0;
    int numchars = datapart->part_n_characters;
    
    for (i = 0; i < numchars; ++i) {
        target[i] = datapart->part_matrix[i + row * numchars];
    }
}

void mfl_copy_from_all_partitions_into_tip_nodedata(mfl_node_t* n, mfl_partition_set_t* partset)
{
    
    
    int i;
    
    assert(n->nodet_tip);
    int rownumber = n->nodet_index;
    
    int numparts = partset->ptset_n_parts;
    assert(n->nodet_num_dat_partitions == numparts);
    
    mfl_datapartition_t* datapart = NULL;
    mfl_nodedata_t* ndata = NULL;
    
    for (i = 0; i < numparts; ++i) {
        
    }
}

void mfl_apply_characters_to_tips(mfl_tree_t* t, mfl_handle_s* handle, mfl_partition_set_t* parts)
{
    int i = 0;
    int numtaxa = handle->n_taxa;
    mfl_nodearray_t nds = t->treet_treenodes;
    
    for (i = 0; i < numtaxa; ++i) {
        // For each partition, copy into the downpass set of the tip
        nds[i]->nodet_dataparts;
    }
}