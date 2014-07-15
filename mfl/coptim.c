/*
 *  coptim.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 2/15/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "morphy.h"

#define MFY_SUBTREE_REINSERTION_LOOP \
for (i = 0; i < nchar; ++i, ++srctemps, ++tgt1apos, ++tgt2apos) { \
    if ( !(*srctemps & (*tgt1apos | *tgt2apos)) ) { \
    
#define MFY_REINSERTION_LOOP_END \
    } \
} \

int mfl_compare (const void * a, const void * b)
{
    return ( *(charstate*)a - *(charstate*)b );
}

int mfl_n_unique_vals_in_array(int *array, int length)
{
    int i, n = 0;
    
    qsort(array, length, sizeof(int), mfl_compare);
    
    for (i = 1; i < length; ++i) {
        if (array[i] > array[i-1]) {
            ++n;
        }
    }
    
    return n;
}

void mfl_insert_matrix_column(int src_width, int tgt_width, charstate *src_column, 
                              charstate *tgt_column, int ntax)
{
    int i;
    
    for (i = 0; i < ntax; ++i) {
        *(tgt_column + i * tgt_width) = *(src_column + i * src_width);
    }
    
}

void mfl_split_matrix(chardata *cdata, charstate *matrix, int ntax, int nchar)
{
    int i=0, j=0;
    int pos_no_na = 0, pos_w_na=0;
    charstate *column_top;
    bool has_na = false;
    
    for (i = 0; i < nchar; ++i) {
        has_na = false;
        
        for (j = 0; j < ntax; ++j) {
            if (matrix[i + nchar *j] == 1) {
                has_na = true;
                break;
            }
        }
        
        if (has_na) {
            column_top = cdata->cd_wgaps + pos_w_na;
            mfl_insert_matrix_column(nchar, cdata->cd_n_wgaps, matrix + i, column_top, ntax);
            ++pos_w_na;
        }
        else {
            column_top = cdata->cd_nogaps + pos_no_na;
            mfl_insert_matrix_column(nchar, cdata->cd_n_nogaps, matrix + i, column_top, ntax);
            ++pos_no_na;
        }
    }
}

chardata *mfl_chardata_for_search(charstate *matrix, int ntax, int nchar)
{
    int i=0, j=0;
    int no_na =0, w_na=0; // Counters for number of characters with applicable/inapplicable values
    
    chardata *cdata = mfl_new_chardata();
    
    // Count inapplicables;
    for (i = 0; i < nchar; ++i) {
        for (j = 0; j < ntax; ++j) {
            if (matrix[i + nchar * j] == 1) {
                ++w_na;
                break;
            }
        }
    }
    
    no_na = nchar - w_na;
    
    cdata->cd_n_nogaps = no_na;
    cdata->cd_n_wgaps = w_na;
    cdata->cd_nogaps = (charstate*)malloc(ntax * no_na * sizeof(charstate));
    if (cdata->cd_nogaps == NULL) {
        dbg_printf("Malloc failure for cd_nogaps in mfl_chardata_for_search\n");
    }
    cdata->cd_wgaps = (charstate*)malloc(ntax * w_na * sizeof(charstate));
    if (cdata->cd_wgaps == NULL) {
        dbg_printf("Malloc failure for cd_wgaps in mfl_chardata_for_search\n");
    }
    
    mfl_split_matrix(cdata, matrix, ntax, nchar);
    
    return cdata;
}

int *mfl_get_character_minchanges(charstate *matrix, int ntax, int nchar/*, int *nwithgaps*/)
{
    /* Loop over all columns in the matrix, count the number of states in each,
     * count the number with gap entries.*/

    int i, j, k;
    charstate temp;
    
    int *minchanges_p = NULL;
    minchanges_p = (int*)malloc(nchar * sizeof(int));
    if (minchanges_p == NULL) {
        dbg_printf("Malloc failure: failed to allocate minchanges_p in coptim.c\n");
    }
    
    int *matrix_colum = NULL;
    /* plus 1 for null terminator */
    matrix_colum = (int*)malloc((ntax + 1) * sizeof(int));
    if (matrix_colum == NULL) {
        dbg_printf("Malloc failure: failed to allocate matrix_column in coptim.c\n");

    }
    
    //memset(nwithgaps, 0, nchar * sizeof(int));
    
    for (i = 0; i < nchar; ++i) {
        
        k = 0;
        memset(matrix_colum, 0, ntax * sizeof(int));
        
        for (j = 0; j < ntax; ++j) {
            if (matrix[i + j * (nchar)] != 1) {
                /*if (nwithgaps) {
                    nwithgaps[i] = 1;
                }
            }
            else {*/
                temp = matrix[i + j * (nchar)];
                if (temp != IS_APPLIC && (temp == (temp & (-(temp))))) {
                    matrix_colum[k] = matrix[i + j * (nchar)];
                    ++k;
                }
            }
        }
        
        matrix_colum[k] = '\0';
        minchanges_p[i] = mfl_n_unique_vals_in_array(matrix_colum, k-1);
    }
    
    free(matrix_colum);
    
    dbg_printf("\nMinimum changes per character:\n");
    for (i = 0; i < nchar; ++i) {
        dbg_printf("%i ", minchanges_p[i]);
    }
    dbg_printf("\n");
    
    return minchanges_p;
}

charstate * mfl_convert_tipdata(char *txtsrc, int ntax, int nchar, bool na_as_missing)
{
    
    /* One of Morphy's most important and unique features will be distinguishing
     * between a missing data entry "?" and a character-inapplicability entry "-".
     * Most existing programs just treat "?" and "-" as identical values (-1).
     * In Morphy, they can be treated differently, such that ancestral states
     * are not reconstructed for taxa which cannot logically have them */
    
    int i, j;
    
    charstate *tipdata = (charstate*) malloc(ntax * nchar * sizeof(charstate));
    
    if (!na_as_missing) {
        dbg_printf("Gap symbol ('-') treated as character inapplicability\n");
    }
    else {
        dbg_printf("Gap symbol ('-') treated as missing data\n");
    }

    
    for (i = 0, j = 0; txtsrc[i]; ++i, ++j) {
        if ((txtsrc[i] - '0') >= 0 && (txtsrc[i] - '0') <= 9) {
            tipdata[j] = 1 << (txtsrc[i] - '0' + 1);
            
        }
        else if (txtsrc[i] == '{' || txtsrc[i] == '(') {
            ++i;
            tipdata[j] = 0;
            while (txtsrc[i] != '}' && txtsrc[i] != ')') {
                if ((txtsrc[i] - '0') >= 0 && (txtsrc[i] - '0') <= 9) {
                    tipdata[j] = tipdata[j] | (1 << (txtsrc[i] - '0' + 1));
                    ++i;
                }
            }
        }
        else if (txtsrc[i] == '?') {
            tipdata[j] = IS_APPLIC;
        }
        else if (txtsrc[i] == '-') {
            if (na_as_missing) {
                tipdata[j] = IS_APPLIC;
            }
            else {
                tipdata[j] = 1;
            }

        }
        else if (txtsrc[i] == '\n') {
            --j;
        }
        else {
            --j;
        }
    }
    
    return tipdata;
}

void mfl_apply_tipdata(tree *currenttree, charstate *tipdata, int ntax, int nchar)
{   
    int i, j;
    
    for (i = 0; i < ntax; ++i) {
        //currenttree->trnodes[i]->apomorphies = &tipdata[i * nchar];
        currenttree->trnodes[i]->tempapos = &tipdata[i * nchar];
        
        if (!currenttree->trnodes[i]->apomorphies) {
            currenttree->trnodes[i]->apomorphies = (charstate*)malloc(nchar * sizeof(charstate));
        }
        for (j = 0; j < nchar; ++j) {
            currenttree->trnodes[i]->apomorphies[j] = currenttree->trnodes[i]->tempapos[j];
        }
    }
}

int mfl_locreopt_cost(node *src, node *tgt1, node *tgt2, int nchar, int diff)
{
    /* Returns cost of inserting subtree src between tgt1 and tgt2 following
     * the algorithms described by Ronquist (1998. Cladistics) and Goloboff 
     * (1993, 1996. Cladistics).*/
    
    int i;
    int cost = 0;
    charstate *srctemps = src->apomorphies;
    charstate *tgt1apos = tgt1->apomorphies;
    charstate *tgt2apos = tgt2->apomorphies;
    
    
    for (i = 0; i < nchar; ++i, ++srctemps, ++tgt1apos, ++tgt2apos) {
        if ( !(*srctemps & (*tgt1apos | *tgt2apos)) ) {
            ++cost;
            if (cost > diff) {
                return cost;
            }
        }
    }
    
    return cost;
}

int mfl_subtr_reinsertion(node *src, node *tgt1, node *tgt2, int nchar)
{
    /* Returns cost of reinserting the subtree src at the original place from
     * which it was clipped (tgt1 and tgt2). This score is used to compute the 
     * difference in length between the two subtrees so that neither the total 
     * tree length nor the length of the individual subtrees needs to be 
     * calculated (Ronquist 1998). */
    
    int i;
    int cost = 0;
    charstate *srctemps = src->apomorphies;
    charstate *tgt1apos = tgt1->apomorphies;
    charstate *tgt2apos = tgt2->apomorphies;

    for (i = 0; i < nchar; ++i, ++srctemps, ++tgt1apos, ++tgt2apos) {
        if ( !(*srctemps & (*tgt1apos | *tgt2apos)) ) {
            ++cost;
        }
    }

    return cost;
}

int mfl_locreopt_cost_inapplicables(node *src, node *tgt1, node *tgt2, int nchar, int diff)
{
    /* Returns cost of inserting subtree src between tgt1 and tgt2 following
     * the algorithms described by Ronquist (1998. Cladistics) and Goloboff
     * (1993, 1996. Cladistics).*/
    
    int i;
    int cost = 0;
    charstate *srctemps = src->apomorphies;
    charstate *tgt1apos = tgt1->apomorphies;
    charstate *tgt2apos = tgt2->apomorphies;
    
    
    for (i = 0; i < nchar; ++i, ++srctemps, ++tgt1apos, ++tgt2apos) {
        if ( !(*srctemps & (*tgt1apos | *tgt2apos)) ) {
            ++cost;
            if (cost > diff) {
                return cost;
            }
        }
    }
    
    return cost;
}


void mfl_reopt_subtr_root(node *n, int nchar)
{
    
    int i;
    charstate *lft_chars = n->next->edge->apomorphies;
    charstate *rt_chars = n->next->next->edge->apomorphies;
    charstate *anc_chars = n->apomorphies;
    
    if (!n->tip) {
        
        for (i = 0; i < nchar; ++i) {
            
            /*if (lft_chars[i] & IS_APPLIC && rt_chars[i] & IS_APPLIC) {
                anc_chars[i] = (lft_chars[i] | rt_chars[i]) & IS_APPLIC;
            }
            else {*/
                anc_chars[i] = lft_chars[i] | rt_chars[i];
            //}
        }
    }
}

void mfl_subtree_count(node *leftdesc, node *rightdesc, node *ancestor, int nchar, int *trlength)
{
    int i;
    charstate lft_chars, rt_chars;
    
    assert(!ancestor->tip);
    
    for (i = 0; i < nchar; ++i) {
        if (leftdesc->tempapos[i] & rightdesc->tempapos[i]) 
        {
            if ((leftdesc->tempapos[i] & 1) && (rightdesc->tempapos[i] & 1)) {
                ancestor->tempapos[i] = (leftdesc->tempapos[i] & rightdesc->tempapos[i]) | 1;
            }
            else {
                ancestor->tempapos[i] = leftdesc->tempapos[i] & rightdesc->tempapos[i];
            }
        }
        else
        {
            lft_chars = leftdesc->tempapos[i];
            rt_chars = rightdesc->tempapos[i];
            ancestor->tempapos[i] = lft_chars | rt_chars;

            if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                ancestor->tempapos[i] = ancestor->tempapos[i] & IS_APPLIC;
                if (trlength) {
                    *trlength = *trlength + 1;
                }
            }
        }
    }
}

int mfl_get_sttreelen(tree *testtree, charstate *tipdata, int ntax, int nchar, int *besttreelen)
{
    int treelen = 0;
    int *treelen_p = &treelen;
    
    mfl_apply_tipdata(testtree, tipdata, ntax, nchar);
    //mfl_subtree_postorder(testtree->root, treelen_p, nchar);
    mfl_count_postorder(testtree->root, treelen_p, nchar, NULL);
    treelen = 0;
    mfl_fitch_preorder(testtree->root, nchar, treelen_p);
    
    return *treelen_p;
}


/*
 * FUNCTIONS FOR WAGNER PARSIMONY 
 */
int mfl_wagner_count(charstate lchar, charstate rchar)
{
    int length = 0;
    
    if (lchar==rchar) {
        return 0;
    }
    
    if (lchar > rchar) {
        while (!(lchar & rchar)) {
            lchar = lchar >> 1;
            ++length;
        }
    }
    else {
        while (!(lchar & rchar)) {
            rchar = rchar >> 1;
            ++length;
        }
    }
    return length;
}
/*
 * End Wagner functions
 */


/*
 *
 * FUNCTIONS FOR FITCH PARSIMONY
 * These functions perform the step counting operations for non-additive 
 * characters Functions with "reopt" in the name are employed during the 
 * reoptimization of subtrees during branch-swapping operations.
 *
 */

void mfl_fitch_prelim(node *leftdesc, node *rightdesc, node *ancestor, int nchar, int *trlength, int *besttreelen)
{
    int i;
    charstate lft_chars, rt_chars;
    
    assert(!ancestor->tip);
    
    for (i = 0; i < nchar; ++i) {
        if (leftdesc->tempapos[i] & rightdesc->tempapos[i]) 
        {
            // **ERROR: This bit requires some considerable work!**
            if ((leftdesc->tempapos[i] & 1) || (rightdesc->tempapos[i] & 1)) {
                ancestor->tempapos[i] = (leftdesc->tempapos[i] & rightdesc->tempapos[i]) | 1;
            }
            else {
                ancestor->tempapos[i] = leftdesc->tempapos[i] & rightdesc->tempapos[i];
            }
        }
        else
        {            
            lft_chars = leftdesc->tempapos[i];
            rt_chars = rightdesc->tempapos[i];
            
            ancestor->tempapos[i] = lft_chars | rt_chars;

            if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                ancestor->tempapos[i] = ancestor->tempapos[i] & IS_APPLIC;
                if (trlength) {
                    *trlength = *trlength + 1;
                }
            }
        }
    }
}

void mfl_fitch_prelim_applicables(node *leftdesc, node *rightdesc, node *ancestor, int nchar, int *trlength, int *besttreelen)
{
    int i;
    charstate lft_chars, rt_chars;
    
    assert(!ancestor->tip);
    
    for (i = 0; i < nchar; ++i) {
        if (leftdesc->tempapos[i] & rightdesc->tempapos[i])
        {
            ancestor->tempapos[i] = leftdesc->tempapos[i] & rightdesc->tempapos[i];
        }
        else
        {
            lft_chars = leftdesc->tempapos[i];
            rt_chars = rightdesc->tempapos[i];
            
            ancestor->tempapos[i] = lft_chars | rt_chars;
            *trlength = *trlength + 1;
        }
    }
}

void mfl_fitch_final(node *n, node *anc, int nchar, int *trlength)
{
    
    assert(!n->tip);
    
    int i;
    charstate *ntemps = n->tempapos;
    charstate *napos = n->apomorphies;
    charstate *ancapos = anc->apomorphies;
    charstate lft_chars, rt_chars, temp;
    charstate *lft_c_ptr, *rt_c_ptr;
    
    
    lft_c_ptr = n->next->edge->tempapos;
    rt_c_ptr = n->next->next->edge->tempapos;
    
    for (i = 0; i < nchar; ++i) {
        
        lft_chars = *(lft_c_ptr + i);
        rt_chars = *(rt_c_ptr + i);
        
        if (ancapos[i]==1) {
            if (ntemps[i] & 1) {
                napos[i] = 1;
            }
            else {
                napos[i] = ntemps[i];
                
                if (trlength) {
                    if (!(lft_chars & rt_chars)) {
                        if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                            *trlength = *trlength + 1;
                        }
                    }
                }
            }
        }
        else {
            
            if ((ntemps[i] & ancapos[i]) == ancapos[i]) 
            {
                
                napos[i] = ntemps[i] & ancapos[i];
                assert(napos[i] != 0);
                
                if (!(lft_chars & rt_chars)) {
                    if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                        if (trlength) {
                            *trlength = *trlength + 1;
                        }
                    }
                }
                
            }
            else {
                
                if ( lft_chars & rt_chars ) { //III
                    //V
                    
                    if ((ancapos[i] & IS_APPLIC) && (((lft_chars | rt_chars) & IS_APPLIC) && (lft_chars | rt_chars) & 1)) {
                        temp = ntemps[i] | (ancapos[i] & ((lft_chars & IS_APPLIC) | (rt_chars & IS_APPLIC)));
                        
                        /* This is potentially dangerous as it could cause doubling of counts*/
                        if (!(ancapos[i] & (lft_chars & rt_chars))) {
                            if (trlength) {
                                *trlength = *trlength + 1;
                            }
                        }
                    }
                    else {
                        temp = (ntemps[i] | (ancapos[i] & (lft_chars | rt_chars)));
                    }
                    napos[i] = temp;
                    
                    
                    assert(napos[i] != 0);
                }
                else {
                    //IV
                    
                    napos[i] = ntemps[i] | ancapos[i];
                    
                    if (ntemps[i] & IS_APPLIC && ancapos[i] & IS_APPLIC) {
                        napos[i] = napos[i] & IS_APPLIC;
                    }
                    
                    assert(napos[i] != 0);

                    if (trlength) {
                        if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                            *trlength = *trlength + 1;
                        }
                    }
                }
                
            }
            
        }

    }
}

void mfl_fitch_final_na(node *n, node *anc, int nchar, int *trlength)
{
    
    assert(!n->tip);
    
    int i;
    charstate *ntemps = n->tempapos;
    charstate *napos = n->apomorphies;
    charstate *ancapos = anc->apomorphies;
    charstate lft_chars, rt_chars, temp;
    charstate *lft_c_ptr, *rt_c_ptr;
    
    
    lft_c_ptr = n->next->edge->tempapos;
    rt_c_ptr = n->next->next->edge->tempapos;
    
    for (i = 0; i < nchar; ++i) {
        
        lft_chars = *(lft_c_ptr + i);
        rt_chars = *(rt_c_ptr + i);
        
        if (ancapos[i]==1) {
            if (ntemps[i] & 1) {
                napos[i] = 1;
            }
            else {
                napos[i] = ntemps[i];
                
                if (trlength) {
                    if (!(lft_chars & rt_chars)) {
                        if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                            *trlength = *trlength + 1;
                        }
                    }
                }
            }
        }
        else {
            
            if ((ntemps[i] & ancapos[i]) == ancapos[i])
            {
                
                napos[i] = ntemps[i] & ancapos[i];
                assert(napos[i] != 0);
                
                if (!(lft_chars & rt_chars)) {
                    if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                        if (trlength) {
                            *trlength = *trlength + 1;
                        }
                    }
                }
                
            }
            else {
                
                if ( lft_chars & rt_chars ) { //III
                    //V
                    
                    if ((ancapos[i] & IS_APPLIC) && (((lft_chars | rt_chars) & IS_APPLIC) && (lft_chars | rt_chars) & 1)) {
                        temp = ntemps[i] | (ancapos[i] & ((lft_chars & IS_APPLIC) | (rt_chars & IS_APPLIC)));
                        
                        /* This is potentially dangerous as it could cause doubling of counts*/
                        if (!(ancapos[i] & (lft_chars & rt_chars))) {
                            if (trlength) {
                                *trlength = *trlength + 1;
                            }
                        }
                    }
                    else {
                        temp = (ntemps[i] | (ancapos[i] & (lft_chars | rt_chars)));
                    }
                    napos[i] = temp;
                    
                    
                    assert(napos[i] != 0);
                }
                else {
                    //IV
                    
                    napos[i] = ntemps[i] | ancapos[i];
                    
                    if (ntemps[i] & IS_APPLIC && ancapos[i] & IS_APPLIC) {
                        napos[i] = napos[i] & IS_APPLIC;
                    }
                    
                    assert(napos[i] != 0);
                    
                    if (trlength) {
                        if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                            *trlength = *trlength + 1;
                        }
                    }
                }
                
            }
            
        }
        
    }
}

void mfl_fitch_final_applicables(node *n, node *anc, int nchar, int *trlength)
{
    
    assert(!n->tip);
    
    int i;
    charstate *ntemps = n->tempapos;
    charstate *napos = n->apomorphies;
    charstate *ancapos = anc->apomorphies;
    charstate lft_chars, rt_chars, temp;
    charstate *lft_c_ptr, *rt_c_ptr;
    
    
    lft_c_ptr = n->next->edge->tempapos;
    rt_c_ptr = n->next->next->edge->tempapos;
    
    for (i = 0; i < nchar; ++i) {
        
        lft_chars = *(lft_c_ptr + i);
        rt_chars = *(rt_c_ptr + i);
        
        if ((ntemps[i] & ancapos[i]) == ancapos[i])
        {
            
            napos[i] = ntemps[i] & ancapos[i];
            assert(napos[i] != 0);
            
        }
        else {
            
            if ( lft_chars & rt_chars ) { //III
                //V
                temp = (ntemps[i] | (ancapos[i] & (lft_chars | rt_chars)));
                napos[i] = temp;
                
                assert(napos[i] != 0);
            }
            else {
                //IV
                napos[i] = ntemps[i] | ancapos[i];
                assert(napos[i] != 0);
            }
            
        }
        
    }
}

void mfl_reopt_fitch(node *leftdesc, node *rightdesc, node *ancestor, int nchar, int *changing)
{
    int i, c;
    charstate lft_chars, rt_chars, temp;
    bool allsame = true;
    
    charstate *ldtemps = leftdesc->tempapos;
    charstate *rdtemps = rightdesc->tempapos;
    charstate *antemps = ancestor->tempapos;
    
    for (c = 0; changing[c]; ++c) {
        
        i = changing[c]-1;
        
        
        if (ldtemps[i] & rdtemps[i]) 
        {
            
            if ((ldtemps[i] & 1) && (rdtemps[i] & 1)) {
                temp = (ldtemps[i] & rdtemps[i]) | 1;
            }
            else {
                temp = ldtemps[i] & rdtemps[i];
            }

            assert(temp != 0);
            if (temp != antemps[i]) {
                antemps[i] = temp;
                allsame = false;
            }
        }
        else
        {
            lft_chars = ldtemps[i];
            rt_chars = rdtemps[i];
            
            temp = lft_chars | rt_chars;
            
            if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                temp = temp & IS_APPLIC;
            }
            
            assert(temp != 0);
            
            if (temp != antemps[i]) {
                antemps[i] = temp;
                allsame = false;
                assert(temp != 0);
            }
        }
    }
    if (allsame) {
        ancestor->success = true;
    }
}

void mfl_reopt_fitch_final(node *n, node *anc, int nchar, int *changing)
{
    int i, c;
    charstate lft_chars, rt_chars;
    charstate temp = 0;
    charstate *ntemps = n->tempapos;
    charstate *napos = n->apomorphies;
    charstate *ancapos = anc->apomorphies;
    
    bool allsame = true;
    
    for (c = 0; changing[c]; ++c) {
        
        i = changing[c]-1;
        
        if (ancapos[i] == 1) {
            if (ntemps[i] & 1) {
                temp = 1;
            }
            else {
                temp = ntemps[i];
            }
            if (temp != napos[i]) {
                napos[i] = temp;
                allsame = false;
            }
        }
        else {
            
            if ((ntemps[i] & ancapos[i]) == ancapos[i]) 
            {
                temp = ntemps[i] & ancapos[i];
                assert(temp != 0);
                if (temp != napos[i]) {
                    napos[i] = temp;
                    allsame = false;
                }
            }
            else {
                lft_chars = n->next->edge->tempapos[i];
                rt_chars = n->next->next->edge->tempapos[i];
                
                if ( lft_chars & rt_chars ) { //III
                    //V
                    //temp = ( ntemps[i] | ( ancapos[i] & (lft_chars | rt_chars)));
                    
                    if ((ancapos[i] & IS_APPLIC) && (((lft_chars | rt_chars) & IS_APPLIC) & 1)) {
                        temp = ntemps[i] | (ancapos[i] & ((lft_chars & IS_APPLIC) | (rt_chars & IS_APPLIC)));
                    }
                    else {
                        temp = (ntemps[i] | (ancapos[i] & (lft_chars | rt_chars)));
                    }
                    
                    assert(temp != 0);
                    if (temp != napos[i]) {
                        napos[i] = temp;
                        allsame = false;
                    }
                }
                else {
                    //IV
                        
                    temp = ntemps[i] | ancapos[i];
                    
                    if (ntemps[i] & IS_APPLIC && ancapos[i] & IS_APPLIC) {
                        temp = temp & IS_APPLIC;
                    }
                    
                    assert(temp != 0);
                    
                    if (temp != napos[i]) {
                        napos[i] = temp;
                        allsame = false;
                    }
                }
            }
        }
    }
        
    if (allsame) {
        n->finished = true;
    }
}

void mfl_tip_reopt(tree *t, int ntax, int nchar)
{
    int i;
    for (i = 0; i < ntax; ++i) {
        mfl_tip_apomorphies(t->trnodes[i], t->trnodes[i]->edge, nchar, NULL);
    }
}

void mfl_set_rootstates(node *n, int nchar, int *trlength)
{
    bool allocdtemps = false;
    
    if (!n->tempapos) {
        n->tempapos = (charstate*)malloc(nchar * sizeof(charstate));
        allocdtemps = true;
        assert(n->tempapos);
    }
    
    mfl_fitch_prelim(n->next->edge, n->next->next->edge, n, nchar, trlength, NULL);
    memcpy(n->apomorphies, n->tempapos, nchar * sizeof(charstate));
    
    /*int i;
    for (i = 0; i < nchar; ++i) {
        assert(n->apomorphies[i] == n->tempapos[i]);
    }*/
    
    if (allocdtemps) {
        free(n->tempapos);
        n->tempapos = NULL;
    }
}

void mfl_reopt_rootstates(node *n, int nchar, int *changing)
{
    int i, c;
    charstate lft_chars, rt_chars, temp;
    // cjd - fix compiler warning
    //bool allsame = true;
    
    charstate *ldtemps = n->next->edge->tempapos;
    charstate *rdtemps = n->next->next->edge->tempapos;
    charstate *antemps = n->apomorphies;
    
    for (c = 0; changing[c]; ++c) {
        i = changing[c]-1;
        if (ldtemps[i] & rdtemps[i]) 
        {
            
            if ((ldtemps[i] & 1) && (rdtemps[i] & 1)) {
                temp = (ldtemps[i] & rdtemps[i]) | 1;
            }
            else {
                temp = ldtemps[i] & rdtemps[i];
            }
            
            
            temp = ldtemps[i] & rdtemps[i];
            assert(temp != 0);
            if (temp != antemps[i]) {
                antemps[i] = temp;
                //allsame = false;
            }
        }
        else
        {
            lft_chars = ldtemps[i];
            rt_chars = rdtemps[i];
            
            temp = lft_chars | rt_chars;
            
            if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                temp = temp & IS_APPLIC;
            }
            
            if (temp != antemps[i]) {
                antemps[i] = temp;
                //allsame = false;
            }
        }
    }
    /*if (allsame) {
     n->finished = true;
     }*/
}

void mfl_tip_apomorphies(node *tip, node *anc, int nchar, int *changing)
{
    /* Reconstructs the tip set if it is polymorphic */
    
    int i, c;
    charstate *tiptemp = tip->tempapos;
    charstate *tipapos = tip->apomorphies;
    charstate *ancapos = anc->apomorphies;
    
    if (changing != NULL) {
        for (c = 0; changing[c]; ++c) {
            i = changing[c]-1;
            
            if (tiptemp[i] != IS_APPLIC) {
                if (tiptemp[i] != (tiptemp[i] & (-(tiptemp[i]))) || tiptemp[i] != 1) {
                    if (tiptemp[i] & ancapos[i]) {
                        tipapos[i] = tiptemp[i] & ancapos[i];
                    }
                    else {
                        tipapos[i] = tiptemp[i];
                    }
                }
                else {
                    tipapos[i] = tiptemp[i];
                }
            }
            else {
                tipapos[i] = ancapos[i];
            }
            
        }
    }
    else {
        for (i = 0; i < nchar; ++i) {
            if (tiptemp[i] != IS_APPLIC) {
                if (tiptemp[i] != (tiptemp[i] & (-(tiptemp[i]))) || tiptemp[i] != 1) {
                    if (tiptemp[i] & ancapos[i]) {
                        tipapos[i] = tiptemp[i] & ancapos[i];
                    }
                    else {
                        tipapos[i] = tiptemp[i];
                    }
                }
                else {
                    tipapos[i] = tiptemp[i];
                }
            }
            else {
                tipapos[i] = ancapos[i];
            }
        }
    }
    
}

/*
 * End Fitch functions
 */

void mfl_subtree_postorder(node *n, int *trlength, int nchar)
{
    node *p;
    
    if (!n->apomorphies) {
        n->apomorphies = (charstate*)malloc(nchar * sizeof(charstate));
        memset(n->apomorphies, 0, nchar * sizeof(charstate));
        if (n->next) {
            mfl_join_apomorphies(n);
        }
    }
    
    if (n->tip) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_subtree_postorder(p->edge, trlength, nchar);
        p = p->next;
    }
    if (!n->tempapos) {
        n->tempapos = (charstate*)malloc(nchar * sizeof(charstate));
    }
    mfl_subtree_count(n->next->edge, n->next->next->edge, n, nchar, trlength);
}

void mfl_count_postorder(node *n, int *trlength, int nchar, int *besttreelen)
{
    node *p;
    
    if (!n->apomorphies) {
        n->apomorphies = (charstate*)malloc(nchar * sizeof(charstate));
        memset(n->apomorphies, 0, nchar * sizeof(charstate));
        if (n->next) {
            mfl_join_apomorphies(n);
        }
    }
    if (!n->tempapos) {
        n->tempapos = (charstate*)malloc(nchar * sizeof(charstate));
        memset(n->tempapos, 0, nchar * sizeof(charstate));
    }
    
    if (n->tip) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_count_postorder(p->edge, trlength, nchar, besttreelen);
        p = p->next;
    }

    mfl_fitch_prelim(n->next->edge, n->next->next->edge, n, nchar, trlength, besttreelen);
    n->nodelen = *trlength;
}

void mfl_partial_downpass(node *n, tree *t, int numnodes, int ntax, int nchar, int *changing)
{
    
    node *p;
    p = n;
    
    mfl_erase_clippath(t, numnodes);
    mfl_definish_tree(t, numnodes);
    mfl_desuccess_tree(t, numnodes);
    //mfl_temproot(t, 0, ntax);

    n->edge->clippath = true;
    
    while (p->edge /*&& !p->tip*/) {
        
        if (p->tocalcroot || !p->edge || p->tip) {

            p->clippath = true;
            
            if (!p->tip) {
                mfl_reopt_fitch(p->next->edge, p->next->next->edge, p, nchar, changing);
            }
            
            if (p->success || p->edge->tip || p->tip) {
                if (p->edge->tip || p->tip) {
                    node *lt = p, *rt = p->edge;
                    mfl_join_nodes(t->trnodes[ntax]->next, lt);
                    mfl_join_nodes(t->trnodes[ntax]->next->next, rt);
                    mfl_reopt_preorder(t->trnodes[ntax], nchar, changing);
                    mfl_join_nodes(lt, rt);
                }
                else {
                    mfl_reopt_fitch_final(p, p->edge, nchar, changing);
                    p->finished = false;
                    mfl_reopt_preorder(p, nchar, changing);
                }

                p->success = false;
                break;
            } 
            p = p->edge;
        }
        
        p = p->next;
    }
    
}

bool mfl_reopt_postorder(node *n, int nchar, int *changing)
{
    
    node *p;
    bool fromclip = false;
    bool allsame = true;
     
    if (n->clip) {
        n->clippath = true;
        return true;
    }
    
    if (n->tip) {
        return false;
    }
    
    p = n->next;
    while (p != n) {
        if (mfl_reopt_postorder(p->edge, nchar, changing)) {
            fromclip = true;
        }
        if (!p->edge->success) {
            allsame = false;
        }
        p = p->next;
    }
    
    if (fromclip) {
        if (!allsame) {
            mfl_reopt_fitch(n->next->edge, n->next->next->edge, n, nchar, changing);
            n->changed = true;
        }
        n->clippath = true;
        return true;
    }
    else {
        return false;
    }

}

void mfl_postorder_allviews(node *n, int *trlength, int nchar, int *besttreelen)
{
    node *p;
    int weight = 0;
    
    if (!n->apomorphies) {
        n->apomorphies = (charstate*)malloc(nchar * sizeof(charstate));
        memset(n->apomorphies, 0, nchar * sizeof(charstate));
        if (n->next) {
            mfl_join_apomorphies(n);
        }
    }
    
    if (n->tip) {// || n->visited) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_postorder_allviews(p->edge, trlength, nchar, besttreelen);
        weight = weight + p->edge->vweight;
        p = p->next;
    }
    
    if (!n->tempapos) {
        n->tempapos = (charstate*)malloc(nchar * sizeof(charstate));
        memset(n->tempapos, 0, nchar * sizeof(charstate));
    }

    mfl_fitch_prelim(n->next->edge, n->next->next->edge, n, nchar, trlength, besttreelen);
    n->vweight = weight;
}

void mfl_reopt_preorder(node *n, int nchar, int *changing)
{
    
    node *dl, *dr;
    
    if (n->tip) {
        return;
    }
    
    if (n->finished && !n->clippath) {
        return;
    }
    
    if (!n->edge || n->isroot) {
        mfl_reopt_rootstates(n, nchar, changing);
    }
    
    dl = n->next->edge;
    dr = n->next->next->edge;
    
    if (!dl->tip) {
        mfl_reopt_fitch_final(dl, n, nchar, changing);
    }
    else if (!n->finished) {
        mfl_tip_apomorphies(dl, n, nchar, changing);
    }
    
    if (!dr->tip) {
        mfl_reopt_fitch_final(dr, n, nchar, changing);
    }
    else if (!n->finished) {
        mfl_tip_apomorphies(dr, n, nchar, changing);
    }
    
    mfl_reopt_preorder(n->next->edge, nchar, changing);
    mfl_reopt_preorder(n->next->next->edge, nchar, changing);
    
}

void mfl_reopt_preorder_ii(node *n, int nchar, int *changing)
{
    
    node *dl, *dr;
    
    if (n->tip) {
        return;
    }
    
    if (n->finished) {
        return;
    }
    
    if (!n->edge || n->isroot) {
        mfl_set_rootstates(n, nchar, NULL);
    }
    
    dl = n->next->edge;
    dr = n->next->next->edge;
    
    if (!dl->tip) {
        mfl_reopt_fitch_final(dl, n, nchar, changing);
    }
    else if (!n->finished) {
        mfl_tip_apomorphies(dl, n, nchar, changing);
    }
    
    if (!dr->tip) {
        mfl_reopt_fitch_final(dr, n, nchar, changing);
    }
    else if (!n->finished) {
        mfl_tip_apomorphies(dr, n, nchar, changing);
    }
    
    mfl_reopt_preorder_ii(n->next->edge, nchar, changing);
    mfl_reopt_preorder_ii(n->next->next->edge, nchar, changing);
    
}

void mfl_fitch_preorder(node *n, int nchar, int *trlength)
{
    node *p, *dl, *dr;
    
    if (n->tip) {// || n->clip) {
        return;
    }
    
    if (!n->edge || n->isroot) {
        mfl_set_rootstates(n, nchar, trlength);
    }
    
    dl = n->next->edge;
    dr = n->next->next->edge;
    
    if (!dl->tip) {
        mfl_fitch_final(dl, n, nchar, trlength);
    }
    else {
        mfl_tip_apomorphies(dl, n, nchar, NULL);
    }
    
    if (!dr->tip) {
        mfl_fitch_final(dr, n, nchar, trlength);
    }
    else {
        mfl_tip_apomorphies(dr, n, nchar, NULL);
    }
    
    p = n->next;
    while (p != n) {
        mfl_fitch_preorder(p->edge, nchar, trlength);
        p = p->next;
    }
}

int mfl_get_subtreelen(node *n, int ntax, int nchar, int *besttreelen)
{
    int treelen = 0;
    int *treelen_p = &treelen;
    
    mfl_count_postorder(n, treelen_p, nchar, besttreelen);
    mfl_fitch_preorder(n, nchar, NULL);
    
    return *treelen_p;
}


int mfl_get_treelen(tree *t, int ntax, int nchar, int *besttreelen)
{
    int treelen = 0;
    int *treelen_p = &treelen;
    
    if (!t->root) {
        mfl_temproot(t, 0, ntax);
        mfl_count_postorder(t->root, treelen_p, nchar, besttreelen);
        mfl_fitch_preorder(t->root, nchar, NULL);
        mfl_undo_temproot(ntax, t);
    }
    else {
        mfl_count_postorder(t->root, treelen_p, nchar, besttreelen);
        *treelen_p = 0;
        mfl_fitch_preorder(t->root, nchar, treelen_p);
    }

    return *treelen_p;
}

int *mfl_get_subtr_changing(node *n, node *up, node *dn, int nchar)
{
    /* Compares preliminary and final states in the source tree in heuristic 
     * searches and lists characters that are expected to require reoptimization*/
    
    int i, j;
    
    charstate *prelims = n->tempapos;
    charstate *finals = n->apomorphies;
    int *changing = (int*)malloc((nchar + 1) * sizeof(int));
    
    j = 0;
    
    for (i = 0; i < nchar; ++i) {
        if (prelims[i] != finals[i]) {
            changing[j] = i + 1;
            ++j;
        }
    }
    
    changing[j] = '\0';
    
    return changing;
}

int *mfl_get_tgt_changing(node *n, node *crownprt, node *rootprt, int nchar)
{
    int i, j;
    
    int *changing = (int*)malloc((nchar + 1) * sizeof(int));
    charstate *crownprelims = crownprt->tempapos;
    charstate *rootfinals = rootprt->apomorphies;
    charstate *nprelims = n->tempapos;
    charstate *nfinals = n->apomorphies;
    
    j = 0;
    for (i = 0; i < nchar; ++i) {
        if (crownprelims[i] != nprelims[i] || rootfinals[i] != nfinals[i]) {
            changing[j] = i + 1;
            ++j;
        }
    }
    
    changing[j] = '\0';
    
    return changing;
}

mfl_changing *mfl_get_changing(node *base, node *subtr, node *crownprt, node *rootprt, int nchar)
{
    int i, j, k;
    
    mfl_changing *changing = (mfl_changing*)malloc(sizeof(mfl_changing));
    changing->srcchanging = (int*)malloc((nchar + 1) * sizeof(int));
    changing->tgtchanging = (int*)malloc((nchar + 1) * sizeof(int));
    
    charstate *prelims = base->tempapos;
    charstate *finals = base->apomorphies;
    
    charstate *crownprelims = crownprt->tempapos;
    charstate *rootfinals = rootprt->apomorphies;
    charstate *nprelims = subtr->tempapos;
    charstate *nfinals = subtr->apomorphies;
    
    j = 0;
    k = 0;
    for (i = 0; i < nchar; ++i) {
        if (prelims[i] != finals[i]) {
            *(changing->srcchanging + j) = i + 1;
            ++j;
        }
        if (crownprelims[i] != nprelims[i] || rootfinals[i] != nfinals[i]) {
            *(changing->tgtchanging + k) = i + 1;
            ++k;
        }
        
    }
    
    changing->srcchanging[j] = '\0';
    changing->tgtchanging[k] = '\0';
    
    return changing;
}

void mfl_wipe_states(node *n, int nchar)
{
    memset(n->apomorphies, 0, nchar * sizeof(charstate));
    memset(n->tempapos, 0, nchar * sizeof(charstate));
}

void mfl_set_calcroot(node *n)
{
    
    node *p;
    
    if (n->tip) {
        return;
    }
    
    p = n->next;
    
    while (p != n) {
        p->tocalcroot = false;
        mfl_set_calcroot(p->edge);
        p = p->next;
    }
    
    n->tocalcroot = true;
    
}

void mfl_copy_originals(node *n, charstate *originals, int nchar)
{
    if (!n->origfinals) {
        n->origfinals = (charstate*)malloc(nchar * sizeof(charstate));
    }
    if (!n->origtemps) {
        n->origtemps = (charstate*)malloc(nchar * sizeof(charstate));
    }
    
    memcpy(n->origfinals, n->apomorphies, nchar * sizeof(charstate));
    memcpy(n->origtemps, n->tempapos, nchar * sizeof(charstate));
}

void mfl_restore_originals(node *n, charstate *originals, int nchar)
{
    if (!n->tip) {
        memcpy(n->tempapos, n->origtemps, nchar * sizeof(charstate));
    }
    memcpy(n->apomorphies, n->origfinals, nchar * sizeof(charstate));
}

void mfl_save_origstates(tree *t, int ntax, int numnodes, int nchar)
{
    int i;
    node *p;
    nodearray tns = t->trnodes;
    
    for (i = 0; i < numnodes; ++i) {
        
        p = tns[i];
        
        if (p->next) {
            while (!p->tocalcroot) {
                p = p->next;
                if (p == tns[i]) {
                    break;
                    dbg_printf("Error: node has no apomorphies array\n");
                }
            }
        }
        
        mfl_copy_originals(p, p->tempapos, nchar);
    }
}

void mfl_restore_origstates(tree *t, int ntax, int numnodes, int nchar)
{
    int i;
    node *p;
    nodearray tns = t->trnodes;
    
    for (i = 0; i < numnodes; ++i) {
        p = tns[i];
        
        if (p->next) {
            while (!p->tocalcroot) {
                p = p->next;
                if (p == tns[i]) {
                    break;
                    dbg_printf("Error: node has no apomorphies array\n");
                }
            }
        }
        mfl_restore_originals(p, p->origfinals, nchar);
    }
}

bool mfl_find_lastchanged(node *n, int nchar, int *changing)
{
    node *p;
    //cjd - fix compiler warning
    //bool found = false;
    
    if (n->tip) {
        return false;
    }
    
    if (n->changed) {
        mfl_reopt_preorder(n, nchar, changing);
        return true;
    }
    
    p = n->next;
    while (p != n) {
        if (mfl_find_lastchanged(p->edge, nchar, changing)) {
            //found = true;
            break;
        }
        p = p->next;
    }
    
    return false;
}

void mfl_allviews_traversal(node *n, tree *t, int ntax, int nchar, int *treelen, int *besttreelen)
{
    
    
    if (n->start) {
        mfl_allviews_traversal(n->edge, t, ntax, nchar, treelen, besttreelen);
        return;
    }
    
    
    if (n->tip || n->edge->tip) {
        
        //dbg_printf("working on node %i\n", n->index);
        t->root = NULL;
        
        mfl_join_nodes(n->edge, t->trnodes[ntax]->next->next);
        mfl_join_nodes(n, t->trnodes[ntax]->next);
        
        t->root = t->trnodes[ntax];
        
        //mfl_reopt_postorder(t->root, nchar);
        
        t->root->visited = 0;
        
        mfl_join_nodes(t->trnodes[ntax]->next->next->edge, n);
        
        if (n->tip) {
            return;
        }
    }
    
    mfl_allviews_traversal(n->next->edge, t, ntax, nchar, treelen, besttreelen);
    mfl_allviews_traversal(n->next->next->edge, t, ntax, nchar, treelen, besttreelen);
}

void mfl_reopt_subtr(node *base, tree *t, int nchar, int numnodes, int *changing)
{
    
    /* For subtree reoptimization only. */
    
    base->isroot = true;
    mfl_desuccess_tree(t, numnodes);
    mfl_definish_tree(t, numnodes);
    mfl_reopt_preorder_ii(base, nchar, changing);
    base->isroot = false;
    
}



void mfl_trav_allviews(node *n, tree *t, int ntax, int nchar, int *changing)
{
    
    /* For subtree reoptimization only. */
    
    mfl_definish_tree(t, 2 * ntax - 1);
    mfl_desuccess_tree(t, 2 * ntax - 1);
    mfl_temproot(t, 0, ntax);
    mfl_erase_clippath(t, 2 * ntax - 1);
    mfl_reopt_postorder(t->root, nchar, changing);
    mfl_find_lastchanged(t->root, nchar, changing);
    mfl_undo_temproot(ntax, t);
    
}

int mfl_all_views(tree *t, int ntax, int nchar, int *besttreelen)
{
    int treelen = 0, /*fptreelen,*/ dptreelen;
    int *treelen_p = &treelen;
    int *dptreelen_p = &dptreelen;
    bool wasrooted = false;
    
    mfl_devisit_tree(t->trnodes, 2 * ntax - 1);
    mfl_definish_tree(t, 2 * ntax - 1);

    *treelen_p = 0;
    *dptreelen_p = 0;
    
    if (!t->root) {
        mfl_temproot(t, 0, ntax);
        wasrooted = true;
    }
    
    mfl_postorder_allviews(t->root, treelen_p, nchar, besttreelen);
    t->root->visited = 0;
    
    //fptreelen = *treelen_p;
    mfl_set_calcroot(t->root);

    t->root->visited = 0;
    
    mfl_fitch_preorder(t->root, nchar, dptreelen_p);

    if (wasrooted) {
        mfl_undo_temproot(ntax, t);
    }
    
    //dbg_printf("postorder length: %i\n", fptreelen);
    //dbg_printf("preorder  length: %i\n", dptreelen);
    return dptreelen;
}
