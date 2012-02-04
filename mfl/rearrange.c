/*
 *  rearrange.c
 *  Morphy
 *
 *  Functions used in tree rearrangements for heuristic searches.
 *
 */

#include "morphy.h"

void mfl_bswap(node *p, node *q)
{
    /* Exchanges the position of two nodes */
    
    node *pbase, *qbase;
    
    pbase = p->outedge;
    qbase = q->outedge;
    
    p->outedge = qbase;
    q->outedge = pbase;
    pbase->outedge = q;
    qbase->outedge = p;
    
}

void mfl_remove_branch(node *n)
{
    node *p, *q, *nb;
    
    nb = n->outedge;
    p = nb->next->outedge;
    q = nb->next->next->outedge;
    nb->next->outedge = NULL;
    nb->next->next->outedge = NULL;
    joinNodes(p, q);
}

void mfl_insert_branch(node *br, node *target, int ntax)
{
    // Inserts a branch with a ring base into another branch
    
    node *br1, *br2, *tdesc, *p;
    
    tdesc = target->outedge;
    
    if (br->tip) {
        br1 = br->outedge->next;
        br2 = br1->next;
    }
    else {
        br1 = mfl_seek_ringnode(br, ntax);
        br2 = mfl_seek_ringnode(br1->next, ntax);
    }

    if (br1->outedge || br2->outedge) {
        printf("Error in branch insertion\n");
        return;
    }
    
    joinNodes(br1, target);
    joinNodes(br2, tdesc);
}

/* Nearest-neighbor interchange (NNI) */

void mfl_nni_traversal(node *n, tree *swapingon, tree **savedtrees, int ntax, 
                       int nchar, int numnodes, long int *current, 
                       charstate *tipdata, bool *undertreelimit, 
                       long int *currentbesttree, bool *foundbettertree)
{
    int trlength = 0;
    node *p;
    
    if (n->start) {
        mfl_nni_traversal(n->outedge, swapingon, savedtrees, ntax, nchar, numnodes, 
                          current, tipdata, undertreelimit, currentbesttree,
                          foundbettertree);
        return;
    }
    
    if (n->tip || !(*undertreelimit)) {
        return;
    }
    
    if (!n->outedge->tip)
    {
        /* Nearly identical steps are conducted twice because all nearest-
         * neighbor interchanges produce two distinct tree topologies */
        
        if (*current + 1 <= TREELIMIT) {
            mfl_bswap(n->next->outedge, n->outedge->next->outedge);
            mfl_root_tree(swapingon, 1, ntax);
            trlength = mfl_get_treelen(swapingon, tipdata, ntax, nchar);
            unroot(ntax, swapingon);
            if (trlength < *currentbesttree) {
                *foundbettertree = true;
                *currentbesttree = trlength;
                swapingon->length = trlength;
                mfl_reinit_treebuffer(savedtrees, swapingon, current, numnodes);
                trlength = 0;
                *current = *current + 1;
                return;
            }
            if (trlength == *currentbesttree) {
                *foundbettertree = false;
                savedtrees[*current] = copytree(swapingon, ntax, numnodes);
                savedtrees[*current]->index = *current;
                savedtrees[*current]->length = trlength;
                trlength = 0;
                *current = *current + 1;
            }
            mfl_bswap(n->next->outedge, n->outedge->next->outedge);
            
            if (*current + 1 <= TREELIMIT){ 
                mfl_bswap(n->next->next->outedge, n->outedge->next->outedge);
                mfl_root_tree(swapingon, 1, ntax);
                trlength = mfl_get_treelen(swapingon, tipdata, ntax, nchar);
                unroot(ntax, swapingon);
                if (trlength < *currentbesttree) {
                    *foundbettertree = true;
                    *currentbesttree = trlength;
                    swapingon->length = trlength;
                    mfl_reinit_treebuffer(savedtrees, swapingon, current, numnodes);
                    trlength = 0;
                    *current = *current + 1;
                    return;
                }
                if (trlength == *currentbesttree) {
                    *foundbettertree = false;
                    savedtrees[*current] = copytree(swapingon, ntax, numnodes);
                    savedtrees[*current]->index = *current;
                    savedtrees[*current]->length = trlength;
                    trlength = 0;
                    *current = *current + 1;
                }
                mfl_bswap(n->next->next->outedge, n->outedge->next->outedge);
            }
        }
        else {
            printf("Hit tree limit\n");
            *undertreelimit = false;
            return;
        }
    }
    
    p = n->next;
    while (p != n) {
        mfl_nni_traversal(p->outedge, swapingon, savedtrees, ntax, nchar, numnodes, 
                          current, tipdata, undertreelimit, currentbesttree,
                          foundbettertree);
        p = p->next;
    }
}


/* Quite a lot about this function will change. For now, it interfaces with the
 * test in main.c but that won't always be the case. */
void mfl_nni_search(int ntax, int nchar, int numnodes, charstate *tipdata, 
                    tree **savedtrees, int starttreelen)
{
    long int i = 0, j, currentbest = starttreelen;
    long int trbufpos = 0;
    long int *nxtintrbuf = &trbufpos, *currentbest_p = &currentbest;
    long int numreps = 100;
    bool undertreelimit = true, *undertreelimit_p = &undertreelimit;
    bool foundbettertree = false, *foundbettertree_p = &foundbettertree;
    
    do {
        for (j = 0; j == 0 || j < *nxtintrbuf; ++j) {
            mfl_nni_traversal(savedtrees[j]->trnodes[0], savedtrees[j], savedtrees, 
                              ntax, nchar, numnodes, nxtintrbuf, tipdata, 
                              undertreelimit_p, currentbest_p, foundbettertree_p);
        }
        ++i;
    } while (i < numreps);
    
    printf("Next in trbuf: %li\n", *nxtintrbuf);
    
    printf("\nThe optimal tree(s) found by nearest neighbor interchange:\n");
    for (i = 0; i < *nxtintrbuf; ++i) {
        printf("Tree %li:\n", i + 1);
        mfl_root_tree(savedtrees[i], 0, ntax);
        printNewick(savedtrees[i]->root);
        printf(";\n");
        printf("Length: %i\n", savedtrees[i]->length);
    }
    printf("\n\n");
    mfl_clear_treebuffer(savedtrees, nxtintrbuf, numnodes);
    
    free(savedtrees);
}

long int mfl_spr_leftotry(int ntax)
{
    long int rleft = 0;
    rleft = 4 * (ntax - 3) * (ntax - 2);
    return rleft;
}

void mfl_regrafting_traversal(node *n, node *subtr, tree *swapingon, 
                              tree **savedtrees, int ntax, int nchar, int numnodes, 
                              long int *current, charstate *tipdata, 
                              bool *undertreelimit, long int *currentbesttree, 
                              bool *foundbettertree, long int *leftotry)
{
    
    if (*foundbettertree) {
        return;
    }
    
    int trlength;
    node *p, *up;
    
    if (!n->initialized || !n->outedge->initialized) {
        up = n->outedge;
        mfl_insert_branch(subtr, up, ntax);
        mfl_temproot(swapingon, 0, ntax);
        trlength = mfl_get_treelen(swapingon, tipdata, ntax, nchar);
        mfl_undo_temproot(ntax, swapingon);
        *leftotry = *leftotry - 1;
        subtr->skip = 1;
        subtr->next->skip = 1;
        subtr->next->next->skip = 1;
        if (trlength < *currentbesttree) {
            *foundbettertree = true;
            swapingon->length = trlength;
            *currentbesttree = trlength;
            mfl_reinit_treebuffer(savedtrees, swapingon, current, numnodes);
            trlength = 0;
            *current = *current + 1;
            *leftotry = mfl_spr_leftotry(ntax);
            return;
        }
        *foundbettertree = false;
        if (trlength == *currentbesttree) {
            printf("saving an equally parsimonious tree of length: %i\n", trlength);
            printf("copying into position: %li\n", *current);
            //subtr->skip = 1;
            //subtr->next->skip = 1;
            //subtr->next->next->skip = 1;
            savedtrees[*current] = copytree(swapingon, ntax, numnodes);
            savedtrees[*current]->index = *current;
            savedtrees[*current]->length = trlength;
            trlength = 0;
            *current = *current + 1;
        }
        joinNodes(n, up);
    }
    
    if (n->start) {
        mfl_regrafting_traversal(n->outedge, subtr, swapingon, savedtrees, ntax,
                                 nchar, numnodes, current, tipdata, 
                                 undertreelimit, currentbesttree, 
                                 foundbettertree, leftotry);
        return;
    }
    
    if (n->tip) {
        return;
    }

    p = n->next;
    while (p != n && !(*foundbettertree)) {
        mfl_regrafting_traversal(p->outedge, subtr, swapingon, savedtrees, ntax,
                                 nchar, numnodes, current, tipdata, 
                                 undertreelimit, currentbesttree, 
                                 foundbettertree, leftotry);
        if (*foundbettertree) {
            return;
        }
        p = p->next;
    }
}

void mfl_pruning_traversal(node *n, tree *swapingon, tree **savedtrees, int ntax, 
                           int nchar, int numnodes, long int *current, 
                           charstate *tipdata, bool *undertreelimit, 
                           long int *currentbesttree, bool *foundbettertree, long int *leftotry)
{
    node *p, *up, *dn, *subtr;
    
    if (n->start) {
        mfl_pruning_traversal(n->outedge, swapingon, savedtrees, ntax,
                              nchar, numnodes, current, tipdata, 
                              undertreelimit, currentbesttree, foundbettertree, 
                              leftotry);
        if (*foundbettertree) {
            return;
        }
        *leftotry = *leftotry - 1;
        //printf("pruning a subtree\n");
        up = n->outedge->next->outedge;
        dn = n->outedge->next->next->outedge;
        subtr = n->outedge->next->next;
        up->initialized = 1;
        dn->initialized = 1;
        joinNodes(up, dn);
        *leftotry = *leftotry - 1;
        swapingon->trnodes[1]->start = true;
        n->start = false;
        mfl_regrafting_traversal(swapingon->trnodes[1], subtr, swapingon, 
                                 savedtrees, ntax, nchar, numnodes, current, 
                                 tipdata, undertreelimit, currentbesttree, 
                                 foundbettertree, leftotry);
        up->initialized = 0;
        dn->initialized = 0;
        swapingon->trnodes[1]->start = false;
        n->start = true;
        if (*foundbettertree) {
            return;
        }
        joinNodes(up, n->outedge->next);
        joinNodes(dn, n->outedge->next->next);
        return;
    }
    
    if (n->tip || *foundbettertree) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_pruning_traversal(p->outedge, swapingon, savedtrees, ntax,
                              nchar, numnodes, current, tipdata, 
                              undertreelimit, currentbesttree, foundbettertree, 
                              leftotry);
        if (*foundbettertree) {
            return;
        }
        if (!p->outedge->skip) {
            *leftotry = *leftotry - 1;
            //printf("pruning a subtree\n");
            up = p->next->outedge;
            dn = p->next->next->outedge;
            subtr = p->next->next;
            up->initialized = 1;
            dn->initialized = 1;
            joinNodes(up, dn);
            *leftotry = *leftotry - 1;
            mfl_regrafting_traversal(swapingon->trnodes[0], subtr, swapingon, 
                                     savedtrees, ntax, nchar, numnodes, current, 
                                     tipdata, undertreelimit, currentbesttree, 
                                     foundbettertree, leftotry);
            up->initialized = 0;
            dn->initialized = 0;
            if (*foundbettertree) {
                return;
            }
            joinNodes(up, p->next);
            joinNodes(dn, p->next->next);
        }
        else {
            printf("Skipping already swapped clade\n");
            //p->outedge->skip = 0;
            //p->outedge->next->skip = 0;
            //p->outedge->next->next->skip = 0;
        }
        p = p->next;
    }   
}

void mfl_spr_search(int ntax, int nchar, int numnodes, charstate *tipdata, 
                    tree **savedtrees, int starttreelen)
{
    long int i = 0, j, currentbest = starttreelen, remaining;
    long int trbufpos = 0;
    long int *nxtintrbuf = &trbufpos, *currentbest_p = &currentbest;
    long int leftotry = 0;
    long int *leftotry_p = &leftotry;
    long int c_in = 0;
    bool undertreelimit = true, *undertreelimit_p = &undertreelimit;
    bool foundbettertree = false, *foundbettertree_p = &foundbettertree;
    bool quit = false;
    bool lastround = false;
    
    *leftotry_p = mfl_spr_leftotry(ntax);
    
    do {
        printf("number of trees in buffer: %li\n", *nxtintrbuf+1);
        remaining = *nxtintrbuf + 1;
        c_in = *nxtintrbuf;
        printf("current IN: %li\n", c_in);
        printf("Looking for shorter trees on all saved trees:\n");
        *foundbettertree_p = false;
        mfl_pruning_traversal(savedtrees[j]->trnodes[0], savedtrees[j], savedtrees, 
                                ntax, nchar, numnodes, nxtintrbuf, tipdata, 
                                undertreelimit_p, currentbest_p, 
                                foundbettertree_p, leftotry_p);
        savedtrees[*nxtintrbuf-1]->swapped = true;
        printf("j = %li\n", j);
        printf("current: %li\n", *nxtintrbuf);                      
        if (*foundbettertree_p) {
            j = 0;
        }
        else {
            ++j;
        }
    } while (/*!quit*/j < 10);
    
    printf("Next in trbuf: %li\n", *nxtintrbuf);
    
    printf("\nThe optimal tree(s) found by subtree pruning and regrafting:\n");
    for (i = 0; i < *nxtintrbuf; ++i) {
        //printf("Tree %li:\n", i + 1);
        printf("TREE str_%li = [&u] ", i+1);
        mfl_root_tree(savedtrees[i], 0, ntax);
        printNewick(savedtrees[i]->root);
        printf(";\n");
        //printf("Length: %i\n", savedtrees[i]->length);
    }
    printf("\n");
    
    mfl_clear_treebuffer(savedtrees, nxtintrbuf, numnodes);
    free(savedtrees);
}

/* Tree bisection and reconnection (TBR) */