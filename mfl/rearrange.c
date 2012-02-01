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

void mfl_nni_traversal(node *n, tree *swapingon, tree **treeset, int ntax, 
                       int nchar, int numnodes, long int *current, 
                       charstate *tipdata, bool *undertreelimit, 
                       long int *currentbesttree, bool *foundbettertree)
{
    int trlength = 0;
    node *p;
    
    if (n->start) {
        mfl_nni_traversal(n->outedge, swapingon, treeset, ntax, nchar, numnodes, 
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
                mfl_reinit_treebuffer(treeset, swapingon, current, numnodes);
                trlength = 0;
                *current = *current + 1;
                return;
            }
            if (trlength == *currentbesttree) {
                *foundbettertree = false;
                treeset[*current] = copytree(swapingon, ntax, numnodes);
                treeset[*current]->index = *current;
                treeset[*current]->length = trlength;
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
                    mfl_reinit_treebuffer(treeset, swapingon, current, numnodes);
                    trlength = 0;
                    *current = *current + 1;
                    return;
                }
                if (trlength == *currentbesttree) {
                    *foundbettertree = false;
                    treeset[*current] = copytree(swapingon, ntax, numnodes);
                    treeset[*current]->index = *current;
                    treeset[*current]->length = trlength;
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
        mfl_nni_traversal(p->outedge, swapingon, treeset, ntax, nchar, numnodes, 
                          current, tipdata, undertreelimit, currentbesttree,
                          foundbettertree);
        p = p->next;
    }
}


/* Quite a lot about this function will change. For now, it interfaces with the
 * test in main.c but that won't always be the case. */
void mfl_nni_search(int ntax, int nchar, int numnodes, charstate *tipdata, 
                    tree **treeset, int starttreelen)
{
    long int i = 0, j, currentbest = starttreelen;
    long int trbufpos = 0;
    long int *nxtintrbuf = &trbufpos, *currentbest_p = &currentbest;
    long int numreps = 100;
    bool undertreelimit = true, *undertreelimit_p = &undertreelimit;
    bool foundbettertree = false, *foundbettertree_p = &foundbettertree;
    
    do {
        for (j = 0; j == 0 || j < *nxtintrbuf; ++j) {
            mfl_nni_traversal(treeset[j]->trnodes[0], treeset[j], treeset, 
                              ntax, nchar, numnodes, nxtintrbuf, tipdata, 
                              undertreelimit_p, currentbest_p, foundbettertree_p);
        }
        ++i;
    } while (i < numreps);
    
    printf("Next in trbuf: %li\n", *nxtintrbuf);
    
    printf("\nThe optimal tree(s) found by nearest neighbor interchange:\n");
    for (i = 0; i < *nxtintrbuf; ++i) {
        printf("Tree %li:\n", i + 1);
        mfl_root_tree(treeset[i], 0, ntax);
        printNewick(treeset[i]->root);
        printf(";\n");
        printf("Length: %i\n", treeset[i]->length);
    }
    printf("\n\n");
    mfl_clear_treebuffer(treeset, nxtintrbuf, numnodes);
    
    free(treeset);
}

/* Subtree pruning and regrafting (SPR) */

/*    
 if (*current + 1 <= TREELIMIT && !n->initialized) {
 up = n->outedge;
 printf("Before the insertion\n");
 //dump_connections(swapingon->trnodes, ntax, numnodes);
 printNewick(swapingon->trnodes[0]);
 printf("\n");
 joinNodes(n, subtr->next->next);
 joinNodes(up, subtr);
 printf("After the insertion\n");
 //dump_connections(swapingon->trnodes, ntax, numnodes);
 printNewick(swapingon->trnodes[0]);
 printf("\n");
 mfl_temproot(swapingon, 0, ntax);
 printf("Rearrangments left: %lu\n", *leftotry);
 ++r;
 trlength = mfl_get_treelen(swapingon, tipdata, ntax, nchar);
 mfl_undo_temproot(ntax, swapingon);
 
 *leftotry = *leftotry - 1;
 
 if (trlength < *currentbesttree) {
 *foundbettertree = true;
 *currentbesttree = trlength;
 swapingon->length = trlength;
 mfl_reinit_treebuffer(treeset, swapingon, current, numnodes);
 //dump_connections(swapingon->trnodes, ntax, numnodes);
 trlength = 0;
 *current = *current + 1;
 *leftotry = mfl_spr_leftotry(ntax);
 printf("Rearrangments left: %li\n", *leftotry);
 return;
 }
 
 *foundbettertree = false;
 
 if (trlength == *currentbesttree) {
 printf("saving an equally parsimonious tree\n");
 treeset[*current] = copytree(swapingon, ntax, numnodes);
 treeset[*current]->index = *current;
 treeset[*current]->length = trlength;
 trlength = 0;
 *current = *current + 1;
 }
 joinNodes(n, up);
 }
 else if ( *current + 1 >= TREELIMIT) {
 printf("Hit tree limit\n");
 *undertreelimit = false;
 return;
 }    
 
 */

long int mfl_spr_leftotry(int ntax)
{
    long int rleft = 0;
    rleft = 4 * (ntax - 3) * (ntax - 2);
    return rleft;
}

void mfl_regrafting_traversal(node *n, node *subtr, tree *swapingon, 
                              tree **treeset, int ntax, int nchar, int numnodes, 
                              long int *current, charstate *tipdata, 
                              bool *undertreelimit, long int *currentbesttree, 
                              bool *foundbettertree, long int *leftotry)
{
    
    if (*foundbettertree) {
        return;
    }
    
    int trlength;
    node *p, *up;
    
    up = n->outedge;
    mfl_insert_branch(subtr, up, ntax);
    printNewick(swapingon->trnodes[0]);
    printf("\n");
    mfl_temproot(swapingon, 0, ntax);
    //printf("Rearrangments left: %lu\n", *leftotry);
    trlength = mfl_get_treelen(swapingon, tipdata, ntax, nchar);
    mfl_undo_temproot(ntax, swapingon);
    
    if (trlength < *currentbesttree) {
        printf("Found shorter tree:\n");
        printNewick(swapingon->trnodes[0]);
        printf("\nBest tree length: %i\n", trlength);
        *foundbettertree = true;
        mfl_reinit_treebuffer(treeset, swapingon, current, numnodes);
        swapingon->length = trlength;
        trlength = 0;
        *current = *current + 1;
        *leftotry = mfl_spr_leftotry(ntax);
        return;
    }
    if (trlength == *currentbesttree) {
        printf("saving an equally parsimonious tree\n");
        printNewick(swapingon->trnodes[0]);
        printf("\n");
        treeset[*current] = copytree(swapingon, ntax, numnodes);
        treeset[*current]->index = *current;
        treeset[*current]->length = trlength;
        trlength = 0;
        *current = *current + 1;
    }    
    //subtr->outedge = NULL;
    //subtr->next->next->outedge = NULL;
    joinNodes(n, up);
    
    if (n->start) {
        mfl_regrafting_traversal(n->outedge, subtr, swapingon, treeset, ntax,
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
        mfl_regrafting_traversal(p->outedge, subtr, swapingon, treeset, ntax,
                                 nchar, numnodes, current, tipdata, 
                                 undertreelimit, currentbesttree, 
                                 foundbettertree, leftotry);
        if (*foundbettertree) {
            return;
        }
        p = p->next;
    }
}

void mfl_pruning_traversal(node *n, tree *swapingon, tree **treeset, int ntax, 
                           int nchar, int numnodes, long int *current, 
                           charstate *tipdata, bool *undertreelimit, 
                           long int *currentbesttree, bool *foundbettertree, long int *leftotry)
{
    node *p, *up, *dn, *subtr;
    
    printf("starting pruning\n");
    
    if (n->start) {
        mfl_pruning_traversal(n->outedge, swapingon, treeset, ntax,
                              nchar, numnodes, current, tipdata, 
                              undertreelimit, currentbesttree, foundbettertree, 
                              leftotry);
        return;
    }
    
    if (n->tip || *foundbettertree) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_pruning_traversal(p->outedge, swapingon, treeset, ntax,
                              nchar, numnodes, current, tipdata, 
                              undertreelimit, currentbesttree, foundbettertree, 
                              leftotry);
        if (*foundbettertree) {
            return;
        }
        printf("pruning a subtree\n");
        up = p->next->outedge;
        dn = p->next->next->outedge;
        dn->outedge = NULL;
        up->outedge = NULL;
        p->next->outedge = NULL;
        p->next->next->outedge = NULL;
        subtr = p->next->next;
        up->initialized = 1;
        dn->initialized = 1;
        joinNodes(up, dn);
        printNewick(swapingon->trnodes[0]);
        printf("\n");
        *leftotry = *leftotry - 1;
        //dump_connections(swapingon->trnodes, ntax, numnodes);
        mfl_regrafting_traversal(swapingon->trnodes[0], subtr, swapingon, 
                                 treeset, ntax, nchar, numnodes, current, 
                                 tipdata, undertreelimit, currentbesttree, 
                                 foundbettertree, leftotry);
        up->initialized = 0;
        dn->initialized = 0;
        if (*foundbettertree) {
            return;
        }
        joinNodes(up, p->next);
        joinNodes(dn, p->next->next);
        p = p->next;
    }   
}

void mfl_spr_search(int ntax, int nchar, int numnodes, charstate *tipdata, 
                    tree **treeset, int starttreelen)
{
    long int i = 0, j, currentbest = starttreelen, remaining;
    long int trbufpos = 0;
    long int *nxtintrbuf = &trbufpos, *currentbest_p = &currentbest;
    long int leftotry = 0;
    long int *leftotry_p = &leftotry;
    long int numreps = 200, c_in = 0;
    bool undertreelimit = true, *undertreelimit_p = &undertreelimit;
    bool foundbettertree = false, *foundbettertree_p = &foundbettertree;
    bool quit = false;
    
    *leftotry_p = mfl_spr_leftotry(ntax);
    
    do {
        printf("number of trees in buffer: %li\n", *nxtintrbuf+1);
        remaining = *nxtintrbuf + 1;
        c_in = *nxtintrbuf;
        printf("current IN: %li\n", c_in);
        for (j = 0; j == 0 || j < *nxtintrbuf; ++j) {
            //printf("%li global rearrangements left.\n", (*nxtintrbuf) * (mfl_spr_leftotry(ntax)));
            printf("Trying tree: %li\n", j+1);
            *foundbettertree_p = false;
            mfl_pruning_traversal(treeset[j]->trnodes[0], treeset[j], treeset, 
                                  ntax, nchar, numnodes, nxtintrbuf, tipdata, 
                                  undertreelimit_p, currentbest_p, 
                                  foundbettertree_p, leftotry_p);
            printf("\nleaving the traversal\n\n");
            remaining = remaining - 1;
            if (!*foundbettertree_p) {
                quit = true;
                j = *nxtintrbuf;
            }
        }
        ++i;
        //printf("current OUT: %li\n", *nxtintrbuf);
    } while (/*!quit*/ i < numreps && *leftotry_p >= 0);
    
    printf("Next in trbuf: %li\n", *nxtintrbuf);
    
    printf("\nThe optimal tree(s) found by subtree pruning and regrafting:\n");
    for (i = 0; i < *nxtintrbuf; ++i) {
        //printf("Tree %li:\n", i + 1);
        printf("TREE str_%li = [&u] ", i+1);
        mfl_root_tree(treeset[i], 0, ntax);
        printNewick(treeset[i]->root);
        printf(";\n");
        printf("Length: %i\n", treeset[i]->length);
    }
    printf("\n");
    
    mfl_clear_treebuffer(treeset, nxtintrbuf, numnodes);
    free(treeset);
}


/* Tree bisection and reconnection (TBR) */