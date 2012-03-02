/*
 *  rearrange.c
 *  Morphy
 *
 *  Functions used in tree rearrangements for heuristic searches.
 *
 */

#include "morphy.h"
#include <time.h>

mfl_searchrec * mfl_create_searchrec(void)
{
    mfl_searchrec *newsearchrec;
    
    return newsearchrec = (mfl_searchrec*)malloc(sizeof(mfl_searchrec));
}

void mfl_destroy_searchrec(mfl_searchrec *searchrec)
{
    free(searchrec);
}

long long int mfl_rearr_num(bool reset)
{
    static long long int rn = 0;

    if (reset) {
        rn = 0;
    }
    
    return ++rn;
}


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
    mfl_join_nodes(p, q);
}

void mfl_insert_branch(node *br, node *target, int ntax)
{
    // Inserts a branch with a ring base into another branch
    
    node *br1, *br2, *tdesc;
    
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
    
    mfl_join_nodes(br1, target);
    mfl_join_nodes(br2, tdesc);
}

/* Nearest-neighbor interchange (NNI) */

void mfl_nni_traversal(node *n, tree *swapingon, tree **savedtrees, int ntax, 
                       int nchar, int numnodes, long int *current, 
                       charstate *tipdata, bool *undertreelimit, 
                       int *currentbesttree, bool *foundbettertree)
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
            trlength = mfl_get_treelen(swapingon, ntax, nchar, currentbesttree);
            mfl_unroot(ntax, swapingon);
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
                savedtrees[*current] = mfl_copytree(swapingon, ntax, numnodes);
                savedtrees[*current]->index = *current;
                savedtrees[*current]->length = trlength;
                trlength = 0;
                *current = *current + 1;
            }
            mfl_bswap(n->next->outedge, n->outedge->next->outedge);
            
            if (*current + 1 <= TREELIMIT){ 
                mfl_bswap(n->next->next->outedge, n->outedge->next->outedge);
                mfl_root_tree(swapingon, 1, ntax);
                trlength = mfl_get_treelen(swapingon, ntax, nchar, currentbesttree);
                mfl_unroot(ntax, swapingon);
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
                    savedtrees[*current] = mfl_copytree(swapingon, ntax, numnodes);
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
    long int i = 0, j = 0; 
    int currentbest = starttreelen;
    long int trbufpos = 0;
    long int *nxtintrbuf = &trbufpos; 
    int *currentbest_p = &currentbest;
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
                              long int *current, 
                              bool *undertreelimit, int *currentbesttree, 
                              bool *foundbettertree, bool *success, long int *leftotry, int diff)
{
    /* Called from within any subtree pruning algorithm used in either SPR or
     * TBR branch swapping. Traverses a binary tree in preorder, inserting the 
     * subtree at each node and (hopefully) skipping a reinsertion at the original
     * site of pruning (flagged by the "visited" boolean values in those nodes) */
    
    if (*success) {
        return;
    }
    
    if (n->start) {
        mfl_regrafting_traversal(n->outedge, subtr, swapingon, savedtrees, ntax,
                                 nchar, numnodes, current, 
                                 undertreelimit, currentbesttree, 
                                 foundbettertree, success, leftotry, diff);
        return;
    }
    
    int trlength = 0;
    node *up;
    
    n->pathid = subtr->next->outedge->stid;
    n->outedge->pathid = subtr->next->outedge->stid;
    
    if ((!n->visited || !n->outedge->visited) /*&& n->skippath != n->stid*/) {
        up = n->outedge;
        trlength = *currentbesttree - diff + mfl_locreopt_cost(subtr->next->outedge, n, up, nchar, diff);
        
        //n->skippath = n->stid;
        //n->outedge->skippath = n->stid;
        
        *leftotry = *leftotry + 1;
        //printf("Left to try: %li\n", *leftotry);
        if (trlength < *currentbesttree) {
            mfl_insert_branch(subtr, up, ntax);
            *foundbettertree = true;
            *success = true;
            swapingon->length = trlength;
            *currentbesttree = trlength;
            mfl_reinit_treebuffer(savedtrees, swapingon, current, numnodes);
            trlength = 0;
            *current = *current + 1;
            //*leftotry = mfl_spr_leftotry(ntax);
            swapingon->bipartitions = mfl_tree_biparts(swapingon, ntax, numnodes);
            return;
        }
        *foundbettertree = false;
        *success = false;
        if (trlength == *currentbesttree) 
        {
            mfl_insert_branch(subtr, up, ntax);
            if (!mfl_compare_alltrees(swapingon, savedtrees, ntax, numnodes, current)) 
            {
                //printf("saving an equally parsimonious tree of length: %i; writing to %li\n", trlength, *current);
                savedtrees[*current] = mfl_copytree(swapingon, ntax, numnodes);
                savedtrees[*current]->index = *current;
                savedtrees[*current]->length = trlength;
                savedtrees[*current]->swapped = false;
                *current = *current + 1;
            }
            else {
                //printf("tree is duplicate\n");
                mfl_join_nodes(n, up);
            }
            trlength = 0;
        }
        mfl_join_nodes(n, up);
    }
    
    if (n->tip) {
        return;
    }

    mfl_regrafting_traversal(n->next->outedge, subtr, swapingon, savedtrees, ntax,
                                 nchar, numnodes, current, 
                                 undertreelimit, currentbesttree, 
                                 foundbettertree, success, leftotry, diff);
    if (*foundbettertree) {
        return;
    }
    mfl_regrafting_traversal(n->next->next->outedge, subtr, swapingon, savedtrees, ntax,
                             nchar, numnodes, current, 
                             undertreelimit, currentbesttree, 
                             foundbettertree, success, leftotry, diff);
    if (*foundbettertree) {
        return;
    }
}

void mfl_pruning_traversal(node *n, tree *swapingon, tree **savedtrees, int ntax, 
                           int nchar, int numnodes, long int *current, bool *undertreelimit, 
                           int *currentbesttree, bool *foundbettertree, bool *success, long int *leftotry)
{

    /* Traverses a binary tree clipping out a subtree in postorder and passing 
     * a pointer to the subtree to mfl_regrafting_traversal. */
    
    node *p, *clipnode, *up, *dn, *subtr;
    int diff = 0;
    int trl = 0, cbest = 0;
    int *trlp = &trl, *cbestp = &cbest;
    
    if (n->start) {
        mfl_pruning_traversal(n->outedge, swapingon, savedtrees, ntax,
                              nchar, numnodes, current, 
                              undertreelimit, currentbesttree, foundbettertree, 
                              success, leftotry);
        return;
    }
    
    if (n->tip || *success) {
        return;
    }
    
    
    clipnode = n->next->next->outedge;
    
    p = n->next;
    while (p != n) {
        mfl_pruning_traversal(p->outedge, swapingon, savedtrees, ntax,
                              nchar, numnodes, current, 
                              undertreelimit, currentbesttree, foundbettertree, 
                              success, leftotry);
        if (*success) {
            return;
        }
        up = p->next->outedge;
        dn = p->next->next->outedge;
        
        subtr = p->next->next;
        
        //if (!subtr->next->outedge->skip) {
            
            //subtr->next->outedge->skip = true;
            
            if (!(up->tip && dn->tip)) {
                up->visited = 1;
                dn->visited = 1;
                
                clipnode->clip = true;
                mfl_join_nodes(up, dn);
                
                // Reoptimize the clipped tree
                //printf("reoptimizing cliptree\n");
                mfl_trav_allviews(swapingon->trnodes[0], swapingon, ntax, nchar, trlp, cbestp);
                
                // Reoptimize the subtree
                //mfl_reopt_postorder(subtr->next->outedge, nchar);
                mfl_reopt_subtr_root(subtr->next->outedge, nchar);
                
                // Determine the cost of local reinsertion
                diff = mfl_subtr_reinsertion(subtr->next->outedge, up, dn, nchar);
                
                //mfl_tip_reopt(swapingon, ntax, nchar);
                
                //*leftotry = *leftotry - 1;
                mfl_regrafting_traversal(swapingon->trnodes[0], subtr, swapingon, 
                                         savedtrees, ntax, nchar, numnodes, current, 
                                         undertreelimit, currentbesttree, 
                                         foundbettertree, success, leftotry, diff);
                
                up->visited = 0;
                dn->visited = 0;
                clipnode->clip = false;
                mfl_devisit_tree(swapingon->trnodes, numnodes);
                if (*success) {
                    return;
                }
                subtr->next->outedge->skip = false;
                mfl_join_nodes(up, p->next);
                mfl_join_nodes(dn, p->next->next);
            }
        //}
        /*else {
            mfl_trav_allviews(swapingon->trnodes[0], swapingon, ntax, nchar, trlp, cbestp);
        }*/
        p = p->next;
        clipnode = n->next->outedge;
    }   
}

/* Tree bisection and reconnection (TBR) */

void mfl_reroot_subtree(node *n, node *subtr, node *up, node *dn, tree *swapingon, 
                        tree **savedtrees, int ntax, int nchar, int numnodes, 
                        long int *current, 
                        bool *undertreelimit, int *currentbesttree, 
                        bool *foundbettertree, bool *success, long int *leftotry, int diff)
{
    
    diff = 0;
    node *nnaybor;
    
    nnaybor = n->outedge;
    
    // Reroot subtr at current place
    mfl_join_nodes(n, subtr->next->outedge->next);
    mfl_join_nodes(n->outedge, subtr->next->outedge->next->next);
    
    // Reoptimize the subtree
    //mfl_reopt_postorder(subtr->next->outedge, nchar);
    mfl_reopt_subtr_root(subtr->next->outedge, nchar);
    
    // Determine the cost of local reinsertion
    diff = mfl_subtr_reinsertion(subtr->next->outedge, up, dn, nchar);
    //printf("cost of reinserting at original place: %i\n", diff);
    
    mfl_regrafting_traversal(swapingon->trnodes[0], subtr, swapingon, savedtrees, ntax,
                             nchar, numnodes, current, 
                             undertreelimit, currentbesttree, 
                             foundbettertree, success, leftotry, diff);
    mfl_join_nodes(n, nnaybor);
    
    if (n->tip) {
        return;
    }
    
    /*mfl_reroot_subtree(n->next->outedge, subtr, up, dn, swapingon, 
                       savedtrees, ntax, nchar, numnodes, current, 
                       undertreelimit, currentbesttree, 
                       foundbettertree, success, leftotry, diff);
    mfl_reroot_subtree(n->next->next->outedge, subtr, up, dn, swapingon, 
                       savedtrees, ntax, nchar, numnodes, current, 
                       undertreelimit, currentbesttree, 
                       foundbettertree, success, leftotry, diff);*/
    
}






void mfl_bisection_traversal(node *n, tree *swapingon, tree **savedtrees, int ntax, 
                             int nchar, int numnodes, long int *current, bool *undertreelimit, 
                             int *currentbesttree, bool *foundbettertree, bool *success, long int *leftotry)
{
    node *p, *q, *clipnode, *up, *dn, *subtr;
    node *lft, *rt;
    int diff = 0;
    int trl = 0, cbest = 0;
    int *trlp = &trl, *cbestp = &cbest;
    
    if (n->start) {
        mfl_bisection_traversal(n->outedge, swapingon, savedtrees, ntax,
                              nchar, numnodes, current, 
                              undertreelimit, currentbesttree, foundbettertree, 
                              success, leftotry);
        return;
    }
    
    if (n->tip || *success) {
        return;
    }
    
    clipnode = n->next->next->outedge;
    
    p = n->next;
    while (p != n) {
        mfl_bisection_traversal(p->outedge, swapingon, savedtrees, ntax,
                              nchar, numnodes, current, 
                              undertreelimit, currentbesttree, foundbettertree, 
                              success, leftotry);
        if (*success) {
            return;
        }
        up = p->next->outedge;
        dn = p->next->next->outedge;
        
        subtr = p->next->next;
        
        if (!(up->tip && dn->tip)) {
            up->visited = 1;
            dn->visited = 1;
            clipnode->clip = true;
            mfl_join_nodes(up, dn);
            
            // Reoptimize the clipped tree
            mfl_trav_allviews(swapingon->trnodes[0], swapingon, ntax, nchar, trlp, cbestp);
            
            if (subtr->next->outedge->tip || (subtr->next->outedge->next->outedge->tip && subtr->next->outedge->next->next->outedge->tip)) {
                // Reoptimize the subtree
                //mfl_reopt_postorder(subtr->next->outedge, nchar);
                mfl_reopt_subtr_root(subtr->next->outedge, nchar);
                
                // Determine the cost of local reinsertion
                diff = mfl_subtr_reinsertion(subtr->next->outedge, up, dn, nchar);
                //printf("cost of reinserting at original place: %i\n", diff);
                
                //mfl_tip_reopt(swapingon, ntax, nchar);
                
                *leftotry = *leftotry - 1;
                mfl_regrafting_traversal(swapingon->trnodes[0], subtr, swapingon, 
                                         savedtrees, ntax, nchar, numnodes, current, 
                                         undertreelimit, currentbesttree, 
                                         foundbettertree, success, leftotry, diff);
            }
            else {
                q = subtr->next->outedge;
                while (!q->tip) {
                    q = q->next;
                    q = q->outedge;
                }
                
                //pop the root
                
                lft = subtr->next->outedge->next->outedge;
                rt = subtr->next->outedge->next->next->outedge;
                mfl_join_nodes(lft, rt);
                
                mfl_reroot_subtree(q->outedge, subtr, up, dn, swapingon, 
                                   savedtrees, ntax, nchar, numnodes, current, 
                                   undertreelimit, currentbesttree, 
                                   foundbettertree, success, leftotry, diff);
                
                //put old root back in
                mfl_join_nodes(lft, subtr->next->outedge->next);
                mfl_join_nodes(rt, subtr->next->outedge->next->next);
            }

            
            up->visited = 0;
            dn->visited = 0;
            clipnode->clip = false;
            mfl_devisit_tree(swapingon->trnodes, numnodes);
            if (*success) {
                return;
            }
            subtr->next->outedge->skip = false;
            mfl_join_nodes(up, p->next);
            mfl_join_nodes(dn, p->next->next);
        }
        p = p->next;
        clipnode = n->next;
    }   
}

void mfl_spr_search(int ntax, int nchar, int numnodes, charstate *tipdata, 
                    tree **savedtrees, int starttreelen)
{
    long int i = 0, j, remaining;
    long int trbufpos = 0;
    long int *nxtintrbuf = &trbufpos; 
    int currentbest = starttreelen, *currentbest_p = &currentbest;
    long int leftotry = 0;
    long int *leftotry_p = &leftotry;
    //long int numreps = 10;
    bool undertreelimit = true, *undertreelimit_p = &undertreelimit;
    bool foundbettertree = false, *foundbettertree_p = &foundbettertree;
    bool quit = false;
    bool success = false;
    bool *success_p = &success;
    
    
    /* A pointer to this whole struct will be used to replace all those 
     * params in the branch-swapping functions */
    mfl_searchrec *searchrec = mfl_create_searchrec();
    searchrec->nextinbuffer = 0;
    searchrec->undertreelimit = true;
    searchrec->bestlength = starttreelen;
    searchrec->foundbettertr = false;
    searchrec->niter_total = 0;
    searchrec->niter_ontree = 0;
    
    double timein = 0;
    double timeout = 0;
    
    *leftotry_p = 0; //mfl_spr_leftotry(ntax);
    
    mfl_devisit_tree(savedtrees[0]->trnodes, numnodes);
    
    *currentbest_p = mfl_all_views(savedtrees[0], ntax, nchar, currentbest_p);
    
    //printf("Length of starting tree again: %i\n", *currentbest_p);
    timein = (double)(clock() / (double)CLOCKS_PER_SEC);
    do {
        //printf("\nSwapping on tree %i\n", j);
        remaining = *nxtintrbuf + 1;
        *foundbettertree_p = false;
        *success_p = false;
        mfl_apply_tipdata(savedtrees[j], tipdata, ntax, nchar);
        mfl_all_views(savedtrees[j], ntax, nchar, currentbest_p);
        mfl_pruning_traversal(savedtrees[j]->trnodes[0], savedtrees[j], savedtrees, 
                              ntax, nchar, numnodes, nxtintrbuf, 
                              undertreelimit_p, currentbest_p, 
         foundbettertree_p, success_p, leftotry_p);
        /*mfl_bisection_traversal(savedtrees[j]->trnodes[0], savedtrees[j], savedtrees, 
                                ntax, nchar, numnodes, nxtintrbuf, 
                                undertreelimit_p, currentbest_p, 
                                foundbettertree_p, success_p, leftotry_p);*/
        if (*foundbettertree_p) {
            //*leftotry_p = 0;
            j = 0;
        }
        else {
            savedtrees[j]->swapped = true;
            //mfl_devisit_tree(savedtrees[j]->trnodes, numnodes);
            ++j;
        }
        
        if (j >= *nxtintrbuf) {
            printf("number of rearrangements tried: %li\n", *leftotry_p);
            quit = true;
        }
        
    } while (!quit);
    
    timeout = (double)(clock() / (double)CLOCKS_PER_SEC);
    
    printf("Total search time: %g\n", timeout - timein);
    
    printf("Next in trbuf: %li\n", *nxtintrbuf);
    
    printf("\nThe optimal tree(s) found by subtree pruning and regrafting:\n");
    for (i = 0; i < *nxtintrbuf; ++i) {
        //printf("Tree %li:\n", i + 1);
        printf("TREE str_%li = [&U] ", i+1);
        mfl_root_tree(savedtrees[i], 0, ntax);
        printNewick(savedtrees[i]->root);
        printf(";\n");
        //printf("Length: %i\n", savedtrees[i]->length);
    }
    printf("\n");
    
    mfl_clear_treebuffer(savedtrees, nxtintrbuf, numnodes);
    free(savedtrees);
    mfl_destroy_searchrec(searchrec);
}

/*mfl_heuristic(mfl_handle_t *mfl_handle)
{
    // Setup heuristic search and call appropriate branch-swappers
}*/