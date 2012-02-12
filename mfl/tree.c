/*
 *  tree.c
 *  Morphy
 *
 *  Miscellaneous functions for trees
 *
 */

#include "morphy.h"

struct node * mfl_seek_internal(int ntax, int numnodes, node **nds)
{
    /* Searches for an unused internal node */
    /* NB: This function needs to return some kind of error msg
     * if no unused nodes are found. */
    
    int i;
    node *unused = NULL;
    node *p;
    bool isUsed = false;
    
    for (i = ntax + 1; i < numnodes; ++i) {
        if (!nds[i]->next /*&& !nds[i]->initialized*/ && !nds[i]->outedge) {
            unused = nds[i];
            i = 2 * ntax;
        }
    }
    
    if (!unused) 
    {
        for (i = ntax + 1; i < numnodes; ++i) 
        {
            isUsed = false;
            if (nds[i]->next) 
            {
                p = nds[i];                
                do
                {
                    if (p->outedge) 
                    {
                        isUsed = true;
                        p = nds[i];
                    }
                    else 
                    {
                        p = p->next;
                    }                
                } while (p != nds[i]);
                
                if (!isUsed) {
                    unused = nds[i];
                }
            }
        }
    }
    
    if (!unused) {
        printf("Error in tree memory allocation\n");
        return unused;
    }
    else {
        //unused->initialized = 1;
        return unused;
    }
}

struct node * mfl_seek_ringnode(node *n, int ntax)
{
    /*-Used by copytree-*/
    node *p;
    bool rootnode = false;
    
    if (n->index == ntax) {
        rootnode = true;
    }
    
    if (!rootnode && !n->outedge) {
        return n;
    } else {
        p = n->next;
        while (p != n) {
            if (!p->outedge) {
                return p;
            }
            p = p->next;
        }
    }
    
    printf("did not find an available node in ring\n");
    return n;
}

void mfl_set_vweight(node *n)
{
    node *p;
    p = n->next;
    while (p != n) {
        n->vweight = n->vweight + p->outedge->vweight;
    }
}

void mfl_close_ring(node *n)
{
    /* Makes sure there isn't a dangling next pointer*/
    node *p;
    
    p = n;
    
    do {
        if (p->next) {
            p = p->next;
        }
        if (!p->next) {
            printf("Report error: dangling next pointer at %p\n", p);
            p->next = n;
            p = p->next; 
        }
    }
    while (p != n);
}

void mfl_as_ring(node *n, int ntax)
{
    
    if (n->next) {
        mfl_close_ring(n);
        deletering(n);
    }
    
    newring(n, ntax);
}

void mfl_as_noring(node *n)
{
    if (n->next) {
        mfl_close_ring(n);
        deletering(n);
    }
}

void mfl_reindex_tree(nodearray nds, int ntax, int numnodes)
{
    int i;
    
    for (i = 0; i < numnodes; ++i)
    {
        /* First half of trnodes are initialized to this */
        if (i < ntax)
        {
            nds[i]->index = i;
        }
        else /* second half are allocated as a newring */
        {
            nds[i]->index = i;
            mfl_set_index(nds[i]);
        }        
    }
}

void mfl_set_ring_to_n(node *n)
{
    node *p;
    
    if (n->next) {
        p = n->next;
        while (p != n) 
        {
            p->index = n->index;
            p->apomorphies = n->apomorphies;
            p->initialized = n->initialized;
            p->skip = n->skip;
            p = p->next;
        }
    }
    mfl_set_order(n);
    
}

void mfl_reset_ring(node *n)
{
    node *p;
    
    n->initialized = 0;
    n->apomorphies = 0;
    n->skip = 0;
    
    if (n->next) {
        p = n->next;
        while (p != n) 
        {
            p->index = n->index;
            p->apomorphies = n->apomorphies;
            p->initialized = n->initialized;
            p->skip = n->skip;
            p = p->next;
        }
    }
    mfl_set_order(n);    
}

void mfl_collapse(node *n1, nodearray nds)
{
    node *an1, *an2, *n2, *tmp;
    
    an1 = n1->outedge;
    an2 = an1->next;
    tmp = an2;
    while (an2->next != an1) {
        an2 = an2->next;
    }
    an2->next = n1->next;
    
    n2 = n1->next;
    while (n2->next != n1) {
        n2 = n2->next;
    }
    n2->next = an1->next;
    an1->next = NULL;
    n1->next = NULL;
    n1->outedge = NULL;
    an1->outedge = NULL;
    free(an1);
    nds[n1->index] = n1;
    mfl_set_order(tmp);
    mfl_set_vweight(tmp);
}

void mfl_arb_resolve(node *n, node **nds, int ntax, int numnodes)
{
    /* Arbitrarily resolves a non-binary node and leaves it binary*/
    
    int c, i, j, ord;
    node *p, *q, *in, *in2;
    
    // Make sure node there is a polytomy, otherwise exit resolve()
    mfl_set_order(n);
    ord = n->order;
    if (ord <= 3 && n->index != ntax) {
        printf("mfl_arb_resolve() called in error on node %i\n", n->index);
        return;
    }
    
    // Find available internal node(s) in nds
    in = mfl_seek_internal(ntax, numnodes, nds);
    in2 = in;
    if (in->next) {
        mfl_as_noring(in);
    }
    if (in->outedge) {
        if (in->outedge->outedge) {
            in->outedge->outedge = NULL;
        }
        in->outedge = NULL;
    }
    
    /* Randomly select branches to join to it. Cycles through the ring (skipping 
     * n itself) ntax % 10 times. This doesn't use a true random number 
     * generator, but should be both sufficiently arbitrary with respect to the
     * tree's topology, but sufficiently deterministic to be repeatable. */
    
    c = ntax % 10;
    
    for (i = 0; i <= (ord - 3); ++i)
    {
        p = n->next;
        for (j = 0; j <= c; ++j) 
        {
            q = p;
            p = p->next;
            if (p == n) 
            {
                q = p;
                p = p->next;
            }
        }
        in2->next = p;
        in2 = in2->next;
        q->next = p->next;
        p->next = NULL;
    }
    
    p->next = allocnode();
    if (p->next->outedge) {
        p->next->outedge = NULL;
    }
    p = p->next;
    p->next = in;
    
    if (ntax % 2) 
    {
        n = n->next;
    }
    else 
    {
        n = n->next->next;
    }
    
    mfl_insert_branch(in, n, ntax);
}

int mfl_determ_order(node *n)
{
    /*determines the number of branchings in a node*/
    int i = 0;
    node *p;
    
    if (n->outedge) {
        i = 1;
    }
    
    if (n->tip || !n->next) {
        return i;
    }
    
    p = n->next;
    while (p != n) 
    {
        if (p->outedge) {
            ++i;
        }
        p = p->next;
    }
    
    return i;
}

void mfl_set_order(node *n)
{
    int ord;
    node *p;
    ord = mfl_determ_order(n);
    
    n->order = ord;
    if (n->next) 
    {
        p = n->next;
        while (p != n) {
            p->order = ord;
            p = p->next;
        }
    }
}

void mfl_clear_order(node *n)
{
    node *p;
    
    p = n->next;
    while (p != n) {
        p->order = 0;
        p = p->next;
    }
    
    n->order = 0;
}

void mfl_set_index(node *n)
{
    node *p;
    
    //printf("setting to index %i\n", n->index);
    
    if (n->next) {
        p = n->next;
        while (p != n) 
        {
            p->index = n->index;
            p = p->next;
        }
    }   
}

void mfl_deinit_tree(tree *t, int numnodes)
{
    int i;
    node *p;
    
    for (i = 0; i < numnodes; ++i) 
    {
        t->trnodes[i]->initialized = 0;
        
        if (t->trnodes[i]->next) 
        {
            p = t->trnodes[i]->next;
            
            do {
                p->initialized = 0;
                if (p->next) 
                {
                    p = p->next;
                }
                else 
                {
                    p = t->trnodes[i]; /* Should prevent a crash in the event of a dangling next pointer*/
                }

            } while (p != t->trnodes[i]);
            
        }
    }
}

void mfl_put_branch_in_ring(node *n, node *rnode)
{
    /* Given a branch (two nodes joined by their outedge pointers), places the 
     * branch in a ring */
    
    node *rnode2 = NULL;
    
    if (rnode->next) {
        rnode2 = rnode->next;
    }
    
    rnode->next = n;
    
    if (rnode2) {
        n->next = rnode2;
    }
    else {
        n->next = rnode;
    }
    
    rnode->order = rnode->order + 1;
    mfl_set_order(rnode);
}


void mfl_temproot(tree *trtoroot, int root, int ntax)
{
    node *lftbr, *rtbr;
    
    lftbr = trtoroot->trnodes[root];
    rtbr = lftbr->outedge;
    
    if (!trtoroot->trnodes[ntax]->next) {
        newring(trtoroot->trnodes[ntax], ntax);
    }
    
    joinNodes(lftbr, trtoroot->trnodes[ntax]->next);
    joinNodes(rtbr, trtoroot->trnodes[ntax]->next->next);
    
    trtoroot->root = trtoroot->trnodes[ntax];
    trtoroot->trnodes[0]->start = false;
}

void mfl_undo_temproot(int ntax, tree *trtounroot)
{
    node *lftbr, *rtbr;
    
    lftbr = trtounroot->trnodes[ntax]->next->outedge;
    rtbr = trtounroot->trnodes[ntax]->next->next->outedge;
    
    joinNodes(lftbr, rtbr);
    
    trtounroot->root = NULL;
    trtounroot->trnodes[0]->start = true;
}


int mfl_tree_enumerator(void)
{
    static long long int treenum = 0;
    
    return ++treenum;
}

void mfl_resize_treebuffer(tree **treebuffer, int *treelimit, int sizeincrease)
{
    tree **newtreebuffer;
    newtreebuffer = (tree **)realloc(treebuffer, (*treelimit + sizeincrease) * sizeof(tree*));
    if (!newtreebuffer) {
        printf("Insufficient memory for new treelimit\n");
        newtreebuffer = treebuffer;
    }
    else {
        *treelimit = *treelimit + sizeincrease;
        printf("Treelimit changed to %i\n", *treelimit);
    }
}

void mfl_clear_treebuffer(tree **treebuffer, long int *numsavedtrees, int numnodes)
{
    int i;
    for (i = 0; i < *numsavedtrees; ++i) {
        if (treebuffer[i]) {
            freetree(treebuffer[i], numnodes);
        }
    }
    *numsavedtrees = 0;
}

void mfl_reinit_treebuffer(tree **treebuffer, tree *newbest, long int *numsavedtrees, int numnodes)
{    
    printf("Found shorter tree. Reinitializing the treebuffer\n");
    printf("Best tree length: %i\n", newbest->length);
    int i;
    for (i = 0; i < *numsavedtrees; ++i) {
        if (treebuffer[i] && treebuffer[i] != newbest) {
            freetree(treebuffer[i], numnodes);
        }
        if (treebuffer[i] == newbest) {
            treebuffer[i] = NULL;
        }
    }
    
    newbest->index = 0;
    newbest->swapped = false;
    treebuffer[0] = newbest;
    *numsavedtrees = 0;
}