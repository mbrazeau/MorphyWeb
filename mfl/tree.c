/*
 *  tree.c
 *  Morphy
 *
 *  Miscellaneous functions for trees
 *
 */

#include "morphy.h"

struct node * seekInternal(int ntax, node **nds)
{
    /* Searches for an unused internal node */
    /* NB: This function needs to return some kind of error msg
     * if no unused nodes are found. */
    
    int i;
    node *unused = NULL;
    node *p;
    bool isUsed = false;
    
    for (i = ntax + 1; nds[i]; ++i) {
        if (!nds[i]->next && !nds[i]->initialized && !nds[i]->outedge) {
            unused = nds[i];
            i = 2 * ntax;
        }
    }
    
    if (!unused) 
    {
        for (i = ntax + 1; nds[i]; ++i) 
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
        unused->initialized = 1;
        return unused;
    }
}

void closeRing(node *n)
{
    /* Makes sure there isn't a dangling next pointer*/
    node *p;
    
    p = n;
    
    do {
        if (p->next) {
            p = p->next;
        }
        if (!p->next) {
            p->next = n;
            p = p->next; 
        }
    }
    while (p != n);
}

void asRing(node *n)
{
    
    if (n->next) {
        closeRing(n);
        deletering(n);
    }
    
    newring(n);
}

void asNoring(node *n)
{
    if (n->next) {
        closeRing(n);
        deletering(n);
    }
}

void collapse(node *n)
{
    /*collapses a branch*/
    
    node *an1, *an2, *b1, *b2;
    node *desc1;
    node *p;
    
    // Set a pointer to the ancestor of n
    an1 = n->outedge;
    // Set a pointer to next node in ring from ancestor of n
    an2 = an1->next;
    
    // Set a pointer to the left-most descendant of n's ring
    desc1 = n->next->outedge;
    
    // Find the segment of the ring that contains the descendant branches
    b1 = n->next->next;
    p = b1;
    while (p->next != n) {
        p = p->next;
    }
    b2 = p;
    
    n->next->next = n;
    n->outedge = NULL;
    n->next->outedge = NULL;
    
    joinNodes(desc1, an1);
    
    an1->next = b1;
    b2->next = an2;
    
}

int mfl_determ_order(node *n)
{
    /*determines the number of branchings in a node*/
    int i = 0;
    node *p;
    
    if (n->outedge) {
        i = 1;
    }
    
    if (n->tip) {
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
    printf("node order: %i\n", i);
    
}

void mfl_set_order(node *n)
{
    int ord;
    node *p;
    ord = mfl_determ_order(n);
    
    n->order = ord;
    if (n->next) {
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
    
    if (n->next) {
        p = n->next;
        while (p != n) {
            p->index = n->index;
            p = p->next;
        }
    }   
}

void deinit_tree(tree *t)
{
    int i;
    node *p;
    
    for (i = 0; t->trnodes[i]; ++i) 
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

void putBranchInRing(node *n, node *rnode)
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

void resolve(node *n, node **nds, int ntax, int numnodes)
{
    /* Arbitrarily resolves a non-binary node and leaves it as binary*/
    
    int i, ord;
    node *p, *in;
    
    // Make sure node there is a polytomy, otherwise exit resolve()
    ord = mfl_determ_order(n);
    n->order = ord;
    if (ord <= 3) {
        return;
    }
    
    // Find available internal node(s) in nds
    in = seekInternal(ntax, nds);
    asNoring(in);
    
    // Randomly select branches to join to it
    
    
    
}