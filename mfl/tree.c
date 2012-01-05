/*
 *  tree.c
 *  Morphy
 *
 *  Miscellaneous functions for trees
 *
 */

#include "morphy.h"
extern int numnodes;

struct node * seekInternal(int ntax, node **nds)
{
    /* Searches for an unused internal node */
    /* NB: This function needs to return some kind of error msg
     * if no unused nodes are found. */
    
    int i;
    node *unused = NULL;
    node *p;
    
    for (i = ntax + 1; i < numnodes; ++i) {
        if (!nds[i]->next) {
            unused = nds[i];
            i = numnodes;
        }
    }
    
    if (!unused) {
        for (i = ntax + 1; i < numnodes; ++i) 
        {
            p = nds[i]->next;
            
            while (p != nds[i]) 
            {
                if (!p->outedge) 
                {
                    unused = p;
                    p = nds[i];
                }
                else 
                {
                    p = p->next;
                }                
            }
        }
    }
    
    if (!unused) {
        printf("Error in tree memory allocation\n");
        return unused;
    }
    else {
        return unused;
    }
}

void closeRing(node *n)
{
    /* Makes sure there isn't a dangling next pointer*/
    node *p;
    
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

int determOrder(node *n)
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
    while (p != n) {
        if (p->outedge) {
            ++i;
        }
        p = p->next;
    }
    
    return i;
    
}

void setOrder(node *n)
{
    int ord;
    node *p;
    ord = determOrder(n);
    
    n->order = ord;
    if (n->next) {
        p = n->next;
        while (p != n) {
            p->order = ord;
            p = p->next;
        }
    }
}

void clearOrder(node *n)
{
    node *p;
    
    p = n->next;
    while (p != n) {
        p->order = 0;
        p = p->next;
    }
    
    n->order = 0;
}

void setIndex(node *n)
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
    setOrder(rnode);
}

void insertBranch(node *br, node *target)
{
    /* Inserts a branch with a ring base into another branch*/
    
    node *p, *bout, *tdesc;
    
    tdesc = target->outedge;
    
    // Find an available node in the ring 
    p = br->next;
    while (p != br) {
        if (!p->outedge) {
            bout = p;
            p = br;
        }
        else {
            p = p->next;
        }
    }
    
    joinNodes(br, target);
    joinNodes(bout, tdesc);
    
}

void resolve(node *n, node **nds, int ntax)
{
    /* Arbitrarily resolves a non-binary node and leaves it as binary*/
    
    int i, ord;
    node *p, *in;
    
    // Make sure node there is a polytomy, otherwise exit resolve()
    ord = determOrder(n);
    n->order = ord;
    if (ord <= 3) {
        return;
    }
    
    // Find available internal node(s) in nds
    in = seekInternal(ntax, nds);
    asNoring(in);
    
    // Randomly select branches to join to it
    
    
    
}