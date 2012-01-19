/*
 *  readnewick.c
 *  Morphy
 *
 *  Functions for reading a tree as a Newick-format string and converting it 
 *  into a tree struct.
 *
 */

#include "morphy.h"

struct node * cpyfromNWK(char *nwktr, int nwklen, int ntax, int numnodes, 
                         int *pos, nodearray nds, bool isRooted)
{
    
    int i, tipnum;
    char tipbuf[10];
    node *n, *nlst, *p;
    
    n = mfl_seek_internal(ntax - 1, numnodes, nds);
    n->initialized = 1;
    nlst = n;
    
    do {
        
        *pos = *pos + 1;
        
        if (nwktr[*pos] == ')')
        {
            nlst->next = n;
            return n;
        }
        
        if (isdigit(nwktr[*pos]))
        {
            for (i = 0; nwktr[*pos] != ',' && nwktr[*pos] != ')'; ++i, ++*pos) 
            {
                tipbuf[i] = nwktr[*pos];
            }
            
            tipbuf[i] = '\0';
        
            tipnum = atoi(tipbuf);
            
            nlst->next = allocnode();
            nlst = nlst->next;
            
            joinNodes(nds[tipnum - 1], nlst);
            
            if (nwktr[*pos] == ')')
            {
                nlst->next = n;
                return n;
            }
            
        }
        
        if (nwktr[*pos] == '(') 
        {
            nlst->next = allocnode();
            nlst = nlst->next;
            p = cpyfromNWK(nwktr, nwklen, ntax, numnodes, pos, nds, isRooted);
            joinNodes(p, nlst);
        }
        
    } while (nwktr[*pos]);
    
    nlst->next = n;
    return n;
}

void NWK_roothandl(char *nwktr, int nwklen, int ntax, int numnodes, 
                   tree *newtree, bool isRooted)
{
    int ctr = 0;
    int *ctrp;
    
    ctrp = &ctr;
    
    newtree->trnodes[ntax] = cpyfromNWK(nwktr, nwklen, ntax, numnodes, ctrp, 
                                        newtree->trnodes, isRooted);
    
    newtree->root = newtree->trnodes[ntax];

    if (!isRooted) {
        unroot(ntax, newtree);
    }
    
    /* debugging print */
    printf("[&");
    if (newtree->root) 
    {
        printf("R] = ");
        printNewick(newtree->root);
    }
    else 
    {
        printf("U] = ");
        printNewick(newtree->trnodes[0]);
    }
    printf("\n");
    /* end debugging print */    

}

struct tree * readNWK (char *nwktr, bool isRooted)
{

    /* Returns a pointer to a tree based on a Newick-format string. The Newick 
     * format does not have to have an arbitrary root to be read. If there is an
     * error in reading the Newick tree, a NULL pointer is returned. */
    
	int i = 0;
    int nwklen;
	int opencount = 0, closecount = 0;	// number of open and closed brackets 
                                        // respectively
	int numnodes_a = 0, numtaxa = 0; // numnodes_a is the actual number of nodes 
                                     // used, numnodes_p is the number possible
	int *position; // Position counter for the nwk string in cpyfromNWK
    bool intaxname = false;
    tree *newtree;

    nwklen = strlen(nwktr) - 1; // Minus 1 for the semicolon ending the string
    
	do {
		if (nwktr[i] == '(') {
			++opencount;
            intaxname = false;
        }
		
		if (nwktr[i] == ',') {
            intaxname = false;
        }
		
		if (nwktr[i] == ')') {
			++closecount;
            intaxname = false;
        }
        
        if (isalnum(nwktr[i]) && !intaxname) {
            ++numtaxa;
            intaxname = true;
        }
		
		++i;
		
	} while (nwktr[i] != ';');

	if (opencount > closecount) 
    {
		printf("Error in Newick format tree: '(' outnumber ')'.\n");
        newtree = NULL;
		return newtree;
	}
	if (opencount < closecount) 
    {
		printf("Error in Newick format tree: ')' outnumber '('.\n");
		newtree = NULL;
		return newtree;
	}
	
	numnodes_a = 2 * numtaxa - 1;
	
	newtree = alloc_noring(numtaxa, numnodes_a);
    mfl_deinit_tree(newtree, numnodes_a);
	
    i = 0;
    position = &i;
    
    NWK_roothandl(nwktr, nwklen, numtaxa, numnodes_a, newtree, isRooted);
    
    return newtree;
}