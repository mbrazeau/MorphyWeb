/*
 *  readnewick.c
 *  Morphy
 *
 *  Functions for reading a tree as a Newick-format string and converting it 
 *  into a tree struct.
 *
 */

#include <ctype.h>
#include "morphy.h"

struct node * cpyfromNWK(char *nwktr, int nwklen, int ntax, int numnodes, int *pos, nodearray nds, bool isRooted)
{
    
    int i, tipnum;
    char tipbuf[10];
    node *n, *nlst, *p;
    
    if (isRooted) 
    {
        n = nds[ntax];
        nlst = n;
        isRooted = false;
    }
    else 
    {
        n = seekInternal(ntax, numnodes, nds);
        nlst = n;
    }
    
    
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

void cpRootedNWK(char *nwktr, int nwklen, int ntax, int numnodes, tree *newtree, bool isRooted)
{
    int ctr = 0;
    int *ctrp;
    
    ctrp = &ctr;
    
    newtree->trnodes[ntax] = cpyfromNWK(nwktr, nwklen, ntax, numnodes, ctrp, newtree->trnodes, isRooted);
    
    newtree->root = newtree->trnodes[ntax];
    
    printNewick(newtree->root);
    printf("\n");
    
}

struct tree * readNWK (char *nwktr, bool isRooted)
{

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
			++numnodes_a;
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
    deinit_tree(newtree);
	
    i = 0;
    position = &i;
    
    cpRootedNWK(nwktr, nwklen, numtaxa, numnodes_a, newtree, isRooted);
    
    return newtree;
}


void testNWKreading(void)
{
    bool isRooted = true;
    
    char *nwktree;
    
    char newickTree[] = "(2,((1,(5,3),4),6));";  
    printf("The newick string: %s;\n", newickTree);
    
    nwktree = newickTree;
    
    readNWK(nwktree, isRooted);
    
}