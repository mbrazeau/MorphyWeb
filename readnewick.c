/*
 *  readnewick.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 10/14/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


void nwk_pass(char *nwk_string, void *any_funct)
{
	/*passes over the newick format tree until it reaches the terminal semicolon and *
	 *calls any function pointed to by void (*funct). */
	
	void (*funct)(char *);
	
	do {
		
		
		++i;
	} while (nwk_string[i] != ';');
	
}

void nwk_one_rootd (char *nwk_string)
{

	int i;
	int opencount = 0, closecount = 0;	// number of open and closed brackets respectively
	int numnodes_a = 0, numnodes_p = 0, numtaxa = 0; //numnodes_a is the actual number of nodes used, numnodes_p is the number possible
	tree *newtree;

	do {
		if (nwk_string[i] == '(') {
			++opencount;
			++numnodes;
		}
		
		if (nwk_string[i] == ',')
			++numtaxa;
		
		if (nwk_string[i] == ')')
			++closecount;
		
		++i;
		
	} while (nwk_string[i] != ';');

	numtaxa = numtaxa + 1;

	if (opencount > closecount) {
		printf("Error in newick format: '(' outnumber ')'.\n");
		return;
	}
	if (opencount < closecount) {
		printf("Error in newick format: ')' outnumber '('.\n");
		return;
	}
	
	numnodes_p = 2 * numtaxa - 1;
	
	newtree = malloc(sizeof(tree));
	if (newtree == NULL) {
		printf("Error in readnewick: Failed to allocate memory for tree.\n");
		return;
	}
	
	newtree.trnodes = malloc(numnodes_p * sizeof(node));
	
	for (i = 0; i < numnodes_p; ++i) {
		newtree.trnodes[i] = malloc(sizeof(node));
	}
	
	
	
}