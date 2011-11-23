/*
 *  random.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 11-07-25.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 *	This file contains code for generating random rooted (pre-specified) and unrooted trees
 *
 */

#include "morphy.h"


extern int ntax;
long long int numnodes;

/*bool isodd (float seedn)
 {
 if (seedn % 2 == 0)
 return false;
 else 
 return true;
 
 }*/

/* Taxon shuffle using array shuffle by Ben Pfaff http://benpfaff.org/writings/clc/shuffle.html 
 * ****NOTE******: Will be modified to employ the GSL random number generator*/

void shuffle(int *taxarray, int length)
{
	int i, j, t;
	
	
    if (length > 1) {
		for (i = 0; i < length - 1; i++) {
			j = i + rand() / (RAND_MAX / (length - i) + 1);
			t = taxarray[j];
			taxarray[j] = taxarray[i];
			taxarray[i] = t;
		}
    }
}

struct tree *randtrunk(tree *newtrunk, node *startn)
{
	/* Generates the 'trunk' of the tree which consists of the network of all
	 * internal nodes in the tree */
	
	
	
	return (newtrunk);
}

struct tree *randrooted (tree *randtree)
{
	/* Returns a random tree with an arbitrary root*/
	
	randtree = randunrooted(randtree);
	randtree->trnodes[0]->start = false;
	rootOnTerminal(randtree, 1);	// The second argument should eventually be replaced by a randomly drawn number between 0 and ntax-1
	
	return (randtree);
}

struct tree *rand_w_root (int root)
{
	/* Returns a random ingroup topology on a given root
	 * and will arbitrarily resolve the outgroup*/
	return;
	
}

struct tree *randunrooted (tree *randtree)
{
	/* Returns a random unrooted tree*/
	
	int i;
	int *taxarray;
	node *p, *q;
	
	taxarray = malloc(ntax * sizeof(int));
	init_taxarray(taxarray);
	
	shuffle(taxarray, ntax);
	
	randtree = alloctree(randtree);
	
	randtree->trnodes[0]->start = true;
	randtree->trnodes[0]->outedge = randtree->trnodes[taxarray[0]];
	
	// Join all the internal nodes (except the root) together
	for (i = 1; i <= (ntax - 2); ++i) {
		p = randtree->trnodes[ntax + i]->next->next;
		q = randtree->trnodes[ntax + i + 1];
		p->outedge = q;
		q->outedge = p;		
	}
	
	// Add all the tips to the appropriate internal nodes
	randtree->trnodes[taxarray[0] - 1]->outedge = randtree->trnodes[ntax + 1];
	randtree->trnodes[ntax + 1]->outedge = randtree->trnodes[taxarray[0] - 1];
	
	for (i = 1; i < ntax - 1; ++i) {
		randtree->trnodes[taxarray[i] - 1]->outedge = randtree->trnodes[ntax + i]->next;
		randtree->trnodes[ntax + i]->next->outedge = randtree->trnodes[taxarray[i] - 1];
	} 
	
	randtree->trnodes[2 * ntax - 2]->next->next->outedge = randtree->trnodes[taxarray[i] - 1];
	randtree->trnodes[taxarray[i] - 1]->outedge = randtree->trnodes[2 * ntax - 2]->next->next;
	
	free(taxarray);
	
	printNewick(randtree->trnodes[0]);
	printf("\n");
	
	return (randtree);
}