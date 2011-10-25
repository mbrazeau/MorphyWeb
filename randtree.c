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
	
	int i;
	int *taxarray;
	node *p, *q;
	
	taxarray = malloc(ntax * sizeof(int));
	init_taxarray(taxarray);
	
	/*for (i = 0; i < ntax; ++i) {
		printf("%i ", taxarray[i]);
	}
	printf("\n");*/
	
	shuffle(taxarray, ntax);
	
	/*for (i = 0; i < ntax; ++i) {
		printf("%i ", taxarray[i]);
	}
	printf("\n");*/
	
	randtree = alloctree(randtree);
	
	randtree->root = randtree->trnodes[ntax]; // Set the pointer to the root.
	
	for (i = 0; i <= (ntax - 2); ++i) {
		p = randtree->trnodes[ntax + i]->next->next;
		q = randtree->trnodes[ntax + i + 1];
		p->outedge = q;
		q->outedge = p;		
	}
	
	for (i = 0; i < ntax - 1; ++i) {
		randtree->trnodes[taxarray[i] - 1]->outedge = randtree->trnodes[ntax + i]->next;
		randtree->trnodes[ntax + i]->next->outedge = randtree->trnodes[taxarray[i] - 1];
	} 
	
	randtree->trnodes[2 * ntax - 2]->next->next->outedge = randtree->trnodes[taxarray[i] - 1];
	randtree->trnodes[taxarray[i] - 1]->outedge = randtree->trnodes[2 * ntax - 2]->next->next;
	
	//printNewick(randtree->root);
	//printf(";\n\n");
	
	free(taxarray);
	
	return (randtree);
}

struct tree *rand_w_root (void)
{
	/* Returns a random ingroup topology on a given root
	 * and will arbitrarily resolve the outgroup*/
	return;
	
}

struct tree *randunrooted (void)
{
	/* Returns a random unrooted tree*/
	
	return;
}