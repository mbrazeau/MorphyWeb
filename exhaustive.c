/*
 *  exhaustive.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 11-07-02.
 *  Copyright 2011. All rights reserved.
 *
 */

#include "morphy.h"
#include <math.h>


extern tree *treeSet;
extern int ntax;

long long int factorial(long long int n)
{
	long long int result;
	
	if (n == 0)
		result = 1;
	else 
		result = n * factorial(n - 1);
	return result;
}

long long int numrooted(int ntaxa)
{
	long long int numerator, denominator, treecombs;
	
	if (ntaxa < 3) {
		printf("Error in ntax: too few taxa\b");		// Should have output to stderr, too
		return 0;
		
	} else {	
		numerator = 2 * ntaxa - 3;
		numerator = factorial(numerator);
		denominator = (long long int)pow(2, ntaxa - 2);		// cast as type long long because pow returns type double
		denominator = denominator * factorial(ntaxa - 2);
		
		treecombs = numerator / denominator;
		return treecombs;
	}
}

long long int nunrooted(int ntaxa)
{
	long long int numerator, denominator, treecombs;
	
	if (ntaxa < 3) {
		printf("Error in ntax: too few taxa\b");		// Should have output to stderr, too
		return 0;
		
	} else {	
		numerator = 2 * ntaxa - 5;
		numerator = factorial(numerator);
		denominator = (long long int)pow(2, ntaxa - 3);
		denominator = denominator * factorial(ntaxa - 3);
		
		treecombs = numerator / denominator;
		return treecombs;
	}
}


int posits_unr(int i)
{
	/* Calculates the number of positions for the ith taxon
	 * in an unrooted tree*/
	
	int result;
	
	result = 2 * i - 5;

	return result;
}

void devisit(node *n)
{
	node *p;
	
	if (n->start) {
		n->visited = 0;
		devisit(n->outedge);
		return;
	}	
	
	if (n->tip && n->outedge->next->outedge->tip && !n->outedge->next->outedge->start) {
		n->visited = 0;
		return;
	}
	
	if (n->tip && !n->start) {
		n->visited = 0;
		return;
	}
	
	n->visited = 0;
	
	p = n->next;
	while (p != n) {
		devisit(p->outedge);
		p = p->next;
		p->visited = 0;
	}	
}

void newbranch(node *desc, node *ancest, tree *basetree, int taxon, int calln)
{	
	/*NB: calln is here replacing an earlier version which used *counter */
	
	node *newtip, *newn;
	
	//basetree->trnodes[taxon - 1] = malloc(sizeof(node));
	//basetree->trnodes[taxon] = malloc(sizeof(node));
	
	newtip = basetree->trnodes[taxon - 1];
	newtip->apomorphies = malloc(3 * sizeof(char));
	if (newtip->apomorphies == NULL) {
		printf("In newbranch(): failed to allocate memory for apomorphies in newtip\n");
		return;
	}
	basetree->trnodes[ntax + calln + 1] = (node *) malloc(sizeof(node));
	if (basetree->trnodes[ntax + calln + 1] == NULL) {
		printf("In newbranch(): failed to allocate memory for new branch node in basetree\n");
		return;
	}
	newn = basetree->trnodes[ntax + calln + 1];
	
	newring(newn);
	
	newn->visited = 1;
	
	ancest->outedge = newn;
	newn->outedge = ancest;
	
	newn->next->outedge = desc;
	desc->outedge = newn->next;
	
	newn->next->next->outedge = newtip;
	newtip->outedge = newn->next->next;
	
	newtip->tip = taxon;
	//*taxon = *taxon + 1;
}

void insert_allp(node *n, tree *origtree, int taxon, int calln, int *counter)
{
	/* Adds a terminal in its next possible location in the sequence. This function is
	 * called iteratively on copies of a starting tree from within allunrooted. 
	 * and keeps track of the number of recursive calls made on itself through a 
	 * pointer, counter. It then executes its traversal until the counter reaches the 
	 * call number, at which point it adds a branch. */
	
	node *p;
	
	if (*counter == calln && !n->outedge->start) {
		newbranch(n, n->outedge, origtree, taxon, calln);
		return;
	}

	if (n->tip && !n->start) {
		return;
	}

	if (n->start) {
		n = n->outedge;
	}
	
	p = n->next;
	
	while (p != n) {
		*counter = *counter + 1;
		if (*counter <= calln) {
			insert_allp(p->outedge, origtree, taxon, calln, counter);
			p = p->next;
		} else {
			return;
		}
	}
}	


void allunrooted(/*tree *treearray, int ntax*/)
{
	
	/* Sets up the only unrooted tertiary tree for three terminal nodes (in this case: 1, 2, and 3)
	 * */
	
	long int i, j, positions, prevpos, count, taxon;
	int *counter = &count;
	long long int ntrees, copies;
	long long int *cpyptr = &copies;
	tree *treearray;
	tree *origtree, *treecopy;
	node *ring1;
	
	printf("Compute all unrooted trees\n");	
	
	ntrees = nunrooted(ntax);
	printf("Number of rooted trees for %i taxa: %g\n", ntax, (double)ntrees);
	

	treearray = (tree*) malloc((ntrees + 1) * sizeof(tree));
	if (treearray == NULL) {
		printf("tree malloc unssuccessful.\n");
	}

	
	for (i = 0; i < ntrees; ++i) {
		treearray[i].trnodes = (nodearray) malloc((2 * ntax - 1) * sizeof(node));
		if (treearray[i].trnodes == NULL) {
			printf("In allunrooted(): failed to allocate memory for trnodes of treearray.\n");
		}
	}
	
	printf("node arrays malloc'd: %li\n", (long)i);
	
	// Create the first unrooted tertiary tree of first three taxa.
	
	for (i = 0; i < ntax; ++i) {
		treearray[0].trnodes[i] = (node *) malloc(sizeof(node));
		if (treearray[0].trnodes[i] == NULL) {
			printf("In allunrooted(): failed to allocate memory for %li node of base tree.\n", (long)(i + 1));
		}
		treearray[0].trnodes[i]->tip = (i + 1);
		treearray[0].trnodes[i]->visited = 0;
	}
	
	treearray[0].trnodes[ntax] = (node *) malloc(sizeof(node)); // Reserved for the eventual root
	if (treearray[0].trnodes[ntax] == NULL) {
		printf("In allunrooted(): unable to allocate memory for root node.\n");
	}
	treearray[0].trnodes[ntax + 1] = (node *) malloc(sizeof(node)); // The node that will join the three taxa
	if (treearray[0].trnodes[ntax + 1] == NULL) {
		printf("In allunrooted(): unable to allocate memory for initial node.\n");
	}

	ring1 = treearray[0].trnodes[ntax + 1];
	
	newring(ring1);
	
	ring1->outedge = treearray[0].trnodes[0];
	ring1->next->outedge = treearray[0].trnodes[1];
	ring1->next->next->outedge = treearray[0].trnodes[2];
	
	treearray[0].trnodes[0]->outedge = ring1;
	treearray[0].trnodes[1]->outedge = ring1->next;
	treearray[0].trnodes[2]->outedge = ring1->next->next;
	
	treearray[0].trnodes[0]->start = true;
	
	copies = 1;
	
	printf("The base star tree:\n");
	printNewick(treearray[0].trnodes[0]);
	printf(";");
	printf("\n");
	
	origtree = &treearray[0];
	prevpos = 1;
	
	int call_number, currtree = 0;
	//int *iterations = &itercount;
	
	tree **offspring;
	
	taxon = 4;
	
	positions = posits_unr(taxon);
	
	do {
		
		for (i = 0; i < prevpos; ++i) {
			origtree = &treearray[i];
			offspring = (tree **) malloc((positions + 1) * sizeof(tree));
			if (offspring == NULL) {
				printf("In allunrooted(): unable to allocate memory for tree copies.\n");
			}
			offspring[0] = origtree;
			*counter = 1;
			
			for (j = 0; j < (positions - 1); ++j) {
				treecopy = &treearray[currtree + j + 1];
				copytree(origtree, treecopy, cpyptr);
				offspring[j + 1] = treecopy;
				
				/*printf("Offspring:");
				printNewick(offspring[j + 1]->trnodes[0]);
				printf("\n");*/
			}
			
			currtree = currtree + j;
			//printf("curr tree: %i\n", currtree);
			call_number = 0;
			
			for (j = 0; j < positions; ++j) {
				++call_number;
				insert_allp(offspring[j]->trnodes[0], origtree, taxon, call_number, counter);
				*counter = 1;
				/**begin debugging code*/
				printf("Debug tree: ");
				printNewick(offspring[j]->trnodes[0]);
				printf(";");
				printf("\n");
				/**end debugging code*/
			}
			
			free(offspring);
		}
		
		++taxon;
		prevpos = positions * prevpos;
		positions = posits_unr(taxon);		
		
		if (taxon == ntax) {
			printf("Adding final taxon\n");
		}
		
	} while (taxon <= ntax);
	
	printf("%i trees printed\n", currtree + 1);
	
	int c;
	
	do {
		printf("Enter c to continue: \n");
		c = getchar();
	} while (c != 'c');
	

	node *start;
	
	nodearray trn_ptr;
	
	//long long int destree = 0;
	
	printf("Destroying trees.\n");
	for (i = 0; i < ntrees; ++i) {
		//printf(">>\n");
		trn_ptr = treearray[i].trnodes;
		start = treearray[i].trnodes[0];
		//printf("Tree to be destroyed: %i\n", (i + 1));
		//printNewick(start);
		//printf("\n");
		//detree(start);
		denode(trn_ptr);
		//++destree;
		//printf("Destroying tree %lli\n", destree);
	}
	
	free(treearray);
	
	//int c;
	
	do {
		printf("Enter q to quit: \n");
		c = getchar();
	} while (c != 'q');
}