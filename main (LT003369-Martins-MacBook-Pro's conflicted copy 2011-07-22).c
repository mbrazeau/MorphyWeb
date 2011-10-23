/*
 *  morphy.h
 *  Morphy
 *
 *  Created by Martin Brazeau on 11-04-26.
 *  Copyright 2011. All rights reserved.
 *  Morphy is provided as is with no warranty of any kind.
 *
 */


#include "morphy.h"
#include <gsl/gsl_rng.h>

int ntax = 9;


void printNewick(node *n)
{	
	node *p;
	
	/* In an unrooted tree, function will be called on a terminal node that has 
	 * its start variable set to 1. */
	
	if (n->start) {
		printf("(%i,", n->tip);
		printNewick(n->outedge);
		printf(")");
		return;
	}	
	
	if (n->tip && n->outedge->next->outedge->tip && !n->outedge->next->outedge->start) {
		printf("%i", n->tip);
		return;
	}
	
	if (n->tip && !n->start) {
		printf("%i", n->tip);
		return;
	}
	
	printf("(");
	
	p = n->next;
	while (p != n) {
		printNewick(p->outedge);
		p = p->next;
		if (p != n) {
			printf(",");
		}
	}	
	printf(")");
}

void fitchdown(node *leftdesc, node *rightdesc, node *ancestor, int *stepcount)
{
    int i, j, k = 0;
	char c;
	
	/* Makes synapset of ancestor the intersection of set1 and set2
	 * if the set is non-null */
    for (i = 0; leftdesc->apomorphies[i]; ++i) {
		for (j = 0; rightdesc->apomorphies[j]; ++j) {
			if (leftdesc->apomorphies[i] == rightdesc->apomorphies[j]) {
				ancestor->apomorphies[k] = rightdesc->apomorphies[j];
				++k;
				ancestor->apomorphies[k]= '\0';
			}
		}
    }
	
	/* Makes synapset of ancestor the union of set1 and set2 if 
	 * the intersection set is null */
	if (k == 0) {
		i = 0; 
		j = 0;
		
		do {
			c = leftdesc->apomorphies[i];
			ancestor->apomorphies[i] = c;
			++i;
		} while ( leftdesc->apomorphies[i] );
		
		while (rightdesc->apomorphies[j] ){
			ancestor->apomorphies[i + j] = rightdesc->apomorphies[j];
			++j;
		} 
		
		ancestor->apomorphies[i + j]= '\0';
		
		*stepcount = *stepcount + 1;
	} 

}

void treelen(node *n, int *stepcount)
{	
	node *p;
	
	if (n->tip) {
		return;
	}
	
	p = n->next;
	while (p != n) {
		treelen(p->outedge, stepcount);
		printf("tip: %i; apomorphy: %s\n", p->outedge->tip, p->outedge->apomorphies);
		p = p->next;
	}
	
	fitchdown(n->next->outedge, n->next->next->outedge, n, stepcount);
	printf("treelen returning\n");
}

void nacount(node *n, int *stepcount)
{	
	/* Will eventually become a function analogous to treelen but that checks to see if a
	 * non-applicable score is used. In which case, it will skip the call to fitchdown, and 
	 * 'wait' to see if the next outgroup is similarly sharing a non-applicable score before 
	 * before generating the synapset for the node. */
	
	node *p;
	
	if (n->tip) {
		return;
	}
	
	p = n->next;
	while (p != n) {
		nacount(p->outedge, stepcount);
		printf("tip: %i; apomorphy: %s\n", p->outedge->tip, p->outedge->apomorphies);
		p = p->next;
	}
	
	fitchdown(n->next->outedge, n->next->next->outedge, n, stepcount);
	printf("treelen returning\n");
}

void newring(node *r1)
{
	/* Generates a new internal node composed of a ring of node structures
	 * joined in a unidirectional manner by their next pointers. Used by any
	 * function that might dynamically add branches at run-time. */
	
	node *r2, *r3;

	r2 = malloc(sizeof(node));
	r3 = malloc(sizeof(node));
	
	r1->next = r2;
	r2->next = r3;
	r3->next = r1;
	
	r1->tip = r2->tip = r3->tip = 0;
	r1->index = r2->index = r3->index = 0;
	r1->start = r2->start = r3->start = false;
	r1->dummy = r2->dummy = r3->dummy = false;
	
	r1->apomorphies = r2->apomorphies = r3->apomorphies = malloc(5 * sizeof(char));
}



void addBranch(node *desc, node *ancest, tree *newtips, int *taxon)
{	
	/* Adds a sister terminal to an existing terminal node in the tree.
	 * Receives pointers of the ancestor and descendant of a branch, pointers to a tree struct
	 * containing an array of nodes and a pointer to an integer which is used to draw the 
	 * correct pointer for the insertion. Mostly useless, but served as an early test for
	 * my understanding of how to dynamically add branches to the tree. */
	
	node *newtip, *newn;
	
	newtips->trnodes[*taxon - 6] = malloc(sizeof(node));
	newtips->trnodes[*taxon] = malloc(sizeof(node));
	
	newtip = newtips->trnodes[*taxon - 6];
	newtip->apomorphies = malloc(3 * sizeof(char));
	newn = newtips->trnodes[*taxon];
	
	newring(newn);
	
	ancest->outedge = newn;
	newn->outedge = ancest;
	
	newn->next->outedge = desc;
	desc->outedge = newn->next;
	
	newn->next->next->outedge = newtip;
	newtip->outedge = newn->next->next;
	
	newtip->tip = *taxon;
	*taxon = *taxon + 1;
}

void addTip(node *n, tree *newtips, int *taxon)
{
	node *p;
	
	if (n->tip) {
		addBranch(n, n->outedge, newtips, taxon);
		return;
	}
	
	p = n->next;
	while (p != n) {
		addTip(p->outedge, newtips, taxon);
		p = p->next;
	}
}
	
void applyData(tree *newtips, char **newtipdata, int ntaxa)
{	
	int i;
	
	for (i = 0; i < ntaxa; ++i) {
		newtips->trnodes[i]->apomorphies = newtipdata[i];
	}
	
}


void growcopy(node *templ, node *target, tree *newtree, int ntax, int *iter)
{	
	/* growcopy traverses the original tree, adding branches to newtree in preorder.
	 * templ is a template node in the original tree, and target is the node in the
	 * new tree to which a new branch and node are to be attached */
	
	node *p;
	
	if (templ->tip) {
		newtree->trnodes[templ->tip - 1]->outedge = target;
		target->outedge = newtree->trnodes[templ->tip - 1];;
		newtree->trnodes[templ->tip - 1]->tip = templ->tip;
		newtree->trnodes[templ->tip - 1]->visited = templ->visited;
		newtree->trnodes[templ->tip - 1]->start = templ->start;
		
		return;
	}
	
	if (templ->outedge != NULL && !templ->outedge->start) {
		
		node *newr1;
		
		newtree->trnodes[ntax + *iter] = malloc(sizeof(node));
		newr1 = newtree->trnodes[ntax + *iter];
		newr1->visited = templ->visited;
		newring(newr1);
		
		newr1->outedge = target;
		target->outedge = newr1;
		
		target = newr1->next;
		
		*iter = *iter + 1;

	} else {
		target = target->next;
	}
	
	p = templ->next;
	while (p != templ) {
		growcopy(p->outedge, target, newtree, ntax, iter);
		p = p->next;
		target = target->next;
	}
}


void copytree(tree *origtree, tree *newtree, int ntax, long long int *counter)
{
	int i, totalnodes, iteration = 1;
	int *iter = &iteration;
	node *base1;
	
	totalnodes = 2 * ntax - 1;
	
	// Allocate memory for the tree and the node array.
	
	if (!newtree) {	
		newtree = malloc(sizeof(tree));
	}
	
	newtree->trnodes = malloc(totalnodes * sizeof(node));
	
	for (i = 0; i < totalnodes; ++i) {
		newtree->trnodes[i] = malloc(sizeof(node));
		if (i >= ntax) {
			newtree->trnodes[i]->tip = 0;
		}
	}
	
	if (origtree->root == NULL) {
		
		/**** If the tree is unrooted ****/
		/* Generates the starting tip node and its corresponding internal
		 * node first. */
		 
		//printf("unrooted copying\n");
		newtree->root == NULL;
		origtree->trnodes[0]->start = true;
		newtree->trnodes[0]->start = origtree->trnodes[0]->start;
		newtree->trnodes[0]->tip = origtree->trnodes[0]->tip;
		newring(newtree->trnodes[ntax + 1]);
		newtree->trnodes[0]->outedge = newtree->trnodes[ntax + 1];
		newtree->trnodes[ntax + 1]->outedge = newtree->trnodes[0];
		growcopy(origtree->trnodes[0]->outedge, newtree->trnodes[ntax + 1], newtree, ntax, iter);
		
		/*printf("copiedtree: ");
		printNewick(newtree->trnodes[0]);
		printf("\n");*/
	} else {
		
		/**** If the tree is rooted ****/
		/* Simply generates the root ring first and begins from there */
		
		base1 = newtree->trnodes[ntax];
		newring(base1);
		newtree->root = newtree->trnodes[ntax];
		growcopy(origtree->root, newtree->root, newtree, ntax, iter);
		
		printf("copiedtree: ");
		printNewick(newtree->root);
		printf("\n");
	}
	
	*counter++;
}


void deroot(tree *rootedtree)
{
	node *proot, *leftdesc, *rightdesc;
	
	proot = rootedtree->root;
	leftdesc = proot->next->outedge;
	rightdesc = proot->next->next->outedge;

	leftdesc->outedge = rightdesc;
	rightdesc->outedge = leftdesc;
	
	rootedtree->root = NULL;
	rootedtree->trnodes[0]->start = true; // Could, in the future, be user-defined.
	
}

void travunrooted(node *n)
{
	node *p;
	
	if (n->tip && !n->start) {
		printf("tip: %i\n", n->tip);
		return;
	}
	
	if (n->start) {
		printf("starttip: %i\n", n->tip);
		travunrooted(n->outedge);
		//n->start = 0;
		return;
	}
	
	p = n->next;
	while (p != n) {
		travunrooted(p->outedge);
		p = p->next;
	}	
	printf("Int node\n");
	
}

int main (void) 
{
	
	int length = 0;
	int *stepcount = &length;

	/* A temporary hard-coded tree */
	/* This is used to develop code for traversing and manipulating trees until
	 * code is written for generating arbitrary trees at run-time. */
	
	node *root;

	 
	/* Terminal nodes */
	node t1, t2, t3, t4, t5;
	
	t1.tip = 1;
	t1.apomorphies = "1";
	t1.start = false;
	t2.tip = 2;
	t2.apomorphies = "0";
	t2.start = false;
	t3.tip = 3;
	t3.apomorphies = "0";
	t3.start = false;
	t4.tip = 4;
	t4.apomorphies = "1";
	t4.start = false;
	t5.tip = 5;	
	t5.apomorphies = "1";
	t5.start = false;
	
	/* Internal nodes */
	node r1, r2, r3;		// The root node
	
	r1.tip = 0;
	r2.tip = 0;
	r3.tip = 0;
	r1.start = r2.start = r3.start = false;
	
	char *datar1 = malloc(5 * sizeof(char));
	r1.apomorphies = datar1;	//Allocate the apomorphies array
	
	node i11, i12, i13;		// The node of t1 and t2
	
	i11.tip = 0;
	i12.tip = 0;
	i13.tip = 0;
	i11.start = i12.start = i13.start = false;
	
	char *datai11 = malloc(5 * sizeof(char));
	i11.apomorphies = datai11;		//Allocate the apomorphies array
	
	node i21, i22, i23;		// The node of t3, t4, and t5
	
	i21.tip = 0;
	i22.tip = 0;
	i23.tip = 0;
	i21.start = i22.start = i23.start = false;
	
	char *datai21 = malloc(5 * sizeof(char));
	i21.apomorphies = datai21;		//Allocate the apomorphies array
	
	node i31, i32, i33;		// The node of t4 and t5
	
	i31.tip = 0;
	i32.tip = 0;
	i33.tip = 0;
	i31.start = i32.start = i33.start = false;
	
	char *datai31 = malloc(5 * sizeof(char));
	i31.apomorphies = datai31;
	//Allocate the apomorphies array
	
	/* Set up the internal node rings */
	
	r1.next = &r2;
	r2.next = &r3;
	r3.next = &r1;
	
	i11.next = &i12;
	i12.next = &i13;
	i13.next = &i11;
	
	i21.next = &i22;
	i22.next = &i23;
	i23.next = &i21;
	
	i31.next = &i32;
	i32.next = &i33;
	i33.next = &i31;
	
	/* Set up the internal nodes to point to each other */
	
	root = &r1;
	r1.outedge = NULL;
	
	r2.outedge = &i11;
	i11.outedge = &r2;
	
	r3.outedge = &i21;
	i21.outedge = &r3;
	
	i23.outedge = &i31;
	i31.outedge = &i23;
	
	/* Add the terminal nodes */
	
	i12.outedge = &t1;
	t1.outedge = &i12;
	
	i13.outedge = &t2;
	t2.outedge = &i13;
	
	i22.outedge = &t3;
	t3.outedge = &i22;
	
	i32.outedge = &t4;
	t4.outedge = &i32;
	
	i33.outedge = &t5;
	t5.outedge = &i33;
	

	/* End of hard-coded tree*/
	
	/* New nodes for tree */
	
	tree newTips;
	newTips.trnodes = malloc(22 * sizeof(node));
	
	int taxon = 6;
	int *ntaxptr = &taxon;
	
	tree *tipset = &newTips;
	
	char *newtipdata[] = { "0", "1", "0", "1", "1"}; 
	
	//announceTips(root);
	printNewick(root);
	printf("\n");
	treelen(root, stepcount);
	printf("Tree length: %i\n", *stepcount);
	printf("root chars: %s\n", r1.apomorphies);
	
	tree originaltree;
	originaltree.root = &r1;
	
	tree *origtreeptr = &originaltree;
	
	tree *trcopyptr1 = NULL;
	
	long long int tempc = 0;
	long long int *tempcounter = &tempc;
	
	copytree(origtreeptr, trcopyptr1, 5, tempcounter);
	
	addTip(root, tipset, ntaxptr);
	//announceTips(root);
	printNewick(root);
	printf("\n");
	
	length = 0;
	
	applyData(tipset, newtipdata, 5);
	treelen(root, stepcount);
	printf("Tree length: %i\n", *stepcount);
	
	origtreeptr = &originaltree;
	
	printNewick(origtreeptr->root);
	printf("\n");
	
	//tree treecopy;
	//tree *trcopyptr2;
	
	//copytree(origtreeptr, trcopyptr2, 10, tempcounter);
	//announceTips(trcopyptr->root);
	//printNewick(trcopyptr->root);
	/*The real business */
	
	//tree *treeset;
	//allrooted(treeset, 10);
	
	tree exhaustive;
	tree *exhaustme = &exhaustive;
	
	allunrooted(exhaustme, ntax);
	
	/*origtreeptr->trnodes[0] = &t1;
	printf("derooting hardcoded tree\n");
	deroot(origtreeptr);
	printNewick(&t1);
	
	t1.start = 1;
	//travunrooted(&t1);
	
	tree *trcopyptr3;
	
	
	copytree(origtreeptr, trcopyptr3, 10, tempcounter);
	printf("dummywait\n");
	//travunrooted(&t1);
	printf("dummywait\n");
	printNewick(&t1);*/
	
	return 0;
	
}
