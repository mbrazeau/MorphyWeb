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

#define MORPHY_NUM_ITERATIONS 15
int ntax = 9;
int outtaxa[MAX_OG_SIZE];
int intaxa[MAX_IG_SIZE];
int maxstates = 5;
int numnodes;
bool OGdefined=false;
nodearray ingroup; 
nodearray outgroup;

void init_taxarray(int *taxarray)
{
	int i;
	
	for (i = 0; i < ntax; ++i) {
		taxarray[i] = i + 1;
	}
}

struct node * allocnode()
{
	node *newNode;
	newNode = (node *)malloc(sizeof(node));
	if (newNode == NULL)
	{
		printf("Error: failed to allocate new node.\n");
	}
	else
	{
		memset(newNode, 0, sizeof(node));
	}
	return newNode;
}

struct tree *alloctree()
{
	int i; //Loop counter
    tree *newtree;
	
	newtree = (tree*)malloc(sizeof(tree));
	
	if (newtree == NULL)
    {
		printf("Error: failed to allocate new tree.\n");
		return (struct tree*) 0;
	}

	newtree->trnodes = (node **)malloc( (numnodes) * sizeof(node*));
	
	for (i = 0; i < (2 * ntax); ++i)
    {
		newtree->trnodes[i] = allocnode();
	}
	
	for (i = 0; i < (2 * ntax); ++i)
    {
        /* First half of trnodes are initialized to this */
        if (i < ntax)
        {
		    newtree->trnodes[i]->tip = i + 1;
		    newtree->trnodes[i]->index = i + 1;
			newtree->trnodes[i]->next = NULL;
        }
        else /* second half are allocated as a newring */
        {
		    newtree->trnodes[i]->index = i;
		    newring(newtree->trnodes[i]);
        }
		
		newtree->trnodes[i]->outedge = NULL;
	}
	
	newtree->root = NULL;
	
	return (newtree);
}

/*
 * freetree - deletes an entire tree, by doing the mirror image
 * of alloctree.
 */
void freetree(tree *newtree)
{
    int i;
    
    /* free all of the trnodes */
	for (i = 0; i < (2 * ntax); ++i)
    {
        if (newtree->trnodes[i]->next)
        {
            deletering(newtree->trnodes[i]);
        }
        free(newtree->trnodes[i]);
    }
    /* free the trnode list */
    free(newtree->trnodes);
    /* free the tree */
    free(newtree);
	
}

void printNewick(node *n)
{	
	/* Prints the tree in Newick format (i.e. using brackets plus commas to separate
	 * equally ranked objects. In an unrooted tree, function will be called on a 
	 * terminal node that has its start variable set to 1. */
	
	node *p;
	
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
	
	r2 = allocnode();
	r3 = allocnode();
		
	r1->next = r2;
	r2->next = r3;
	r3->next = r1;
	
	r1->tip = r2->tip = r3->tip = 0;
	r1->index = r2->index = r3->index;
	r1->start = r2->start = r3->start = false;
	r1->dummy = r2->dummy = r3->dummy = false;
	
	r1->apomorphies = r2->apomorphies = r3->apomorphies = (char*) malloc(maxstates * sizeof(char));
}

/*
 * deletering - deletes a node (ring) by doing the mirror image of newring.
 */
void deletering(node *r1)
{
    node *r2, *r3;
    r2 = r1->next;
    r3 = r2->next;
    /* Note: apomorphies is only allocated once, r2 and r3 both point to the r1 apomorphies. */
    free(r1->apomorphies);
    free(r2);
    free(r3);
}

void addBranch(node *desc, node *ancest, tree *newtips, int *taxon)
{	
	/* Adds a sister terminal to an existing terminal node in the tree.
	 * Receives pointers of the ancestor and descendant of a branch, pointers to a tree struct
	 * containing an array of nodes and a pointer to an integer which is used to draw the 
	 * correct pointer for the insertion. Mostly useless, but served as an early test for
	 * my understanding of how to dynamically add branches to the tree. */
	
	node *newtip, *newn;
	
	newtips->trnodes[*taxon - 6] = (node *) malloc(sizeof(node));
	if (newtips->trnodes[*taxon - 6] == NULL) {
		printf("malloc failed in addBranch newtips->trnodes[*taxon - 6]\n");
	}
	newtips->trnodes[*taxon] = (node *) malloc(sizeof(node));
	if (newtips->trnodes[*taxon] == NULL) {
		printf("malloc failed in addBranch newtips->trnodes[*taxon]\n");
	}
	
	newtip = newtips->trnodes[*taxon - 6];
	newtip->apomorphies = malloc(3 * sizeof(char));
	if (newtip->apomorphies == NULL) {
		printf("malloc failed in addBranch newtip->apomorphies\n");
	}
	
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

void applyData(tree *currenttree, char **tipdata, int ntaxa, int *start)
{	
	int i;
	
	for (i = 0; i < ntaxa; ++i) {
		currenttree->trnodes[i]->apomorphies = tipdata[i + *start];
	}
	
	*start = *start + i;	// Maybe i - 1
}


void growcopy(node *templ, node *target, tree *newtree, int *iter)
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
		
		//newtree->trnodes[ntax + *iter] = (node *) malloc(sizeof(node));
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
		growcopy(p->outedge, target, newtree, iter);
		p = p->next;
		target = target->next;
	}
}


struct tree * copytree(tree *origtree /*, long long int *counter*/)
{
	int iteration = 1;
	int *iter = &iteration;
	node *base1;
	tree *treecopy;
	
	treecopy = alloctree();
	
	if (origtree->root == NULL) {
		
		/**** If the tree is unrooted ****/
		/* Generates the starting tip node and its corresponding internal
		 * node first. */
		
		treecopy->root == NULL;
		origtree->trnodes[0]->start = true;
		treecopy->trnodes[0]->start = origtree->trnodes[0]->start;
		treecopy->trnodes[0]->tip = origtree->trnodes[0]->tip;
		newring(treecopy->trnodes[ntax + 1]);
		treecopy->trnodes[0]->outedge = treecopy->trnodes[ntax + 1];
		treecopy->trnodes[ntax + 1]->outedge = treecopy->trnodes[0];
		growcopy(origtree->trnodes[0]->outedge, treecopy->trnodes[ntax + 1], treecopy, iter);
		
		/*printf("copiedtree: ");
		 printNewick(newtree->trnodes[0]);
		 printf("\n");*/
	} else {
		
		/**** If the tree is rooted ****/
		/* Simply generates the root ring first and begins from there */
		
		base1 = treecopy->trnodes[ntax];
		newring(base1);
		treecopy->root = treecopy->trnodes[ntax];
		growcopy(origtree->root, treecopy->root, treecopy, iter);
		
		printf("copiedtree: ");
		printNewick(treecopy->root);
		printf("\n");
	}
	
	//*counter++;
	
	return treecopy;
}



void point_bottom(node *n, node **nodes, int *counter)
{
	/* Re-sets the pointers to the internal nodes in the tree's
	 * node array to point to the 'bottom' (rootward) node in the ring*/
	
	node *p;
	
	if (n->tip) {
		return;
	}
	
	if (n->outedge) {
		nodes[*counter] = n;
		*counter = *counter + 1;
	}	
	
	p = n->next;
	while (p != n) {
		point_bottom(p->outedge, nodes, counter);
		p = p->next;		
	}
	
}

void rootOnTerminal(tree *trtoroot, int root)
{
	
	/*Roots the tree between a terminal (leaf) and an internal node*/
	
	int counter = ntax + 1;
	int *count_ptr;
	node *nodeptr, *r2, *r3;
	
	nodeptr = trtoroot->trnodes[root]->outedge;
	r2 = trtoroot->trnodes[ntax]->next;
	r3 = trtoroot->trnodes[ntax]->next->next;
	
	trtoroot->trnodes[root]->outedge = r2;
	r2->outedge = trtoroot->trnodes[root];
	
	r3->outedge = nodeptr;
	nodeptr->outedge = r3;
	
	trtoroot->root = trtoroot->trnodes[ntax];
	trtoroot->trnodes[ntax]->outedge = NULL;
	
	count_ptr = &counter;
	
	point_bottom(trtoroot->root, trtoroot->trnodes, count_ptr);
	
}

void unroot(tree *rootedtree)
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

void pauseit(void)
{
	int c;
	
	do {
		printf("Enter c to continue: \n");
		c = getchar();
	} while (c != 'c');
	c = getchar();
}

void rand_tree (void) 
{
	int i; //Loop counter
	
	/*generate a random tree*/
	
	tree **randtrees;
	randtrees = (tree **) malloc(MORPHY_NUM_ITERATIONS * sizeof(tree*));
	if (randtrees == NULL) {
		printf("Error in main(): failed malloc for randtrees\n");
		return;
	}
	
	for (i = 0; i < MORPHY_NUM_ITERATIONS; ++i) {
		randtrees[i] = randrooted();
		printNewick(randtrees[i]->root);
		printf(";\n");
	}
	
	for (i = 0; i < MORPHY_NUM_ITERATIONS; ++i) {
        freetree(randtrees[i]);
	}
	
	free(randtrees);

}

void numberOfNodes(void)
{
	numnodes = 2 * ntax - 1;
}

int main(void)
{
	numberOfNodes();
    /*rand_tree();
	pauseit();
    rand_tree();
	pauseit();*/
	
	tree *originaltree;
	tree *copiedtree;
	
	originaltree = randrooted();
	printNewick(originaltree->root);
	printf("\n");
	
	copiedtree = copytree(originaltree);
	printNewick(copiedtree->root);
	printf("\n");
	
	return 0;
}
