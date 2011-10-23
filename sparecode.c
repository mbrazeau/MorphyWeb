/*
 *  sparecode.c
 *  Morphy
 *	
 *	This file is for spare pieces of code, test code, debugging code 
 *	or other code that is not essential to the functioning of the 
 *	putative release version of this program
 * 
 */


/* The general tree traversal function */

void traverse(node *n)
{	
	node *p;
	
	if (n->tip) {
		return;
	}
	
	p = n->next;
	while (p != n) {
		traverse(p->outedge);
		p = p->next;
	}
}


/* Traverses a tree and declares when it has reached a tip*/

void announceTips(node *n)
{	
	node *p;
	
	printf("calling announceTips()\n");
	
	if (n->tip) {
		printf("Reached tip %i\n", n->tip);
		return;
	}
	
	p = n->next;
	while (p != n) {
		announceTips(p->outedge);
		p = p->next;
	}	
	printf("announceTips() returning\n");
}


/* Counts the number of states in an apomorphy array */
 
void numstates(node *n)
{	
	node *p;
	
	printf("calling numstates\n");
	if (n->tip) {
		n->numstates = strlen((char*)n->apomorphies);
		printf("numstates: %i\n", n->numstates);
		return;
		
	}
	
	p = n->next;
	while (p != n) {
		numstates(p->outedge);
		n->numstates = strlen((char*)n->apomorphies);
		printf("numstates: %i\n", n->numstates);
		p = p->next;
	}
}

/*Traverses unrooted tree */

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


void addBranch(node *desc, node *ancest, tree *newtips, int *taxon)
{	
	/* Adds a sister terminal to an existing terminal node in the tree.
	 * Receives pointers of the ancestor and descendant of a branch, pointers to a tree struct
	 * containing an array of nodes and a pointer to an integer which is used to draw the 
	 * correct pointer for the insertion. Mostly useless, but served as an early test for
	 * my understanding of how to dynamically add branches to the tree. */
	
	node *newtip, *newn;
	
	newtips->trnodes[*taxon - 6] = malloc(sizeof(node));
	if (newtips->trnodes[*taxon - 6] == NULL) {
		printf("malloc failed in addBranch newtips->trnodes[*taxon - 6]\n");
	}
	newtips->trnodes[*taxon] = malloc(sizeof(node));
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