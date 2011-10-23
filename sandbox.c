/*
 *  sandbox.c
 *  Morphy
 *
 *  A place for sketching out code
 *
 */

/*includes so predictive text works*/
#include "morphy.h"
/*end include*/

void addBranch(node *desc, node *ancest, tree *newtips, int *taxon)
{	
	node *newtip, *newn1, *newn2, *newn3;
	
	newtips->trnodes[*taxon - ntax] = malloc(sizeof(node));
	newtips->trnodes[*taxon] = malloc(sizeof(node));
	
	newtip = newtips->trnodes[*taxon - ntax];
	newtip->apomorphies = malloc(3 * sizeof(char));
	newn1 = newtips->trnodes[*taxon];
	//newn1->apomorphies = malloc(5 * sizeof(char));
	
	newn2 = malloc(sizeof(node));
	newn3 = malloc(sizeof(node));
	
	newn1->next = newn2;
	newn2->next = newn3;
	newn3->next = newn1;
	
	newn2->apomorphies = newn3->apomorphies = newn1->apomorphies = malloc(5 * sizeof(char));
	
	ancest->outedge = newn1;
	newn1->outedge = ancest;
	
	newn1->next->outedge = desc;
	desc->outedge = newn1->next;
	
	newn1->next->next->outedge = newtip;
	newtip->outedge = newn1->next->next;
	
	newtip->tip = *taxon;
	newn1->tip = newn2->tip = newn3->tip = 0;
	*taxon = *taxon + 1;
}



void copytree(tree *origtree, tree *newtree, int ntax, long long int counter)
{
	int i, totalnodes;
	
	totalnodes = 2 * ntax - 1;
	
	for (i = 0; i < totalnodes; ++i) {
		newtree->trnodes[i] = malloc(sizeof(node));
	}
}


void newring(node *r1)
{
	node *r2, *r3;
	
	r1 = malloc(sizeof(node));
	r2 = malloc(sizeof(node));
	r3 = malloc(sizeof(node));
	
	r1->next = r2;
	r2->next = r3;
	r3->next = r1;
	
	r1->tip = r2->tip = r3->tip = 0;
	r1->index = r2->index = r3->index = 0;
	r1->visited = r2->visited = r3->visited = 0;
	r1->start = r2->start = r3->start = false;
	r1->dummy = r2->dummy = r3->dummy = false;
	
}

void insert_allp(node *n, tree *origtree, tree *treecopy, long long int *copies, int taxon)
{
	node *p;
	
	if (n->tip && !n->start && !n->visited) {
		n->visited = 1;
		copytree(origtree, treecopy, ntax, copies);
		newbranch(n, n->outedge, origtree, taxon);
		insert_allp(n->outedge, origtree, treecopy, copies, taxon);
		return;
	}
	
	if (n->start && !n->visited) {
		n->visited = 1;
		n->outedge->visited = 1;
		copytree(origtree, treecopy, ntax, copies);
		newbranch(n->outedge, n, origtree, taxon);
		return;
	}
	
	if (n->start && n->visited) {
		insert_allp(n->outedge, origtree, treecopy, copies, taxon);
		return;
	}
	
	if (!n->start && !n->tip && n->visited) {
		return;
	}
	
	p = n->next;
	n->visited = 1;
	while (p != n) {
		if (!p->outedge->visited) {
			p->visited = 1;
			p->outedge->visited = 1;
			copytree(origtree, treecopy, ntax, copies);
			newbranch(p->outedge, p, origtree, taxon);
			return;
		}
		
		insert_allp(p->outedge, origtree, treecopy, copies, taxon);
		p = p->next;
	}	
	
}


void insert_allp(node *n, tree *origtree, tree *treecopy, long long int *copies, int taxon)
{
	node *p;
	
	if (n->tip && !n->start) {
		n->visited = 1;
		return;
	}
	
	if (n->start && !n->visited) {
		n->visited = 1;
		newbranch(n, n->outedge, origtree, taxon);
		origtree = (origtree + 1);
		insert_allp(n->outedge, origtree, treecopy, copies, taxon);
		return;
	}
	
	if (!n->tip && n->outedge->visited) {
		n->visited = 1;
		insert_allp(n->next, origtree, treecopy, copies, taxon);
		return;
	}
	
	if (!n->visited) {
		newbranch(n, n->outedge, origtree, taxon);
		return;
	}
	
	p = n->next;
	while (p != n) {
		insert_allp(p->outedge, origtree, treecopy, copies, taxon);
		p = p->next;
	}
	
}	


void insert_allp(node *n, tree *origtree, tree *treecopy, long long int *copies, int taxon)
{
	node *p;
	
	if (n->tip && !n->start && !n->visited) {
		n->visited = 1;
		copytree(origtree, treecopy, ntax, copies);
		newbranch(n, n->outedge, origtree, taxon);
		return;
	}
	
	if (n->start && !n->outedge->visited) {
		n->visited = 1;
		n->outedge->visited = 1;
		copytree(origtree, treecopy, ntax, copies);
		newbranch(n->outedge, n, origtree, taxon);
		return;
	}
	
	if (n->start && n->outedge->visited) {
		insert_allp(n->outedge, origtree, treecopy, copies, taxon);
		return;
	}
	
	if (n->tip && n->visited && !n->start) {
		return;
	}
	
	if (!n->tip && !n->outedge->visited) {
		//n->visited = 1;
		n->outedge->visited = 1;
		copytree(origtree, treecopy, ntax, copies);
		newbranch(n->outedge, n, origtree, taxon);
		return;
	}
	
	p = n->next;
	while (p != n) {
		if (!p->outedge->visited) {
			insert_allp(p->outedge, origtree, treecopy, copies, taxon);
			return;
		}
		//insert_allp(p->outedge, origtree, treecopy, copies, taxon);
		p = p->next;
	}
	
}

void insert_allp(node *n, tree *origtree, int taxon, int calln, int *counter)
{
	node *p;
	
	if (n->tip && !n->start) {
		return;
	}
	
	if (counter == calln) {
		newbranch(n, n->outedge, origtree, taxon);
		*counter++;
		return;
	}
	
	p = n->next;
	while (p != n) {
		insert_allp(p->outedge, origtree, taxon, calln, counter);
		p = p->next;
	}
	
}	

// The add all algorithm

		
int taxadded = 3;
int taxon = 4;
int prevpos = 1;

int itercount = 1;
int *iterations = &itercount;

positions = posits_unr(taxon);

origtree = &treearray[i + positions + /*****something goes in here*****/];

do {
	
	for (i = 0; i < prevpos; ++i) {
		
		origtree = &treearray[i];
	
		for (j = 0; j < positions; ++j) {
			treecopy = &treearray[i + *iterations];
			copytree(origtree, treecopy, ntax, cpyptr);
			insert_allp(origtree->trnodes[0], origtree, taxon, call_number, counter);
			
			*counter = 1;
			*iterations = *iterations + 1;
			
			//put newick printing debug code here
		}
		
	}
	
	++taxon;

	positions = posits_unr(taxon);		
	
} while (taxon <= ntax);



/***guff***/
for	(i = 0; i < positions; ++i) {
	origtree = &treearray[i + positions + /*****something goes in here*****/];
	insert_allp(origtree->trnodes[0], origtree, taxon, call_number, counter);
	*counter = 1;
	printNewick(origtree->trnodes[0]);
	printf("\n");
}

/**begin debugging code*/
printNewick(origtree->trnodes[0]);
printf("\n");
/**end debugging code*/


/**
 * Begin the random number generator */

int randrange(int lowbound, int upbound)
{
	
    int i, n, rn;
    unsigned long int m;
    gsl_rng *r;
    const gsl_rng_type *T;    
	
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
	

	
    torn = malloc(n * sizeof(unsigned long int));
	
    m = 10;
	
    for (i = 0; i <= n; ++i) {
		torn[i] = gsl_rng_uniform_int(r, m);
		printf("%lu\n", (torn[i]) + 1);
    }
    
    return 0;
}

/* End of random number generator*/

/*A function that partitions the ingroup and outgroup*/

void inout_partition(tree t, node *ingroup, node *outgroup, int *outtaxa)
{
	int i, j, ogsize, igsize;
	int *intaxa;
	
	for (i = 0; outtaxa[i]; ++i) {
		outgroup[i] = t.trnodes[outtaxa[i]];
	}
	
	ogsize = i - 1;
	
	igsize = ntax - ogsize;
	intaxa = malloc(igsize * sizeof(int));
	
	for (i = 0; i < igsize; ++i) {
		for (j = 0; j < ogsize; ++j) {
			if (i != outtaxa[j]) {
				intaxa[i] = i + j;
			}
		}
	}
	
	for (i = 0; i < igsize; ++i) {
		ingroup[i] = t.trnodes[intaxa[i]];
	}
	
	free(intaxa);
	
}

/*A function that allocates a tree and all of its nodes for a given ntax*/

void alloctree(tree *newtree)
{
	newtree = (tree *)malloc(sizeof(tree));

	newtree->trnodes = (nodearray)malloc( (2 * ntax - 1) * sizeof(node *));
										 
}


/*an attempt to use the gls rng to shuffle the taxon list*/

void shuffle(int *taxarray, int length)

{
		int i, j, t;
		gsl_rng *r;
		const gsl_rng_type * T;
		
		gsl_rng_env_setup();
		
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
		
		unsigned long int gsl_RNG_MAX = gsl_rng_max (r);
		
		if (length > 1) {
			for (i = 0; i < length - 1; i++) {
				j = i + gsl_rng_get(r) / (gsl_RNG_MAX / (length - i) + 1);
				t = taxarray[j];
				taxarray[j] = taxarray[i];
				taxarray[i] = t;
			}
		}
	}