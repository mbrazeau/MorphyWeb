/*
 *  scraps.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 11-07-02.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */



node *tempnode;
tempnode = malloc(sizeof(node));


for (i = 0; i < ntaxa; ++i) {
	//newtips->trnodes[0] = malloc(sizeof(node));
	tempnode = newtips->trnodes[0];
	//tempnode->apomorphies = malloc(3 * sizeof(char));
	tempnode->apomorphies = newtipdata[0];


	
	origtree = &treearray[0];
	treecopy = &treearray[1];
	insert_allp(origtree->trnodes[0], origtree, treecopy, cpyptr, 4/*taxon*/);
	printf("first addseq iteration:");	
	printNewick(treearray[0].trnodes[0]);
	printf(";");
	printf("\n");
	printNewick(treearray[1].trnodes[0]);
	printf(";");
	printf("\n");
	
	
	origtree = &treearray[1];
	treecopy = &treearray[2];
	insert_allp(origtree->trnodes[0], origtree, treecopy, cpyptr, 4/*taxon*/);
	printf("second addseq iteration:");
	printNewick(treearray[1].trnodes[0]);
	printf(";");
	printf("\n");
	printNewick(treearray[2].trnodes[0]);
	printf(";");
	printf("\n");
	
	origtree = &treearray[2];
	treecopy = &treearray[3];
	insert_allp(origtree->trnodes[0], origtree, treecopy, cpyptr, 4/*taxon*/);
	printf("third addseq iteration:");
	printNewick(treearray[2].trnodes[0]);
	printf(";");
	printf("\n");
	printNewick(treearray[3].trnodes[0]);
	printf(";");
	printf("\n");
	
	printf("devisiting\n");
	devisit(treearray[0].trnodes[0]);
	
	origtree = &treearray[0];
	treecopy = &treearray[4];
	insert_allp(origtree->trnodes[0], origtree, treecopy, cpyptr, 5/*taxon*/);
	printf("fourth addseq iteration:");
	printNewick(treearray[0].trnodes[0]);
	printf(";");
	printf("\n");
	printNewick(treearray[4].trnodes[0]);
	printf(";");
	printf("\n");
	
	origtree = &treearray[4];
	treecopy = &treearray[5];
	insert_allp(origtree->trnodes[0], origtree, treecopy, cpyptr, 5/*taxon*/);
	printf("fifth addseq iteration:");
	printNewick(treearray[4].trnodes[0]);
	printf(";");
	printf("\n");
	printNewick(treearray[5].trnodes[0]);
	printf(";");
	printf("\n");
	
	origtree = &treearray[5];
	treecopy = &treearray[6];
	insert_allp(origtree->trnodes[0], origtree, treecopy, cpyptr, 5/*taxon*/);
	printf("sixth addseq iteration:");
	printNewick(treearray[5].trnodes[0]);
	printf(";");
	printf("\n");
	printNewick(treearray[6].trnodes[0]);
	printf(";");
	printf("\n");
	
	origtree = &treearray[6];
	treecopy = &treearray[7];
	insert_allp(origtree->trnodes[0], origtree, treecopy, cpyptr, 5/*taxon*/);
	printf("seventh addseq iteration:");
	printNewick(treearray[6].trnodes[0]);
	printf(";");
	printf("\n");
	printNewick(treearray[7].trnodes[0]);
	printf(";");
	printf("\n");
	
	origtree = &treearray[7];
	treecopy = &treearray[8];
	insert_allp(origtree->trnodes[0], origtree, treecopy, cpyptr, 5/*taxon*/);
	printf("eight addseq iteration:");
	printNewick(treearray[7].trnodes[0]);
	printf(";");
	printf("\n");
	printNewick(treearray[8].trnodes[0]);
	printf(";");
	printf("\n");

}

/* from the earlier attempt to generate all trees manually*/


 for (i = 0; i < 2; ++i) {
 treecopy = &treearray[i + 1];
 copytree(origtree, treecopy, ntax, cpyptr);
 
 }
 
 taxon = 4;
 prevpos = 0;
 positions = posits_unr(taxon);
 
 for (i = 0; i < positions; ++i) {
 origtree = &treearray[i + prevpos];
 insert_allp(origtree->trnodes[0], origtree, taxon, (i + 1), counter);
 *counter = 1;
 printNewick(origtree->trnodes[0]);
 printf("\n");
 }
 prevpos = positions;
 ++taxon;
 
 positions = posits_unr(taxon);
 
 for (i = 0; i < positions; ++i) {
 treecopy = &treearray[i + 3];
 copytree(origtree, treecopy, ntax, cpyptr);
 }
 
 origtree = &treearray[1];
 
 for (i = 0; i < positions; ++i) {
 treecopy = &treearray[i + 6];
 copytree(origtree, treecopy, ntax, cpyptr);
 }
 
 origtree = &treearray[2];
 
 for (i = 0; i < positions; ++i) {
 treecopy = &treearray[i + 11];
 copytree(origtree, treecopy, ntax, cpyptr);
 }
 
 *counter = 1;
 
 for (i = 0; i < positions; ++i) {
 origtree = &treearray[i];
 insert_allp(origtree->trnodes[0], origtree, taxon, (i + 1), counter);
 *counter = 1;
 printNewick(origtree->trnodes[0]);
 printf("\n");
 }
 
 *counter = 1;
 
 for (i = 0; i < positions; ++i) {
 origtree = &treearray[i + 5];
 insert_allp(origtree->trnodes[0], origtree, taxon, (i + 1), counter);
 *counter = 1;
 printNewick(origtree->trnodes[0]);
 printf("\n");
 }
 
 *counter = 1;
 
 for (i = 0; i < positions; ++i) {
 origtree = &treearray[i + 10];
 insert_allp(origtree->trnodes[0], origtree, taxon, (i + 1), counter);
 *counter = 1;
 printNewick(origtree->trnodes[0]);
 printf("\n");
 }



do {
	
	for (i = 0; i <= prevpos; ++i) {
		
		origtree = &treearray[i];
		
		for (j = 0; j < positions; ++j) {
			++call_number;
			treecopy = &treearray[i + *iterations];
			copytree(origtree, treecopy, ntax, cpyptr);
			insert_allp(origtree->trnodes[0], origtree, taxon, (j + 1), counter);
			origtree = treecopy;
			*counter = 1;
			*iterations = *iterations + 1;
			
			/**begin debugging code*/
			printf("Debug tree: ");
			printNewick(origtree->trnodes[0]);
			printf("\n");
			/**end debugging code*/
		}
		
	}
	
	++taxon;
	prevpos = positions;
	positions = posits_unr(taxon);		
	
} while (taxon <= ntax);

/*-------OLD MAIN FUNCTION--------------*/

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


/*-------END OLD MAIN-------------------*/