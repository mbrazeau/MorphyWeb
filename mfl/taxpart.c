/*
 *  taxpart.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 10/15/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "morphy.h"

void wipe_Og(int outtaxa[], nodearray outgroup)
{
	int i; // Loop counter
	
	for (i = 0; outtaxa[i]; ++i) {
		outtaxa[i] = '\0';
        if (outgroup)
        {
            free(outgroup);
        }
	}
}

void wipe_Ig(int intaxa[], nodearray ingroup)
{
	int i; // Loop counter
	
	for (i = 0; intaxa[i]; ++i) {
		intaxa[i] = '\0';
        if (ingroup)
        {
            free(ingroup);
        }
	}
}

void defOutgroup(int ntax, int outtaxa[], nodearray outgroup, int intaxa[], 
                 nodearray ingroup, bool *OGdefined)
{

	int i = 0, j = 0, k = 0; // Loop counters
	int OGsize, IGsize, taxnum; // Size of the outgroup, ingroup, and a temporary taxon number holder
	char c, input[100];
	char buffer[5];
	bool in_og = false;

	dbg_printf("Specify outgroup taxa:");
	
	/*takes input from the user and stores it as a string called input*/
	do {
		c = getchar();
		input[i] = c;
		++i;
	} while (c != '\n');
	
	input[i - 1] = ';';
	input[i] = '\0';
	
	/*parses the string in input and stores the integers in outtaxa*/
	for (i = 0; input[i]; ++i) {
		if (input[i] >= '0' && input[i] <= '9') {
			buffer[j] = input[i];
			++j;
			buffer[j] = '\0';
		} 		
		if (input[i] == ' ' || input[i] == ';') {
			taxnum = atoi(buffer);
			if (taxnum > ntax) {
				dbg_printf("Error: taxon %i not in dataset.\n", taxnum);
				dbg_printf("No taxon %i added to outgroup.\n", taxnum);
			} else if (taxnum <= ntax) {
				outtaxa[k] = taxnum;
				++k;
			}
			
			j = 0;
		}
	}
	
	outtaxa[k] = '\0';
	
	OGsize = k;
	
	IGsize = ntax - OGsize;
	
	if (IGsize < 3) {
		dbg_printf("Error: resulting ingroup is trivial\n");
		for (i = 0; outtaxa[i]; ++i) {
			outtaxa[i] = '\0';
		}
		dbg_printf("Taxa not transferred to outgroup\n");
		return;
	}
	
	*OGdefined = true;
	
	
	/*takes the complement of ntax and outtaxa to give the values in inttaxa, the ingroup terminals*/
	k = 0;
	
	for (i = 1; i <= ntax; ++i) {		
		for (j = 0; j <= OGsize && in_og == false; ++j) {
			if (i == outtaxa[j]) {
				in_og = true;
			}
		}
		if (in_og == false) {
			intaxa[k] = i;
			++k;
		}
		in_og = false;
	}
	
	intaxa[k] = '\0';
	
	/*feedback to the user that verifies the contents of outtaxa and intaxa*/
	
	dbg_printf("Taxa transferred to outgroup: ");
	for (i = 0; outtaxa[i]; ++i)
		dbg_printf("%i ", outtaxa[i]);
	dbg_printf("\n");
	
	dbg_printf("Taxa belonging to ingroup: ");
	for (i = 0; intaxa[i]; ++i)
		dbg_printf("%i ", intaxa[i]);
	dbg_printf("\n");
	
}

/* Add a sequential outgroup */

/* Add a polychotomous outgroup */

/* Add a constrained outgroup topology */
