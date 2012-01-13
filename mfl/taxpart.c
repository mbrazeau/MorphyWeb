/*
 *  taxpart.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 10/15/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "morphy.h"

extern int outtaxa[MAX_OG_SIZE];
extern int intaxa[MAX_IG_SIZE];
extern bool OGdefined;
extern nodearray ingroup; 
extern nodearray outgroup;


int strToInt (char string[])
{
    int i = 0, intValue, result = 0;
	
    while (string[i] >= '0' && string[i] <= '9')
    {
		intValue = string[i] - '0';
		result = result * 10 + intValue;
		++i;
    }
	
	return result;
}

void wipe_Og(void)
{
	int i; // Loop counter
	
	for (i = 0; outtaxa[i]; ++i) {
		outtaxa[i] = '\0';
	}
}

void wipe_Ig(void)
{
	int i; // Loop counter
	
	for (i = 0; intaxa[i]; ++i) {
		intaxa[i] = '\0';
	}
}

void defOutgroup(int ntax)
{

	int i = 0, j = 0, k = 0; // Loop counters
	int OGsize, IGsize, taxnum; // Size of the outgroup, ingroup, and a temporary taxon number holder
	char c, input[100];
	char buffer[5];
	bool in_og = false;

	printf("Specify outgroup taxa:");
	
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
			taxnum = strToInt(buffer);
			if (taxnum > ntax) {
				printf("Error: taxon %i not in dataset.\n", taxnum);
				printf("No taxon %i added to outgroup.\n", taxnum);
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
		printf("Error: resulting ingroup is trivial\n");
		for (i = 0; outtaxa[i]; ++i) {
			outtaxa[i] = '\0';
		}
		printf("Taxa not transferred to outgroup\n");
		return;
	}
	
	OGdefined = true;
	
	
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
	
	printf("Taxa transferred to outgroup: ");
	for (i = 0; outtaxa[i]; ++i)
		printf("%i ", outtaxa[i]);
	printf("\n");
	
	printf("Taxa belonging to ingroup: ");
	for (i = 0; intaxa[i]; ++i)
		printf("%i ", intaxa[i]);
	printf("\n");
	
}