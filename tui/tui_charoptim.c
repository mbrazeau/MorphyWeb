//
//  tui_charoptim.c
//  Morphy
//
//  Created by mbrazeau on 06/06/2016.
//
//

#include "morphy.h"
#include "tuimfy.h"



void tui_test_lengths(tui_test_topol_t* ttopol)
{
    if (ttopol->tt_expected_length == ttopol->tt_measured_length) {
        dbg_ppass("expected and calcuated lengths matched");
    }
    else {
        if (ttopol->tt_measured_length > ttopol->tt_expected_length) {
            dbg_pfail("calculated length exceeds expected length");
        }
        else {
            dbg_pfail("calculated length is shorter than expected");
        }
    }
}

void tui_test_character_optimisation(void)
{
 
    int _12taxmax = 2;
    int _16taxmax = 2;
    
    char* trees12tax[_12taxmax];
    trees12tax[0] = (char*)"UNTITLED = [&R] ((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));";
    trees12tax[1] = (char*)"UNTITLED = [&R] ((((((1,2),3),4),5),6),(7,((8,9),(10,(11,12)))));";
    
    char* trees16tax[_16taxmax];
    trees16tax[0] = (char*)"Tree1 = [&R] ((((1,2),((3,4),5)),(6,7)),((((8,9),10),(11,(12,13))),(14,(15,16))));";
    trees16tax[1] = (char*)"Tree2 = [&R] ((((12,2),(7,16)),(5,4)),(((6,(9,10)),((((14,1),(13,11)),8),15)),3));";
    
    tui_test_topol_t* topoldat[4];
    
    
}