/*
 *  mfl_brswap.c
 *
 *  THE MORPHY FUNCTION LIBRARY
 *  A library for phylogenetic analysis with emphasis on parsimony and
 *  morphology (but someday other methods)
 *
 *  Copyright (C) 2016  by Martin D. Brazeau, Thomas Guillerme,
 *  and Chris Desjardins
 *
 *  Created by Martin Brazeau and Chris Desjardins on 1/26/12.
 *  Some data structs, routines and ideas derived from:
 *      - The PHYLIP package by Joe Felsenstein
 *          <http://evolution.genetics.washington.edu/phylip.html>
 *      - MrBayes by John Huelsenbeck and Fredrik Ronquist
 *          <http://mrbayes.sourceforge.net/>
 *
 *  Any bugs, errors, inefficiences and general amateurish handling are our own
 *  and most likely the responsibility of MDB. We make no guarantees.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include "morphy.h"

void mfl_temp_nodeweight_set(mfl_node_t* n)
{
    mfl_node_t *p = NULL;
    int weightcount = 0;
    
    if (n->nodet_tip) {
        n->nodet_isbottom = true;
        return;
    }
    
    p = n->nodet_next;
    do {
        mfl_temp_nodeweight_set(p->nodet_edge);
        weightcount += p->nodet_edge->nodet_weight;
        p = p->nodet_next;
    } while (p != n);
    
    n->nodet_weight = weightcount;
    n->nodet_isbottom = true;
    
    return;
}


void mfl_save_topology(mfl_tree_t* t, mfl_treebuffer_t* trbuf, mfl_searchrec_t* searchrec)
{
    mfl_tree_t* topolcopy = mfl_copy_tree_topology(t);
    mfl_append_tree_to_treebuffer(topolcopy, trbuf, searchrec->sr_handle_ptr);
}

bool mfl_check_node_wt(int weight, mfl_nodeweights_t* ndweightlist)
{
    assert(weight < ndweightlist->nwl_max);
    
    if (ndweightlist->nwl_weights_list[weight-1]) {
        dbg_printf("This weight has been seen before\n");
        return true;
    }
    else {
        return false;
    }
}


void mfl_set_node_wt(bool value, int weight, mfl_nodeweights_t* ndweightlist)
{
    assert(weight < ndweightlist->nwl_max);
    
    ndweightlist->nwl_weights_list[weight-1] = value;
}


mfl_nodeweights_t* mfl_create_node_wt_list(int num_taxa)
{
    
    int i = 0;
    int num_weights = 0;
    mfl_nodeweights_t* newnwl = (mfl_nodeweights_t*)mfl_malloc(sizeof(mfl_nodeweights_t), 0, __FXN_NAME__);
    
    num_weights = num_taxa - 1; // The number of useful weights in branch breaking
    
    newnwl->nwl_max = num_weights;
    
    newnwl->nwl_weights_list = (bool*)mfl_malloc(num_weights * sizeof(bool), 0, __FXN_NAME__);
    
    return newnwl;
}


inline mfl_node_t* mfl_clip_branch(mfl_node_t* n, mfl_cliprec_t* cliprec)
{
    mfl_node_t* retnode = NULL;
    
    cliprec->src1 = n->nodet_edge->nodet_next;
    cliprec->src2 = cliprec->src1->nodet_next;
    cliprec->tgt1 = cliprec->src1->nodet_edge;
    cliprec->tgt2 = cliprec->src2->nodet_edge;
    
    if (cliprec->src1->nodet_isbottom) {
        retnode = cliprec->src1;
    } else {
        retnode = cliprec->src2;
    }
    
    mfl_disconnect_node(cliprec->src1);
    mfl_disconnect_node(cliprec->src2);
    
    mfl_join_node_edges(cliprec->tgt1, cliprec->tgt2);
    
    return retnode;
}


inline void mfl_restore_branching(mfl_cliprec_t* cliprec)
{
    mfl_join_node_edges(cliprec->src1, cliprec->tgt1);
    mfl_join_node_edges(cliprec->src2, cliprec->tgt2);
}


inline void mfl_temp_rebranching(mfl_node_t* src, mfl_node_t* tgt, mfl_cliprec_t* regraft)
{
    regraft->tgt1 = tgt;
    regraft->tgt2 = tgt->nodet_edge;
    regraft->src1 = src;
    
    if (src->nodet_next->nodet_edge) {
        regraft->src2 = src->nodet_next->nodet_next;
    }
    else {
        regraft->src2 = src->nodet_next;
        assert(!regraft->src2->nodet_edge);
    }
    
    mfl_join_node_edges(regraft->src1, regraft->tgt1);
    mfl_join_node_edges(regraft->src2, regraft->tgt2);
}


inline void mfl_undo_temp_rebranching(mfl_cliprec_t* regraft)
{
    regraft->src1->nodet_edge = NULL;
    regraft->src2->nodet_edge = NULL;
    mfl_join_node_edges(regraft->tgt1, regraft->tgt2);
}


mfl_node_t* mfl_get_src_root(mfl_node_t *n)
{
    mfl_node_t* p = n->nodet_next;
    
    while (!p->nodet_edge) {
        p = p->nodet_next;
    }
    
    return p->nodet_edge;
}

void mfl_subtree_reroot_traversal(mfl_node_t* subtr, mfl_node_t* src, mfl_node_t* tgt, mfl_cliprec_t* cliprec, mfl_searchrec_t* searchrec, bool neighbour_rule)
{

    mfl_node_t *p;
    mfl_cliprec_t reroot;
    
    // Reroot shit here
    mfl_temp_rebranching(cliprec->src1, subtr, &reroot);
    
    // Perform regrafting
    mfl_regraft_subtree(src, tgt, searchrec, neighbour_rule);
    
    // Undo reroot here
    mfl_undo_temp_rebranching(&reroot);
    
    if (subtr->nodet_tip) {
        return;
    }
    
    p = subtr->nodet_next;
    
    do {
        mfl_subtree_reroot_traversal(p->nodet_edge, src, tgt, cliprec, searchrec, neighbour_rule);
        p = p->nodet_next;
    } while (p != subtr);
    
}


void mfl_reroot_subtree(mfl_node_t* src, mfl_node_t* tgt, mfl_searchrec_t* searchrec, bool neighbour_rule)
{
    
    mfl_cliprec_t cliprec;
    mfl_node_t *p = NULL;
    mfl_node_t *subtr = NULL;
    
    
    mfl_regraft_subtree(src, tgt, searchrec, neighbour_rule);

    p = mfl_get_src_root(src);
    
    if (!p->nodet_tip) {
        
        mfl_clip_branch(p->nodet_edge, &cliprec); // Remove the root of the subtree from the subtree
        
        if (!cliprec.tgt1->nodet_tip) {
            if (neighbour_rule) {
                mfl_subtree_reroot_traversal(cliprec.tgt1->nodet_next->nodet_edge, src, tgt, &cliprec, searchrec, false);
                mfl_subtree_reroot_traversal(cliprec.tgt1->nodet_next->nodet_next->nodet_edge, src, tgt, &cliprec, searchrec, false);
            } else {
                mfl_subtree_reroot_traversal(cliprec.tgt1, src, tgt, &cliprec, searchrec, neighbour_rule);
            }
        }
        
        if (!cliprec.tgt2->nodet_tip) {
    
            if (neighbour_rule) {
                mfl_subtree_reroot_traversal(cliprec.tgt2->nodet_next->nodet_edge, src, tgt, &cliprec, searchrec, false);
                mfl_subtree_reroot_traversal(cliprec.tgt2->nodet_next->nodet_next->nodet_edge, src, tgt, &cliprec, searchrec, false);
            } else {
                mfl_subtree_reroot_traversal(cliprec.tgt2, src, tgt, &cliprec, searchrec, neighbour_rule);
            }
        }
        
        mfl_restore_branching(&cliprec);
    }
}


void mfl_regrafting_traversal(mfl_node_t* tgt, mfl_node_t* src, mfl_searchrec_t* searchrec)
{
    // Traverse the target tree, attempting the reinsertion
    
    mfl_node_t* p = tgt->nodet_next;
    mfl_cliprec_t regraft;
    
    // Insert branch
    mfl_temp_rebranching(src, tgt, &regraft);
    
    // Copy the tree, append it to the buffer
    mfl_save_topology(searchrec->sr_swaping_on, searchrec->sr_treebuffer, searchrec);
    
    // Count the number of rearrangements
    mfl_undo_temp_rebranching(&regraft);
    
    
    ++searchrec->sr_rearrangement_counter;
    
    if (tgt->nodet_tip) {
        return;
    }
    
    if (p->nodet_edge) {
        mfl_regrafting_traversal(p->nodet_edge, src, searchrec);
    }
    p = p->nodet_next;
    
    if (p->nodet_edge) {
        mfl_regrafting_traversal(p->nodet_edge, src, searchrec);
    }
    
    
    return;
}


void mfl_regraft_subtree(mfl_node_t* src, mfl_node_t* tgt, mfl_searchrec_t* searchrec, bool neighbor_rule)
{
    mfl_node_t* tgt_opp = tgt;
    
    mfl_node_t* p = NULL;
    mfl_node_t* q = NULL;
    
    do {
        
        if (!tgt_opp->nodet_tip) {
            
            p = tgt_opp->nodet_next;
            
            
            if (neighbor_rule) {
                
                do {
                    
                    q = p->nodet_edge;
                    
                    if (!q->nodet_tip) {
                        
                        q = q->nodet_next;
                        
                        do {
                            mfl_regrafting_traversal(q->nodet_edge, src, searchrec);
                            q = q->nodet_next;
                        } while (q != p->nodet_edge);
                        
                    }
                    
                    p = p->nodet_next;
                    
                } while (p != tgt_opp);
                
            } else {
                
                do {
                    mfl_regrafting_traversal(p->nodet_edge, src, searchrec);
                    p = p->nodet_next;
                } while (p != tgt_opp);
                
            }
            
        }
        
        tgt_opp = tgt_opp->nodet_edge;
        
    } while (tgt_opp != tgt);
    
}


void mfl_pruning_traversal(mfl_node_t* n, mfl_searchrec_t* searchrec)
{
    mfl_node_t* p = NULL;
    mfl_node_t* q = NULL;
    mfl_cliprec_t cliprec;
    bool neighbour_rule = false;
    
    if (n->nodet_tip) {
        return;
    }
    
    q = n;
    
    do {
        
        q = q->nodet_next;
        
        if (q == n) {
            
            p = mfl_clip_branch(n->nodet_edge , &cliprec);
            neighbour_rule = false;
            
        }
        else {
            mfl_pruning_traversal(q->nodet_edge, searchrec);
            p = mfl_clip_branch(q->nodet_edge, &cliprec);
            neighbour_rule = true;
        }
    
        if (searchrec->sr_bswaptype == MFL_BST_SPR) {
            // SPR
            mfl_regraft_subtree(p, cliprec.tgt1, searchrec, neighbour_rule);
        }
        else if (searchrec->sr_bswaptype == MFL_BST_TBR) {
            // TBR
            if (q != n) {
                mfl_reroot_subtree(p, cliprec.tgt1, searchrec, neighbour_rule);
            }
            else if (n->nodet_edge->nodet_tip) {
                mfl_reroot_subtree(p, cliprec.tgt1, searchrec, neighbour_rule);
            }
            else if (n->nodet_edge->nodet_next->nodet_edge->nodet_tip && n->nodet_edge->nodet_next->nodet_next->nodet_edge->nodet_tip) {
                mfl_reroot_subtree(p, cliprec.tgt1, searchrec, neighbour_rule);
            }
        }
        else {
            dbg_eprintf("nice try, asshole");
        }
            
        mfl_restore_branching(&cliprec);
        
    } while (q != n);
    
    return;
}


void mfl_destroy_node_wt_list(mfl_nodeweights_t* nodewtlst)
{
    free(nodewtlst->nwl_weights_list);
    free(nodewtlst);
}

void tui_spr_test_environment(void)
{
    // NOTE: This is a tui function and should be moved there.
    
    int temp_max_trees = 50000;
    
    mfl_handle_s* testhandle = mfl_t2s(mfl_create_handle());
    
    char* cliptesttree = NULL;
    
    cliptesttree = "tree1=[&R] (1,(2,(((((((((((((((((((((3,39),12),(11,(53,64))),30),(42,62)),48),(25,32)),74),21),((((((6,61),76),17),67),(8,45)),((((((((((14,22),38),(16,18)),((37,58),75)),(59,73)),15),26),68),(51,56)),36))),((((13,(40,((46,55),54))),49),((((((29,34),(33,63)),72),57),65),35)),23)),70),44),27),(31,43)),(((9,((19,41),(20,28))),24),47)),71),((4,10),69)),((50,78),52)),7),(((5,66),77),60))));";
    //cliptesttree = "tree1=[&R] (1,((2,(79,(80,((81,((82,88),(86,87))),(83,(84,85)))))),(((((((((((((((((((((3,39),12),(11,(53,64))),30),(42,62)),48),(25,32)),74),21),((((((6,61),76),17),67),(8,45)),((((((((((14,22),38),(16,18)),((37,58),75)),(59,73)),15),26),68),(51,56)),36))),((((13,(40,((46,55),54))),49),((((((29,34),(33,63)),72),57),65),35)),23)),70),44),27),(31,43)),(((9,((19,41),(20,28))),24),47)),71),((4,10),69)),((50,78),52)),7),(((5,(66,((89,90),((91,((92,98),(96,97))),(93,(94,(95,(99,100)))))))),77),60))));";
    //cliptesttree = "temp_examp6=[&R] (((((1,4),5),3),2),(6,(7,((8,9),(10,(11,12))))));";
    //cliptesttree = "tree_857 = [&R] (1,((4,(2,((6,(7,((8,9),(10,(11,12))))),3))),5));";
    //cliptesttree = "tree_838 = [&R] ((((12,(((((3,2),6),7),(8,9)),10)),11),1),(4,5));";
    cliptesttree = "tree_795 = [&R] (1,(((((7,11),(10,((2,((4,8),6)),3))),5),9),12));";
    cliptesttree = "temp_examp6=[&R] (1,(4,(5,(3,(2,(6,(7,(8,(9,(10,(11,12)))))))))));";
    cliptesttree = "temp_examp6=[&R] (1,(4,(5,(3,(2,(6,(7,(8,(9,(10,(11,12)))))))))));";
    cliptesttree = "temp_examp6=[&R] (1,(4,(5,(3,(2,(6,(7,(8,(9,(10,(11,12)))))))))));";
    cliptesttree = "temp_examp6=[&R] (1,(2,(((7,(6,(5,(4,3)))),(((12,11),10),9)),8)));";
    //cliptesttree = "temp_examp6=[&R] ((1,2),(3,4));";
    //cliptesttree = "temp_examp6=[&R] ((1,2),3);";
    //cliptesttree = "temp_examp6=[&R] ((1,2),(3,(4,5)));";
    //cliptesttree = "temp_examp6=[&R] (1,(2,(3,(4,(5,6)))));";
    //cliptesttree = "temp_examp6=[&R] (5,(4,(3,(2,1))));";
    //char* cliptesttree = "temp_examp6=[&R] ((1,(2,(6,7))),(3,(4,5)));";
    //cliptesttree = "temp_examp6=[&R] ((1,((2,8),(6,7))),(3,(4,5)));";
    //char* cliptesttree = "UNTITLED = [&R] ((((((((1,2),3),((((4,5),(((((6,7),8),(9,10)),(11,12)),((13,(14,15)),(16,17)))),18),((19,20),21))),((((((22,23),(24,25)),((((((26,27),28),29),30),(31,(32,33))),((34,35),(((36,37),((38,39),40)),41)))),(42,((43,((44,45),46)),(47,(48,(49,50)))))),51),((((52,53),54),((55,56),57)),(58,((59,60),61))))),(((62,((63,64),(65,66))),(67,(68,(69,(70,71))))),(((72,73),((74,(((75,76),77),(78,79))),80)),81))),((((82,83),(((84,85),86),(87,(88,89)))),(90,91)),((92,93),((94,(95,((96,97),98))),(99,100))))),(((((((101,102),(103,(((104,105),106),107))),(((108,109),(110,111)),(((112,(113,114)),(115,116)),(117,(118,119))))),(((120,((121,((122,(((123,((124,125),126)),127),(128,129))),(130,(131,132)))),(((133,((134,135),(136,((137,138),139)))),(((140,141),142),(((143,144),145),((146,147),148)))),(149,(((150,151),152),153))))),((154,155),156)),((157,(158,(159,160))),((((161,162),163),(164,165)),((166,167),168))))),(((((169,(170,171)),172),((173,174),175)),(((176,177),(178,(179,((180,181),(182,183))))),184)),(((185,186),((187,(188,189)),190)),((((((191,192),193),(((194,195),196),(197,198))),(199,((200,201),(202,203)))),204),(((205,206),207),(208,(209,210))))))),(((((211,212),(((213,214),(215,216)),(((217,218),((219,220),221)),222))),((223,224),(((225,226),(227,228)),(((229,(230,((231,232),233))),(((234,235),(236,237)),238)),(((239,240),((241,242),243)),(244,245)))))),(((((246,247),248),249),250),251)),((((((252,253),(254,255)),((256,((257,258),259)),(260,261))),(((262,263),(264,((265,266),(267,268)))),((((269,(270,271)),((272,273),(((274,275),276),277))),((278,(279,(280,281))),(((282,283),(284,285)),286))),(287,288)))),(((289,((290,((291,292),293)),(294,295))),(296,(297,((298,299),(300,301))))),(((302,303),304),((305,((306,307),((308,(309,310)),(311,312)))),(((313,(314,315)),316),317))))),((318,319),(((320,321),(322,323)),((324,325),(326,327))))))),(328,((((329,330),331),332),(333,334))))),(((((335,336),((337,((338,339),(340,341))),((342,343),344))),(((((345,346),347),348),349),(350,(351,352)))),(((353,354),((((355,356),357),((358,(359,360)),(361,362))),(((363,((364,365),366)),(((367,368),369),(370,371))),((((372,(373,374)),(375,376)),(((377,(378,379)),(380,(381,(382,383)))),384)),(((385,386),387),((388,389),390)))))),((((((391,392),(393,((394,395),((396,397),398)))),(399,400)),(((401,(((402,403),404),(405,(406,407)))),((((408,409),(410,(411,(412,413)))),414),415)),(((416,((417,(418,419)),420)),(((421,422),423),(424,425))),(((426,427),428),((429,430),431))))),(((432,433),((434,435),((436,437),438))),((439,(440,441)),(442,(443,(((444,445),(446,447)),448)))))),((449,450),(451,452))))),(((((453,((454,455),456)),(457,458)),((((((459,460),461),((462,(((463,(464,465)),(466,467)),468)),(469,(470,(471,472))))),((473,(((474,475),476),477)),((478,(((479,480),481),482)),((((483,484),485),(486,487)),(((488,489),(490,491)),((492,493),((494,495),(496,497)))))))),(((498,(499,((500,501),502))),(((503,504),(505,506)),(((507,508),509),(((510,511),(512,513)),(514,515))))),((516,517),((518,519),(520,((521,522),523)))))),((524,(525,526)),((((527,((528,(529,530)),531)),532),(((533,(534,(535,(536,537)))),(((((538,539),540),541),((542,543),(((544,(545,546)),547),548))),549)),((((550,551),552),(553,(554,555))),((556,557),558)))),(((((559,560),561),(562,(563,564))),((565,566),567)),((568,569),(570,(571,(572,573))))))))),((((((574,((575,(576,577)),578)),(((579,580),581),582)),((583,584),(585,(586,587)))),(((((588,((589,(590,(591,592))),(593,(594,(595,596))))),(597,(598,599))),(((600,601),((((602,603),604),605),(606,607))),((608,(609,610)),611))),(((612,(613,614)),615),(616,(((617,(618,619)),(620,(621,(622,623)))),(((624,625),626),627))))),((628,(629,(630,631))),632))),((633,((634,635),(636,637))),(((((638,((639,640),641)),(((642,643),644),(((645,(646,647)),(648,649)),650))),(651,(((652,653),654),(655,656)))),(657,((658,659),660))),(661,662)))),((((((((663,664),665),(666,667)),(((668,(669,670)),(671,672)),(673,674))),(675,676)),((((((677,678),679),((680,((681,(682,683)),((684,685),686))),((687,(688,689)),(690,691)))),((692,693),694)),(((695,696),((697,(698,699)),(700,(701,702)))),(703,((704,(705,706)),707)))),708)),(((((((709,710),711),(((712,(713,714)),715),((716,(717,718)),(719,(720,721))))),(((722,(723,724)),((725,726),727)),728)),(((729,(730,(((731,(732,733)),734),735))),((736,((737,(738,739)),((740,741),742))),((743,(744,745)),((746,747),748)))),((((749,(750,751)),((752,753),(754,(755,(((756,757),(758,759)),760))))),761),(762,763)))),((((764,765),(766,767)),((((768,769),(770,771)),(772,773)),((774,775),776))),(((777,778),((779,((780,781),782)),(783,784))),((785,((786,787),((788,789),790))),((((791,(792,793)),794),((795,796),797)),(798,799)))))),((((((800,(801,802)),803),804),(805,((806,807),(808,809)))),((810,811),((812,((813,814),(815,(816,817)))),818))),(((819,820),((821,(822,(823,824))),(((825,826),827),(828,829)))),((((830,831),(832,(833,834))),((835,836),((837,(838,839)),(840,841)))),(842,843)))))),((844,845),((((846,((847,(848,849)),(((850,(851,852)),(853,854)),((855,((856,857),858)),859)))),(((860,((861,862),(((863,(864,(865,866))),867),(868,869)))),(((870,871),(872,(873,(874,(875,876))))),(877,878))),(879,880))),((((881,882),(883,(884,885))),((886,887),(888,889))),890)),((891,(892,893)),(((894,895),896),(((897,(898,899)),900),901)))))))),((((((902,903),904),905),(906,((907,(908,909)),((910,911),(912,(913,914)))))),(((((915,916),917),918),((919,(((920,(921,922)),(923,924)),925)),(926,(927,928)))),(929,(930,931)))),((((932,(933,(934,935))),936),((937,938),939)),((((940,941),(((942,943),((944,945),946)),(947,(948,((949,950),951))))),((((952,953),((954,(955,((956,957),(958,959)))),960)),(961,962)),(((963,964),(965,966)),((((967,968),(969,970)),(971,972)),973)))),((((974,975),(976,977)),(((978,979),((980,981),(982,983))),((984,(985,986)),((987,((988,989),(990,991))),(992,993))))),((994,995),(996,((997,(998,999)),1000))))))))));";
    
    
    mfl_tree_t* testree = mfl_convert_newick_to_mfl_tree_t(cliptesttree, 0);
    
    char* treedraw = mfl_drawtree(testree);
    dbg_printf("%s", treedraw);
    
    testhandle->n_taxa = testree->treet_num_taxa;
    
    mfl_searchrec_t* searchrec = mfl_create_searchrec(testhandle);
    
    mfl_unroot_tree(testree);
    
    if (testree->treet_root) {
        searchrec->sr_swap_entry = testree->treet_root;
    } else {
        searchrec->sr_swap_entry = testree->treet_start;
    }
    
    searchrec->sr_treebuffer = mfl_alloc_treebuffer(temp_max_trees);
    searchrec->sr_swaping_on = testree;
    searchrec->sr_handle_ptr = testhandle;
    searchrec->sr_nodweights = mfl_create_node_wt_list(testree->treet_num_taxa);
    searchrec->sr_bswaptype = MFL_BST_TBR;
    
    mfl_append_tree_to_treebuffer(testree, searchrec->sr_treebuffer, testhandle);
    
    time_t before = 0;
    time_t after = 0;
    
    before = clock();
    mfl_pruning_traversal(searchrec->sr_swap_entry, searchrec);
    after = clock();
    
    dbg_printf("\nNumber of rearrangements attempted for %i taxa: %lli\n", testree->treet_num_taxa, searchrec->sr_rearrangement_counter);
//    dbg_printf("\nAnd here are the trees:\n");
    dbg_printf("Time: %f\n\n", (float)(after-before)/CLOCKS_PER_SEC);

    char* newicktr = NULL;
    int i = 0;
    for (i = 0; i < searchrec->sr_treebuffer->tb_num_trees; ++i) {
        newicktr = mfl_convert_mfl_tree_t_to_newick(searchrec->sr_treebuffer->tb_savedtrees[i], 0);
        dbg_printf("TREE tree_%i = ", (i+1));
        dbg_printf("%s\n", newicktr);
        free(newicktr);
    }
    
    mfl_free_tree(testree);

}



bool mfl_heuristic_search(mfl_handle_s *mfl_handle)
{
    /* Eventually, this will parse the information in the handle and set up the 
     * heuristic search. For now, it just returns zero because we're not even 
     * close to that yet. */
    return 0;
}