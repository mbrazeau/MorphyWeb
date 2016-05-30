//
//  mfl_compare.c
//  Morphy
//
//  Created by mbrazeau on 05/05/2016.
//
//

#include "morphy.h"

void mfl_sort_bipartition_set(mfl_partition_set_t* biparts)
{
    
}


int mfl_binary_search_bipartition(mfl_bitset_t* bipartition, mfl_partition_set_t* targetlist)
{
    return 0;
}

int mfl_compare_all_bipartitions(mfl_partition_set_t* partset1, mfl_partition_set_t* partset2)
{
    
}


mfl_partition_set_t* mfl_create_bipartition_set(int num_taxa)
{
    mfl_partition_set_t* newbipartset = (mfl_partition_set_t*)malloc(sizeof(mfl_partition_set_t));
    if (!newbipartset) {
        dbg_eprintf("unable to allocate memory for new bipartition set");
        return NULL;
    }
    else {
        memset(newbipartset, 0, sizeof(mfl_partition_set_t));
    }
    
}

void mfl_destroy_bipartition_set(mfl_partition_set_t* bptset)
{
    
}


/* This is mostly a temporary function, as it is likely that the setting of 
 * bipartitions will be handled simultaneously by other tree traversals. 
 */
void mfl_set_bipartitions(mfl_node_t* n)
{
    assert(n->nodet_bipart->bts_bitfields);
    
    mfl_node_t* p = NULL;
    
    if (n->nodet_tip) {
        return;
    }
    
    p = n->nodet_next;
    do {
        mfl_set_bipartitions(p->nodet_edge);
        mfl_bts_OR(p->nodet_edge->nodet_bipart, n->nodet_bipart, n->nodet_bipart);
        p = p->nodet_next;
    } while (p != n);
    
}


//Allocate the memory to the edgetable
mfl_edgetable_t* mfl_initiate_edgetable_t(int num_tips)
{
    int i = 0;
    
    //Intialise the table
    mfl_edgetable_t* new_edgetable = NULL;
    
    //malloc the edge table
    new_edgetable = (mfl_edgetable_t*)mfl_malloc(sizeof(mfl_edgetable_t), 0);
    
    //add the number of entries
    new_edgetable->numentries = ( 2 * num_tips - 1); //TG: if tree is unrooted numentries should be 2*num_tips - 1 (this generates a leak of 1 integer).
    
    //malloc the size of the edge table
    new_edgetable->edgetable = (int*)mfl_malloc(new_edgetable->numentries * sizeof(int), 0);
    
//    for (i = 0; i < new_edgetable->numentries; ++i) {
//        new_edgetable->edgetable[i] = i;
//    }
    
    return new_edgetable;
}

//Destroys a mfl_edgetable_t object
void mfl_destroy_edgetable(mfl_edgetable_t* edgetable)
{
    if (edgetable->edgetable) {
        free(edgetable->edgetable);
    }
    free(edgetable);
}

void mfl_edgetable_traversal(mfl_node_t* start, mfl_edgetable_t* edgetable, int* node_counter, int* counter)
{
    
    mfl_node_t* node = start;
    
    if(node->nodet_tip != 0) {
        // When node is a tip, always increment tip counter and print the tip.
        edgetable->edgetable[*counter] = node->nodet_tip;
    } else {
        
        // Increment node counter only when next is not a tip
        if(node->nodet_next->nodet_tip ) {
            ++*node_counter;
        }
        
        edgetable->edgetable[*counter] = *node_counter;
    }
    
    ++*counter;
    node = node->nodet_next;
    
    do {
        mfl_edgetable_traversal(node->nodet_edge, edgetable, node_counter, counter);
        
        node = node->nodet_next;
        
    } while (node != start);
    
    
}

void mfl_get_edgetable(mfl_edgetable_t* edgetable, mfl_tree_t* tree)
{
    // variables before the traversal
    int counter = 0;
    int node_counter = tree->treet_num_taxa;     // node count starts at number of tips + 1
    
    mfl_edgetable_traversal(tree->treet_start, edgetable, &counter, &node_counter);
    
}


// Function idea
bool mfl_compare_edge_tables(mfl_edgetable_t* t1, mfl_edgetable_t* t2)
{
    if (t1->numentries != t2->numentries) {
        return false;
    }
    
    if (!memcmp(t1, t2, t1->numentries * sizeof(int))) {
        return true; // memcmp returns 0 if there's a match
    }
    else {
        return false;
    }
}


void tui_test_edgetables(void)
{
    char* cliptesttree = NULL;
    
    //cliptesttree = "temp_examp6=[&R] ((1,2),(3,4));";
    cliptesttree = "temp_examp6=[&R] ((1,2),3);";
    //cliptesttree = "temp_examp6=[&R] ((1,2),(3,(4,5)));";
    //cliptesttree = "temp_examp6=[&R] (1,(2,(3,(4,(5,6)))));";
    //cliptesttree = "temp_examp6=[&R] (5,(4,(3,(2,1))));";
    //cliptesttree = "temp_examp6=[&R] ((1,(2,(6,7))),(3,(4,5)));";
    
    mfl_tree_t* testree = mfl_convert_newick_to_mfl_tree_t(cliptesttree, 0);
    
    // Allocate the table
    mfl_edgetable_t* test_edgetable = mfl_initiate_edgetable_t(testree->treet_num_taxa);
    
    test_edgetable->edgetable[0] = 8;
    test_edgetable->edgetable[1] = 7;
    test_edgetable->edgetable[3] = 9;
    test_edgetable->edgetable[4] = 1;
    test_edgetable->edgetable[5] = 6;
    test_edgetable->edgetable[6] = 5;
    test_edgetable->edgetable[7] = 1;
    test_edgetable->edgetable[8] = 2;
    
    int i = 0;
    int j = 9;
    
    for (i = 0; i < j; ++i) {
        dbg_printf("%i is %i\n", i, test_edgetable->edgetable[i]);
    }
    
    //mfl_get_edgetable(test_edgetable, testree);
    
    // Destroy the table
//    mfl_destroy_edgetable(test_edgetable);
}
