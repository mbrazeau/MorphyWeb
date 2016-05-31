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
mfl_edgetable_t* mfl_initiate_edgetable_t(int num_tips, bool is_rooted)
{
    //Intialise the table
    mfl_edgetable_t* new_edgetable = NULL;
    
    //malloc the edge table
    new_edgetable = (mfl_edgetable_t*)mfl_malloc(sizeof(mfl_edgetable_t), 0);
    
    //add the number of entries
    if(!is_rooted) {
        new_edgetable->numentries = ( 2 * num_tips - 3); // The last node is ignored
    } else {
        new_edgetable->numentries = ( 2 * num_tips - 2); // The last node is ignored
    }
    
    //malloc the size of the edge table
    new_edgetable->edgetable = (int*)mfl_malloc(new_edgetable->numentries * sizeof(int), 0);
    
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

//Sets the edge reference in a node ring (arbitrarilly on the bottom node)
void mfl_set_edge_ref_in_ring(mfl_node_t* node, int reference)
{
    mfl_node_t* node_entry = node; //TG: minor suggestion here is to put this variable into the while loop. The idea is because this function will be called a lot, it'll save looking for the bottom 1/3rd of the time.
    
    
    if(mfl_check_node_is_bottom(node)) {
        node->nodet_edge_ref = reference;
    } else {
        // Search for the bottom node
        while (!mfl_check_node_is_bottom(node)) {
            node = node->nodet_next;
            //Stop if back to node entry!
            assert(node != node_entry); //TG: thought using an assert here might be good.
        }
        node->nodet_edge_ref = reference;
    }
}

//Getting the node edge ref value node ring (stored in bottom node!)
int mfl_get_edge_ref_from_ring(mfl_node_t* node)
{
    int node_reference = NULL;
    mfl_node_t* node_entry = node; //TG: minor suggestion here is to put this variable into the while loop. The idea is because this function will be called a lot, it'll save looking for the bottom 1/3rd of the time.
    
    if(mfl_check_node_is_bottom(node)) {
        node_reference = node->nodet_edge_ref;
        return node_reference;
    } else {
        // Search for the bottom node
        while (!mfl_check_node_is_bottom(node)) {
            node = node->nodet_next;
            //Stop if back to node entry!
            assert(node != node_entry); //TG: thought using an assert here might be good.
        }
        node_reference = node->nodet_edge_ref;
        return node_reference;
    }
}


void mfl_add_nodesref_traversal(mfl_node_t* start, int* node_counter)
{
    mfl_node_t *node = start;
    if (node->nodet_tip != 0) {
        return;
    }
    node = node->nodet_next;
    
    do {
        
        if(!node->nodet_edge) {
            mfl_add_nodesref_traversal(node->nodet_next->nodet_edge, node_counter);
        } else {
            mfl_add_nodesref_traversal(node->nodet_edge, node_counter);
        }
        
        if(mfl_get_edge_ref_from_ring(node) == 0) {
            ++*node_counter;
            mfl_set_edge_ref_in_ring(node, *node_counter);
        }
        
        node = node->nodet_next;
    } while (node != start);

}

void mfl_get_edgetable(mfl_edgetable_t* edgetable, mfl_tree_t* tree)
{
    // variables before the traversal
    int num_taxa = tree->treet_num_taxa;//TODO: Should be num taxa active!
    int node_counter = num_taxa;
    int node_reference = 0;
    int counter = 0;
    mfl_node_t* current_node = NULL;

    // declare the first node with an entry at tip 1
    current_node = tree->treet_treenodes[counter];
    
    // set the edge reference in the node ring
    mfl_set_edge_ref_in_ring(current_node->nodet_edge, node_counter);
    
    // store this value in the edgetable
    edgetable->edgetable[counter] = node_counter;
    
    //increment the counter
    ++counter;
    
    for (counter; counter < edgetable->numentries; ++counter) {
        // Get the next node in the array
        current_node = tree->treet_treenodes[counter];
        // Get the node ring reference
        if(!current_node->nodet_tip) {
            node_reference = mfl_get_edge_ref_from_ring(current_node);
        } else {
            node_reference = mfl_get_edge_ref_from_ring(current_node->nodet_edge);
        }
        
        if(node_reference != 0){
            edgetable->edgetable[counter] = node_reference;
        } else {
            // name the nodes in a post order traversal
            if(!current_node->nodet_tip) {
                mfl_add_nodesref_traversal(current_node, &node_counter);
            } else {
                mfl_add_nodesref_traversal(current_node->nodet_edge, &node_counter);
            }
            // Get the node ring reference and store in edgetable
            if(!current_node->nodet_tip) {
                node_reference = mfl_get_edge_ref_from_ring(current_node);
            } else {
                node_reference = mfl_get_edge_ref_from_ring(current_node->nodet_edge);
            }
            edgetable->edgetable[counter] = node_reference;
        }
    }
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

void tui_print_edgetable(mfl_edgetable_t* edgetable)
{
    int i = 0;
    dbg_printf("Tip/node connects to tip/node\n");
    
    for(i = 0; i < edgetable->numentries; ++i) {
        dbg_printf("%i connects to %i\n", i, edgetable->edgetable[i]);
    }
}



void tui_test_edgetables(void)
{
    char* cliptesttree = NULL;
    
    cliptesttree = "temp_examp6=[&R] ((1,2),(3,4));";
    //cliptesttree = "temp_examp6=[&R] ((1,2),3);";
    //cliptesttree = "temp_examp6=[&R] ((1,2),(3,(4,5)));";
    //cliptesttree = "temp_examp6=[&R] (1,(2,(3,(4,(5,6)))));";
    //cliptesttree = "temp_examp6=[&R] (5,(4,(3,(2,1))));";
    //cliptesttree = "temp_examp6=[&R] ((1,(2,(6,7))),(3,(4,5)));";
    
    mfl_edgetable_t* test_edgetable = NULL;
    mfl_tree_t* testree = mfl_convert_newick_to_mfl_tree_t(cliptesttree, 0);
    
    if(!testree->treet_root) {
        mfl_assign_bottom_node(testree->treet_start);
        test_edgetable = mfl_initiate_edgetable_t(testree->treet_num_taxa, 0);
    } else {
        mfl_assign_bottom_node(testree->treet_root);
        test_edgetable = mfl_initiate_edgetable_t(testree->treet_num_taxa, 1);
    }
    
    mfl_get_edgetable(test_edgetable, testree);

    tui_print_edgetable(test_edgetable);
//    
//    free(test_edgetable);
    // Destroy the table
//    mfl_destroy_edgetable(test_edgetable);
}
