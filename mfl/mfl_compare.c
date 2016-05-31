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
    
    new_edgetable->num_tips = num_tips;
    
    //add the number of entries
    new_edgetable->numentries = ( 2 * num_tips - 2);
    if(!is_rooted) {
        new_edgetable->num_nodes = num_tips - 2;
    } else {
        new_edgetable->num_nodes = num_tips - 1;
    }
    
    //malloc the size of the edge table
    new_edgetable->edgetable = (int*)mfl_malloc(new_edgetable->numentries * sizeof(int), 0);
    
    //malloc the tips/nodes addresses
    new_edgetable->edge_tips = (mfl_node_t**)mfl_malloc(new_edgetable->num_tips * sizeof(mfl_node_t*), 0);
    new_edgetable->edge_nodes = (mfl_node_t**)mfl_malloc(new_edgetable->num_nodes * sizeof(mfl_node_t*), 0);

    return new_edgetable;
}

// Add the links to the tips
void mfl_get_edgetable_tips(mfl_edgetable_t* edgetable, mfl_tree_t* tree)
{
    int i = 0;
    
    for (i = 0; i < edgetable->num_tips; ++i) {
        edgetable->edge_tips[i] = tree->treet_treenodes[i];
    }
}

//Destroys a mfl_edgetable_t object
void mfl_destroy_edgetable(mfl_edgetable_t* edgetable)
{
    if (edgetable->edgetable) {
        free(edgetable->edgetable);
    }
    if (edgetable->edge_tips) {
        free(edgetable->edge_tips);
    }
    if (edgetable->edge_nodes) {
        free(edgetable->edge_nodes);
    }
    free(edgetable);
}

/*!
 @description finds the bottom node in a node ring.
 @param node a mfl_node_t pointer to a node in a node ring.
 @return the pointer to the bottom node
 */
mfl_node_t* mfl_find_bottom_node_in_ring(mfl_node_t* node)
{
    mfl_node_t* node_entry = node;
    
    if(mfl_check_node_is_bottom(node)) {
        return node;
    } else {
        while (!mfl_check_node_is_bottom(node)) {
            node = node->nodet_next;
            
            assert(node != node_entry); //Break if goes into an infinite loop!
        }
        return node;
    }
}

/*!
 @description adds a nodet_edge_ref to a bottom node.
 @param node a mfl_node_t pointer to a bottom node in a node ring.
 @param reference an int to be used as the node edge reference
 */
void mfl_set_edge_ref_in_ring(mfl_node_t* node, int reference)
{
    mfl_node_t* node_bottom = NULL;
    node_bottom = mfl_find_bottom_node_in_ring(node);
    node_bottom->nodet_edge_ref = reference;
}

/*!
 @description get the nodet_edge_ref of a bottom node.
 @param node a mfl_node_t pointer to a bottom node in a node ring.
 @return the nodet_edge_ref (int).
 */
int mfl_get_edge_ref_from_ring(mfl_node_t* node)
{
    int node_reference = NULL;
    mfl_node_t* node_bottom = NULL;
    node_bottom = mfl_find_bottom_node_in_ring(node);
    node_reference = node_bottom->nodet_edge_ref;
    return node_reference;
}


void mfl_add_nodesref_traversal(mfl_node_t* start, int* node_counter, mfl_edgetable_t* edgetable)
{
    int num_tips = edgetable->num_tips;
    mfl_node_t* node = start;
    if (node->nodet_tip != 0) {
        return;
    }
    node = node->nodet_next;
    
    do {
        
        if(!node->nodet_edge) {
            mfl_add_nodesref_traversal(node->nodet_next->nodet_edge, node_counter, edgetable);
        } else {
            mfl_add_nodesref_traversal(node->nodet_edge, node_counter, edgetable);
        }
        
        if(mfl_get_edge_ref_from_ring(node) == 0) {
            ++*node_counter;
            mfl_set_edge_ref_in_ring(node, *node_counter);
            edgetable->edge_nodes[*node_counter - edgetable->num_tips] = mfl_find_bottom_node_in_ring(node);
        }
        
        node = node->nodet_next;
    } while (node != start);

}

void mfl_get_edgetable(mfl_edgetable_t* edgetable, mfl_tree_t* tree)
{
    // variables before the traversal
    mfl_node_t* dummy_test_node = NULL;
    int node_counter = edgetable->num_tips;
    int node_reference = 0;
    int counter = 0;
    mfl_node_t* current_node = NULL;
    
    // Set the tips in the edgetable
    mfl_get_edgetable_tips(edgetable, tree);
    
    // Do the first tip
    current_node = edgetable->edge_tips[counter];
    mfl_set_edge_ref_in_ring(current_node->nodet_edge, node_counter);
    edgetable->edgetable[counter] = node_counter;
    edgetable->edge_nodes[counter] = mfl_find_bottom_node_in_ring(current_node->nodet_edge);
    ++counter;
    
    // Loop through the other tips
    for (counter; counter < edgetable->num_tips; ++ counter) {
        current_node = edgetable->edge_tips[counter];
        node_reference = mfl_get_edge_ref_from_ring(current_node->nodet_edge);
        
        if(node_reference != 0) {
            edgetable->edgetable[counter] = node_reference;
        } else {
            mfl_add_nodesref_traversal(current_node->nodet_edge, &node_counter, edgetable);
            node_reference = mfl_get_edge_ref_from_ring(current_node->nodet_edge);
            edgetable->edgetable[counter] = node_reference;
        }
    }
    
    //TODO: There's a weird behaviour of the edgetable->edge_tips that gets filled past his capacity
    //Maybe link to the general problem of mallocing node_t (mallocs always more!)
    

    // Then loop through the nodes
    for (counter; counter < edgetable->numentries ; ++ counter) {
        current_node = edgetable->edge_nodes[counter-edgetable->num_tips];
        // Condition if node is the root
        if(!current_node->nodet_edge) {
            node_reference = mfl_get_edge_ref_from_ring(current_node);
        } else {
            // Condition if node is not the root
            if(!current_node->nodet_edge->nodet_tip) {
                node_reference = mfl_get_edge_ref_from_ring(current_node->nodet_edge);
            } else {
                // Conditions if tree is unrooted
                if(!current_node->nodet_next->nodet_edge->nodet_tip) {
                    node_reference = mfl_get_edge_ref_from_ring(current_node->nodet_next->nodet_edge);
                } else {
                    node_reference = mfl_get_edge_ref_from_ring(current_node->nodet_next->nodet_next->nodet_edge);
                }
            }
        }
        edgetable->edgetable[counter] = node_reference;
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
    
    //cliptesttree = "temp_examp6=[&R] ((1,2),(3,4));";
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
