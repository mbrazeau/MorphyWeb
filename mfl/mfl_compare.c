//
//  mfl_compare.c
//  Morphy
//
//  Created by mbrazeau on 05/05/2016.
//
//

#include "morphy.h"


// Create an empty bipartition table
mfl_bipartition_table* mfl_initialise_bipartition_table(void) {
    
    //Intialise the table
    mfl_bipartition_table* new_bipartition_table = NULL;
    
    //malloc the edge table
    new_bipartition_table = (mfl_bipartition_table*)mfl_malloc(sizeof(mfl_bipartition_table), 0);
    
    return new_bipartition_table;
}

// Destroy a bipartition table
void mfl_destroy_bipartition_table(mfl_bipartition_table* bipartition_table)
{
    if (bipartition_table->bipartition_occurence_counter) {
        bipartition_table->bipartition_occurence_counter = NULL;
        free(bipartition_table->bipartition_occurence_counter);
    }
    if (bipartition_table->bipartitions) {
        bipartition_table->bipartitions = NULL;
        free(bipartition_table->bipartitions);
    }
    free(bipartition_table);
}

// Appends the malloc for a bipartition table
void mfl_append_malloc_bipartition_table(mfl_bipartition_table* bipartition_table)
{
    //Append memory
    bipartition_table->bipartition_occurence_counter = (int*)realloc(bipartition_table->bipartition_occurence_counter, bipartition_table->number_of_bipartitions+1 * sizeof(int));
    bipartition_table->bipartitions = (int*)realloc(bipartition_table->bipartitions, bipartition_table->number_of_bipartitions+1 * sizeof(int));
    //Set new values to 0
    bipartition_table->bipartition_occurence_counter[bipartition_table->number_of_bipartitions+1] = 0;
    bipartition_table->bipartitions[bipartition_table->number_of_bipartitions+1] = 0;
}

/*!
 @description Gets the bipartition bitfield value.
 @param node a mfl_node_t pointer to a node in a node ring.
 @return and int that is the bit value of the bipartition
 */
int mfl_get_node_bipartition(mfl_node_t* node)
{
    // Return some bipartition integer (using bitwise business)
    int bipartition = 0;
    int i = 0;
    int j = 0;
    int bit_position = 0;
    mfl_bitfield_t bitfield = 1;
    
    //TODO: not sure if the i is mandatory + j seems to loop way to far
    for (i = 0; i < node->nodet_bipart->bts_nfields; ++i) {
        for (j = 0; j < MORPHY_BTS_IN_BITSET; ++j) {
            bit_position = node->nodet_bipart->bts_bitfields[i] & (bitfield << j);
            if (bit_position) {
                bipartition ^= (-bit_position ^ bipartition) & (1 << j);
            }
        }
    }
    dbg_printf("Bipartition = %i\n", bipartition);
    return bipartition;
}

//Matching a bipartition value in the list of bipartitions
int mfl_match_bipartition(int bipartition, mfl_bipartition_table* bipartition_table)
{
    int i = 0;
    
    //Loop through the recorded biparititions
    for (i = 0; i < bipartition_table->number_of_bipartitions; ++i){
        if(bipartition == bipartition_table->bipartitions[i]) {
            //Return the bipartition number
            return i;
        }
    }
    
    return -1;
}

//Traversal for getting all the bipartitions
void mfl_get_bipartition_traversal(mfl_node_t* node, mfl_bipartition_table* bipartition_table)
{
    int current_bipartition = 0;
    int current_bipartition_position = -1; // Initialised to be -1 (no position; c.f. 0 that is the first position)
    mfl_node_t* start = NULL;
    
    if (node->nodet_tip) {
        return;
    }
    
    start = node->nodet_next;
    
    do {
        mfl_get_bipartition_traversal(start->nodet_edge, bipartition_table);
        
        start = start->nodet_next;
        
    } while (start != node);

    current_bipartition = mfl_get_node_bipartition(node);
    //Get the current bipartition position
    current_bipartition_position = mfl_match_bipartition(current_bipartition, bipartition_table);
    //Increment the biparitition table
    if(current_bipartition_position != -1){
        // Increment the occurence of this bipartition
        ++bipartition_table->bipartition_occurence_counter[current_bipartition_position];
    } else {
        // Append the bipartition table size
        mfl_append_malloc_bipartition_table(bipartition_table);
        // Add the bipartition to the table
        bipartition_table->bipartitions[bipartition_table->number_of_bipartitions] = current_bipartition;
        // Increment the occurence of this bipartition
        ++bipartition_table->bipartition_occurence_counter[bipartition_table->number_of_bipartitions];
        // Increment the total number of bipartitions
        ++bipartition_table->number_of_bipartitions;
    }
    
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
    if(!is_rooted) {
        new_edgetable->numentries = ( 2 * num_tips - 2);
        new_edgetable->num_nodes = num_tips - 2;
    } else {
        new_edgetable->numentries = ( 2 * num_tips - 1);
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
    }
    
    do {
        node = node->nodet_next;
        assert(node != node_entry);
    } while (!mfl_check_node_is_bottom(node));
    
    return node;
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

void mfl_get_edgetable(mfl_edgetable_t* edgetable, mfl_tree_t* tree)
{
    // variables before the traversal
    int node_counter = edgetable->num_tips;
    assert(node_counter == tree->treet_num_taxa);
    int node_reference = 0;
    int counter = 0;
    mfl_node_t* current_node = NULL;
    mfl_node_t* target_node = NULL;
    
    // Set the tips in the edgetable
    mfl_get_edgetable_tips(edgetable, tree);

    // Go through the tips
    for (counter = 0; counter < edgetable->num_tips; ++counter) {
        current_node = edgetable->edge_tips[counter];
        node_reference = mfl_get_edge_ref_from_ring(current_node->nodet_edge);
        
        if(node_reference != 0) {
            edgetable->edgetable[counter] = node_reference;
        } else {
            // Set the node number linking to the tip
            mfl_set_edge_ref_in_ring(current_node->nodet_edge, node_counter);
            edgetable->edgetable[counter] = node_counter;
            edgetable->edge_nodes[node_counter - edgetable->num_tips] = mfl_find_bottom_node_in_ring(current_node->nodet_edge);
            ++node_counter;
        }
    }
    
    assert(counter);
    // Go through the nodes
    for (; counter < edgetable->numentries ; ++ counter) {
        current_node = edgetable->edge_nodes[counter - edgetable->num_tips];
        
        // Select the node to target (usually node_edge but not if root is current node or for some nodes in unrooted trees)
        if(!current_node->nodet_edge) {
            // Edge is null, root points to itself
            target_node = current_node;
        } else {
            // TODO: This can be done as a loop; safer for non-binary nodes
            if(!current_node->nodet_edge->nodet_tip) {
                target_node = current_node->nodet_edge;
            } else {
                if(!current_node->nodet_next->nodet_edge->nodet_tip) {
                    target_node = current_node->nodet_next->nodet_edge;
                } else {
                    target_node = current_node->nodet_next->nodet_next->nodet_edge;
                }
            }
        }
        
        node_reference = mfl_get_edge_ref_from_ring(target_node);
        
        
        if(node_reference != 0) {
            edgetable->edgetable[counter] = node_reference;
        } else {
            // Set the node number linking tot the node
            mfl_set_edge_ref_in_ring(target_node, node_counter);
            edgetable->edgetable[counter] = node_counter;
            edgetable->edge_nodes[node_counter - edgetable->num_tips] = mfl_find_bottom_node_in_ring(target_node);
            ++node_counter;
        }
    }
}


// Comparing edgetables
bool mfl_compare_edge_tables(mfl_edgetable_t* t1, mfl_edgetable_t* t2)
{
    if (t1->numentries != t2->numentries) {
        return false;
    }
    
    if (!memcmp(t1->edgetable, t2->edgetable, t1->numentries * sizeof(int))) {
        return true; // memcmp returns 0 if there's a match
    }
    else {
        return false;
    }
}
