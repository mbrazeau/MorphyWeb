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

void mfl_get_edgetable(mfl_edgetable_t* edgetable, mfl_tree_t* tree)
{
    // variables before the traversal
    int node_counter = edgetable->num_tips;
    int node_reference = 0;
    int counter = 0;
    mfl_node_t* current_node = NULL;
    mfl_node_t* target_node = NULL;
    
    // Set the tips in the edgetable
    mfl_get_edgetable_tips(edgetable, tree);

    // Go through the tips
    for (counter; counter < edgetable->num_tips; ++counter) {
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
    
    // Go through the nodes
    for (counter; counter < edgetable->numentries ; ++ counter) {
        current_node = edgetable->edge_nodes[counter - edgetable->num_tips];
        
        // Select the node to target (usually node_edge but not if root is current node or for some nodes in unrooted trees)
        if(!current_node->nodet_edge) {
            // Edge is null, root points to itself
            target_node = current_node;
        } else {
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
    
    //cliptesttree = "temp_examp6=[&U] ((1,2),(3,4));";
    //cliptesttree = "temp_examp6=[&U] ((1,2),(3,(4,5)));";
    //cliptesttree = "temp_examp6=[&U] (1,(2,(3,(4,(5,6)))));";
    //cliptesttree = "temp_examp6=[&U] (5,(4,(3,(2,1))));";
    //cliptesttree = "temp_examp6=[&U] ((1,(2,(6,7))),(3,(4,5)));";
    //cliptesttree = "equal_test=[&U] ((1,(2,3)), (4,(5,6)));";
    //cliptesttree = "equal_test=[&U] ((4,(5,6)), (1,(2,3)));";
    //cliptesttree = "equal_test=[&U] (2, ((4,7), ((1,(3,5)), (8,(6,9)))));";
    //cliptesttree = "equal_test=[&U] (2, ((4,7), ((8,(6,9)), (1,(3,5)))));";
    //cliptesttree = "equal_test=[&U] (2, (((8,(6,9)), (1,(3,5))), (4,7)));";
    cliptesttree = "tree1=[&U] (1,(2,(((((((((((((((((((((3,39),12),(11,(53,64))),30),(42,62)),48),(25,32)),74),21),((((((6,61),76),17),67),(8,45)),((((((((((14,22),38),(16,18)),((37,58),75)),(59,73)),15),26),68),(51,56)),36))),((((13,(40,((46,55),54))),49),((((((29,34),(33,63)),72),57),65),35)),23)),70),44),27),(31,43)),(((9,((19,41),(20,28))),24),47)),71),((4,10),69)),((50,78),52)),7),(((5,66),77),60))));";
    
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
    free(test_edgetable);
    // Destroy the table
//    mfl_destroy_edgetable(test_edgetable);
}
