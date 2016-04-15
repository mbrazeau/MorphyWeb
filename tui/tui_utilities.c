//
//  tui_utilities.c
//  Morphy
//
//  Created by mbrazeau on 12/04/2016.
//
//

/*  In this file:
 *  General-purpose testing and output functions for Morphy.
 *
 */

#include "morphy.h"
#include "tuimfy.h"


void tui_print_node_data(mfl_node_t* p, const char *calling_fxn)
{
    dbg_printf("Data for node @: %p\n", p);
    
    dbg_printf("\tp->nodet_edge: %p\n", p->nodet_edge);
    dbg_printf("\tp->nodet_next: %p\n", p->nodet_next);
    dbg_printf("\tp->nodet_tip: %i\n", p->nodet_tip);
    dbg_printf("\tp->nodet_index: %i\n", p->nodet_index);
    dbg_printf("\tp->nodet_isingroup %i\n", p->nodet_isingroup);
    dbg_printf("\tp->nodet_isbottom %i\n", p->nodet_isbottom);
    dbg_printf("\tp->nodet_downpass_visited: %i\n", p->nodet_downpass_visited);
    dbg_printf("\tp->nodet_uppass_visited: %i\n", p->nodet_uppass_visited);
    //dbg_printf("\t\n");
    //dbg_printf("\t\n");
    dbg_printf("\tp->nodet_isroot: %i\n", p->nodet_isroot);
    //dbg_printf("\t\n");
    dbg_printf("\t\t(%s() called by: %s() )\n", __FXN_NAME__, calling_fxn);
    
}

int tui_tip_check(mfl_node_t* n, const char* calling_fxn, const int* verbose)
{
    if (verbose) {
        tui_print_node_data(n, __FXN_NAME__);
    }
    
    if (!n->nodet_tip) {
        dbg_eprintf("queried node is not a tip\n");
        dbg_pcall(calling_fxn);
        return 1;
    }
    else {
        /* Perform other checks on tip values.*/
        return 0;
    }
    
}


int tui_check_node_is_root(mfl_node_t* p, const int* verbose, const char* calling_fxn)
{
    
    bool allchecks = true;
    
    if (p->nodet_isroot) {
        dbg_printf("Node %p called by %s() is declared as root\n\n", p, calling_fxn);
    }
    else {
        dbg_printf("Node %p called by %s() is NOT declared as root\n\n", p, calling_fxn);
        allchecks = false;
    }
    
    if (*verbose) {
        dbg_printf("From %s()\n", __FXN_NAME__);
        tui_print_node_data(p, __FXN_NAME__);
        dbg_printf("\n\n");
    }
    
    if (allchecks) {
        return 0;
    }
    else {
        return 1;
    }
}


void* tui_check_binary_traversal(mfl_node_t *p, const int* verbose, const char* calling_fxn)
{
    void* ret = NULL;
    void* err = NULL;
    mfl_node_t* q = p;
    
    if (p->nodet_tip) {
        if (p->nodet_next) {
            dbg_eprintf("terminal node has valid nodet_next pointer\n\n");
            dbg_pcall(calling_fxn);
        }
        else {
            tui_tip_check(p, __FXN_NAME__, verbose);
        }
        
        return 0;
    }
    
    do {
        err = NULL;
        q = q->nodet_next;
        if (q->nodet_edge) {
            err = tui_check_binary_traversal(q->nodet_edge, verbose, __FXN_NAME__);
        }
        else {
            tui_check_node_is_root(q, verbose, __FXN_NAME__);
        }
        if (err) {
            ret = err;
        }
    } while (q != p);
    
    if (!ret) {
        if (mfl_node_is_n_ary(p, 2)) {
            return p;
        }
        else {
            return 0;
        }
    }
    else {
        return ret;
    }
    
}


/**
 # int tui_count_num_taxa(mfl_tree_t *t)
 Counts the number of taxa, but also has a built-in error check.
 It will loop over the taxa array counting the number of taxa that have the
 tip value set. However, it will also attempt to discover whether there are 
 unexpected non-tip nodes inside the node array.
 @Param t pointer to mfl_tree_t
 @Returns int number of taxa
 */
int tui_count_num_taxa(mfl_tree_t *t)
{

    int i = 0;
    int num_taxa = 0;
    bool intipssubarray = true;
    
    do {
        if (t->treet_treenodes[i]->nodet_tip) {
            if (intipssubarray) {
                ++num_taxa;
            }
            else {
                dbg_eprintf("unexpected non-tip node in the terminal sub-array");
                return 0;
            }
        }
        else {
            intipssubarray = false;
        }
    } while (t->treet_treenodes[i]);
    
    return num_taxa;
}


/**
 # tui_check_tree_for_dangling_pointers(mfl_tree_t* t, int num_nodes, int *verbose)
 Scans the tree checking for dangling pointers and other pointer errors. These would
 lead to FATAL ERRORS in processes on the tree.
 @param mfl_tree_t*
 @param int
 @param int*
 @return int error value, 0 for no error
 */
int tui_check_tree_for_dangling_pointers(mfl_tree_t* t, int num_nodes, int *verbose)
{
    int i = 0;
 
    int err = 0;
    
    for (i = 0; i < num_nodes; ++i) {
        
        if (!t->treet_treenodes[i]->nodet_edge) {
            if (!t->treet_treenodes[i]->nodet_isroot) {
                dbg_printf("ERROR in input tree: dangling pointer at: %p\n ", t->treet_treenodes[i]);
                err = 1;
                if (*verbose) {
                    dbg_printf("Error found in: ");
                    tui_print_node_data(t->treet_treenodes[i], __FXN_NAME__);
                }
            }
        }
        else if (i < t->treet_num_taxa) {
            if (t->treet_treenodes[i]->nodet_edge->nodet_tip) {
                dbg_printf("ERROR in %s(): terminals connected without ancestral internal node\n", __FXN_NAME__);
                err = 1;
                
            }
            if (t->treet_treenodes[i]->nodet_next) {
                dbg_printf("ERROR in %s(): terminal node has unexpected non-NULL value at nodet_next\n", __FXN_NAME__);
                err = 1;
            }
            if (t->treet_treenodes[i]->nodet_tip != (i-1)) {
                dbg_printf("ERROR in %s(): unexpected tip number reassignment at: %p\n",__FXN_NAME__, t->treet_treenodes[i]);
                if (*verbose) {
                    tui_print_node_data(t->treet_treenodes[i], __FXN_NAME__);
                }
                err = 1;
            }
        }
    }
    
    return err;
}


/**
 #int tui_check_for_anastomosis(mfl_tree_t* t, int num_nodes, int *verbose)
 Checks for multiple edge pointers pointing to the same node. Such a situation
 could cause a FATAL ERROR in any tree operation routine.
 @param mfl_tree_t* t, the subject tree
 @param int num_nodes, the number of nodes in the whole structure
 @param int* verbose, declaring verbosity
 @returns int the number of errors detected
 */
int tui_check_for_anastomosis(mfl_tree_t* t, int num_nodes, int *verbose)
{
    
    int i = 0;
    int j = 0;
    int count = 0;
    int err = 0;
    
    mfl_nodearray_t testnodes = (mfl_nodearray_t)malloc((num_nodes+1)*sizeof(mfl_node_t*));
    
    memcpy(testnodes, t->treet_treenodes, (num_nodes+1)*sizeof(mfl_nodearray_t*)); // CHECK THAT THIS IS CORRECT SIZING
    
    for (i = 0; i < num_nodes; ++i) {
        count = 0;
        
        for (j = 0; j < num_nodes; ++j) {
            if (testnodes[i] == t->treet_treenodes[j]) {
                ++count;
            }
        }
        if (count > 1) {
            dbg_printf("ERROR detected by %s: mfl_node_t %p interacts with more than one edge\n\n",__FXN_NAME__, testnodes[i]);
            if (*verbose) {
                tui_print_node_data(testnodes[i], __FXN_NAME__);
            }
            ++err;
        }
    }
    
    free(testnodes);
    
    return err;
}

/**
 ## int tui_check_broken_tree(mfl_tree_t *t, int *verbose)
 Attempts to determine if a tree is broken by checking for: dangling pointers 
 intended for nodes; invalid pointers to nodes; cyclicity and anastomosis.
 
 @param mfl_tree_t* tree to be checked.
 @param int* yes/no value for verbose output
 @returns int
 
 */
int tui_check_broken_tree(mfl_tree_t *t, int *verbose)
{
    /* 
     *  What defines a broken tree?
     *  Anastomosis.
     *  Cyclicity.
     *  Pointers that shouldn't be
     *  Dangling pointers that shouldn't be.
     *
     */
    
    int err = 0;
    int num_nodes = 0;
    int num_taxa = 0;
    
    if (!t->treet_num_taxa) {
        dbg_eprintf("tree does not contain treet_num_taxa value!");
        dbg_eprintf(". . . Calculating estimated num_taxa from tips.\n");
        if (!(num_taxa = tui_count_num_taxa(t))) {
            dbg_eprintf("tree has no terminals or has invalid terminal sub-array");
            dbg_printf("Your goddamned tree at %p is broken.\n", t);
            return -1;
        }
    }
    else {
        num_taxa = t->treet_num_taxa;
    }
    
    
    num_nodes = mfl_calculate_number_of_nodes_to_allocate(t->treet_num_taxa);
    
    // Check for dangling pointers and tip node misbehaviour
    err = tui_check_tree_for_dangling_pointers(t, num_nodes, verbose);
    if (err) {
        dbg_pfail("\nYour goddamned tree is broken.");
        dbg_printf("\tThe goddamned broken tree at %p\n", t);
    }
    else {
        dbg_ppass("input tree connections verified OK");
    }
    
    // Checking anastomosis. Each node record should be accessed by no more than one other edge.
    
    err = tui_check_for_anastomosis(t, num_nodes, verbose);
    
    // Checking cyclicity. Each node in a ring should form a closed cycle and only point to nodes
    //      intended to be internal nodes via their nodet_next pointer. Thus, they should have their
    //      nodet_tip value set to 0 and they should not be found in the array between [0 and num_taxa)
    
    // Pointers that should be. Harder to define, but they should point to valid memory or to NULL.
    
    
    
    return err;
}


int tui_check_all_binary(mfl_tree_t *querytree, const int *verbose)
{
    void* ptr = NULL;
    
    mfl_node_t *start_ptr = NULL;
    
    if (ptr = tui_check_binary_traversal(start_ptr, verbose, __FXN_NAME__) ) {
        dbg_eprintf(__FXN_NAME__, "non-binary node detected at: ");
        dbg_printf("%p", ptr);
    }
}


// Printing a newick string from a mfl_tree_t
void tui_print_newick_recursive(mfl_node_t *start)
{
    
    // Set the node to start
    mfl_node_t *node = start;
    
    if (node->nodet_tip != 0) {
        // print a tip
        dbg_printf("%i", node->nodet_tip);
        
        return;
        
    } else {
        // open a clade
        dbg_printf("(");
    }
    
    // Move to the next node
    node = node->nodet_next;
    
    do {
        if (node != start->nodet_next) {
            dbg_printf(",");
        }
        // Go through the next nodes and print the next tip
        tui_print_newick_recursive(node->nodet_edge);
        
        // close a clade if the next edge is not a tip
        if(node->nodet_edge->nodet_tip == 0) {
            dbg_printf(")");
        }

        // Move to the next node
        node = node->nodet_next;
    
    // Stop once reached back the start tree
    } while (node != start);
    
    return;
}

void tui_print_newick(mfl_node_t *start)
{
    // Running the recursive loop
    tui_print_newick_recursive(start);
    // Closing the tree
    dbg_printf(");");
}
