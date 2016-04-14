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


int tui_check_broken_tree(mfl_tree_t *t)
{
    /* 
     *  What defines a broken tree?
     *  Anastomosis.
     *  Cyclicity.
     *  Pointers that shouldn't be
     *  Dangling pointers that shouldn't be.
     *
     */
    
    // Checking anastomosis. Each node record should be accessed by no more than one other edge.
    
    // Checking cyclicity. Each node in a ring should form a closed cycle and only point to nodes
    //      intended to be internal nodes via their nodet_next pointer. Thus, they should have their
    //      nodet_tip value set to 0 and they should not be found in the array between [0 and num_taxa)
    
    // Pointers that should be. Harder to define, but they should point to valid memory or to NULL.
    
    // Dangling pointers. Related to the above.
    
    return 0;
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
    tui_print_newick_recursive(mfl_node_t *start);
    // Closing the tree
    dbg_printf(");");
}
