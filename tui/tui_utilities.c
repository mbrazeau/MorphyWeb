//
//  tui_utilities.c
//  Morphy
//
//  Created by mbrazeau on 12/04/2016.
//
//

#include "morphy.h"

void tui_print_node_data(mfl_node_t* p, const char* calling_fxn);
int tui_tip_check(mfl_node_t* n, const char* calling_fxn, const int* verbose);
void* tui_check_binary_traversal(mfl_node_t *p, const int* verbose, const char* calling_fxn);
int tui_check_node_is_root(mfl_node_t* p, const int* verbose, const char* calling_fxn);


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

int tui_check_all_binary(mfl_tree_t *querytree, const int *verbose)
{
    void* ptr = NULL;
    
    mfl_node_t *start_ptr = NULL;
    
    if (ptr = tui_check_binary_traversal(start_ptr, verbose, __FXN_NAME__) ) {
        dbg_eprintf(__FXN_NAME__, "non-binary node detected at: ");
        dbg_printf("%p", ptr);
    }
}