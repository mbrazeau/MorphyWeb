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


void mfl_tree_check_traversal(mfl_node_t *node)
{
    /* 
     *** UPDATE FOR MORE EXTENSIVE CHECKING ***
     */
    mfl_node_t *p = NULL;
    int i = 0;
    
    if (node->nodet_tip) {
        dbg_printf("tip: %i\n", node->nodet_tip);
        return;
    }
    
    p = node->nodet_next;
    
    do {
        mfl_tree_check_traversal(p->nodet_edge);
        p = p->nodet_next;
    } while (p != node);
    
    dbg_printf("Node %i: ", node->nodet_index);
    
    p = node->nodet_next;
    
    do {
        ++i;
        dbg_printf("child %i: %i ", i, p->nodet_edge->nodet_index);
        p = p->nodet_next;
    } while (p != node);
    dbg_printf("\n");
    
    return;
}

int tui_print_warning(const char *msg, const char *fxn_name, tui_testrec* testrec)
{
    /*
        Although not yet implemented, this function can be updated to take the warning message and print it to 
        screen. The problem is passing variable numbers of arguments for format conversions
     */
    
    if (testrec->tr_num_warn == ULLONG_MAX) {
        
        if (testrec->tr_verbosity) {
            
            dbg_printf("WARNING: number of warnings has overflowed the buffer at %llu\n", testrec->tr_num_warn);
            
            testrec->tr_warn_overflow = testrec->tr_num_warn;
            testrec->tr_num_warn = 1;
        }
    }
    else {
        ++testrec->tr_num_warn;
    }
    
    return 0;
}

void tui_conclude_test(const char *testname, tui_testrec *testrec)
{
    
    dbg_printf("%s concluded with %llu errors and %llu warnings\n\n", testname, testrec->tr_num_err, testrec->tr_num_warn);
    
    if (testrec->tr_err_overflow) {
        dbg_printf("\t Number of errors overflowed the buffer: %llu additional errors (just to rub it in)\n", testrec->tr_err_overflow);
    }
    if (testrec->tr_warn_overflow) {
        dbg_printf("\t Number of warnings overflowed the buffer: %llu additional warnings\n", testrec->tr_err_overflow);
    }
    
    if (testrec->tr_verbosity == TUI_SILENT) {
        dbg_printf("\tAt least %llu warnings were suppressed. To view these, re-run the test in VERBOSE mode", testrec->tr_num_warn);
    }
    dbg_printf("\n\n");
    
}


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
    
    q = p;
    do {
        err = NULL;
        
        if (q->nodet_edge) {
            err = tui_check_binary_traversal(q->nodet_edge, verbose, __FXN_NAME__);
        }
        if (err) {
            ret = err;
        }
        q = q->nodet_next;
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
 Counts the number of terminals, but also has a built-in error check.
 It will loop over the taxa array counting the number of taxa that have the
 tip value set. However, it will also attempt to discover whether there are 
 unexpected non-tip nodes inside the node array.
 @Param t (mfl_tree_t*) an input tree.
 @Returns Number of taxa (int)
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


bool tui_check_reciprocal_edge(mfl_node_t *n, mfl_tree_t* t, int num_nodes, int *verbose)
{
    int i = 0;
    bool backlink = 0;
    
    for (i = 0; i < num_nodes; ++i) {
        if (t->treet_treenodes[i]->nodet_edge == n) {
            backlink = true;
        }
    }
    
    return backlink;
}

/*!
 Scans the tree checking for dangling pointers and other pointer errors. These would
 lead to FATAL ERRORS in processes on the tree.
 @param t (mfl_tree_t*)
 @param num_nodes (int)
 @param verbose (int*)
 @return int error value, 0 for no error
 */
int tui_check_tree_for_connection_errors(mfl_tree_t* t, int num_nodes, int *verbose)
{
    int i = 0;
 
    int err = 0;
    
    for (i = 0; i < num_nodes; ++i) {
        
        if (!t->treet_treenodes[i]->nodet_edge) {
            if (t->treet_treenodes[i] != t->treet_root) {
                
                if (*verbose) {
                    dbg_printf("WARNING possible dangling edge pointer at: %p\n ", t->treet_treenodes[i]);
                    if (*verbose) {
                        dbg_printf("Error found in: ");
                        tui_print_node_data(t->treet_treenodes[i], __FXN_NAME__);
                    }
                }
                
                if (!tui_check_reciprocal_edge(t->treet_treenodes[i], t, num_nodes, verbose)) {
                    if (*verbose) {
                        dbg_printf("However, no reciprocal backlink found.\n\n");
                    }
                    err = 0;
                }
                else {
                    dbg_printf("ERROR in %s(): node %p points to valid node with no reciprocal link\n", __FXN_NAME__, t->treet_treenodes[i]);
                    if (*verbose) {
                        dbg_printf("Error found in: ");
                        tui_print_node_data(t->treet_treenodes[i], __FXN_NAME__);
                    }
                    err = 1;
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
            if (t->treet_treenodes[i]->nodet_tip != (i+1)) {
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


/*!
 Checks for multiple edge pointers pointing to the same node. Such a situation
 could cause a FATAL ERROR in any tree operation routine.
 @param t (mfl_tree_t*) the subject tree
 @param num_nodes (int) the number of nodes in the whole structure
 @param verbose (int*) declaring verbosity
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
            if (testnodes[i] == t->treet_treenodes[j]->nodet_edge) {
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


/*mfl_tree_t* tui_create_tree_for_testing(int num_taxa, bool israndom, char* input_tree, bool isconstraint)
{
    
    
    
}*/


int tui_check_all_node_ring_circularity(const mfl_tree_t *t, int num_nodes, int *verbose)
{
    int i = 0;
    int err = 0;
    int count = 0;
    mfl_node_t** p = NULL;
    mfl_nodearray_t nds = t->treet_treenodes;
    
    for (i = t->treet_num_taxa; i < num_nodes; ++i) {
        
        // Check for valid next pointer.
        if (!nds[i]->nodet_next) {
            if (nds[i]->nodet_edge) {
                err = 1;
                dbg_eprintf("node with valid nodet_edge has no valid nodet_next pointer");
                dbg_printf("\tNode @: %p", nds[i]);
                if (*verbose) {
                    tui_print_node_data(nds[i], __FXN_NAME__);
                }
            }
        }
        
        // Check for mid-branch internal cycle ("galls")
        if (nds[i]->nodet_next) {
            if (nds[i]->nodet_next->nodet_next) {
                if (nds[i]->nodet_next->nodet_next == nds[i]) {
                    err = 1;
                    dbg_eprintf("Theres a gall in the branch");
                    dbg_printf("\tNode @: %p", nds[i]);
                    if (*verbose) {
                        tui_print_node_data(nds[i], __FXN_NAME__);
                    }
                }
            }
        }

        
        // Check against nodet_next pointer to non-internal node.
        if (nds[i]->nodet_next) {
            if (nds[i]->nodet_next->nodet_tip) {
                err = 1;
                dbg_eprintf("nodet_next pointer in internal node points to terminal node");
                dbg_printf("\tNode @: %p", nds[i]);
                if (*verbose) {
                    tui_print_node_data(nds[i], __FXN_NAME__);
                }
            }
        }
        
        // Check against infinite loops
        // This is accomplished by counting the number of times a node is accessed by another nodet_next pointer.
        
        count = 0;
        
        p = &nds[t->treet_num_taxa];
        
        do {
            if ((**p).nodet_next == nds[i]) {
                ++count;
            }
            ++p;
        } while (*p < nds[num_nodes-1]);
        
        if (count > 1) {
            err = 1;
            dbg_eprintf("node accessed by more than one nodet_next pointer in array. Condition would lead to infinite loop");
            dbg_printf("\tNode @: %p\n", nds[i]);
            if (*verbose) {
                tui_print_node_data(nds[i], __FXN_NAME__);
            }
        }
        
    }
    
    
    return err;
}


/*!
 Attempts to determine if a tree is broken by checking for: dangling pointers
 intended for nodes; invalid pointers to nodes; cyclicity and anastomosis.
 
 @param t (mfl_tree_t*) tree to be checked.
 @param verbose (int*) yes/no value for verbose output
 @returns int
 */
int tui_check_broken_tree(mfl_tree_t *t, int *verbose)
{
    dbg_printf("BEGIN TEST %s()\n", __FXN_NAME__);
    
    int err = 0;
    int ret = 0;
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
    err = tui_check_tree_for_connection_errors(t, num_nodes, verbose);
    
    // Checking anastomosis. Each node record should be accessed by no more than one other edge.
    
    ret = tui_check_for_anastomosis(t, num_nodes, verbose);
    
    if (ret) {
        err = ret;
    }
    ret = 0;
    
    
    // Checking cyclicity. Each node in a ring should form a closed cycle and only point to nodes
    //      intended to be internal nodes via their nodet_next pointer. Thus, they should have their
    //      nodet_tip value set to 0 and they should not be found in the array between [0 and num_taxa)
    
    ret = tui_check_all_node_ring_circularity(t, num_nodes, verbose);
    
    if (ret) {
        err = ret;
    }
    ret = 0;
    
    
    if (err) {
        dbg_pfail("\nYour goddamned tree is broken.");
        dbg_printf("\tThe goddamned tree at %p. That's the one that's broken.\n", t);
    }
    else {
        dbg_ppass("input tree connections verified OK");
    }
    
    dbg_printf("\nEND broken tree test\n\n");
    return err;
}


int tui_check_all_binary(mfl_tree_t *querytree, const int *verbose)
{
    void* ptr = NULL;
    
    mfl_node_t *start_ptr = NULL;
    
    if ((ptr = tui_check_binary_traversal(start_ptr, verbose, __FXN_NAME__)) ) {
        dbg_eprintf("non-binary node detected at: ");
        dbg_printf("%p", ptr);
        
        return 1;
    }
    else {
        return 0;
    }
}


// Printing a newick string from a mfl_tree_t
void tui_print_newick_recursive(mfl_node_t *n)
{
    
    // Set the node to start
    mfl_node_t *node = n;
    
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
        if (node != n->nodet_next) {
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
    } while (node != n);
    
    return;
}


void tui_print_newick(mfl_node_t *start)
{
    // Running the recursive loop
    tui_print_newick_recursive(start);
    // Closing the tree
    dbg_printf(");");
}


void tui_test_checktree_(void)
{
    
    char temp_example1[] = "temp_examp2=[&R] (2,(6,((3,4),(5,1))));";
    char temp_example2[] = "temp_examp4=[&U] (2,((3,4),(5,20),1));";
    
    dbg_printf("\n\nTesting an assembled tree for breaks\n\n");
    mfl_tree_t* testtr = mfl_convert_newick_to_mfl_tree_t(temp_example1, 0);
    
    
    /* Cycle a ring */
    //testtr->treet_treenodes[6]->nodet_next->nodet_next = testtr->treet_treenodes[6]->nodet_next;
    
    int verbose = TUI_SILENT;
    tui_check_broken_tree(testtr, &verbose);
    
    mfl_free_tree(testtr);
    
    /* Do it again with non-consecutive tip numbers. */
    testtr = mfl_convert_newick_to_mfl_tree_t(temp_example2, 0);

    tui_check_broken_tree(testtr, &verbose);
    
    mfl_free_tree(testtr);

}

void tui_print_edgetable(mfl_edgetable_t* edgetable)
{
    int i = 0;
    dbg_printf("Tip/node connects to tip/node\n");
    
    for(i = 0; i < edgetable->numentries; ++i) {
        dbg_printf("%i connects to %i\n", i, edgetable->edgetable[i]);
    }
}

//Converts all bitfield into a series of integers
void tui_print_bitset(mfl_bitset_t* bitset)
{
    int i = 0;
    dbg_printf("%llu", bitset->bts_bitfields[i]);
    ++i;
    for(i = 1; i < bitset->bts_nfields; ++i){
        dbg_printf(", ");
        dbg_printf("%llu", bitset->bts_bitfields[i]);
    }
}
//Prints the biparition table
void tui_print_bipartition_tables(mfl_bipartition_table* bipartition_table)
{
    int i = 0;
    dbg_printf("Total number of biparitions = %i\n", bipartition_table->bipartitions_list->bsl_num_sets);
    for(i = 0; i < bipartition_table->bipartitions_list->bsl_num_sets; ++i){
        dbg_printf("Biparition ");
        tui_print_bitset(bipartition_table->bipartitions_list->bsl_bitsets[i]);
        dbg_printf(" - counted %i times\n", bipartition_table->bipartition_occurence_counter[i]);
    }
}

// Run a matrix
void tui_runmatrix(char matrix[], int expect)
{
    int num_taxa = 12;//78
    int num_chars = 1;//236;
    int num_og_tax = 0;
    char* testnewick;
    
    testnewick = (char*)"UNTITLED = [&R] ((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));";
    mfl_tree_t* testtree = mfl_convert_newick_to_mfl_tree_t(testnewick, num_taxa);
    
    mfl_handle_s* handle = mfl_t2s(mfl_create_handle());
    
    // Setup the handle:
    handle->n_taxa = num_taxa;
    handle->n_chars = num_chars;
    handle->n_outgroup_taxa = num_og_tax;
    handle->n_to_hold = 3;
    handle->addseq_type = MFL_AST_ASIS;
    
    handle->input_data = matrix;
    mfl_partition_set_t* dataparts = mfl_generate_search_data(handle);
    
    mfl_setup_input_tree_with_node_data(testtree, dataparts);
    tui_check_broken_tree(testtree, false);
    mfl_fullpass_tree_optimisation(testtree, dataparts);
    
    dbg_printf("\n===================\n");
    dbg_printf("Martix: %s\n", matrix);
    dbg_printf("Is of length: %i\n", testtree->treet_parsimonylength);
    dbg_printf("Should be %i.", expect);
    dbg_printf("\n===================\n");
    //
    mfl_free_tree(testtree);
}
