/*
 *  mfl.h
 *  Morphy
 *
 *  Created by Martin Brazeau and Chris Desjardins on 1/26/12.
 *  Copyright 2012. All rights reserved.
 *
 */

typedef enum
{
    MFL_PT_NUM_TAX,
    MFL_PT_NUM_NODES,
    MFL_PT_NUM_CHAR,
    MFL_PT_NUM_ITERATIONS,
    MFL_PT_NUM_TREES,
    MFL_PT_BRANCH_SWAP_TYPE,
    MFL_PT_RATCHET_SEARCH,
    MFL_PT_TIP_DATA,
    MFL_PT_ADD_SEQUENCE_TYPE,
} mfl_param_t;

mfl_handle_t* mfl_create_handle();
void mfl_destroy_handle(mfl_handle_t* mfl_handle);

/* MDB: Now that I see some of this new code, I think I now get what the 
 * problem with these functions is. Indeed, the lists of included and excluded 
 * characters and taxa are effectively parameters in any of the analyses. If I
 * understand this correctly, a command to change the taxon/character inclusion
 * would just result in a new mfl_handle with all the appropriate set parameters
 */

#ifdef MFL_EXCLUDES
mfl_set_excludes(mfl_handle_t*, taxa_to_exclude, chars_to_exclude); // Optionally remove some taxa or characters from the analysis.
mfl_set_includes(mfl_handle_t*, taxa_to_include, chars_to_include); // Optionally reinstate removed taxa or characters to the analysis
#endif


#ifdef VERSION_2_0
mfl_define_outgroup_taxa(list_of_taxon_names OR int tip_numbers); // Partition the taxon set into ingroup and outgroup
mfl_constrain_outgroup_topology(); // user constrains the search to always include a particular subtree topology
#endif


bool mfl_set_parameter(mfl_handle_t* mfl_handle, mfl_param_t param_type, void *param_data);

/* 
 * CJD: Need to come up with return codes for these, I assume it is a list of trees... 
 * but is that really what we want? The user is going to do a heuristic search (for example)
 * and then they will print the results to the screen, or to a file... so maybe these functions
 * just return a string with the data formatted ready for output? OR maybe they just return a count
 * of the number of trees found? and we then have another function like:
 *
string mfl_get_trees(mfl_handle_t* mfl_handle);
 * 
 * Which pulls all of the trees out of the mfl_handle and formats them into a string ready for output
 * to the user, then the interface can just directly write the string to the screen, or to a file, etc...
 * What other use cases are there?
 *
 * MDB: Perhaps they return a struct mfl_results*. Contains the relevant data
 * such as the length of the shortest tree(s), number of trees save, number of 
 * rearrangments tried (where relevant), a pointer to the list of saved
 * trees, etc. We need to save the tree structs because the user may choose to 
 * compute consensus trees and only save those (especially if > 15 trees were 
 * saved) or they'll want compute full ancestral states reconstructions. These will
 * require the original tree structs in memory.
 *
 */
mfl_heuristic       (mfl_handle_t* mfl_handle);
mfl_exhaustive      (mfl_handle_t* mfl_handle);
mfl_branchandbound  (mfl_handle_t* mfl_handle);
mfl_bootstrap       (mfl_handle_t* mfl_handle);
mfl_decay           (mfl_handle_t* mfl_handle);
mfl_jackknife       (mfl_handle_t* mfl_handle);


#ifdef VERSION_1_5
mfl_safe_taxonomic_reduction(list_of_ordered_characters);
mfl_ratchet(n_chars_to_perturb, reweightorjackknife);   // A type of super-fast heuristic search, but is subordinate to mfl_heuristic
mfl_consensus(tree **savedtrees, bool strict, bool majority_rule, int majority_rule_optional_cutoff, int adams);
mfl_collapse_zerolength(during_search, if_min_zero_length, if_max_zero_length);
#endif

/* 
 * CJD: I think the resize treebuffer is not an external function, if the user changes the number of trees to store, 
 * then the mfl_set_param will be called with MFL_PT_NUM_TREES, and that will invoke the resize_treebuffer automagically.
 * 
 * The clear_treebuffer, I can see as being an external function, however I am not sure if that is the interface to it.
 * I am not sure why the params numsavedtrees, and numnodes are needed, could that data be stored in the tree datastruct?
 * Or in mfl_handle? in which case the mfl_clear_treebuffer function just takes mfl_handle_t* mfl_handle as a parameter.
 * Additionally if we are going to make tree a part of the API, then it should really be called mfl_tree... but I am not 100%
 * sure that we need to export the tree datastruct...
 *
 * MDB: Right. These are probably both 'behind the scenes'. Not sure how
 * numsavedtrees could be stored usefully in the tree struct. More likely in 
 * mfl_handle (or even in the results struct I suggested).
 */
/*-Already written-*/
void mfl_resize_treebuffer(tree **treebuffer, int *treelimit, int sizeincrease); // Called if the user gives a command to store more than the current/default number of trees
void mfl_clear_treebuffer(tree **treebuffer, long int *numsavedtrees, int numnodes); // Called when the user gives a command to clear all trees in memory
