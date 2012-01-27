/*
 *  mfl.h
 *  Morphy
 *
 *  Created by Martin Brazeau and Chris Desjardins on 1/26/12.
 *  Copyright 2012. All rights reserved.
 *
 */

/* 
 * The underlying model here is that a datafile has been opened, parsed, and 
 * read. mfl_analysis_environment has a series of variables within it set 
 * according to the file data or some values calculated based on that. From 
 * here, it's possible to call any number of the subsequent functions relating 
 * to data analysis and tree searching. Some of the parameters would come from
 * command-line input. Others would be values available within the scope of 
 * mfl_analysis_environment, set either by the input file data or by functions
 *
 */

mfl_analysis_environment(data_from_input_file); // When morphy has data from the reader, this starts up otherwise, the user gets an error msg that there's no active data file
mfl_set_excludes(taxa_to_exclude, chars_to_exclude); // Optionally remove some taxa or characters from the analysis.
mfl_set_includes(taxa_to_include, chars_to_include); // Optionally reinstate removed taxa or characters to the analysis
mfl_define_outgroup_taxa(list_of_taxon_names OR int tip_numbers); // Partition the taxon set into ingroup and outgroup
mfl_constrain_outgroup_topology(); // user constrains the search to always include a particular subtree topology
mfl_set_optimization_type(int *list_of_characters, char *optimization_type);
mfl_heuristic(int ntax, int numnodes, int nchar, int *list_of_included taxa, int *list_of_included_characters, int *tipdata, addsequence_type, number_of_iterations, branchswapping_type, trees_to_hold_at_each_step, is_ratchet_search);
mfl_exhaustive(int ntax, int numnodes, int nchar, int *list_of_included taxa, int *list_of_included_characters, int *tipdata);
mfl_branchandbound(int ntax, int numnodes, int nchar, int *list_of_included taxa, int *list_of_included_characters, int *tipdata);
mfl_bootstrap(int ntax, int numnodes, int nchar, int *list_of_included taxa, int *list_of_included_characters, int *tipdata);
mfl_jackknife(int ntax, int numnodes, int nchar, int *list_of_included taxa, int *list_of_included_characters, int *tipdata);
mfl_safe_taxonomic_reduction(list_of_ordered_characters);
mfl_ratchet(n_chars_to_perturb, reweightorjackknife);   // A type of super-fast heuristic search, but is subordinate to mfl_heuristic
mfl_consensus(tree **savedtrees, bool strict, bool majority_rule, int majority_rule_optional_cutoff, int adams);
mfl_collapse(during_search, if_min_zero_length, if_max_zero_length);
mfl_set_treelimit(int max_number_of_trees_in_mem, tree **savedtrees);
mfl_save_trees_to_file(tree **savedtrees); // Depends on there being trees in memory
mfl_show_all_trees(tree **savedtrees); // Prints out all the trees in memory, if any
mfl_show_data(); // Prints the data matrix on screen

/*-Already written-*/
void mfl_resize_treebuffer(tree **treebuffer, int *treelimit, int sizeincrease); // Called if the user gives a command to store more than the current/default number of trees
void mfl_clear_treebuffer(tree **treebuffer, long int *numsavedtrees, int numnodes); // Called when the user gives a command to clear all trees in memory