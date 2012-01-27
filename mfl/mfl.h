/*
 *  mfl.h
 *  Morphy
 *
 *  Created by Martin Brazeau and Chris Desjardins on 1/26/12.
 *  Copyright 2012. All rights reserved.
 *
 */

mfl_analysis_environment();
mfl_set_excludes(taxa_to_exclude, chars_to_exclude); // Optionally remove some taxa or characters from the analysis.
mfl_set_includes(taxa_to_include, chars_to_include); // Optionally reinstate removed taxa or characters to the analysis
mfl_define_outgroup_taxa(list_of_taxon_names OR int tip_numbers);
mfl_constrain_outgroup_topology();
mfl_set_optimization_type(int *list_of_characters, char *optimization_type);
mfl_heuristic(ntax, numnodes, nchar, addsequence_type, number_of_iterations, branchswapping_type, trees_to_hold_at_each_step, is_ratchet_search);
mfl_exhaustive();
mfl_branchandbound();
mfl_bootstrap();
mfl_jackknife();
mfl_safe_taxonomic_reduction();
mfl_ratchet(n_chars_to_perturb, reweightorjackknife);   // A type of super-fast heuristic search
mfl_consensus();
