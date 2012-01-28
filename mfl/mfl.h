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

/* 
 * The excludes and includes can be handled by the library NCL... 
 * As such I dont think these functions are necessary in mfl.
 * In other words the data I provide to mfl will already have 
 * all the exclusions excluded... :D
 *
 * MDB: Not necessarily. For instance, I might run an analysis. Find that 
 * there's crappy resolution in the strict consensus tree. But an adams 
 * consensus might show me that that's because one taxon is 'jumping around'. I 
 * might therefore decide to re-run the analysis with a smaller set of the taxa
 * in the datafile. I wouldn't want to have to go back to the datafile and
 * re-write it just to do the experiment.
 *
 * This could actually be under mfl_set_param, actually.
 *
 */
mfl_set_excludes(taxa_to_exclude, chars_to_exclude); // Optionally remove some taxa or characters from the analysis.
mfl_set_includes(taxa_to_include, chars_to_include); // Optionally reinstate removed taxa or characters to the analysis

/*
 * What is the command the user will enter in the Morphy_commands.txt that defines the outgroup taxa?
 *
 * MDB: Haven't added those yet. The list in the txt file wasn't incomplete.
 * Would be something like "outgroup" and a list of taxon names or tip numbers.
 * So it might be "outgroup 3 4 10" or "outgroup Polyodon Acipenser"
 */
mfl_define_outgroup_taxa(list_of_taxon_names OR int tip_numbers); // Partition the taxon set into ingroup and outgroup
/*
 * Again, what command is used for this? 
 * 
 * MDB: This might be something like the above. But the option would be a partial
 * Newick string. "ogtopol (4,((3,10),(Ingroup))". This would create a 
 * constraint tree in memory that would constrain the search to trees that only
 * have the predefined subtree in it. 
 *
 */
mfl_constrain_outgroup_topology(); // user constrains the search to always include a particular subtree topology
/*
 * I propose a function like mfl_set_parameter(mfl_handle*, param_type, param_data*) in leiu of functions like this, 
 * the reasons are: it will make the API a bit smaller, and make upgrade paths easier to handle, in the event that
 * this interface is changed or removed, or new parameters are added..
 *
 * MDB: makes sense. 
 */

mfl_set_parameter(mfl_handle*, param_type, param_data*);

/*
 * This interface is not workable: too many params, a lot of these params
 * should go into a mfl_set_parameter function which modifies values in a data struct, 
 * that datastruct should then be passed into this function. The datastruct should
 * be created with some kind of function like mfl_create_handle() which allocs the memory
 * for the datastruct, and also sets reasonable defaults for all params.
 *
 * Again the datastruct should not be defined in this interface, just a reference to it.
 * End users should not know or care about the guts of the datastruct, we just pass around
 * pointers to it and the mfl uses and abuses it however it sees fit...
 *
 * MDB: I see. I think I just put those in because I wanted to list some minimum criteria
 * that I though each of the operations would need. The programs like Morphy that I've 
 * used in the past usually required a list of similar parameters (minus ntax, numnodes,
 * nchar). For instance, I just kind of dumped out the stuff I would normally
 * type in at the command line in one of those programs. However, if we're
 * using the mfl_set_parameter() function then, as you said, this stuff could
 * stay out of the function prototypes. As long as it's 
 *
 */
mfl_heuristic(int ntax, int numnodes, int nchar, int *list_of_included taxa, int *list_of_included_characters, int *tipdata, addsequence_type, number_of_iterations, branchswapping_type, trees_to_hold_at_each_step, is_ratchet_search);
mfl_exhaustive(int ntax, int numnodes, int nchar, int *list_of_included taxa, int *list_of_included_characters, int *tipdata);
mfl_branchandbound(int ntax, int numnodes, int nchar, int *list_of_included taxa, int *list_of_included_characters, int *tipdata);
mfl_bootstrap(int ntax, int numnodes, int nchar, int *list_of_included taxa, int *list_of_included_characters, int *tipdata);
mfl_decay(int ntax, int numnodes, int nchar, int *list_of_included taxa, int *list_of_included_characters, int *tipdata);
mfl_jackknife(int ntax, int numnodes, int nchar, int *list_of_included taxa, int *list_of_included_characters, int *tipdata);

mfl_safe_taxonomic_reduction(list_of_ordered_characters);
mfl_ratchet(n_chars_to_perturb, reweightorjackknife);   // A type of super-fast heuristic search, but is subordinate to mfl_heuristic
mfl_consensus(tree **savedtrees, bool strict, bool majority_rule, int majority_rule_optional_cutoff, int adams);
mfl_collapse_zerolength(during_search, if_min_zero_length, if_max_zero_length);
/*
 * I am not sure what this function does, why does it take tree ** as a parameter?
 * If it actually does set the treelimit, then again I would propose to make this part
 * of the mfl_set_parameter function
 *
 * MDB: sounds fine.
 *
 */
mfl_set_treelimit(int max_number_of_trees_in_mem, tree **savedtrees);

/*
 * This should be handled by the I/O module (i.e. ncl) not the mfl library,
 * mfl should be able to return a list of trees in a format that can easily be stored
 * to a file, but the actual file I/O does not belong in mfl, same for the show functions...
 * mfl does not interact with the user, it is an API that allows higher level code
 * to interact with the user. If we made a GUI for morphy, then these functions would
 * be useless... What you want is an API that provides the data in such a way that
 * the application can present it to the end user using whatever I/O methods it has
 * available, in a command line program that might be cout/printf, in a GUI it will 
 * be something completely different, and in an opengl 3d visualization application
 * it will be a whole different set of I/O functions again...
 *
 * MDB: Okay. This is all becoming clearer.
 */

/*-Already written-*/
void mfl_resize_treebuffer(tree **treebuffer, int *treelimit, int sizeincrease); // Called if the user gives a command to store more than the current/default number of trees
void mfl_clear_treebuffer(tree **treebuffer, long int *numsavedtrees, int numnodes); // Called when the user gives a command to clear all trees in memory
