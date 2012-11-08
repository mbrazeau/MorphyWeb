/*
 *  mfl.h
 *  Morphy
 *
 *  Created by Martin Brazeau and Chris Desjardins on 1/26/12.
 *  Copyright 2012. All rights reserved.
 *
 */

#pragma once

#include <string>
#include <vector>
#include <stdexcept>

using namespace std;

typedef enum
{
    /* MFL_PT_NUM_TAX (int) Number of taxa in the dataset, and the maximum number of
     * leaves in any tree used by Morphy.*/
    MFL_PT_NUM_TAX,
    
    /* Deleting MFL_PT_NUM_NODES. The number of nodes can be calculated by
     * a library routine */
    //MFL_PT_NUM_NODES,
    
    /* MFL_PT_NUM_CHAR (int) the number of columns ('transformation series') in 
     * the character matrix. Note that polymorphisms (i.e. {01}) only count as a 
     * single cell. Thus it's important to take this value from NCHAR in the 
     * Nexus file and then check for errors */
    MFL_PT_NUM_CHAR,
    
    /* MFL_PT_SEARCH_TYPE (see mfl_search_t enum for parameter) 
     * the type of search procedure used, whether it is exhaustive, 
     * branch-and-bound, or heuristic. The type of search would be a command by
     * itself OR would be a command option when either a bootstrap (BTS) or 
     * jackknife (JK) search is conducted (since these analyses are just 
     * iterative searches, and thus require the user to set a type of search
     * algorithm. */
    MFL_PT_SEARCH_TYPE, 
    
    /* MFL_PT_NUM_ITERATIONS (int) the number of iterations in a in a heuristic 
     * search that uses a random addition sequence. */
    MFL_PT_NUM_ITERATIONS,
    
    /* MFL_PT_TREELIMIT (int) maximum number of trees allowed to be stored in
     * memory. This should start as some default value (100 to 500 is certainly
     * safe). The user may choose a different value prior to the search and 
     * and also choose to allow the program to automatically increase the limit
     * (or simply have no limit). It should also be possible to change this
     * value in the middle of a search. For instance, if the user has not 
     * instructed the program to use automatically resize the tree list, then
     * if the analysis runs out of space (i.e. needs a bigger array to save more
     * trees, then the user would be prompted to increase the tree limit, and 
     * to choose whether they would like that limit to automatically increase */
    MFL_PT_TREELIMIT,
    
    /* MFL_PT_BRANCH_SWAP_TYPE (see mfl_branch_swap_t enum for parameter) There are three main 
     * types of branch swapping algorithm use ONLY in a heuristic search: 
     * nearest-neighbor interchange (NNI), subtree pruning and regrafting (SPR), 
     * and tree bisection and reconnection (TBR). One of these options would be 
     * chosen when the user chooses to initiate a heuristic search. The default 
     * value should be TBR, but that is complicated and so I've not yet written 
     * it. */
    MFL_PT_BRANCH_SWAP_TYPE,
    
    /* MFL_PT_RATCHET_SEARCH (bool) this is an option for speeding up the 
     * heuristic search. It involves initial modifications to the data matrix
     * before searching, followed by a search on the 'normal', unperturbed data.
     * It would be selected when the user enters commands for heuristic search */
    MFL_PT_RATCHET_SEARCH,
    
    /* MFL_PT_TIP_DATA (char*). The sequence of characters in the data matrix. 
     * These get converted to the charstate int32_t by functions in the Morphy
     * library immediately before a search is initiated. The reason why we do it
     * this way is because the way the data gets converted will depend on some 
     * other options set by the user. */
    MFL_PT_INPUT_DATA,
    
    /* MFL_PT_ADD_SEQUENCE_TYPE (see mfl_add_sequence_t enum for parameter). 
     * Selects the manner in which branches are added during the generation of
     * starting trees. */
    MFL_PT_ADD_SEQUENCE_TYPE,
    
    /* MFL_PT_COLLAPSE (bool)  
     * collapse zero-length branches option is in effect during search. */
    MFL_PT_COLLAPSE,
    
    /* MFL_PT_COLLAP_AT (see mfl_set_collapse_at_t enum for parameter)
     * sets whether to collapse branches with minimum 
     * length 0, with maximum length 0, or have two incident nodes with equal 
     * apomorphies reconstruction sets. Default is to collapse nodes if their 
     * MAXIMUM length is 0 */
    MFL_PT_COLLAP_AT,
    
    /* MFL_PT_GAP (mfl_gap_t - missing data is default) sets whether gap symbol ('-') 
     * will be treated as inapplicability or as missing data (i.e. equivalent to '?' score). */
    MFL_PT_GAP,
    
} mfl_param_t;

typedef enum
{
    MFL_RT_NUM_REARRANGMENTS, // (long, or long long) number of swaps made in heuristic search
    MFL_RT_NUM_SAVED_TREES,   // (long int) number of trees saved during search
    MFL_RT_SHORTEST_TREE_LEN, // (int) number of steps of shortest tree found in search
    MFL_RT_SEARCH_TIME,       // (time_t) amount of time taken for the search
    MFL_RT_ISLAND_COUNT,      // (int) number of islands
    MFL_RT_ISLAND_SIZE,       // (int) island size
    MFL_RT_ISLAND_LENGTH,     // (int) island length
} mfl_resultant_data_t;

/* 
 * mfl_search_t - used as a parameter when the MFL_PT_SEARCH_TYPE is being set
 */
typedef enum
{
    MFL_ST_EXHAUSTIVE,
    MFL_ST_BRANCH_BOUND,
    MFL_ST_HEURISTIC,

    MFL_ST_MAX
} mfl_search_t;

/*
 * mfl_branch_swap_t - used as a parameter when the MFL_PT_BRANCH_SWAP_TYPE is being set
 * CJD: you might want to make these more english like, or at least add comments to what each
 * one is...
 */
typedef enum
{
    MFL_BST_TBR, // Tree bisection and reconnection: prunes a subtree, and performs all re-rooting + re-insertion operations
    MFL_BST_SPR, // Subtree pruning and regrafting: prunes a subtree and reinserts it at all locations in the target tree
    MFL_BST_NNI, // Nearest-neighbor interchanges: swaps two sub-trees separated by an internal node

    MFL_BST_MAX
} mfl_branch_swap_t;

/*
 * mfl_add_sequence_t - used as a parameter when the MFL_PT_ADD_SEQUENCE_TYPE is being set
 */
typedef enum
{
    MFL_AST_SIMPLE,
    MFL_AST_RANDOM,
    MFL_AST_ASIS,
    MFL_AST_CLOSEST,

    MFL_AST_MAX
} mfl_add_sequence_t;

/*
 * mfl_set_collapse_at_t - use as a parameter when the MFL_PT_COLLAP_AT is being set
 */
typedef enum
{
    MFL_SC_MAX_LEN,
    MFL_SC_MIN_LEN,
    MFL_SC_EQUAL_RECONSTRUCTION_SETS,

    MFL_SC_MAX
} mfl_set_collapse_at_t;

/*
 * mfl_gap_t - use as a parameter when the  is being set MFL_PT_GAP
 */
typedef enum
{
    MFL_GAP_INAPPLICABLE,
    MFL_GAP_MISSING_DATA,

    MFL_GAP_MAX
} mfl_gap_t;

typedef void* mfl_handle_t;

mfl_handle_t mfl_create_handle();
void mfl_destroy_handle(mfl_handle_t mfl_handle);
void mfl_erase_trees_newick(mfl_handle_t mfl_handle);
void mfl_destroy_resultant_data(mfl_handle_t mfl_handle);

/* MDB: Now that I see some of this new code, I think I now get what the 
 * problem with these functions is. Indeed, the lists of included and excluded 
 * characters and taxa are effectively parameters in any of the analyses. If I
 * understand this correctly, a command to change the taxon/character inclusion
 * would just result in a new mfl_handle with all the appropriate set parameters
 */

#ifdef MFL_EXCLUDES
mfl_set_excludes(mfl_handle_t, taxa_to_exclude, chars_to_exclude); // Optionally remove some taxa or characters from the analysis.
mfl_set_includes(mfl_handle_t, taxa_to_include, chars_to_include); // Optionally reinstate removed taxa or characters to the analysis
#endif


#ifdef VERSION_2_0
mfl_define_outgroup_taxa(list_of_taxon_names OR int tip_numbers); // Partition the taxon set into ingroup and outgroup
mfl_constrain_outgroup_topology(); // user constrains the search to always include a particular subtree topology
#endif

/*
** CJD: Here is an example call:
** bool ret = mfl_set_parameter(mfl_handle, MFL_PT_BRANCH_SWAP_TYPE, (void*)MFL_BST_TBR);
*/
bool mfl_set_parameter(mfl_handle_t mfl_handle, mfl_param_t param_type, void *param_data);
void* mfl_get_parameter(mfl_handle_t mfl_handle, mfl_param_t param_type);

/* 
 * CJD: Need to come up with return codes for these, I assume it is a list of trees... 
 * but is that really what we want? The user is going to do a heuristic search (for example)
 * and then they will print the results to the screen, or to a file... so maybe these functions
 * just return a string with the data formatted ready for output? OR maybe they just return a count
 * of the number of trees found? and we then have another function like:
 *
string mfl_get_trees(mfl_handle_t mfl_handle);
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
 * We can save all of the necessary data in the mfl handle, that is not an issue
 * The issue is: what does the end user actually need to see with his two eyes?
 * What needs to be saved to a file? What format? All of the details of datastructures
 * should be hidden from the end user, if the user wants to compute consensus trees
 * then they will be able to do a mfl_compute_consensus(mfl_handle_t) and that should
 * return some data that is ready for direct I/O. OR there should be extra functions
 * that pulls the results out of the handle and formats them in such a way that it can
 * be cleanly presented to the end user, or stored to disk... Maybe there are different
 * functions to present different formats... or maybe there are parameters for different
 * formats... but just returning some complex data struct and expecting random programmers
 * not to screw it up is asking for trouble.
 *
 * MDB: Ignoring the struct I suggested, I think all of the data I listed above 
 * needs to be saved to file when any resultant trees are saved to file. 
 * The user may choose to save all of their results trees, or just the
 * consensus tree(s). Either way, they would be saved as Newick strings
 * in a Nexus-formatted file that just has a TREES block and (optionally) a TAXA 
 * block if the taxon names are not inlined in the Newick trees (i.e. if only 
 * leaf numbers are used in the Newick strings). As a default for tree searches, 
 * the output to screen should be fairly limited. It should report the length of 
 * the shortest tree found, the number of trees saved and the number of 
 * rearrangements tried (also, number of replicates if random addition sequence 
 * was used). So, the search ends, and the user might get output like this:
 *
 *      Heuristic search completed:
 *          rearrangments tried: 2998374571630 
 *          number of trees saved: 299
 *          shortest tree length: 180
 *
 * The main thing the user is going to want to save are the trees. The trees are 
 * saved as Newick-formatted strings, as described above. The file also has to 
 * contain all of the above information, plus the search parameters.
 *
 * Trees should be printed to screen only when the user asks or when a consensus
 * tree is computed. At present, the only function we have for printing to 
 * screen is printNewick (in main.c, but that'll get moved obviously). Nobody
 * wants to read those. Ideally, we'd print them out as text trees like:
 *
 * taxon_1_name --\
 *                +--\
 * taxon_3_name --/  |
 *                   +-  
 * taxon_2_name -----/
 *
 * Printing either Newick trees or a tree like this requires (for each saved 
 * tree) a traversal on the original tree struct. So, whatever does the printing
 * to screen needs that data OR receives the char string(s) ready for printing. 
 *
 * So, this looks like multiple function calls to the handle.
 *
 * CJD: Yes, this all sounds good, the basic idea is that the mfl API does the 
 * heavy lifting and provides APIs that can be easily used by higher level code,
 * and the APIs just return ready to use information. So I see things like the 
 * following, assuming all of data is of type int, which might not be the case:
 */
int mfl_get_resultant_data(mfl_handle_t mfl_handle, mfl_resultant_data_t resultant_data, int param);
/* 
 * note: the string datatype is a real C++ datatype, like char* 
 * but much easier to use, no malloc/free necessary, and it is 
 * always big enough to store whatever string you put in it... like magic...
 * it is part of the std c++ library, which you can now use in mfl 
 * because it is compiled with the c++ compiler :D tho it requires the line:
 using namespace std;
 * before the string datatype can be used in the code, otherwise you must
 * refer to it as std::string
 *
 * MDB: I don't know if it matters at all at this stage, but for an eventual
 * Mac GUI, these would probably have to have an NSstring return value. In fact,
 * for any C++ data types we return, we should probably prepare to have an
 * analogous function that returns the analogous Objective-C type. Then probably
 * have some conditional compilation setup. It'll save us a few headaches
 * down the road, I think.
 *
 * CJD: Actually you can use the STL strings in Objective-C and then "it should just work"
 * on all platforms. p.s. there is a long history with the phrase "it should just work"
 * in the world of computer science: CS nerds love to say it, and managers hate to hear
 * it... because when does it ever just work?
 *
 * Also making the search functions return bool, tho that might not be what you want?
 * I am not sure, but something needs to go there to get this file to compile.
 */
char** mfl_get_saved_trees_newick(mfl_handle_t mfl_handle);
vector<string> mfl_get_saved_trees       (mfl_handle_t mfl_handle);


/* The functions below return the data on the number (count), size, and length of 
 * each of the islands of optimality found by the search. They can be used to 
 * print the data about the optimality islands to screen. For instance, after 
 * completion of a heuristic search with random addition sequences, an interface 
 * display the following:
 *
 *  Island     Size     Length
 *  ---------------------------
 *       1     3348        318
 *       2     1842        318
 *       3      542        319
 *       4       12        323
 * 
 * To do this, call loop for all i < the value returned by mfl_get_island_count(),
 * where island_number = i.
 */
/* use mfl_get_resultant_data instead
int mfl_get_island_count(mfl_handle_t mfl_handle);
int mfl_get_island_size(mfl_handle_t mfl_handle, int island_number);
int mfl_get_island_length(mfl_handle_t mfl_handle, int island_number);
*/

bool mfl_heuristic           (mfl_handle_t mfl_handle);
bool mfl_exhaustive          (mfl_handle_t mfl_handle);
bool mfl_branchandbound      (mfl_handle_t mfl_handle);
bool mfl_consensus           (mfl_handle_t mfl_handle);
bool mfl_collapse_zerolength (mfl_handle_t mfl_handle);

#ifdef VERSION_1_5
mfl_bootstrap       (mfl_handle_t mfl_handle);
mfl_decay           (mfl_handle_t mfl_handle);
mfl_jackknife       (mfl_handle_t mfl_handle);
mfl_safe_taxonomic_reduction(list_of_ordered_characters);
mfl_ratchet(n_chars_to_perturb, reweightorjackknife);   // A type of super-fast heuristic search, but is subordinate to mfl_heuristic
#endif

class mfl_exception : public runtime_error
{
    public:
        mfl_exception(mfl_handle_t mfl_handle, string const& msg);
        virtual const char * what() const throw();
    protected:
        mfl_handle_t m_mfl_handle;
};
