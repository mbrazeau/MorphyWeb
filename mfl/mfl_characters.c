/*
 *
 * mfl_characters.c
 *
 * Functions for processing the character input. These functions assign
 * the character types to their relevant partitions as mfl_charstate type
 * vectors.
 *
 */

#include "morphy.h"

void mfl_move_in_nexus_multistate(char **col);
mfl_charstate* mfl_allocate_nodal_character_set(int num_characters);
void mfl_free_nodedata(mfl_nodedata_t *olddata);
mfl_nodedata_t* mfl_alloc_datapart(void);

mfl_nodedata_t* mfl_alloc_datapart(void)
{
    mfl_nodedata_t* newdatapart = NULL;
    
    newdatapart = (mfl_nodedata_t*)malloc(sizeof(mfl_nodedata_t));
    
    if (!newdatapart) {
        dbg_printf("ERROR in mfl_alloc_datapart(): unable to allocate memory for new data partition.\n");
        return NULL;
    }
    else {
        memset(newdatapart, 0, sizeof(mfl_nodedata_t));
    }
    
    return newdatapart;
}


void mfl_free_nodedata(mfl_nodedata_t *olddata)
{
    if (olddata->nd_prelim_set) {
        free(olddata->nd_prelim_set);
        olddata->nd_prelim_set = NULL;
    }
    if (olddata->nd_final_set) {
        free(olddata->nd_final_set);
        olddata->nd_final_set = NULL;
    }
    if (olddata->nd_subtree_prelim_set) {
        free(olddata->nd_subtree_prelim_set);
        olddata->nd_subtree_prelim_set = NULL;
    }
    if (olddata->nd_subtree_final_set) {
        free(olddata->nd_subtree_final_set);
        olddata->nd_subtree_final_set = NULL;
    }
    
    // Depending on how cost matrices are handled, they might be freed here, too.
    
    free(olddata);
}


mfl_charstate* mfl_allocate_nodal_character_set(int num_characters)
{
    mfl_charstate *newcharset = NULL;
    
    newcharset = (mfl_charstate*)malloc(num_characters * sizeof(mfl_charstate));
    if (!newcharset) {
        dbg_printf("ERROR in mfl_allocate_nodal_character_set(): unable to allocate memory for character data.\n");
        return NULL;
    }
    else {
        memset(newcharset, 0, num_characters * sizeof(mfl_charstate));
    }
    
    return newcharset;
}


// Processing the input datamatrix which is received as a char*
    // Need NTAX and NCHAR from the associated Nexus file. These should live in the handle.
    // Need to determine the number of partitions to create
        // This depends on:
            // How gap characters are treated
            // How many optimisation method types are used

/* Processing input */

// Initialisation function

mfl_charstate mfl_convert_gap_character(mfl_optimisation_t opt_method)
{
    if (opt_method == MFL_GAP_INAPPLICABLE || opt_method == MFL_GAP_NEWSTATE) {
        //
        return 1;
    }
    else {
        return ~0;
    }
}

char *mfl_move_past_eq_sign(char *input)
{
    /* Move through a line and return the position of the first character after
     * an '=' sign; return original input position if no such symbol found. */
    char *eq_ptr = NULL;
    eq_ptr = input;
    do {
        if (*eq_ptr == '=') {
            ++eq_ptr;
            return eq_ptr;
        }
        ++eq_ptr;
    } while (*eq_ptr);
    
    return input;
}

int mfl_check_state_number_support(char *datatype_list)
{
    
    /* Checks if the number of input states is supported by the default maximum
     * number of states. Returns the */
    
    int num_symbols = 0;
    char *dtype_ptr = NULL;
    
    // In case the input contains the tokens for Format and DataType as used in
    // Nexus files, this will move past the equals sign that begins the
    // assignment.
    dtype_ptr = mfl_move_past_eq_sign(datatype_list);
    
    while (dtype_ptr) {
        if (!isspace(*dtype_ptr)) {
            ++num_symbols;
        }
        ++dtype_ptr;
    }
    
    assert(num_symbols != 0);
    
    if (num_symbols > MORPHY_MAX_STATE_NUMBER) {
        dbg_printf("ERROR in mfl_check_number_support(): number of input states exceeds maximum allowed (%i).\n", MORPHY_MAX_STATE_NUMBER);
        dbg_printf("State list will be truncated\n");
        return MORPHY_MAX_STATE_NUMBER;
    }
    else {
        return num_symbols;
    }
    
}

void mfl_set_datatype_converter_from_nexus(char* datype_converter, char* datatype_list, int num_states)
{
    /* Sets an new char array to be a space-less vector of symbols use as 
     * character states in a Nexus-derived "Format Symbols" command. A lookup
     * into this will return the required shift value. */
    int i = 0;
    char *dtype_ptr = NULL;
    
    dtype_ptr = datatype_list;
    
    while (dtype_ptr) {
        if (!isspace(*dtype_ptr)) {
            datype_converter[i] = *dtype_ptr;
            ++i;
        }
    }
}

int mfl_convert_nexus_symbol_to_shift_value(char in, char *datype_converter)
{
    int i = 0;
    
    do {
        if (in == datype_converter[i]) {
            break;
        }
        ++i;
    }
    while (datype_converter[i]);
    
    return i + MORPHY_SPECIAL_STATE_PAD;
}

int mfl_convert_alphanum_to_shift_value(char in)
{
    /* Takes an alphanumeric value on the scale 0-9ABC...xyz and returns an int
     * representing the number of shifts needed for setting its corresponding
     * bit in an mfl_charstate. This function is used if there is no pre-
     * specified list of state symbols. In effect, Morphy will convert them to a
     * character state corresponding to their order as a default. Will not work
     * well with more than 62 states since it only uses letters and numbers. */
    
    if (isdigit(in)) {
        return in - '0' + MORPHY_SPECIAL_STATE_PAD;
    }
    else if (isalpha(in))
    {
        if (isupper(in)) {
            return in - 'A' + 11 + MORPHY_SPECIAL_STATE_PAD;
        }
        else {
            return in - 'a' + 11 + 26 + MORPHY_SPECIAL_STATE_PAD;
        }
    }
    else {
        dbg_printf("ERROR in mfl_convert_alphanum_to_shift_value(): input is not alphanumeric\n");
        return 0;
    }
}

void mfl_move_in_nexus_multistate(char **col)
{
    if (**col != '(' && **col != '{') {
        dbg_printf("WARNING in mfl_move_in_nexus_multistate(): *col is not start of multistate entry\n");
        return;
    }
    
    if (**col == '(') {
        do {
            ++(*col);
        } while (**col != ')');
    }
    else {
        do {
            ++(*col);
        } while (**col != '}');
    }
}

mfl_charstate mfl_convert_nexus_multistate(char *xstates, char *datype_converter)
{
    u_int64_t inbit = 0;
    int shift_val;
    mfl_charstate xstate_bits = 0;
    
    do {
        inbit = 1;
        shift_val = mfl_convert_nexus_symbol_to_shift_value(*xstates, datype_converter);
        inbit = inbit << shift_val;
        ++xstates;
    } while (xstates && *xstates != ')' && *xstates != '}');
    
    return xstate_bits;
}

int mfl_get_numstates_from_matrix(char *matrix)
{
    /* When no state symbols are specified and default reading is in effect, 
     * the number of unique symbols used can be counted directly from the matrix
     * */
}

mfl_charstate mfl_convert_single_symbol_categorical(char symbol)
{
    int shift_value;
    u_int64_t bit_setter;
    mfl_charstate charstate;

    

}

mfl_charstate mfl_convert_state_symbol_DNA(char state_symbol)
{
    
}

int mfl_calculate_data_partitions_required(mfl_handle_t)
{
    
}

mfl_datapartition_t* mfl_convert_and_partition_input_data(mfl_handle_t* mfl_handle)
{
    
    int *list_of_fitch_characters = NULL;
    int *list_of_wagner_characters = NULL;
    int *list_of_dollo_characters = NULL;
    int *list_of_irreversible_characters = NULL;
    
    
}