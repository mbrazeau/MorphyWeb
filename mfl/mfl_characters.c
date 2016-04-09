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


mfl_datapartition_t* mfl_convert_and_partition_input_data(mfl_handle_t* mfl_handle)
{
    
    int *list_of_fitch_characters = NULL;
    int *list_of_wagner_characters = NULL;
    int *list_of_dollo_characters = NULL;
    int *list_of_irreversible_characters = NULL;
    
    
}


bool mfl_is_nexus_stop_position(char a)
{
    bool is_stop = false;
    
    if (a) {
        if (isspace(a) || a == '-' || a == ';' || a == ',') {
            is_stop = true;
        }
    }
    else {
        is_stop = true;
    }
    
    return is_stop;
}


void mfl_skip_spaces(char **current)
{
    if (**current == ' ') {
        do {
            ++(*current);
        } while (isspace(**current));
    }
}

int mfl_read_nexus_type_int(char **current)
{
    int nexus_int = 0;
    
    // Read char until whitespace, a dash '-', semicolon, or end of string
    do {
        nexus_int = 10 * nexus_int;
        nexus_int = nexus_int + (**(current) - '0');
        ++(*current);
    } while (!mfl_is_nexus_stop_position(**current));
    
    mfl_skip_spaces(current);
    dbg_printf("Returning int %i from a Nexus sub-command\n", nexus_int);
    
    return nexus_int;
}


void mfl_set_include_value(int vectornum, bool includeval, bool* includes)
{
    /* In any set of bool include array, this sets the entry's value to the
     * desired true/false value */
    
    
    if (includeval) {
        includes[vectornum-1] = true;
    }
    else {
        includes[vectornum-1] = false;
    }
}


void mfl_set_include_range(int first, int last, bool includeval, bool* includes)
{
    int i = 0;
    
    for (i = first; i <= last; ++i) {
        mfl_set_include_value(i, includeval, includes);
    }
}


void mfl_move_current_to_digit(char** current)
{
    if (!isdigit(**current)) {
        do {
            ++(*current);
        } while (!isdigit(**current));
    }
}


void mfl_set_inclusion_list(bool* includes, bool includeval, int listmax, char *subcommand)
{
    int first = 0;
    int last = 0; // first and last for identifier ranges (e.g. taxa 1-4 or characters 11-17)
    char *current = NULL;
    current = subcommand;
    
    if (!isdigit(*current)) {
        mfl_move_current_to_digit(&current);
    }
    
    do {
        first = 0;
        last = 0;
        if (isdigit(*current)) {
            // Read token as either a number or value range.
            first = mfl_read_nexus_type_int(&current);
        }
        else {
            ++current;
        }
        
        if (*current == '-') {
            mfl_move_current_to_digit(&current);
            last = mfl_read_nexus_type_int(&current);
        }
        
        if (first) {
            if (last) {
                if (last > listmax) {
                    last = listmax;
                    
                    dbg_printf("WARNING in mfl_set_inclusion_list(): Proposed range is outside of list maximum. Proposed range will be truncated\n");
                }
                mfl_set_include_range(first, last, includeval, includes);
                
            }
            else {
                if (first <= listmax) {
                    mfl_set_include_value(first, includeval, includes);
                }
                else {
                    dbg_printf("WARNING in mfl_set_inclusion_list(): value will attempt to index outside of range. Ignorning proposed index\n");
                }
            }
        }
    } while (*current && *current != ';');
}


bool* mfl_alloc_character_inclusion_list(int num_chars)
{
    bool *new_inclusion_list = NULL;
    
    new_inclusion_list = (bool*)malloc(num_chars *sizeof(bool));
    if (!new_inclusion_list) {
        dbg_printf("ERROR in mfl_alloc_character_inclusion_list(): unable to allocate memory\n");
    }
    else {
        memset(new_inclusion_list, MORPHY_DEFAULT_CHARACTER_INCLUDE, num_chars * sizeof(bool));
    }
    
    return new_inclusion_list;
}


void mfl_free_inclusion_list(bool *inclist)
{
    free(inclist);
}


bool* mfl_read_nexus_exset_subcmd(char *subcommand, int Nexus_NCHARS)
{
    bool* includelist = NULL;
    
    includelist = mfl_alloc_character_inclusion_list(Nexus_NCHARS);
    subcommand = mfl_move_past_eq_sign(subcommand);
    mfl_set_inclusion_list(includelist, false, Nexus_NCHARS, subcommand);
    
    return includelist;
}


int mfl_calculate_data_partitions_required(mfl_handle_t)
{
    
}


void tui_test_character_including()
{
    int i = 0;
    int num_chars = 20;
    char subcmd1[] = "ExSet * Exclude= 1-5 8 17;";
    char subcmd2[] = "Exclude = 1-5, 8 17;";
    char subcmd3[] = "18-51 100";
    char subcmd4[] = "10-15 2-6";
    
    char *subcmd = NULL;
    subcmd = subcmd3;
    dbg_printf("Processing the following subcommand: %s\n\n", subcmd);
    
    bool* includelist = mfl_read_nexus_exset_subcmd(subcmd, num_chars);
    
    dbg_printf("\nIncluding characters:\n");
    for (i = 0; i < num_chars; ++i) {
        if (includelist[i]) {
            dbg_printf("%i ", i + 1);
        }
    }
    dbg_printf("\n");
    
    dbg_printf("\nExcluding characters:\n");
    for (i = 0; i < num_chars; ++i) {
        if (!includelist[i]) {
            dbg_printf("%i ", i + 1);
        }
    }
    dbg_printf("\n");
    
    mfl_free_inclusion_list(includelist);
}