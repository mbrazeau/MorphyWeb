#include <sys/stat.h>
#include "nui.h"

/* 
 * Use the preprocessor to define a derived class for each command,
 * the derived class implements what some people would call a functionoid.
 * A functionoid is basically a pure virtual function in a base class
 * that is overloaded a bunch of times in derived classes, with each
 * derived class implementing different functionality. In this case
 * the different functionality is all of the commands... The overloaded
 * functioniod MenuFunction simply calls a function in the CNexusUserInterface
 * class, that function is different for each derived class and is always takes
 * the name of the "type" parameter (Ex: Spacer, OpenNexusFile, SaveFile, etc...)
 *
 * The behavior of the MenuFunction is defined to return true if the
 * program is to continue to operate, and false if the program should
 * terminate. Possibly due to errors, or if the user selects the quit 
 * command.
 *
 * Note: The preprocessor ## concatenation operator is used to create
 * the class name based on the value of the "type" macro parameter
 */
#define NEW_COMMAND_DEFINE(type) \
    class type : public CNexusMenuBase \
    { \
    public:\
        type(const char * strCommand, const char * strHelpText) : CNexusMenuBase(strCommand, strHelpText){}\
        bool MenuFunction(CNexusUserInterface *pNexusUserInterface)\
        {\
            return pNexusUserInterface->f##type();\
        }\
    };

/*
 * The followind actually defines the derived class for each command on the main menu
 */
NEW_COMMAND_DEFINE(CNexusMenuSpacer         )
NEW_COMMAND_DEFINE(CNexusMenuOpenNexusFile  )
NEW_COMMAND_DEFINE(CNexusMenuSaveFile       )
NEW_COMMAND_DEFINE(CNexusMenuCloseNexusFile )

NEW_COMMAND_DEFINE(CNexusMenuHelp           )
NEW_COMMAND_DEFINE(CNexusMenuQuit           )
NEW_COMMAND_DEFINE(CNexusMenuAbout          )
NEW_COMMAND_DEFINE(CNexusMenuCommandLog     )
NEW_COMMAND_DEFINE(CNexusMenuStatus         )
NEW_COMMAND_DEFINE(CNexusMenuChdir          )

NEW_COMMAND_DEFINE(CNexusMenuExclude        )
NEW_COMMAND_DEFINE(CNexusMenuInclude        )
NEW_COMMAND_DEFINE(CNexusMenuOutgroup       )
NEW_COMMAND_DEFINE(CNexusMenuIngroup        )
NEW_COMMAND_DEFINE(CNexusMenuChar           )
NEW_COMMAND_DEFINE(CNexusMenuSet            )

NEW_COMMAND_DEFINE(CNexusMenuHeuristicSearch)
NEW_COMMAND_DEFINE(CNexusMenuExhaust        )
NEW_COMMAND_DEFINE(CNexusMenuBNB            )
NEW_COMMAND_DEFINE(CNexusMenuBootstrap      )
NEW_COMMAND_DEFINE(CNexusMenuJackknife      )
NEW_COMMAND_DEFINE(CNexusMenuSTR            )

NEW_COMMAND_DEFINE(CNexusMenuConsens        )
NEW_COMMAND_DEFINE(CNexusMenuCollapse       )
NEW_COMMAND_DEFINE(CNexusMenuReport         )

/*
 * The followind actually defines the derived class for each command on the set menu
 */
NEW_COMMAND_DEFINE(CNexusMenuSearchType     )
NEW_COMMAND_DEFINE(CNexusMenuBranchSwapType )
NEW_COMMAND_DEFINE(CNexusMenuAddSeqType     )
NEW_COMMAND_DEFINE(CNexusMenuCollapseAt     )
NEW_COMMAND_DEFINE(CNexusMenuCollapseZero   )
NEW_COMMAND_DEFINE(CNexusMenuNumIterations  )
NEW_COMMAND_DEFINE(CNexusMenuTreeLimit      )
NEW_COMMAND_DEFINE(CNexusMenuRatchetSearch  )

NEW_COMMAND_DEFINE(CNexusMenuMainMenu       )

/*
 * The UI constructor puts the menu together, it stores each
 * menu option in a stl vector.
 */
CNexusUserInterface::CNexusUserInterface()
{
    m_mflHandle = NULL;
    m_pNexusParse = NULL;
    m_strCwd = "./";
    m_pMainMenu = new CNexusMenuData("Main Menu");
    if (!m_pMainMenu)
    {
        throw "Unable to allocate memory for main menu";
    }
    m_pMainMenu->AddMenuItem(new CNexusMenuSpacer      (NULL, "File"));
    m_pMainMenu->AddMenuItem(new CNexusMenuOpenNexusFile   ("O", "Open a nexus file"));
    m_pMainMenu->AddMenuItem(new CNexusMenuCloseNexusFile  ("C", "Close a nexus file"));
    m_pMainMenu->AddMenuItem(new CNexusMenuSaveFile        ("S", "Save according to the options"));

    m_pMainMenu->AddMenuItem(new CNexusMenuSpacer      (NULL, "Program"));
    m_pMainMenu->AddMenuItem(new CNexusMenuCommandLog      ("LOG" ,"Toggles record of commands, variable states, etc"));
    m_pMainMenu->AddMenuItem(new CNexusMenuStatus          ("STAT","Prints status of all current settings, eg. logmode on/off"));
    m_pMainMenu->AddMenuItem(new CNexusMenuChdir           ("CD"  ,"change working directory"));
    m_pMainMenu->AddMenuItem(new CNexusMenuHelp            ("H"   , "Help"));
    m_pMainMenu->AddMenuItem(new CNexusMenuAbout           ("A"   , "About"));
    m_pMainMenu->AddMenuItem(new CNexusMenuQuit            ("Q"   , "Quit"));

    m_pMainMenu->AddMenuItem(new CNexusMenuSpacer      (NULL, "Data"));
    m_pMainMenu->AddMenuItem(new CNexusMenuExclude         ("EXC" , "Exclude taxa or characters"));
    m_pMainMenu->AddMenuItem(new CNexusMenuInclude         ("INC" , "Include excluded taxa or characters"));
    m_pMainMenu->AddMenuItem(new CNexusMenuOutgroup        ("OUTG", "Assign taxa to outgroup"));
    m_pMainMenu->AddMenuItem(new CNexusMenuIngroup         ("ING" , "Return taxa from outgroup to ingrou"));
    m_pMainMenu->AddMenuItem(new CNexusMenuChar            ("CHAR", "Modify a character's type"));
    m_pMainMenu->AddMenuItem(new CNexusMenuSet             ("SET" , "Set morphy configuration parameters"));

    m_pMainMenu->AddMenuItem(new CNexusMenuSpacer      (NULL, "Analysis"));
    m_pMainMenu->AddMenuItem(new CNexusMenuHeuristicSearch ("HS" , "Begin a heuristic search"));
    m_pMainMenu->AddMenuItem(new CNexusMenuExhaust         ("EXS", "Begin an exhaustive search"));
    m_pMainMenu->AddMenuItem(new CNexusMenuBNB             ("BNB", "Begin a branch-and-bound search"));
    m_pMainMenu->AddMenuItem(new CNexusMenuBootstrap       ("BTS", "Begin a bootstrap analysis"));
    m_pMainMenu->AddMenuItem(new CNexusMenuJackknife       ("JK" , "Begin a jackknife analysis"));
    m_pMainMenu->AddMenuItem(new CNexusMenuSTR             ("STR", "Perform a safe taxonomic reduction"));

    m_pMainMenu->AddMenuItem(new CNexusMenuSpacer      (NULL, "Results"));
    m_pMainMenu->AddMenuItem(new CNexusMenuConsens         ("CONSENS" , "Compute consensus tree for trees in memory"));
    m_pMainMenu->AddMenuItem(new CNexusMenuCollapse        ("COLLAPSE", "Collapse zero-length branches, condense the tree set"));
    m_pMainMenu->AddMenuItem(new CNexusMenuReport          ("REPORT"  , "Print a report about the current open nexus file"));

    m_pSetMenu = new CNexusMenuData("Set parameters");
    if (!m_pSetMenu)
    {
        throw "Unable to allocate memory for set menu";
    }
    m_pSetMenu->AddMenuItem(new CNexusMenuSearchType       ("SEARCHTYPE", "Set the search type for JK and BTS searches"));
    m_pSetMenu->AddMenuItem(new CNexusMenuBranchSwapType   ("BRANCHSWAP", "Set branch swap type for heuristic searches"));
    m_pSetMenu->AddMenuItem(new CNexusMenuAddSeqType       ("ADDSEQ"    , "Selects the manner in which branches are added during the generation of starting trees"));
    m_pSetMenu->AddMenuItem(new CNexusMenuCollapseAt       ("COLLAPSEAT", "Configure when to collapse nodes"));
    m_pSetMenu->AddMenuItem(new CNexusMenuCollapseZero     ("COLLAPSEZERO", "Enable collapsing of zero length branches during search"));
    m_pSetMenu->AddMenuItem(new CNexusMenuNumIterations    ("NUMITE"    , "Set the number of iterations for a heuristic search"));
    m_pSetMenu->AddMenuItem(new CNexusMenuTreeLimit        ("TREELIMIT" , "Set the maximum number of trees allowed to be stored in memory"));
    m_pSetMenu->AddMenuItem(new CNexusMenuRatchetSearch    ("RATCHET"   , "Set the ratchet search parameter"));

    m_pSetMenu->AddMenuItem(new CNexusMenuHelp             ("H"         , "Help"));
    m_pSetMenu->AddMenuItem(new CNexusMenuMainMenu         ("Q"         , "Return to main menu"));
}

void CNexusUserInterface::ChangeMenu(CNexusMenuData *pMenu)
{
    m_pMenu = pMenu;
    fCNexusMenuHelp(false);
}


/*
 * The UI destructor deletes all memory.
 */
CNexusUserInterface::~CNexusUserInterface()
{
    delete m_pMainMenu;
    delete m_pSetMenu;
    fCNexusMenuCloseNexusFile(false);
    if (m_fCommandLog)
    {
        m_fCommandLog.close();
    }
    DestroyHandle();
}

/*
 * Anytime the user is to give input, the input must come through this function.
 */
void CNexusUserInterface::GetUserInput(string strPrompt, string *strInput)
{
    cout<<strPrompt;
    cin>>*strInput;
    if (cin.fail())
    {
        *strInput = "q";
        throw "Stardard input failure";
    }
    if (m_fCommandLog)
    {
        m_fCommandLog<<*strInput<<endl;
    }
}

/*
 * Actually read input from the user, and issue the selected commands
 */
void CNexusUserInterface::DoMenu()
{
    string strInput;
    fCNexusMenuAbout(false);
    ChangeMenu(m_pMainMenu);
    try
    {
        do
        {
            strInput.clear();
            GetUserInput(m_pMenu->GetPrompt() ,&strInput);
        } while (m_pMenu->RunSelection(strInput, this));
    }
    catch (const char *e)
    {
        cout<<endl<<"Error: "<<e<<endl;
    }
}

bool CNexusUserInterface::SetMorphyOpenParams()
{
    stringstream ss;
    int nTax = m_pNexusParse->m_cTaxa->GetNTax();
    int nChar = m_pNexusParse->m_cChars->GetNCharTotal();
    int i;

    mfl_set_parameter(m_mflHandle, MFL_PT_NUM_TAX, (void*)nTax);
    mfl_set_parameter(m_mflHandle, MFL_PT_NUM_CHAR, (void*)nChar);
    for (i = 0; i < nTax; i++)
    {
        m_pNexusParse->m_cChars->WriteStatesForTaxonAsNexus(ss, i, 0, nChar);
    }
    mfl_set_parameter(m_mflHandle, MFL_PT_INPUT_DATA, (void*)ss.str().c_str());
    return true;
}

/*
 * Run the open file command, just prompt the user for input
 * and attempt to read the nexus file
 */
bool CNexusUserInterface::fCNexusMenuOpenNexusFile()
{
    string strFilename;

    if (!m_pNexusParse)
    {
        GetUserInput(" Enter filename: " + m_strCwd, &strFilename);
        strFilename = m_strCwd + strFilename;
        m_pNexusParse = new CNexusParse();
        if (m_pNexusParse)
        {
            CreateHandle();
            if (m_pNexusParse->ReadNexusFile(&strFilename, NULL))
            {
                SetMorphyOpenParams();
                cout<<" "<<strFilename<<" open successfully"<<endl;
            }
            else
            {
                fCNexusMenuCloseNexusFile(false);
                cout<<" Error: Unable to read "<<strFilename<<endl;
            }
        }
        else
        {
            cout<<" Error: Unable to open "<<strFilename<<endl;
        }
    }
    else
    {
        cout<<" Error: Nexus file "<<m_pNexusParse->m_cNexusReader->GetInFileName()<<" is already open"<<endl;
    }
    return true;
}

void CNexusUserInterface::DestroyHandle()
{
    if (m_mflHandle)
    {
        mfl_destroy_handle(m_mflHandle);
        m_mflHandle = NULL;
    }
}

void CNexusUserInterface::CreateHandle()
{
    DestroyHandle();
    m_mflHandle = mfl_create_handle();
}

bool CNexusUserInterface::fCNexusMenuSaveFile       ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

/*
 * Close a file opened with the OpenNexusFile command
 * bVerbose is defaulted to = true
 */
bool CNexusUserInterface::fCNexusMenuCloseNexusFile(bool bVerbose)
{
    if (m_pNexusParse)
    {
        delete m_pNexusParse;
        m_pNexusParse = NULL;
        if (bVerbose)
        {
            cout<<" Successfully closed file..."<<endl;
        }
    }
    else if (bVerbose)
    {
        cout<<" Error: No file is currently open"<<endl;
    }
    DestroyHandle();

    return true;
}

/*
 * Print the menu
 */
bool CNexusUserInterface::fCNexusMenuHelp           (bool bForceShowMenu)
{
    m_pMenu->Help(bForceShowMenu);
    return true;
}

bool CNexusUserInterface::fCNexusMenuQuit           ()
{
    cout<<"Goodbye!"<<endl;
    return false;
}

bool CNexusUserInterface::fCNexusMenuAbout          (bool bShowBuildTime)
{
    cout<<"Morphy NUI Version: "<<NUI_MAJOR_VERSION<<"."<<NUI_MINOR_VERSION<<endl;
    cout<<"Copyright 2012 (C) Martin Brazeau and Chris Desjardins. All rights reserved."<<endl;
    cout<<"This program uses the NCL by Paul O. Lewis."<<endl;
    if (bShowBuildTime)
    {
        cout<<"Build time: "<<__DATE__<<" "<<__TIME__<<endl;
    }
    return true;
}

bool CNexusUserInterface::fCNexusMenuCommandLog     ()
{
    string strFilename;
    
    if (!m_fCommandLog)
    {
        GetUserInput(" Enter log filename: " + m_strCwd, &strFilename);
        strFilename = m_strCwd + strFilename;
       
        m_fCommandLog.open(strFilename.c_str());
        if (m_fCommandLog)
        {
            cout<<" Successfully opened '"<<strFilename<<"'"<<endl;
        }
        else
        {
            cout<<" Error: Unable to open '"<<strFilename<<"'"<<endl;
        }
    }
    else
    {
        cout<<" Error: Log file "<<m_fCommandLog.GetFileName()<<" is already open"<<endl;
    }
    return true;
}

bool CNexusUserInterface::fCNexusMenuStatus         ()
{
    string strNexusFile = "File not open";

    if ((m_pNexusParse) && (m_pNexusParse->m_cNexusReader))
    {
        strNexusFile = m_pNexusParse->m_cNexusReader->GetInFileName();
    }
    cout<<endl;
    cout<<"Nexus file: "<<strNexusFile<<endl;
    cout<<"Command log file: "<<m_fCommandLog.GetFileName()<<endl;
    cout<<"Working directory: "<<m_strCwd<<endl;
    cout<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuChdir          ()
{
    string strCwd;
    struct stat st;

    GetUserInput(" Enter new working directory: ", &strCwd);
    if (strCwd[strCwd.length() - 1] != '/')
    {
        strCwd += "/";
    }
    if ((stat(strCwd.c_str(), &st) == 0) && (S_ISDIR(st.st_mode)))
    {
        cout<<" Setting working directory to: "<<strCwd<<endl;
        m_strCwd = strCwd;
    }
    else
    {
        cout<<" Error: Invalid directory '"<<strCwd<<"'"<<endl;
    }
    return true;
}

bool CNexusUserInterface::fCNexusMenuExclude        ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuInclude        ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuOutgroup       ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuIngroup        ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuChar           ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuSet            ()
{
    ChangeMenu(m_pSetMenu);
    return true;
}

bool CNexusUserInterface::fCNexusMenuHeuristicSearch()
{
    return mfl_heuristic(m_mflHandle);
}

bool CNexusUserInterface::fCNexusMenuExhaust        ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuBNB            ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuBootstrap      ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuJackknife      ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuSTR            ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuConsens        ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuCollapse       ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuReport()
{
    if (m_pNexusParse)
    {
        m_pNexusParse->Report();
    }
    else
    {
        cout<<"No file is currently open"<<endl;
    }
    return true;
}

bool CNexusUserInterface::fCNexusMenuSearchType     ()
{
    string strInput;
    mfl_search_t search_type = MFL_ST_MAX;
    GetUserInput(" Enter search type\n  1) Exhaustive\n  2) Branch Bound\n  3) Heuristic\n # ", &strInput);
    if (strInput.compare("1") == 0)
    {
        search_type = MFL_ST_EXHAUSTIVE;
    }
    else if (strInput.compare("2") == 0)
    {
        search_type = MFL_ST_BRANCH_BOUND;
    }
    else if (strInput.compare("3") == 0)
    {
        search_type = MFL_ST_HEURISTIC;
    }
    mfl_set_parameter(m_mflHandle, MFL_PT_SEARCH_TYPE, (void*)search_type);
    return true;
}

bool CNexusUserInterface::fCNexusMenuBranchSwapType ()
{
    string strInput;
    mfl_branch_swap_t swap_type = MFL_BST_MAX;
    GetUserInput(" Enter branch swap type\n  1) Tree bisection\n  2) Subtree pruning\n  3) Nearist neighbor\n # ", &strInput);
    if (strInput.compare("1") == 0)
    {
        swap_type = MFL_BST_TBR;
    }
    else if (strInput.compare("2") == 0)
    {
        swap_type = MFL_BST_SPR;
    }
    else if (strInput.compare("3") == 0)
    {
        swap_type = MFL_BST_NNI;
    }
    mfl_set_parameter(m_mflHandle, MFL_PT_BRANCH_SWAP_TYPE, (void*)swap_type);
    return true;
}

bool CNexusUserInterface::fCNexusMenuAddSeqType     ()
{
    string strInput;
    mfl_add_sequence_t add_seq_type = MFL_AST_MAX;
    GetUserInput(" Enter add seq type\n  1) Simple\n  2) Random\n  3) As is\n  4) Closest\n # ", &strInput);
    if (strInput.compare("1") == 0)
    {
        add_seq_type = MFL_AST_SIMPLE;
    }
    else if (strInput.compare("2") == 0)
    {
        add_seq_type = MFL_AST_RANDOM;
    }
    else if (strInput.compare("3") == 0)
    {
        add_seq_type = MFL_AST_ASIS;
    }
    else if (strInput.compare("4") == 0)
    {
        add_seq_type = MFL_AST_CLOSEST;
    }
    mfl_set_parameter(m_mflHandle, MFL_PT_ADD_SEQUENCE_TYPE, (void*)add_seq_type);
    return true;
}

bool CNexusUserInterface::fCNexusMenuCollapseAt           ()
{
    string strInput;
    mfl_set_collapse_at_t collapse_at = MFL_SC_MAX;
    GetUserInput(" Enter when to collapse branches\n  1) When max length is 0\n  2) When min length is 0\n  3) When two incident nodes with equal apomorphies reconstruction sets\n # ", &strInput);
    if (strInput.compare("1") == 0)
    {
        collapse_at = MFL_SC_MAX_LEN;
    }
    else if (strInput.compare("2") == 0)
    {
        collapse_at = MFL_SC_MIN_LEN;
    }
    else if (strInput.compare("3") == 0)
    {
        collapse_at = MFL_SC_EQUAL_RECONSTRUCTION_SETS;
    }

    mfl_set_parameter(m_mflHandle, MFL_PT_COLLAP_AT, (void*)collapse_at);
    return true;
}

bool CNexusUserInterface::fCNexusMenuCollapseZero         ()
{
    string strInput;
    bool collapse_zero = false;
    GetUserInput(" Select if zero length branches should be collapsed during search\n  1) Collapse\n  2) Do not collapse\n # ", &strInput);
    if (strInput.compare("1") == 0)
    {
        collapse_zero = true;
    }
    mfl_set_parameter(m_mflHandle, MFL_PT_COLLAPSE, (void*)collapse_zero);
    return true;
}

bool CNexusUserInterface::fCNexusMenuNumIterations        ()
{
    string strInput;
    int num_iterations;
    stringstream ss;
    GetUserInput(" Enter number of iterations# ", &strInput);
    ss<<strInput.c_str();
    ss>>num_iterations;
    mfl_set_parameter(m_mflHandle, MFL_PT_NUM_ITERATIONS, (void*)num_iterations);
    return true;
}

bool CNexusUserInterface::fCNexusMenuTreeLimit        ()
{
    string strInput;
    int tree_limit;
    stringstream ss;
    GetUserInput(" Enter number of trees to store# ", &strInput);
    ss<<strInput.c_str();
    ss>>tree_limit;
    mfl_set_parameter(m_mflHandle, MFL_PT_TREELIMIT, (void*)tree_limit);
    return true;
}

bool CNexusUserInterface::fCNexusMenuRatchetSearch        ()
{
    string strInput;
    bool ratchet = false;
    GetUserInput(" Select how the ratchet search parameter should be set\n  1) Enabled\n  2) Disabled\n # ", &strInput);
    if (strInput.compare("1") == 0)
    {
        ratchet = true;
    }
    mfl_set_parameter(m_mflHandle, MFL_PT_RATCHET_SEARCH, (void*)ratchet);
    return true;
}

bool CNexusUserInterface::fCNexusMenuMainMenu       ()
{
    ChangeMenu(m_pMainMenu);
    return true;
}

int main(int argc, char *argv[])
{
    CNexusUserInterface ui;
    ui.DoMenu();
    
    return 0;
}

