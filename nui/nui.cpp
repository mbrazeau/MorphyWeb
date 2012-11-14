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
        type(const char * strCommand, const char * strHelpText, vector<string> assignments) : CNexusMenuBase(strCommand, strHelpText, assignments){}\
        type(const char * strCommand, const char * strHelpText, vector<int>    assignments) : CNexusMenuBase(strCommand, strHelpText, assignments){}\
        bool MenuFunction(CNexusUserInterface *pNexusUserInterface, string *value)\
        {\
            return pNexusUserInterface->f##type(value);\
        }\
    };

#define MAKE_STR_VECTOR(slist)  vector<string>(slist, slist + (sizeof(slist) / sizeof(slist[0])))
#define MAKE_INT_VECTOR(slist)  vector<int>(slist, slist + (sizeof(slist) / sizeof(slist[0])))

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

NEW_COMMAND_DEFINE(CNexusMenuHeuristicSearch)
NEW_COMMAND_DEFINE(CNexusMenuExhaust        )
NEW_COMMAND_DEFINE(CNexusMenuBNB            )
NEW_COMMAND_DEFINE(CNexusMenuBootstrap      )
NEW_COMMAND_DEFINE(CNexusMenuJackknife      )
NEW_COMMAND_DEFINE(CNexusMenuSTR            )

NEW_COMMAND_DEFINE(CNexusMenuConsens        )
//NEW_COMMAND_DEFINE(CNexusMenuCollapse       )
NEW_COMMAND_DEFINE(CNexusMenuReport         )

/*
 * The following actually defines the derived class for each command on the set menu
 */
NEW_COMMAND_DEFINE(CNexusMenuSearchType     )
NEW_COMMAND_DEFINE(CNexusMenuBranchSwapType )
NEW_COMMAND_DEFINE(CNexusMenuAddSeqType     )
NEW_COMMAND_DEFINE(CNexusMenuCollapseAt     )
NEW_COMMAND_DEFINE(CNexusMenuCollapseZero   )
NEW_COMMAND_DEFINE(CNexusMenuNumIterations  )
NEW_COMMAND_DEFINE(CNexusMenuTreeLimit      )
NEW_COMMAND_DEFINE(CNexusMenuRatchetSearch  )
NEW_COMMAND_DEFINE(CNexusMenuGap            )

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
    m_pMainMenu = new CNexusMenuData("Main Menu\nEnter selection#");
    if (!m_pMainMenu)
    {
        throw "Unable to allocate memory for main menu";
    }
    m_pMainMenu->AddMenuItem(new CNexusMenuSpacer      (NULL, "File"));
    m_pMainMenu->AddMenuItem(new CNexusMenuOpenNexusFile   ("open"          , "Open a nexus file"));
    m_pMainMenu->AddMenuItem(new CNexusMenuCloseNexusFile  ("close"         , "Close a nexus file"));
    m_pMainMenu->AddMenuItem(new CNexusMenuSaveFile        ("save"          , "Save according to the options"));

    m_pMainMenu->AddMenuItem(new CNexusMenuSpacer      (NULL, "Program"));
    m_pMainMenu->AddMenuItem(new CNexusMenuCommandLog      ("log"           ,"Toggles record of commands, variable states, etc"));
    m_pMainMenu->AddMenuItem(new CNexusMenuStatus          ("status"        ,"Prints status of all current settings, eg. logmode on/off"));
    m_pMainMenu->AddMenuItem(new CNexusMenuChdir           ("cd"            ,"change working directory"));
    m_pMainMenu->AddMenuItem(new CNexusMenuHelp            ("help"          , "Help"));
    m_pMainMenu->AddMenuItem(new CNexusMenuAbout           ("about"         , "About"));
    m_pMainMenu->AddMenuItem(new CNexusMenuQuit            ("quit"          , "Quit"));

    m_pMainMenu->AddMenuItem(new CNexusMenuSpacer      (NULL, "Data"));
    m_pMainMenu->AddMenuItem(new CNexusMenuExclude         ("exclude"       , "Exclude taxa or characters"));
    m_pMainMenu->AddMenuItem(new CNexusMenuInclude         ("include"       , "Include excluded taxa or characters"));
    m_pMainMenu->AddMenuItem(new CNexusMenuOutgroup        ("outGroup"      , "Assign taxa to outgroup"));
    m_pMainMenu->AddMenuItem(new CNexusMenuIngroup         ("inGroup"       , "Return taxa from outgroup to ingroup"));
    m_pMainMenu->AddMenuItem(new CNexusMenuChar            ("char"          , "Modify a character's type"));

    m_pMainMenu->AddMenuItem(new CNexusMenuSpacer      (NULL, "Analysis"));
    m_pMainMenu->AddMenuItem(new CNexusMenuHeuristicSearch ("heuristic"     , "Begin a heuristic search"));
    m_pMainMenu->AddMenuItem(new CNexusMenuExhaust         ("exhaustive"    , "Begin an exhaustive search"));
    m_pMainMenu->AddMenuItem(new CNexusMenuBNB             ("branchbound"   , "Begin a branch-and-bound search"));
    m_pMainMenu->AddMenuItem(new CNexusMenuBootstrap       ("bootStrap"     , "Begin a bootstrap analysis"));
    m_pMainMenu->AddMenuItem(new CNexusMenuJackknife       ("jackKnife"     , "Begin a jackknife analysis"));
    m_pMainMenu->AddMenuItem(new CNexusMenuSTR             ("reduction"     , "Perform a safe taxonomic reduction"));

    m_pMainMenu->AddMenuItem(new CNexusMenuSpacer      (NULL, "Results"));
    m_pMainMenu->AddMenuItem(new CNexusMenuConsens         ("consensus"     , "Compute consensus tree for trees in memory"));
    //m_pMainMenu->AddMenuItem(new CNexusMenuCollapse        ("collapse"      , "Collapse zero-length branches, condense the tree set"));
    m_pMainMenu->AddMenuItem(new CNexusMenuReport          ("report"        , "Print a report about the current open nexus file"));

    const char* SearchType    [] = {"Exhaustive", "BranchBound", "Heuristic"};
    const char* BranchSwapType[] = {"TreeBisection", "SubtreePruning", "NearistNeighbor"};
    const char* AddSeqType    [] = {"Simple", "Random", "AsIs", "Closest"};
    const char* CollapseAt    [] = {"MaxIs0", "MinIs0", "Equal"};
    const char* CollapseZero  [] = {"Collapse", "NoCollapse"};
    const int   NumIterations [] = {0, 10000000};
    const int   TreeLimit     [] = {0, 10000000};
    const char* RatchetSearch [] = {"Enable", "Disable"};
    const char* Gap           [] = {"Inapplicable", "Missing"};

    m_pMainMenu->AddMenuItem(new CNexusMenuSpacer      (NULL, "Parameters"));
    m_pMainMenu->AddMenuItem(new CNexusMenuSearchType       ("searchType"    , "Set the search type for JK and BTS searches", MAKE_STR_VECTOR(SearchType)));
    m_pMainMenu->AddMenuItem(new CNexusMenuBranchSwapType   ("branchSwap"    , "Set branch swap type for heuristic searches", MAKE_STR_VECTOR(BranchSwapType)));
    m_pMainMenu->AddMenuItem(new CNexusMenuAddSeqType       ("addSeq"        , "Selects the manner in which branches are added during the generation of starting trees", MAKE_STR_VECTOR(AddSeqType)));
    m_pMainMenu->AddMenuItem(new CNexusMenuCollapseAt       ("collapseAt"    , "Configure when to collapse nodes", MAKE_STR_VECTOR(CollapseAt)));
    m_pMainMenu->AddMenuItem(new CNexusMenuCollapseZero     ("collapseZero"  , "Enable collapsing of zero length branches during search", MAKE_STR_VECTOR(CollapseZero)));
    m_pMainMenu->AddMenuItem(new CNexusMenuNumIterations    ("numite"        , "Set the number of iterations for a heuristic search", MAKE_INT_VECTOR(NumIterations)));
    m_pMainMenu->AddMenuItem(new CNexusMenuTreeLimit        ("treeLimit"     , "Set the maximum number of trees allowed to be stored in memory", MAKE_INT_VECTOR(TreeLimit)));
    m_pMainMenu->AddMenuItem(new CNexusMenuRatchetSearch    ("ratchet"       , "Set the ratchet search parameter", MAKE_STR_VECTOR(RatchetSearch)));
    m_pMainMenu->AddMenuItem(new CNexusMenuGap              ("gap"           , "Set whether gap symbol ('-') will be treated as inapplicability or as missing data", MAKE_STR_VECTOR(Gap)));

    m_ioCommands = new CEditLineHist("nui1234567890", &m_fCommandLog);
    m_ioInputFiles = new CEditLineHist("nui12345678901234567890", &m_fCommandLog);
    m_ioLogFiles = new CEditLineHist("nui1234567890", &m_fCommandLog);
    m_ioWorkingDir = new CEditLineHist("nui1234567890", &m_fCommandLog);
    m_ioNumericSubCommands = new CEditLineHist("nui1234567890", &m_fCommandLog);
}

void CNexusUserInterface::ChangeMenu(CNexusMenuData *pMenu)
{
    m_pMenu = pMenu;
    fCNexusMenuHelp(NULL, false);
}

void CNexusUserInterface::Delete(CEditLineHist *pMem)
{
    if (pMem)
    {
        delete (pMem);
    }
}

/*
 * The UI destructor deletes all memory.
 */
CNexusUserInterface::~CNexusUserInterface()
{
    delete m_pMainMenu;
    fCNexusMenuCloseNexusFile(NULL, false);
    if (m_fCommandLog)
    {
        m_fCommandLog.close();
    }
    DestroyHandle();
    Delete(m_ioCommands);
    Delete(m_ioInputFiles);
    Delete(m_ioLogFiles);
    Delete(m_ioWorkingDir);
    Delete(m_ioNumericSubCommands);
}

/*
 * Actually read input from the user, and issue the selected commands
 */
void CNexusUserInterface::DoMenu()
{
    bool ret;
    string strInput;
    fCNexusMenuAbout(NULL, false);
    ChangeMenu(m_pMainMenu);
    do
    {
        try
        {
            strInput.clear();
            m_ioCommands->GetUserInput(m_pMenu->GetPrompt(), &strInput);
            ret = m_pMenu->RunSelection(strInput, this);
        }
        catch (const char *e)
        {
            cout<<endl<<"Error: "<<e<<endl;
        }
        catch (ifstream::failure &e)
        {
            cout<<endl<<"Error: "<<e.what()<<endl;
        }
        catch (NxsException &e)
        {
            cout<<endl<<"Error: "<<e.nxs_what()<<endl;
        }
        catch (const std::exception &e)
        {
            cout<<endl<<"Error: "<<e.what()<<endl;
        }
        catch (...)
        {
            cout<<endl<<"Error: Unknown exception"<<endl;
        }
    } while (ret);
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
bool CNexusUserInterface::fCNexusMenuOpenNexusFile(string *value)
{
    string strFilename;

    if (!m_pNexusParse)
    {
        m_ioInputFiles->GetUserInput(" Enter filename: " + m_strCwd, &strFilename);
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
                fCNexusMenuCloseNexusFile(NULL, false);
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

bool CNexusUserInterface::SaveTranslateTable(myofstream &fSave)
{
    int it;
    int n_taxa = m_pNexusParse->m_cTaxa->GetNumTaxonLabels();
    for (it = 0; it < n_taxa; it++) 
    {
        fSave<<"\t\t"<<it + 1<<" ";
        if (m_pNexusParse->m_cTaxa->NeedsQuotes(it)) 
        {
            fSave<<"'"<<m_pNexusParse->m_cTaxa->GetTaxonLabel(it)<<"'";
        }
        else 
        {
            fSave<<m_pNexusParse->m_cTaxa->GetTaxonLabel(it);
        }
        if (it < n_taxa-1) 
        {
            fSave<<",";
        }
        fSave<<endl;
    }
    return true;
}

bool CNexusUserInterface::SaveNewickStrings(myofstream &fSave)
{
    int it;
    char **newicks = mfl_get_saved_trees_newick(m_mflHandle);
    for (it = 0; newicks[it]; it++)
    {
        fSave<<"\t\tTREE Morphy_"<<it+1<<" = "<< newicks[it]<<endl;
    }

    return true;
}

bool CNexusUserInterface::fCNexusMenuSaveFile       (string *value)
{
    string strFilename;
    myofstream fSave;

    m_ioLogFiles->GetUserInput(" Enter save filename: " + m_strCwd, &strFilename);
    strFilename = m_strCwd + strFilename;
    fSave.open(strFilename.c_str());
    if (fSave)
    {
        fSave<<"#NEXUS"<<endl;
        fSave<<"BEGIN TREES;"<<endl<<"\tTRANSLATE"<<endl;

        SaveTranslateTable(fSave);

        fSave<<"\t\t;"<<endl<<endl;

        SaveNewickStrings(fSave);

        fSave<<"END;";
        fSave.close();
        cout<<" Successfully opened '"<<strFilename<<"'"<<endl;
    }
    else
    {
        cout<<" Error: Unable to open '"<<strFilename<<"'"<<endl;
    }
    return true;
}

/*
 * Close a file opened with the OpenNexusFile command
 * bVerbose is defaulted to = true
 */
bool CNexusUserInterface::fCNexusMenuCloseNexusFile(string *value, bool bVerbose)
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
bool CNexusUserInterface::fCNexusMenuHelp           (string *value, bool bForceShowMenu)
{
    m_pMenu->Help(bForceShowMenu);
    return true;
}

bool CNexusUserInterface::fCNexusMenuQuit           (string *value)
{
    cout<<"Goodbye!"<<endl;
    return false;
}

bool CNexusUserInterface::fCNexusMenuAbout          (string *value, bool bShowBuildTime)
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

bool CNexusUserInterface::fCNexusMenuCommandLog     (string *value)
{
    string strFilename;
    
    if (!m_fCommandLog)
    {
        m_ioLogFiles->GetUserInput(" Enter log filename: " + m_strCwd, &strFilename);
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

bool CNexusUserInterface::fCNexusMenuStatus         (string *value)
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

bool CNexusUserInterface::fCNexusMenuChdir          (string *value)
{
    string strCwd;
    struct stat st;

    m_ioWorkingDir->GetUserInput(" Enter new working directory: ", &strCwd);
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

bool CNexusUserInterface::fCNexusMenuExclude        (string *value)
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuInclude        (string *value)
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuOutgroup       (string *value)
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuIngroup        (string *value)
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuChar           (string *value)
{
    cout<<"Not implemented"<<endl;
    return true;
}

void CNexusUserInterface::PrintIslandData()
{
    mfl_add_sequence_t addseq_type = (mfl_add_sequence_t)(long int)mfl_get_parameter(m_mflHandle, MFL_PT_ADD_SEQUENCE_TYPE);
    if (addseq_type == MFL_AST_RANDOM)
    {
        int count = mfl_get_resultant_data(m_mflHandle, MFL_RT_ISLAND_COUNT,  0);
        int size;
        int len;
        cout<<left<<setw(10)<<"Island"<<setw(10)<<"Size"<<setw(10)<<"Length"<<endl;
        for (int i = 0; i < count; i++)
        {
            size = mfl_get_resultant_data(m_mflHandle, MFL_RT_ISLAND_SIZE ,  i);
            len  = mfl_get_resultant_data(m_mflHandle, MFL_RT_ISLAND_LENGTH, i);
            cout<<setw(10)<<i<<setw(10)<<size<<setw(10)<<len<<endl;
        }
    }
}

void CNexusUserInterface::PrintHsearchData()
{
    cout <<endl<<endl<<"Heuristic search completed" << endl;
    //int rearr = mfl_get_resultant_data(m_mflHandle, MFL_RT_NUM_REARRANGMENTS, 0);
    int savedtr = mfl_get_resultant_data(m_mflHandle, MFL_RT_NUM_SAVED_TREES, 0);
    int bestlen = mfl_get_resultant_data(m_mflHandle, MFL_RT_SHORTEST_TREE_LEN, 0);
    /* We don't want to report rearrangements until we can set user-defined 
     * seeds. Otherwise, the number may vary between machines.*/
    //cout << "    Number of rearrangements tried: "<<rearr<<endl;
    cout << "    Number of trees saved: "<<savedtr<<endl;
    cout << "    Length of shortest tree: "<<bestlen<<endl<<endl;
}

bool CNexusUserInterface::fCNexusMenuHeuristicSearch(string *value)
{
    bool ret = mfl_heuristic(m_mflHandle);
    PrintHsearchData();
    PrintIslandData();
    return ret;
}

bool CNexusUserInterface::fCNexusMenuExhaust        (string *value)
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuBNB            (string *value)
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuBootstrap      (string *value)
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuJackknife      (string *value)
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuSTR            (string *value)
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuConsens        (string *value)
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuCollapse       (string *value)
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::fCNexusMenuReport(string *value)
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

bool CNexusUserInterface::fCNexusMenuSearchType     (string *value)
{
    string strInput;
    mfl_search_t search_type = MFL_ST_MAX;
    m_ioNumericSubCommands->GetUserInput(" Enter search type\n  1) Exhaustive\n  2) Branch Bound\n  3) Heuristic\n # ", &strInput);
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

bool CNexusUserInterface::fCNexusMenuBranchSwapType (string *value)
{
    string strInput;
    mfl_branch_swap_t swap_type = MFL_BST_MAX;
    m_ioNumericSubCommands->GetUserInput(" Enter branch swap type\n  1) Tree bisection\n  2) Subtree pruning\n  3) Nearist neighbor\n # ", &strInput);
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

bool CNexusUserInterface::fCNexusMenuAddSeqType     (string *value)
{
    string strInput;
    mfl_add_sequence_t add_seq_type = MFL_AST_MAX;
    m_ioNumericSubCommands->GetUserInput(" Enter add seq type\n  1) Simple\n  2) Random\n  3) As is\n  4) Closest\n # ", &strInput);
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

bool CNexusUserInterface::fCNexusMenuCollapseAt           (string *value)
{
    string strInput;
    mfl_set_collapse_at_t collapse_at = MFL_SC_MAX;
    m_ioNumericSubCommands->GetUserInput(" Enter when to collapse branches\n  1) When max length is 0\n  2) When min length is 0\n  3) When two incident nodes with equal apomorphies reconstruction sets\n # ", &strInput);
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

bool CNexusUserInterface::fCNexusMenuCollapseZero         (string *value)
{
    string strInput;
    bool collapse_zero = false;
    m_ioNumericSubCommands->GetUserInput(" Select if zero length branches should be collapsed during search\n  1) Collapse\n  2) Do not collapse\n # ", &strInput);
    if (strInput.compare("1") == 0)
    {
        collapse_zero = true;
    }
    mfl_set_parameter(m_mflHandle, MFL_PT_COLLAPSE, (void*)collapse_zero);
    return true;
}

bool CNexusUserInterface::fCNexusMenuNumIterations        (string *value)
{
    string strInput;
    int num_iterations;
    stringstream ss;
    m_ioNumericSubCommands->GetUserInput(" Enter number of iterations# ", &strInput);
    ss<<strInput.c_str();
    ss>>num_iterations;
    mfl_set_parameter(m_mflHandle, MFL_PT_NUM_ITERATIONS, (void*)num_iterations);
    return true;
}

bool CNexusUserInterface::fCNexusMenuTreeLimit        (string *value)
{
    string strInput;
    int tree_limit;
    stringstream ss;
    m_ioNumericSubCommands->GetUserInput(" Enter number of trees to store# ", &strInput);
    ss<<strInput.c_str();
    ss>>tree_limit;
    mfl_set_parameter(m_mflHandle, MFL_PT_TREELIMIT, (void*)tree_limit);
    return true;
}

bool CNexusUserInterface::fCNexusMenuRatchetSearch        (string *value)
{
    string strInput;
    bool ratchet = false;
    m_ioNumericSubCommands->GetUserInput(" Select how the ratchet search parameter should be set\n  1) Enabled\n  2) Disabled\n # ", &strInput);
    if (strInput.compare("1") == 0)
    {
        ratchet = true;
    }
    mfl_set_parameter(m_mflHandle, MFL_PT_RATCHET_SEARCH, (void*)ratchet);
    return true;
}

bool CNexusUserInterface::fCNexusMenuGap                  (string *value)
{
    string strInput;
    bool gap_missing = false;
    m_ioNumericSubCommands->GetUserInput(" Select if the gap symbol ('-') will be treated as inapplicability or as missing data\n  1) Inapplicable\n  2) Missing data\n # ", &strInput);
    if (strInput.compare("1") == 0)
    {
        gap_missing = MFL_GAP_INAPPLICABLE;
    }
    else if (strInput.compare("2") == 0)
    {
        gap_missing = MFL_GAP_MISSING_DATA;
    }
    mfl_set_parameter(m_mflHandle, MFL_PT_GAP, (void*)gap_missing);
    return true;
}

bool CNexusUserInterface::fCNexusMenuMainMenu       (string *value)
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

