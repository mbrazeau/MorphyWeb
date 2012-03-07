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
    class CNexusMenu##type : public CNexusMenuBase \
    { \
    public:\
        CNexusMenu##type(const char * strCommand, const char * strHelpText) : CNexusMenuBase(strCommand, strHelpText){}\
        bool MenuFunction(CNexusUserInterface *pNexusUserInterface)\
        {\
            return pNexusUserInterface->type();\
        }\
    };

/*
 * The followind actually defines the derived class for each command
 */
NEW_COMMAND_DEFINE(Spacer         )
NEW_COMMAND_DEFINE(OpenNexusFile  )
NEW_COMMAND_DEFINE(SaveFile       )
NEW_COMMAND_DEFINE(CloseNexusFile )

NEW_COMMAND_DEFINE(Help           )
NEW_COMMAND_DEFINE(Quit           )
NEW_COMMAND_DEFINE(About          )
NEW_COMMAND_DEFINE(CommandLog     )
NEW_COMMAND_DEFINE(Status         )
NEW_COMMAND_DEFINE(Chdir          )

NEW_COMMAND_DEFINE(Exclude        )
NEW_COMMAND_DEFINE(Include        )
NEW_COMMAND_DEFINE(Outgroup       )
NEW_COMMAND_DEFINE(Ingroup        )
NEW_COMMAND_DEFINE(Char           )

NEW_COMMAND_DEFINE(HeuristicSearch)
NEW_COMMAND_DEFINE(Exhaust        )
NEW_COMMAND_DEFINE(BNB            )
NEW_COMMAND_DEFINE(Bootstrap      )
NEW_COMMAND_DEFINE(Jackknife      )
NEW_COMMAND_DEFINE(STR            )

NEW_COMMAND_DEFINE(Consens        )
NEW_COMMAND_DEFINE(Collapse       )
NEW_COMMAND_DEFINE(Report         )
/*
 * The UI constructor puts the menu together, it stores each
 * menu option in a stl vector.
 */
CNexusUserInterface::CNexusUserInterface()
{
    m_pNexusParse = NULL;
    m_strCwd = "./";

    m_vMenu.push_back(new CNexusMenuSpacer      (NULL, "File"));
    m_vMenu.push_back(new CNexusMenuOpenNexusFile   ("O", "Open a nexus file"));
    m_vMenu.push_back(new CNexusMenuCloseNexusFile  ("C", "Close a nexus file"));
    m_vMenu.push_back(new CNexusMenuSaveFile        ("S", "Save according to the options"));

    m_vMenu.push_back(new CNexusMenuSpacer      (NULL, "Program"));
    m_vMenu.push_back(new CNexusMenuCommandLog      ("LOG" ,"Toggles record of commands, variable states, etc"));
    m_vMenu.push_back(new CNexusMenuStatus          ("STAT","Prints status of all current settings, eg. logmode on/off"));
    m_vMenu.push_back(new CNexusMenuChdir           ("CD"  ,"change working directory"));
    m_vMenu.push_back(new CNexusMenuHelp            ("H"   , "Help"));
    m_vMenu.push_back(new CNexusMenuAbout           ("A"   , "About"));
    m_vMenu.push_back(new CNexusMenuQuit            ("Q"   , "Quit"));

    m_vMenu.push_back(new CNexusMenuSpacer      (NULL, "Data"));
    m_vMenu.push_back(new CNexusMenuExclude         ("EXC" , "Exclude taxa or characters"));
    m_vMenu.push_back(new CNexusMenuInclude         ("INC" , "Include excluded taxa or characters"));
    m_vMenu.push_back(new CNexusMenuOutgroup        ("OUTG", "Assign taxa to outgroup"));
    m_vMenu.push_back(new CNexusMenuIngroup         ("ING" , "Return taxa from outgroup to ingrou"));
    m_vMenu.push_back(new CNexusMenuChar            ("CHAR", "Modify a character's type"));

    m_vMenu.push_back(new CNexusMenuSpacer      (NULL, "Analysis"));
    m_vMenu.push_back(new CNexusMenuHeuristicSearch ("HS" , "Begin a heuristic search"));
    m_vMenu.push_back(new CNexusMenuExhaust         ("EXS", "Begin an exhaustive search"));
    m_vMenu.push_back(new CNexusMenuBNB             ("BNB", "Begin a branch-and-bound search"));
    m_vMenu.push_back(new CNexusMenuBootstrap       ("BTS", "Begin a bootstrap analysis"));
    m_vMenu.push_back(new CNexusMenuJackknife       ("JK" , "Begin a jackknife analysis"));
    m_vMenu.push_back(new CNexusMenuSTR             ("STR", "Perform a safe taxonomic reduction"));

    m_vMenu.push_back(new CNexusMenuSpacer      (NULL, "Results"));
    m_vMenu.push_back(new CNexusMenuConsens         ("CONSENS" , "Compute consensus tree for trees in memory"));
    m_vMenu.push_back(new CNexusMenuCollapse        ("COLLAPSE", "Collapse zero-length branches, condense the tree set"));
    m_vMenu.push_back(new CNexusMenuReport          ("REPORT"  , "Print a report about the current open nexus file"));
}

/*
 * The UI destructor deletes all memory.
 */
CNexusUserInterface::~CNexusUserInterface()
{
    vector<CNexusMenuBase*>::iterator it;
    CNexusMenuBase* pMenu;

    for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
    {
        pMenu = *it;
        if (pMenu)
        {
            delete(pMenu);
        }
    }
    m_vMenu.clear();
    CloseNexusFile(false);
    if (m_fCommandLog)
    {
        m_fCommandLog.close();
    }
}

/*
 * Anytime the user is to give input, the input must come through this function.
 */
void CNexusUserInterface::GetUserInput(string strPrompt, string *strInput)
{
    cout<<strPrompt;
    cin>>*strInput;
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
    About();
    Help();
    do
    {
        GetUserInput("Enter selection# " ,&strInput);
    } while (RunSelection(strInput));
}

/*
 * Loop through the menu vector and look for the command entered by the user.
 * If the command is found then run the command by calling the overloaded 
 * functionoid MenuFunction() with the current CNexusUserInterface instance
 * as the parameter.
 */
bool CNexusUserInterface::RunSelection(string strInput)
{
    vector<CNexusMenuBase*>::iterator it;
    CNexusMenuBase* pMenu;

    for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
    {
        pMenu = *it;
        if ((pMenu) && (pMenu->IsSelection(strInput)))
        {
            bool bRet;
            cout<<endl;
            bRet = pMenu->MenuFunction(this);
            cout<<endl;
            return bRet;
        }
    }
    cout<<" Unknown command: "<<strInput<<endl;
    return true;
}

/*
 * Run the open file command, just prompt the user for input
 * and attempt to read the nexus file
 */
bool CNexusUserInterface::OpenNexusFile()
{
    string strFilename;

    if (!m_pNexusParse)
    {
        GetUserInput(" Enter filename: " + m_strCwd, &strFilename);
        strFilename = m_strCwd + strFilename;
        m_pNexusParse = new CNexusParse();
        if (m_pNexusParse)
        {
            if (m_pNexusParse->ReadNexusFile(&strFilename, NULL))
            {
                cout<<" "<<strFilename<<" open successfully"<<endl;
            }
            else
            {
                CloseNexusFile(false);
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

bool CNexusUserInterface::SaveFile       ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

/*
 * Close a file opened with the OpenNexusFile command
 * bVerbose is defaulted to = true
 */
bool CNexusUserInterface::CloseNexusFile(bool bVerbose)
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
    return true;
}

/*
 * Print the menu
 */
bool CNexusUserInterface::Help           ()
{
    vector<CNexusMenuBase*>::iterator it;
    CNexusMenuBase* pMenu;

    for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
    {
        pMenu = *it;
        if (pMenu)
        {
            cout<<pMenu->GetMenuOutput()<<endl;
        }
    }
    cout<<endl;
    return true;
}

bool CNexusUserInterface::Quit           ()
{
    return false;
}

bool CNexusUserInterface::About          ()
{
    cout<<"Morphy NUI Version: "<<NUI_MAJOR_VERSION<<"."<<NUI_MINOR_VERSION<<endl;
    cout<<"Copyright 2012 (C) Martin Brazeau and Chris Desjardins. All rights reserved."<<endl;
    cout<<"This program uses the NCL by Paul O. Lewis."<<endl;
    cout<<"Build time: "<<__DATE__<<" "<<__TIME__<<endl;
    return true;
}

bool CNexusUserInterface::CommandLog     ()
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

bool CNexusUserInterface::Status         ()
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

bool CNexusUserInterface::Chdir          ()
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

bool CNexusUserInterface::Exclude        ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::Include        ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::Outgroup       ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::Ingroup        ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::Char           ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::HeuristicSearch()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::Exhaust        ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::BNB            ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::Bootstrap      ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::Jackknife      ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::STR            ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::Consens        ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::Collapse       ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::Report()
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

int main(int argc, char *argv[])
{
    CNexusUserInterface ui;
    ui.DoMenu();
    
    return 0;
}

