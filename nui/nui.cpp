#include "nui.h"

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

NEW_COMMAND_DEFINE(Spacer         )
NEW_COMMAND_DEFINE(OpenFile       )
NEW_COMMAND_DEFINE(SaveFile       )
NEW_COMMAND_DEFINE(CloseFile      )

NEW_COMMAND_DEFINE(Help           )
NEW_COMMAND_DEFINE(Quit           )
NEW_COMMAND_DEFINE(About          )
NEW_COMMAND_DEFINE(Log            )
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

CNexusUserInterface::CNexusUserInterface()
{
    m_pNexusParse = NULL;

    m_vMenu.push_back(new CNexusMenuSpacer      (NULL, "File"));
    m_vMenu.push_back(new CNexusMenuOpenFile        ("O", "Open a file"));
    m_vMenu.push_back(new CNexusMenuCloseFile       ("C", "Close a file"));
    m_vMenu.push_back(new CNexusMenuSaveFile        ("S", "Save according to the options"));

    m_vMenu.push_back(new CNexusMenuSpacer      (NULL, "Program"));
    m_vMenu.push_back(new CNexusMenuLog             ("LOG" ,"Toggles record of commands, variable states, etc"));
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
    m_vMenu.push_back(new CNexusMenuReport          ("REPORT", "Print a report about the current open file"));
}

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
    if (m_pNexusParse)
    {
        delete m_pNexusParse;
        m_pNexusParse = NULL;
    }    
}

void CNexusUserInterface::DoMenu()
{
    string strInput;
    Help();
    do
    {
        cout<<"Enter selection# ";
        cin>>strInput;
    } while (RunSelection(strInput));
}

bool CNexusUserInterface::RunSelection(string strInput)
{
    vector<CNexusMenuBase*>::iterator it;
    CNexusMenuBase* pMenu;

    for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
    {
        pMenu = *it;
        if ((pMenu) && (pMenu->IsSelection(strInput)))
        {
            return pMenu->MenuFunction(this);
        }
    }
    cout<<"Unknown command: "<<strInput<<endl;
    return true;
}

bool CNexusUserInterface::OpenFile()
{
    string strFilename;

    cout<<" Enter filename: ";
    cin>>strFilename;
    m_pNexusParse = new CNexusParse(&strFilename, NULL);
    if (m_pNexusParse)
    {
        if (m_pNexusParse->ReadNexusFile())
        {
            cout<<endl<<" "<<strFilename<<" open successfully"<<endl<<endl;
        }
        else
        {
            delete m_pNexusParse;
            m_pNexusParse = NULL;
            cout<<"Error: Unable to read "<<strFilename<<endl;
        }
    }
    else
    {
        cout<<"Error: Unable to open "<<strFilename<<endl;
    }
    return true;
}

bool CNexusUserInterface::SaveFile       ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::CloseFile      ()
{
    if (m_pNexusParse)
    {
        delete m_pNexusParse;
        m_pNexusParse = NULL;
        cout<<endl<<"Successfully closed open file..."<<endl<<endl;
    }
    else
    {
        cout<<endl<<"No file is currently open"<<endl<<endl;
    }
    return true;
}

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
    return true;
}

bool CNexusUserInterface::Quit           ()
{
    return false;
}

bool CNexusUserInterface::About          ()
{
    cout<<endl<<"Copyright 2012 (C) Martin Brazeau, and Chris Desjardins. All rights reserved."<<endl;
    cout<<"This program uses the NCL by Paul O. Lewis."<<endl<<endl;
    return true;
}

bool CNexusUserInterface::Log            ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::Status         ()
{
    cout<<"Not implemented"<<endl;
    return true;
}

bool CNexusUserInterface::Chdir          ()
{
    cout<<"Not implemented"<<endl;
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
        cout<<endl<<"No file is currently open"<<endl<<endl;
    }
    return true;
}

int main(int argc, char *argv[])
{
    CNexusUserInterface ui;
    ui.DoMenu();
    
    return 0;
}

