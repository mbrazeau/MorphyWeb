#include "nui.h"

bool CNexusMenuSpacer::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    return true;
}

bool CNexusMenuOpenFile::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"OPEN FILE STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuSaveFile::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"SAVE FILE STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuCloseFile::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"CLOSE FILE STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuHelp::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    pNexusUserInterface->PrintMenu();
    return true;
}

bool CNexusMenuQuit::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"Thanks for visiting..."<<endl;
    return false;
}

bool CNexusMenuAbout::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"ABOUT STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuLog::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"LOG STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuStatus::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"STATUS STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuChdir::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"CHDIR STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuExclude::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"EXCLUDE STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuInclude::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"INCLUDE STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuOutgroup::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"OUTGROUP STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuIngroup::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"INGROUP STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuChar::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"CHAR STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

                         
bool CNexusMenuHeuristicSearch::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"HEURISTIC SEARCH STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuExhaust::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"EXHAUST STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuBNB::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"BNB STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuBootstrap::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"BOOTSTRAP STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuJackknife::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"JACKKNIFE STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuSTR::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"STR STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

                         
bool CNexusMenuConsens::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"CONSENS STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}

bool CNexusMenuCollapse::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"COLLAPSE STILL NEEDS TO BE IMPLEMENTED"<<endl;
    return true;
}


CNexusUserInterface::CNexusUserInterface()
{
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
}

void CNexusUserInterface::PrintMenu()
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
}

void CNexusUserInterface::DoMenu()
{
    string strInput;
    PrintMenu();
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

int main(int argc, char *argv[])
{
    CNexusUserInterface ui;
    ui.DoMenu();
    /*
    CNexusParse cNexusParse(argv[1], argv[2]);
    cNexusParse.ReadNexusFile();
    cNexusParse.Report();
    */
    return 0;
}

