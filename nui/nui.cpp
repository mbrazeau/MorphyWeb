#include "nui.h"

bool CNexusMenuOpenFile::MenuFunction(CNexusUserInterface *pNexusUserInterface)
{
    cout<<"OPEN FILE STILL NEEDS TO BE IMPLEMENTED"<<endl;
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

CNexusUserInterface::CNexusUserInterface()
{
    m_vMenu.push_back(new CNexusMenuOpenFile  ("O", "Open a file"));
    m_vMenu.push_back(new CNexusMenuCloseFile ("C", "Close a file"));
    m_vMenu.push_back(new CNexusMenuHelp      ("H", "Help"));
    m_vMenu.push_back(new CNexusMenuQuit      ("Q", "Quit"));
}

CNexusUserInterface::~CNexusUserInterface()
{
    vector<CNexusMenuBase*>::iterator it;
    CNexusMenuBase* pMenu;

    for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
    {
        pMenu = *it;
        delete(pMenu);
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
        cout<<pMenu->GetMenuOutput()<<endl;
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
        if (pMenu->IsSelection(strInput))
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

