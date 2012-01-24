#include "nui.h"

CNexusUserInterface::CNexusUserInterface()
{
    m_vMenu.push_back(new CNexusMenuOpenFile("O", "Open a file"));
    m_vMenu.push_back(new CNexusMenuCloseFile("C", "Close a file"));
}

CNexusUserInterface::~CNexusUserInterface()
{
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

void CNexusUserInterface::DoMenu(bool bPrintMenu)
{
    string strInput;
    if (bPrintMenu)
    {
        PrintMenu();
    }
    cout<<"Enter selection# ";
    cin>>strInput;
    RunSelection(strInput);
}

void CNexusUserInterface::RunSelection(string strInput)
{
    vector<CNexusMenuBase*>::iterator it;
    CNexusMenuBase* pMenu;

    for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
    {
        pMenu = *it;
        if (pMenu->IsSelection(strInput))
        {
            pMenu->MenuFunction();
        }
    }
}

int main(int argc, char *argv[])
{
    CNexusUserInterface ui;
    ui.DoMenu(true);
    /*
    CNexusParse cNexusParse(argv[1], argv[2]);
    cNexusParse.ReadNexusFile();
    cNexusParse.Report();
    */
    return 0;
}

