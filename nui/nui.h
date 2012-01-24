#pragma once

#include <vector>
#include "NexusParse.h"


class CNexusMenuBase
{
public:
    CNexusMenuBase(const char * strCommand, const char * strHelpText)
    {
        m_strCommand = strCommand;
        m_strHelpText = strHelpText;
    };
    virtual int MenuFunction() = 0;
    string GetMenuOutput()
    {
        return m_strCommand + ") " + m_strHelpText;
    }
    bool IsSelection(string strInput)
    {
        return strInput == m_strCommand;
    }
private:
    string m_strCommand;
    string m_strHelpText;
};

class CNexusMenuOpenFile : public CNexusMenuBase
{
public:
    CNexusMenuOpenFile(const char * strCommand, const char * strHelpText) : CNexusMenuBase(strCommand, strHelpText){};
    int MenuFunction()
    {
        cout<<"OPEN FILE"<<endl;
        return 0;
    };
};

class CNexusMenuCloseFile : public CNexusMenuBase
{
public:
    CNexusMenuCloseFile(const char * strCommand, const char * strHelpText) : CNexusMenuBase(strCommand, strHelpText){};
    int MenuFunction()
    {
        cout<<"CLOSE FILE"<<endl;
        return 0;
    };
};

class CNexusUserInterface
{
public:
    typedef  int (CNexusUserInterface::*TMenuAction)();
    CNexusUserInterface();
    ~CNexusUserInterface();

    void PrintMenu();
    void DoMenu(bool bPrintMenu);
    void RunSelection(string strInput);

private:
    vector <CNexusMenuBase*> m_vMenu;
};


