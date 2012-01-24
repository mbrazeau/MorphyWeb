#pragma once

#include <vector>
#include "NexusParse.h"

class CNexusUserInterface;

class CNexusMenuBase
{
public:
    CNexusMenuBase(const char * strCommand, const char * strHelpText)
    {
        m_strCommand = strCommand;
        m_strHelpText = strHelpText;
        transform(m_strCommand.begin(), m_strCommand.end(), m_strCommand.begin(), ::toupper);
    }

    string GetMenuOutput()
    {
        return m_strCommand + ") " + m_strHelpText;
    }

    bool IsSelection(string strInput)
    {
        transform(strInput.begin(), strInput.end(), strInput.begin(), ::toupper);        
        return strInput == m_strCommand;
    }

    virtual bool MenuFunction(CNexusUserInterface *pNexusUserInterface) = 0;
private:
    string m_strCommand;
    string m_strHelpText;
};

#define NEW_COMMAND_DEFINE(type) \
    class type : public CNexusMenuBase \
    { \
    public:\
        type(const char * strCommand, const char * strHelpText) : CNexusMenuBase(strCommand, strHelpText){}\
        bool MenuFunction(CNexusUserInterface *pNexusUserInterface);\
    };

NEW_COMMAND_DEFINE(CNexusMenuOpenFile)
NEW_COMMAND_DEFINE(CNexusMenuCloseFile)
NEW_COMMAND_DEFINE(CNexusMenuHelp)
NEW_COMMAND_DEFINE(CNexusMenuQuit)

class CNexusUserInterface
{
public:
    typedef  int (CNexusUserInterface::*TMenuAction)();
    CNexusUserInterface();
    ~CNexusUserInterface();

    void PrintMenu();
    void DoMenu();
    bool RunSelection(string strInput);

private:
    vector <CNexusMenuBase*> m_vMenu;
};


