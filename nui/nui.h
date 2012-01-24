#pragma once

#include <vector>
#include "NexusParse.h"

class CNexusUserInterface;

class CNexusMenuBase
{
public:
    CNexusMenuBase(const char * strCommand, const char * strHelpText)
    {
        if (strCommand)
        {
            m_strCommand = strCommand;
        }
        m_strHelpText = strHelpText;
        transform(m_strCommand.begin(), m_strCommand.end(), m_strCommand.begin(), ::toupper);
    }

    string GetMenuOutput()
    {
        if (m_strCommand.length() > 0)
        {
            return " " + m_strCommand + ") " + m_strHelpText;
        }
        else
        {
            return "\n=== " + m_strHelpText;
        }
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

NEW_COMMAND_DEFINE(CNexusMenuSpacer         )
                                            
NEW_COMMAND_DEFINE(CNexusMenuOpenFile       )
NEW_COMMAND_DEFINE(CNexusMenuSaveFile       )
NEW_COMMAND_DEFINE(CNexusMenuCloseFile      )
                                            
NEW_COMMAND_DEFINE(CNexusMenuHelp           )
NEW_COMMAND_DEFINE(CNexusMenuQuit           )
NEW_COMMAND_DEFINE(CNexusMenuAbout          )
NEW_COMMAND_DEFINE(CNexusMenuLog            )
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
NEW_COMMAND_DEFINE(CNexusMenuCollapse       )

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


