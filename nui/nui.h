#pragma once

#include <vector>
#include "NexusParse.h"

#define NUI_MAJOR_VERSION   0
#define NUI_MINOR_VERSION   1

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

class CNexusUserInterface
{
public:
    typedef  int (CNexusUserInterface::*TMenuAction)();
    CNexusUserInterface();
    ~CNexusUserInterface();

    void DoMenu();
    bool RunSelection(string strInput);

    bool Spacer() 
    {
        return false;
    }
    bool OpenFile       ();
    bool SaveFile       ();
    bool CloseFile      (bool bVerbose = true);
               
    bool Help           ();
    bool Quit           ();
    bool About          ();
    bool Log            ();
    bool Status         ();
    bool Chdir          ();
               
    bool Exclude        ();
    bool Include        ();
    bool Outgroup       ();
    bool Ingroup        ();
    bool Char           ();
               
    bool HeuristicSearch();
    bool Exhaust        ();
    bool BNB            ();
    bool Bootstrap      ();
    bool Jackknife      ();
    bool STR            ();
               
    bool Consens        ();
    bool Collapse       ();
    bool Report         ();


private:
    vector <CNexusMenuBase*> m_vMenu;

    CNexusParse *m_pNexusParse;
};


