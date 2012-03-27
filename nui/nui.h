#pragma once

#include <vector>
#include "NexusParse.h"
#include "mfl.h"

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
        if (strInput.length() > 0)
        {
            transform(strInput.begin(), strInput.end(), strInput.begin(), ::toupper);        
            return strInput == m_strCommand;
        }
        return false;
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
        throw "Cannot invoke the spacer";
        return true;
    }
    bool OpenNexusFile  ();
    bool SaveFile       ();
    bool CloseNexusFile (bool bVerbose = true);
               
    bool Help           ();
    bool Quit           ();
    bool About          (bool bShowBuildTime = true);
    bool CommandLog     ();
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

protected:
    void DestroyHandle();
    void CreateHandle();
    bool SetMorphyOpenParams();

private:
    void GetUserInput(string strPrompt, string *strInput);
    vector <CNexusMenuBase*> m_vMenu;
    CNexusParse *m_pNexusParse;
    string m_strCwd;

    myofstream m_fCommandLog;
    mfl_handle_t m_mflHandle;
};


