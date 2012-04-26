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

    bool fCNexusMenuSpacer() 
    {
        throw "Cannot invoke the spacer";
        return true;
    }
    bool fCNexusMenuOpenNexusFile  ();
    bool fCNexusMenuSaveFile       ();
    bool fCNexusMenuCloseNexusFile (bool bVerbose = true);
               
    bool fCNexusMenuHelp           ();
    bool fCNexusMenuQuit           ();
    bool fCNexusMenuAbout          (bool bShowBuildTime = true);
    bool fCNexusMenuCommandLog     ();
    bool fCNexusMenuStatus         ();
    bool fCNexusMenuChdir          ();
               
    bool fCNexusMenuExclude        ();
    bool fCNexusMenuInclude        ();
    bool fCNexusMenuOutgroup       ();
    bool fCNexusMenuIngroup        ();
    bool fCNexusMenuChar           ();
    bool fCNexusMenuSearchType     ();
               
    bool fCNexusMenuHeuristicSearch();
    bool fCNexusMenuExhaust        ();
    bool fCNexusMenuBNB            ();
    bool fCNexusMenuBootstrap      ();
    bool fCNexusMenuJackknife      ();
    bool fCNexusMenuSTR            ();
               
    bool fCNexusMenuConsens        ();
    bool fCNexusMenuCollapse       ();
    bool fCNexusMenuReport         ();

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


