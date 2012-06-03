#pragma once

#include <vector>
#include "NexusParse.h"
#include "mfl.h"
#include "menu.h"
#include "EditLineHist.h"

#define NUI_MAJOR_VERSION   0
#define NUI_MINOR_VERSION   1

    
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
               
    bool fCNexusMenuHelp           (bool bForceShowMenu = true);
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
    bool fCNexusMenuSet            ();
               
    bool fCNexusMenuHeuristicSearch();
    bool fCNexusMenuExhaust        ();
    bool fCNexusMenuBNB            ();
    bool fCNexusMenuBootstrap      ();
    bool fCNexusMenuJackknife      ();
    bool fCNexusMenuSTR            ();
               
    bool fCNexusMenuConsens        ();
    bool fCNexusMenuCollapse       ();
    bool fCNexusMenuReport         ();

    bool fCNexusMenuSearchType     ();
    bool fCNexusMenuBranchSwapType ();
    bool fCNexusMenuAddSeqType     ();
    bool fCNexusMenuCollapseAt     ();
    bool fCNexusMenuCollapseZero   ();
    bool fCNexusMenuNumIterations  ();
    bool fCNexusMenuTreeLimit      ();
    bool fCNexusMenuRatchetSearch  ();
    bool fCNexusMenuGap            ();

    bool fCNexusMenuMainMenu       ();

protected:
    void DestroyHandle();
    void CreateHandle();
    bool SetMorphyOpenParams();
    void ChangeMenu(CNexusMenuData *pMenu);
    void Delete(CEditLineHist *pMem);

private:
    CNexusMenuData *m_pMainMenu;
    CNexusMenuData *m_pSetMenu;
    CNexusMenuData *m_pMenu;

    CNexusParse *m_pNexusParse;
    string m_strCwd;

    myofstream m_fCommandLog;
    mfl_handle_t m_mflHandle;

    CEditLineHist *m_ioWorkingDir;
    CEditLineHist *m_ioLogs;
    CEditLineHist *m_ioFiles;
    CEditLineHist *m_ioCommands;
    CEditLineHist *m_ioNumericSubCommands;
};


