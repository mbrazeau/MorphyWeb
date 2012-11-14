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

    bool fCNexusMenuSpacer(string *value = NULL) 
    {
        throw "Cannot invoke the spacer";
        return true;
    }
    bool fCNexusMenuOpenNexusFile  (string *value = NULL);
    bool fCNexusMenuSaveFile       (string *value = NULL);
    bool fCNexusMenuCloseNexusFile (string *value = NULL, bool bVerbose = true);
               
    bool fCNexusMenuHelp           (string *value = NULL, bool bForceShowMenu = true);
    bool fCNexusMenuQuit           (string *value = NULL);
    bool fCNexusMenuAbout          (string *value = NULL, bool bShowBuildTime = true);
    bool fCNexusMenuCommandLog     (string *value = NULL);
    bool fCNexusMenuSave           (string *value = NULL);
    bool fCNexusMenuStatus         (string *value = NULL);
    bool fCNexusMenuChdir          (string *value = NULL);
               
    bool fCNexusMenuExclude        (string *value = NULL);
    bool fCNexusMenuInclude        (string *value = NULL);
    bool fCNexusMenuOutgroup       (string *value = NULL);
    bool fCNexusMenuIngroup        (string *value = NULL);
    bool fCNexusMenuChar           (string *value = NULL);
    bool fCNexusMenuSet            (string *value = NULL);
               
    bool fCNexusMenuHeuristicSearch(string *value = NULL);
    bool fCNexusMenuExhaust        (string *value = NULL);
    bool fCNexusMenuBNB            (string *value = NULL);
    bool fCNexusMenuBootstrap      (string *value = NULL);
    bool fCNexusMenuJackknife      (string *value = NULL);
    bool fCNexusMenuSTR            (string *value = NULL);
               
    bool fCNexusMenuConsens        (string *value = NULL);
    bool fCNexusMenuCollapse       (string *value = NULL);
    bool fCNexusMenuReport         (string *value = NULL);

    bool fCNexusMenuSearchType     (string *value = NULL);
    bool fCNexusMenuBranchSwapType (string *value = NULL);
    bool fCNexusMenuAddSeqType     (string *value = NULL);
    bool fCNexusMenuCollapseAt     (string *value = NULL);
    bool fCNexusMenuCollapseZero   (string *value = NULL);
    bool fCNexusMenuNumIterations  (string *value = NULL);
    bool fCNexusMenuTreeLimit      (string *value = NULL);
    bool fCNexusMenuRatchetSearch  (string *value = NULL);
    bool fCNexusMenuGap            (string *value = NULL);

    bool fCNexusMenuMainMenu       (string *value = NULL);

protected:
    void DestroyHandle();
    void CreateHandle();
    bool SetMorphyOpenParams();
    void ChangeMenu(CNexusMenuData *pMenu);
    void Delete(CEditLineHist *pMem);
    bool SaveTranslateTable(myofstream &fSave);
    bool SaveNewickStrings(myofstream &fSave);
    void PrintIslandData();
    void PrintHsearchData();

private:
    CNexusMenuData *m_pMainMenu;
    CNexusMenuData *m_pMenu;

    CNexusParse *m_pNexusParse;
    string m_strCwd;

    myofstream m_fCommandLog;
    mfl_handle_t m_mflHandle;

    CEditLineHist *m_ioWorkingDir;
    CEditLineHist *m_ioLogFiles;
    CEditLineHist *m_ioInputFiles;
    CEditLineHist *m_ioSaveFiles;
    CEditLineHist *m_ioCommands;
    CEditLineHist *m_ioNumericSubCommands;
};


