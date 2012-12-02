#pragma once

#include <vector>
#include "NexusParse.h"
#include "mfl.h"
#include "EditLineHist.h"
#include "NexusMenuData.h"

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

    bool fCNexusMenuSpacer(string *value = NULL, int nMappedVal = -1) 
    {
        throw "Cannot invoke the spacer";
        return true;
    }
    bool fCNexusMenuOpenNexusFile  (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuSaveFile       (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuCloseNexusFile (string *value = NULL, int nMappedVal = -1, bool bVerbose = true);
               
    bool fCNexusMenuHelp           (string *value = NULL, int nMappedVal = -1, bool bForceShowMenu = true);
    bool fCNexusMenuQuit           (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuAbout          (string *value = NULL, int nMappedVal = -1, bool bShowBuildTime = true);
    bool fCNexusMenuCommandLog     (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuSave           (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuStatus         (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuChdir          (string *value = NULL, int nMappedVal = -1);
               
    bool fCNexusMenuExclude        (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuInclude        (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuOutgroup       (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuIngroup        (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuChar           (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuSet            (string *value = NULL, int nMappedVal = -1);
               
    bool fCNexusMenuHeuristicSearch(string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuExhaust        (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuBNB            (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuBootstrap      (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuJackknife      (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuSTR            (string *value = NULL, int nMappedVal = -1);

    bool fCNexusMenuConsens        (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuCollapse       (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuReport         (string *value = NULL, int nMappedVal = -1);

    bool fCNexusMenuSearchType     (string *value = NULL, mfl_search_t nMappedVal = MFL_ST_MAX);
    bool fCNexusMenuBranchSwapType (string *value = NULL, mfl_branch_swap_t nMappedVal = MFL_BST_MAX);
    bool fCNexusMenuAddSeqType     (string *value = NULL, mfl_add_sequence_t nMappedVal = MFL_AST_MAX);
    bool fCNexusMenuCollapseAt     (string *value = NULL, mfl_set_collapse_at_t nMappedVal = MFL_SC_MAX);
    bool fCNexusMenuCollapseZero   (string *value = NULL, bool nMappedVal = false);
    bool fCNexusMenuNumIterations  (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuTreeLimit      (string *value = NULL, int nMappedVal = -1);
    bool fCNexusMenuRatchetSearch  (string *value = NULL, bool nMappedVal = false);
    bool fCNexusMenuGap            (string *value = NULL, mfl_gap_t nMappedVal = MFL_GAP_MAX);

    bool fCNexusMenuMainMenu       (string *value = NULL, int nMappedVal = -1);

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

    void ConfigMenuSearchType();
    void ConfigMenuBranchSwapType();
    void ConfigMenuAddSeqType();
    void ConfigMenuCollapseAt();
    void ConfigMenuCollapseZero();
    void ConfigMenuRatchetSearch();
    void ConfigMenuGap();

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


