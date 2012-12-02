#pragma once

#include "histedit.h"
#include "myostream.h"

const static char *prompt(EditLine *el)
{
    return " ";
}

class CEditLineHist
{
public:
    typedef char    *(*el_pfunc_t)(EditLine *);
    CEditLineHist(const char *pgmName, myofstream *fCommandLog)
    {
        m_pEditLine = el_init(pgmName, stdin, stdout, stderr);
        if (!m_pEditLine)
        {
            throw "Unable to initialize line editor";
        }
        el_set(m_pEditLine, EL_PROMPT, &prompt);
        el_set(m_pEditLine, EL_EDITOR, "emacs");

        m_pHistory = history_init();
        if (!m_pHistory)
        {
            throw "Unable to initialize history";
        }
        history(m_pHistory, &m_histEvent, H_SETSIZE, 800);
        el_set(m_pEditLine, EL_HIST, history, m_pHistory);
        m_pfCommandLog = fCommandLog;
    }

    ~CEditLineHist()
    {
        if (m_pHistory)
        {
            history_end(m_pHistory);
        }
        if (m_pEditLine)
        {
            el_end(m_pEditLine);
        }
    }

    void GetUserInput(string strPrompt, string *strInput)
    {
        const char *pLine;
        int nCount;
        cout<<strPrompt;
        pLine = el_gets(m_pEditLine, &nCount);
        if (nCount > 0)
        {
            history(m_pHistory, &m_histEvent, H_ENTER, pLine);
            *strInput = pLine;
            strInput->erase(std::remove(strInput->begin(), strInput->end(), '\n'), strInput->end());
            if (*m_pfCommandLog)
            {
                (*m_pfCommandLog)<<*strInput<<endl;
            }            
        }
        else
        {
            *strInput = "q";
            throw "Stardard input failure";            
        }
    }
protected:

    EditLine *m_pEditLine;
    History *m_pHistory;
    HistEvent m_histEvent;
    static string m_strPrompt;
    myofstream *m_pfCommandLog;

private:

};
