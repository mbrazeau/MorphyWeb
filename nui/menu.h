#pragma once

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
    virtual ~CNexusMenuBase()
    {
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
protected:
    string m_strCommand;
    string m_strHelpText;
};

class CNexusMenuData
{
public:
    CNexusMenuData(string strMenuTitle)
    {
        bMenuShown = false;
        m_strMenuTitle = strMenuTitle;
    }
    ~CNexusMenuData()
    {
        vector<CNexusMenuBase*>::iterator it;
        CNexusMenuBase* pMenuItem;
       
        for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
        {
            pMenuItem = *it;
            if (pMenuItem)
            {
                delete(pMenuItem);
            }
        }
        m_vMenu.clear();
    }

    void AddMenuItem(CNexusMenuBase* pMenuItem)
    {
        m_vMenu.push_back(pMenuItem);
    }
    
    bool RunSelection(string strInput, CNexusUserInterface *pNexusUserInterface)
    {
        vector<CNexusMenuBase*>::iterator it;
        CNexusMenuBase* pMenuItem;
        
        for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
        {
            pMenuItem = *it;
            if ((pMenuItem) && (pMenuItem->IsSelection(strInput)))
            {
                bool bRet = true;
                cout<<endl;
                try
                {
                    bRet = pMenuItem->MenuFunction(pNexusUserInterface);
                }
                catch (const char *e)
                {
                    cout<<"NUI Error: "<<e<<endl;
                }
                catch (mfl_exception e)
                {
                    cout<<"MFL Error: "<<e.what()<<endl;
                }
                cout<<endl;
                return bRet;
            }
        }
        cout<<" Unknown command: "<<strInput<<endl;
        return true;
    }

    bool Help           (bool bForceShowMenu)
    {
        if ((bMenuShown == false) || (bForceShowMenu == true))
        {
            vector<CNexusMenuBase*>::iterator it;
            CNexusMenuBase* pMenu;
           
            for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
            {
                pMenu = *it;
                if (pMenu)
                {
                    cout<<pMenu->GetMenuOutput()<<endl;
                }
            }
            cout<<endl;
            bMenuShown = true;
        }
        return true;
    }

    string GetPrompt()
    {
        return m_strMenuTitle;
    }

protected:
    bool bMenuShown;
    string m_strMenuTitle;
    vector <CNexusMenuBase*> m_vMenu;
};
