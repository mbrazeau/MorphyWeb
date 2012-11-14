#pragma once

class CNexusUserInterface;

class CNexusMenuBase
{
public:
    CNexusMenuBase(const char * strCommand, const char * strHelpText)
    {
        if (strCommand)
        {
            m_strCommandNormal = strCommand;
            m_strCommand       = strCommand;
        }
        m_strHelpText = strHelpText;
        transform(m_strCommand.begin(), m_strCommand.end(), m_strCommand.begin(), ::tolower);
    }
    virtual ~CNexusMenuBase()
    {
    }

    string GetMenuOutput()
    {
        if (m_strCommand.length() > 0)
        {
            return " " + m_strCommandNormal + ") " + m_strHelpText;
        }
        else
        {
            return "\n=== " + m_strHelpText;
        }
    }

    bool IsSelection(string strInput)
    {
        int index;
        if (strInput.length() > 0)
        {
            transform(strInput.begin(), strInput.end(), strInput.begin(), ::tolower);
            index = m_strCommand.find(strInput);
            return (index == 0);
        }
        return false;
    }

    virtual bool MenuFunction(CNexusUserInterface *pNexusUserInterface) = 0;
protected:
    string m_strCommandNormal;
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
        vector<CNexusMenuBase*> pMenuItems;
        
        pMenuItems = GetMenuSelection(strInput);
        if (pMenuItems.size() == 1)
        {
            bool bRet = true;
            cout<<endl;
            try
            {
                bRet = pMenuItems[0]->MenuFunction(pNexusUserInterface);
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
        PrintError(strInput, pMenuItems);

        return true;
    }

    bool Help(bool bForceShowMenu)
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

    void PrintError(string strInput, vector<CNexusMenuBase*> pMenuItems)
    {
        vector<CNexusMenuBase*>::iterator it;
        cout<<" Unknown command: '"<<strInput<<"'"<<endl<<endl;
        CNexusMenuBase* pMenuItem;
        for (it = pMenuItems.begin(); it < pMenuItems.end(); it++)
        {
            pMenuItem = *it;
            cout<<"  "<<pMenuItem->GetMenuOutput()<<endl;
        }
        cout<<endl;
    }

    vector<CNexusMenuBase*> GetMenuSelection(string strInput)
    {
        vector<CNexusMenuBase*>::iterator it;
        CNexusMenuBase* pMenuItem;
        vector<CNexusMenuBase*> ret;
        for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
        {
            pMenuItem = *it;
            if ((pMenuItem) && (pMenuItem->IsSelection(strInput)))
            {
                ret.push_back(pMenuItem);
            }
        }
        return ret;
    }

    bool bMenuShown;
    string m_strMenuTitle;
    vector <CNexusMenuBase*> m_vMenu;
};
