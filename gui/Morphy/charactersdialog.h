#ifndef CHARACTERSDIALOG_H
#define CHARACTERSDIALOG_H

#include <QDialog>

namespace Ui {
class charactersDialog;
}

class charactersDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit charactersDialog(QWidget *parent = 0);
    ~charactersDialog();
    
private:
    Ui::charactersDialog *ui;
};

#endif // CHARACTERSDIALOG_H
