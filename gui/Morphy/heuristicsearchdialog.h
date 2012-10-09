#ifndef HEURISTICSEARCHDIALOG_H
#define HEURISTICSEARCHDIALOG_H

#include <QDialog>

namespace Ui {
class heuristicSearchDialog;
}

class heuristicSearchDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit heuristicSearchDialog(QWidget *parent = 0);
    ~heuristicSearchDialog();
    
private slots:
    void on_buttonBox_accepted();

private:
    Ui::heuristicSearchDialog *ui;
};

#endif // HEURISTICSEARCHDIALOG_H
