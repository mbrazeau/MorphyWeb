#ifndef TAXADIALOG_H
#define TAXADIALOG_H

#include <QDialog>

namespace Ui {
class taxadialog;
}

class taxadialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit taxadialog(QWidget *parent = 0);
    ~taxadialog();
    
private:
    Ui::taxadialog *ui;
};

#endif // TAXADIALOG_H
