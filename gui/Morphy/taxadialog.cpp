#include "taxadialog.h"
#include "ui_taxadialog.h"

taxadialog::taxadialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::taxadialog)
{
    ui->setupUi(this);
}

taxadialog::~taxadialog()
{
    delete ui;
}
