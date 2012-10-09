#include "charactersdialog.h"
#include "ui_charactersdialog.h"

charactersDialog::charactersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::charactersDialog)
{
    ui->setupUi(this);
}

charactersDialog::~charactersDialog()
{
    delete ui;
}
