#include "heuristicsearchdialog.h"
#include "ui_heuristicsearchdialog.h"
#include "mainwindow.h"
#include <iostream>

heuristicSearchDialog::heuristicSearchDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::heuristicSearchDialog)
{
    ui->setupUi(this);
}

heuristicSearchDialog::~heuristicSearchDialog()
{
    delete ui;
}

void heuristicSearchDialog::on_buttonBox_accepted()
{
    // Want to print "Initiating heuristic search (and then describe search params)"
    //std::cout << "Initiating heuristic search ... [NOT REALLY!]" << std::endl;

}
