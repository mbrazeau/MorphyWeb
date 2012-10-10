#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "heuristicsearchdialog.h"
#include "charactersdialog.h"
#include "taxadialog.h"
#include "QDialog"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionHeuristic_Search_triggered()
{
    heuristicSearchDialog hsearchd;
    if (hsearchd.exec())
    {
        ui->mwTextDisplay->appendPlainText("\nInitiating Heuristic Search...\n");
    }
}

void MainWindow::on_actionSettsTaxaAndOutgroup_triggered()
{
    taxadialog taxd;
    taxd.exec();
}

void MainWindow::on_actionSettsCharacters_triggered()
{
    charactersDialog chard;
    chard.exec();
}


void MainWindow::on_toolButton_2_clicked()
{
    ui->mwTextDisplay->clear();
}
