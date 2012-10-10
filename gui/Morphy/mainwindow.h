#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
class QAction;
class QActionGroup;
class QLabel;
class QMenu;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_actionHeuristic_Search_triggered();

    void on_actionSettsTaxaAndOutgroup_triggered();

    void on_actionSettsCharacters_triggered();

    void on_toolButton_2_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
