/********************************************************************************
** Form generated from reading UI file 'hdrauto.ui'
**
** Created by: Qt User Interface Compiler version 5.5.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_HDRAUTO_H
#define UI_HDRAUTO_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_HDRAutoClass
{
public:
    QAction *actionFile;
    QAction *action_2;
    QAction *actionOpen_File;
    QWidget *centralWidget;
    QMenuBar *menuBar;
    QMenu *menuFile;

    void setupUi(QMainWindow *HDRAutoClass)
    {
        if (HDRAutoClass->objectName().isEmpty())
            HDRAutoClass->setObjectName(QStringLiteral("HDRAutoClass"));
        HDRAutoClass->resize(600, 400);
        actionFile = new QAction(HDRAutoClass);
        actionFile->setObjectName(QStringLiteral("actionFile"));
        action_2 = new QAction(HDRAutoClass);
        action_2->setObjectName(QStringLiteral("action_2"));
        actionOpen_File = new QAction(HDRAutoClass);
        actionOpen_File->setObjectName(QStringLiteral("actionOpen_File"));
        centralWidget = new QWidget(HDRAutoClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        HDRAutoClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(HDRAutoClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 600, 22));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QStringLiteral("menuFile"));
        HDRAutoClass->setMenuBar(menuBar);

        menuBar->addAction(menuFile->menuAction());
        menuFile->addAction(actionOpen_File);

        retranslateUi(HDRAutoClass);

        QMetaObject::connectSlotsByName(HDRAutoClass);
    } // setupUi

    void retranslateUi(QMainWindow *HDRAutoClass)
    {
        HDRAutoClass->setWindowTitle(QApplication::translate("HDRAutoClass", "HDRAuto", 0));
        actionFile->setText(QApplication::translate("HDRAutoClass", "File", 0));
        action_2->setText(QApplication::translate("HDRAutoClass", "\351\226\213\345\225\237\346\252\224\346\241\210", 0));
        actionOpen_File->setText(QApplication::translate("HDRAutoClass", "Open File", 0));
        menuFile->setTitle(QApplication::translate("HDRAutoClass", "File", 0));
    } // retranslateUi

};

namespace Ui {
    class HDRAutoClass: public Ui_HDRAutoClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_HDRAUTO_H
