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
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_HDRAutoClass
{
public:
    QAction *actionFile;
    QAction *action_2;
    QAction *actionOpen_File;
    QWidget *centralWidget;
    QLabel *MaskTitle;
    QLabel *MaskPic;
    QPushButton *StartButton;
    QPushButton *SetMask;
    QLabel *ThresholdTitle1;
    QSlider *Slider1;
    QLabel *ThresholdText1;
    QSlider *Slider2;
    QLabel *ThresholdTitle2;
    QLabel *ThresholdText2;
    QMenuBar *menuBar;
    QMenu *menuFile;

    void setupUi(QMainWindow *HDRAutoClass)
    {
        if (HDRAutoClass->objectName().isEmpty())
            HDRAutoClass->setObjectName(QStringLiteral("HDRAutoClass"));
        HDRAutoClass->resize(742, 585);
        QPalette palette;
        QBrush brush(QColor(255, 255, 255, 255));
        brush.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::WindowText, brush);
        QBrush brush1(QColor(0, 0, 0, 255));
        brush1.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Button, brush1);
        palette.setBrush(QPalette::Active, QPalette::Light, brush1);
        palette.setBrush(QPalette::Active, QPalette::Midlight, brush1);
        palette.setBrush(QPalette::Active, QPalette::Dark, brush1);
        palette.setBrush(QPalette::Active, QPalette::Mid, brush1);
        palette.setBrush(QPalette::Active, QPalette::Text, brush);
        palette.setBrush(QPalette::Active, QPalette::BrightText, brush);
        palette.setBrush(QPalette::Active, QPalette::ButtonText, brush);
        palette.setBrush(QPalette::Active, QPalette::Base, brush1);
        palette.setBrush(QPalette::Active, QPalette::Window, brush1);
        palette.setBrush(QPalette::Active, QPalette::Shadow, brush1);
        palette.setBrush(QPalette::Active, QPalette::AlternateBase, brush1);
        QBrush brush2(QColor(255, 255, 220, 255));
        brush2.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::ToolTipBase, brush2);
        palette.setBrush(QPalette::Active, QPalette::ToolTipText, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::WindowText, brush);
        palette.setBrush(QPalette::Inactive, QPalette::Button, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Light, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Midlight, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Dark, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Mid, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Text, brush);
        palette.setBrush(QPalette::Inactive, QPalette::BrightText, brush);
        palette.setBrush(QPalette::Inactive, QPalette::ButtonText, brush);
        palette.setBrush(QPalette::Inactive, QPalette::Base, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Window, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Shadow, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::AlternateBase, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::ToolTipBase, brush2);
        palette.setBrush(QPalette::Inactive, QPalette::ToolTipText, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::WindowText, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Button, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Light, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Midlight, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Dark, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Mid, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Text, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::BrightText, brush);
        palette.setBrush(QPalette::Disabled, QPalette::ButtonText, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Window, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Shadow, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::AlternateBase, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::ToolTipBase, brush2);
        palette.setBrush(QPalette::Disabled, QPalette::ToolTipText, brush1);
        HDRAutoClass->setPalette(palette);
        actionFile = new QAction(HDRAutoClass);
        actionFile->setObjectName(QStringLiteral("actionFile"));
        action_2 = new QAction(HDRAutoClass);
        action_2->setObjectName(QStringLiteral("action_2"));
        actionOpen_File = new QAction(HDRAutoClass);
        actionOpen_File->setObjectName(QStringLiteral("actionOpen_File"));
        centralWidget = new QWidget(HDRAutoClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        MaskTitle = new QLabel(centralWidget);
        MaskTitle->setObjectName(QStringLiteral("MaskTitle"));
        MaskTitle->setGeometry(QRect(10, 10, 111, 41));
        QFont font;
        font.setFamily(QStringLiteral("Arial"));
        font.setPointSize(20);
        MaskTitle->setFont(font);
        MaskPic = new QLabel(centralWidget);
        MaskPic->setObjectName(QStringLiteral("MaskPic"));
        MaskPic->setGeometry(QRect(20, 50, 500, 500));
        StartButton = new QPushButton(centralWidget);
        StartButton->setObjectName(QStringLiteral("StartButton"));
        StartButton->setEnabled(false);
        StartButton->setGeometry(QRect(560, 50, 151, 41));
        SetMask = new QPushButton(centralWidget);
        SetMask->setObjectName(QStringLiteral("SetMask"));
        SetMask->setEnabled(false);
        SetMask->setGeometry(QRect(560, 160, 151, 41));
        ThresholdTitle1 = new QLabel(centralWidget);
        ThresholdTitle1->setObjectName(QStringLiteral("ThresholdTitle1"));
        ThresholdTitle1->setGeometry(QRect(570, 260, 121, 61));
        QFont font1;
        font1.setFamily(QStringLiteral("Arial"));
        font1.setPointSize(16);
        ThresholdTitle1->setFont(font1);
        Slider1 = new QSlider(centralWidget);
        Slider1->setObjectName(QStringLiteral("Slider1"));
        Slider1->setEnabled(false);
        Slider1->setGeometry(QRect(550, 320, 160, 19));
        Slider1->setMaximum(255);
        Slider1->setValue(7);
        Slider1->setOrientation(Qt::Horizontal);
        ThresholdText1 = new QLabel(centralWidget);
        ThresholdText1->setObjectName(QStringLiteral("ThresholdText1"));
        ThresholdText1->setGeometry(QRect(580, 340, 121, 61));
        ThresholdText1->setFont(font1);
        Slider2 = new QSlider(centralWidget);
        Slider2->setObjectName(QStringLiteral("Slider2"));
        Slider2->setEnabled(false);
        Slider2->setGeometry(QRect(550, 450, 160, 19));
        Slider2->setMaximum(255);
        Slider2->setValue(127);
        Slider2->setOrientation(Qt::Horizontal);
        ThresholdTitle2 = new QLabel(centralWidget);
        ThresholdTitle2->setObjectName(QStringLiteral("ThresholdTitle2"));
        ThresholdTitle2->setGeometry(QRect(570, 390, 121, 61));
        ThresholdTitle2->setFont(font1);
        ThresholdText2 = new QLabel(centralWidget);
        ThresholdText2->setObjectName(QStringLiteral("ThresholdText2"));
        ThresholdText2->setGeometry(QRect(580, 470, 121, 61));
        ThresholdText2->setFont(font1);
        HDRAutoClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(HDRAutoClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 742, 22));
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
        MaskTitle->setText(QApplication::translate("HDRAutoClass", "Mask\357\274\232", 0));
        MaskPic->setText(QString());
        StartButton->setText(QApplication::translate("HDRAutoClass", "Start", 0));
        SetMask->setText(QApplication::translate("HDRAutoClass", "Set Mask", 0));
        ThresholdTitle1->setText(QApplication::translate("HDRAutoClass", "Threshold 1", 0));
        ThresholdText1->setText(QApplication::translate("HDRAutoClass", "7", 0));
        ThresholdTitle2->setText(QApplication::translate("HDRAutoClass", "Threshold 2", 0));
        ThresholdText2->setText(QApplication::translate("HDRAutoClass", "127", 0));
        menuFile->setTitle(QApplication::translate("HDRAutoClass", "File", 0));
    } // retranslateUi

};

namespace Ui {
    class HDRAutoClass: public Ui_HDRAutoClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_HDRAUTO_H
