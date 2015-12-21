#pragma once
#ifndef HDRAUTO_H
#define HDRAUTO_H


#include <QtWidgets/QMainWindow>
#include "ui_hdrauto.h"

#include <math.h>
#include <QDebug>
#include <QFile>
#include <QVector>
#include <QFileDialog>
#include <QTextStream>
#include <QTime>
#include <QStringList>
#include <QMessageBox>

#include <Eigen/Sparse>
template<class T> class Image;

class HDRAuto : public QMainWindow
{
	Q_OBJECT

public:
	HDRAuto(QWidget *parent = 0);
	~HDRAuto();

private:
	Ui::HDRAutoClass ui;

	QTime					CountTime;															//�p��ɶ�

	QFileDialog				*Dialog;
	QVector<QImage *>		ImageStack;															//��Ҧ����Ϧs�i��
	QVector<double>			ShutterStack;														//��Ҧ����n���ɶ��s�_��
	

	// Cut Image to Fit
	bool					ReadImageAndCut(QString, int, QTextStream *, Image<unsigned char> *, double *);//����QT Ū�� �A�h��ϵ��ŤU��

	// HDR �����Ʊ�
	void					DoHDRMain(QString, Image<unsigned char> *, double*, int);			
	void					HDR_image(Image<double> &, Image<unsigned char> *, int, double *, double *, double **);
	void					Gsolve(double **, unsigned char **, double *, int, double *, int, int);


	const bool				DebugMode = false;													//DebugMode �O�_�n�}��									
private slots:
	void					OpenFileEvent();													//�}���ɮת��ƥ�
};

#endif // HDRAUTO_H
