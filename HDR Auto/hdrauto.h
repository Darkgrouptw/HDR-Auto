#pragma once
#ifndef HDRAUTO_H
#define HDRAUTO_H


#include <QtWidgets/QMainWindow>
#include "ui_hdrauto.h"

#include <math.h>
#include <QDebug>
#include <QFile>
#include <QVector>
#include <QVector2D>
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
	void					DoHDRMain(QString, Image<unsigned char> *, Image<double> &, double*, int);	
	bool					Check_In_Area(int, int);											// �P�_�諸��դ��A�b�۾����d���
	bool					IsUseAblePoint(Image<unsigned char>*, int, int, int);				// �P�_��ӡA�Ĥ@�ӬO��쪺�I�b�̫G�����i�W���A�O���O�¦⪺�A�ĤG�ӬO�P�_�b�̷t�����i�W���A�O���O�ժ�

	void					HDR_image(Image<double> &, Image<unsigned char> *, int, double *, double *, double **);
	void					Gsolve(double **, unsigned char **, double *, int, double *, int, int);

	int						Radius;																//�b�|
	int						CenterX;															//���߮y��
	int						CenterY;

	// ��X���̬O���O���B����k
	void					FindMaskArea(QString, Image<double> &, Image<unsigned char> *, int);
	void					FillGrayColor(Image<double> *, int, int, int);						// ��Ƕ��ȶ�i�h
	void					ImgToGray(Image<double> *, Image<unsigned char> *, int);
	void					ImageBinarization(Image<double> *, double);

	// 2�Ȥ�
	bool					IsInQueue(QVector<QVector2D>, int, int);							// �T�w��i�Ӫ��F��A�S������
	bool					IsSameColor(Image<double> *, QVector2D, QVector2D);					// �P�_����I���C��@���@��
	void					Grouping(Image<double> *);											//���s��
	const double			threshold = 0.5;													// ���{�ɭȡA�b�y�G�Ȥƪ��ɭԡA���䪺�I�ƥحn�O => �`�@ x threshold

	const bool				DebugMode = true;													//DebugMode �O�_�n�}��									
private slots:
	void					OpenFileEvent();													//�}���ɮת��ƥ�
};

#endif // HDRAUTO_H
