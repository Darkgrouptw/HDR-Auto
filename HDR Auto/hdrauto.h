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

	QTime					CountTime;															//計算時間

	QFileDialog				*Dialog;
	QVector<QImage *>		ImageStack;															//把所有的圖存進來
	QVector<double>			ShutterStack;														//把所有的曝光時間存起來
	

	// Cut Image to Fit
	bool					ReadImageAndCut(QString, int, QTextStream *, Image<unsigned char> *, double *);//先用QT 讀檔 再去把圖裁剪下來

	// HDR 做的事情
	void					DoHDRMain(QString, Image<unsigned char> *, double*, int);			
	void					HDR_image(Image<double> &, Image<unsigned char> *, int, double *, double *, double **);
	void					Gsolve(double **, unsigned char **, double *, int, double *, int, int);


	const bool				DebugMode = false;													//DebugMode 是否要開啟									
private slots:
	void					OpenFileEvent();													//開啟檔案的事件
};

#endif // HDRAUTO_H
