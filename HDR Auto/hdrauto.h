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

	QTime					CountTime;															//計算時間

	QFileDialog				*Dialog;
	QVector<QImage *>		ImageStack;															//把所有的圖存進來
	QVector<double>			ShutterStack;														//把所有的曝光時間存起來
	

	// Cut Image to Fit
	bool					ReadImageAndCut(QString, int, QTextStream *, Image<unsigned char> *, double *);//先用QT 讀檔 再去把圖裁剪下來
	// HDR 做的事情
	void					DoHDRMain(QString, Image<unsigned char> *, Image<double> &, double*, int);	
	bool					Check_In_Area(int, int);											// 判斷選的典試不適在相機的範圍裡
	bool					IsUseAblePoint(Image<unsigned char>*, int, int, int);				// 判斷兩個，第一個是選到的點在最亮的那張上面，是不是黑色的，第二個是判斷在最暗的那張上面，是不是白的

	void					HDR_image(Image<double> &, Image<unsigned char> *, int, double *, double *, double **);
	void					Gsolve(double **, unsigned char **, double *, int, double *, int, int);

	int						Radius;																//半徑
	int						CenterX;															//中心座標
	int						CenterY;

	// 找出哪裡是不是有遮住的算法
	void					FindMaskArea(QString, Image<double> &, Image<unsigned char> *, int);
	void					FillGrayColor(Image<double> *, int, int, int);						// 把灰階值填進去
	void					ImgToGray(Image<double> *, Image<unsigned char> *, int);
	void					ImageBinarization(Image<double> *, double);

	// 2值化
	bool					IsInQueue(QVector<QVector2D>, int, int);							// 確定抓進來的東西，沒有重複
	bool					IsSameColor(Image<double> *, QVector2D, QVector2D);					// 判斷兩個點的顏色一不一樣
	void					Grouping(Image<double> *);											//分群化
	const double			threshold = 0.5;													// 取臨界值，在座二值化的時候，左邊的點數目要是 => 總共 x threshold

	const bool				DebugMode = true;													//DebugMode 是否要開啟									
private slots:
	void					OpenFileEvent();													//開啟檔案的事件
};

#endif // HDRAUTO_H
