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
	bool					ReadImage(QString, int, QTextStream *, Image<unsigned char> *, double *);//讀檔近來
	// HDR 做的事情
	void					DoHDRMain(QString, Image<unsigned char> *, Image<double> &, double*, int);	
	bool					Check_In_Area(int, int);											// 判斷選的典試不適在相機的範圍裡
	bool					IsUseAblePoint(Image<unsigned char>*, int, int, int);				// 判斷兩個，第一個是選到的點在最亮的那張上面，是不是黑色的，第二個是判斷在最暗的那張上面，是不是白的

	void					HDR_image(Image<double> &, Image<unsigned char> *, int, double *, double *, double **);
	void					Gsolve(double **, unsigned char **, double *, int, double *, int, int);

	const int				Radius = 726;														//半徑
	int						CenterX;															//中心座標
	int						CenterY;

	// 找出哪裡是不是有遮住的算法
	void					CutImageToCube(Image<double> &);									//把邊界剪裁到 1453 x 1453
	void					FindMaskArea(QString, Image<double> &, Image<unsigned char> *, int, Image<double> *);
	void					FillGrayColor(Image<double> *, int, int, int);						// 把灰階值填進去
	void					ImgToGray(Image<double> *, Image<unsigned char> *, int);
	void					ImageBinarization(Image<double> *, double);

	// 2值化
	bool					IsInQueue(QVector<QVector2D> *, int, int);							// 確定抓進來的東西，沒有重複
	bool					IsSameColor(Image<double> *, QVector2D, QVector2D);					// 判斷兩個點的顏色一不一樣
	void					Grouping(Image<double> *);											// 分群化
	bool					CheckConditionInQueue(QVector<QVector2D> *, int **, int, int, Image<double> *);
	const double			threshold = 0.5;													// 取臨界值，在座二值化的時候，左邊的點數目要是 => 總共 x threshold
	
	// Texture Synthesis
	bool					BounaryCheck(int, int, int, int);									// 確定有沒有這個值
	void					TextureSynthesis(Image<double> &,Image<double> * ,QString);			// 有Mask只要填顏色而已~~

	// 把145個值算出來填顏色
	int						ChangeCoordinate_ToPatchIndex(int, int);							// 給仰角跟角度，傳出一個patch index
	int						CountForPatchIndex(int, int);										// 給兩個點判斷說是在哪一個patch
	void					RenderHDR_ToResult(Image<double> &, Image<double> *&);				// 把結果顯示在 Render

	//位移的部分
	const int				StartXPos = 564;
	const int				StartYPos = 155;
	const int				CubeLength = 1453;

	const double			ErrorArea = 0.001;
	
	
	const bool				DebugMode = true;													//DebugMode 是否要開啟									
private slots:
	void					OpenFileEvent();													//開啟檔案的事件
};

#endif // HDRAUTO_H
