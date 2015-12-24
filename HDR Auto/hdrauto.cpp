#include "hdrauto.h"
#include "Image.hpp"

HDRAuto::HDRAuto(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);

	//Connect 事件
	connect(ui.actionOpen_File, SIGNAL(triggered()), this, SLOT(OpenFileEvent()));
}
HDRAuto::~HDRAuto()
{

}

bool HDRAuto::ReadImageAndCut(QString path, int nImage, QTextStream *ss, Image<unsigned char> *imgs, double *B_shutter)
{

	qDebug() << "========== ReadImage And CutImage ==========";
	CountTime.start();

	//ReadImage
	QVector<QImage> ImagePack;
	QString r;

	//if (!QDir(path + "Cut").exists())
		//QDir().mkdir(path + "Cut");

	for (int i = 0; i < nImage; i++)
	{
		(*ss) >> r;
		(*ss) >> B_shutter[i];
		QImage temp(path + r);

		if (temp.width() == 2592 && temp.height() == 1728)
		{
			imgs[i].readImage((path + r).toLocal8Bit().data());
			imgs[i].sRGBTolsRGB_255();
			B_shutter[i] = log(1. / B_shutter[i]);
		}
		else
		{
			qDebug("ReadImageAndCut Image Size is not 2592 x 1728");
			return false;
		}
	}

	CenterX = imgs[0].width / 2;
	CenterY = imgs[0].height / 2;

	//Delete Temp File
	/*if (!DebugMode)
		if (!QDir(path + "Cut").removeRecursively())
		qDebug() << "Remove Temp Directory error";*/

	qDebug() << "ReadImageAndCut => " << CountTime.elapsed() / 1000.0 << " s";
	return true;
}

void HDRAuto::Gsolve(double **x, unsigned char **z, double B_shutter[], int l_size, double w[], int pixelNum, int nImage)
{
	int n = 256;
	Eigen::SparseMatrix<double, Eigen::RowMajor> A1(pixelNum*nImage + n + 1, n + pixelNum);
	Eigen::VectorXd B = Eigen::RowVectorXd::Zero(pixelNum*nImage + n + 1);

	A1.reserve(Eigen::VectorXi::Constant(pixelNum*nImage + n + 1, 3));

	int k = 0;
	for (int i = 0; i < pixelNum; i++)
	{
		for (int j = 0; j < nImage; j++)
		{
			double wij = w[z[i][j]];

			A1.insert(k, z[i][j]) = wij;
			A1.insert(k, n + i) = -wij;

			B[k] = wij * B_shutter[j];
			k++;
		}
	}

	A1.insert(k, 128) = 1.0;
	k++;

	for (int i = 0; i < n - 2; i++)
	{
		A1.insert(k, i) = l_size*w[i + 1];
		A1.insert(k, i + 1) = -2 * l_size*w[i + 1];
		A1.insert(k, i + 2) = l_size*w[i + 1];
		k++;
	}
	A1.makeCompressed();


	Eigen::SimplicialLLT< Eigen::SparseMatrix<double>, Eigen::RowMajor > solver;

	Eigen::VectorXd ATB = A1.transpose() * B;
	Eigen::SparseMatrix<double, Eigen::RowMajor> ATA = A1.transpose() * A1;
	solver.compute(ATA);

	Eigen::VectorXd _x = solver.solve(ATB);

	*x = new double[n + pixelNum];
	double *test = _x.data();
	std::copy(_x.data(), _x.data() + n + pixelNum, *x);
}
void HDRAuto::HDR_image(Image<double> &EnergyMap, Image<unsigned char> *img, int filesize, double *wt, double *B_shutter, double **G)
{
	if (img == NULL || wt == NULL || B_shutter == NULL || G == NULL){
		printf("input not enough: %s == NULL\n",
			img == NULL ? "img" : wt == NULL ? "wt" : B_shutter == NULL ? "B_shutter" : "G");
		//Image<double> d;
		return;
	}

	int img_col = img[0].height;
	int img_row = img[0].width;

	EnergyMap.init(img_col, img_row, 1);
	//printf("init HDR image: %d, %d\n", img_col, img_row);

	for (int i = 0; i < img_col; i++)
	{
		for (int j = 0; j < img_row; j++)
		{
			//printf("i=%d, j=%d\n", i, j);
			double t1 = 0, t2 = 0;
			for (int now = 0; now < filesize; now++)
			{
				int z1 = img[now].R[i][j];
				t1 = t1 + wt[z1] * (G[0][z1] - B_shutter[now]);
				t2 = t2 + wt[z1];
			}
			//printf("R: t1=%lg, t2=%lg, exp=%lg\n", t1, t2, exp(t1/t2));
			EnergyMap.R[i][j] = exp(t1 / t2);
			if (t2 == 0 && t1 == 0) EnergyMap.R[i][j] = 0;

			t1 = 0, t2 = 0;
			for (int now = 0; now < filesize; now++)
			{
				int z2 = img[now].G[i][j];
				t1 = t1 + wt[z2] * (G[1][z2] - B_shutter[now]);
				t2 = t2 + wt[z2];
			}
			//printf("G: t1=%lg, t2=%lg, exp=%lg\n", t1, t2, exp(t1/t2));
			EnergyMap.G[i][j] = exp(t1 / t2);
			if (t2 == 0 && t1 == 0) EnergyMap.G[i][j] = 0;

			t1 = 0, t2 = 0;
			for (int now = 0; now < filesize; now++)
			{
				int z3 = img[now].B[i][j];
				t1 = t1 + wt[z3] * (G[2][z3] - B_shutter[now]);
				t2 = t2 + wt[z3];
			}
			//printf("B: t1=%lg, t2=%lg, exp=%lg\n", t1, t2, exp(t1/t2));
			EnergyMap.B[i][j] = exp(t1 / t2);
			if (t2 == 0 && t1 == 0) EnergyMap.B[i][j] = 0;
		}
	}
}
bool HDRAuto::Check_In_Area(int tx, int ty)
{
	if (Radius >= sqrt((CenterX - tx) * (CenterX - tx) + (CenterY - ty) * (CenterY - ty)))
		return true;
	return false;
}
bool HDRAuto::IsUseAblePoint(Image<unsigned char>* imgs, int nImage, int tx, int ty)
{
	if (imgs[0].R[tx][ty] == 0 && imgs[0].G[tx][ty] == 0 && imgs[0].B[tx][ty] == 0)
		return false;
	else if (imgs[nImage - 1].R[tx][ty] == 255 && imgs[nImage - 1].G[tx][ty] == 255 && imgs[nImage - 1].B[tx][ty] == 255)
		return false;
	return true;
}
void HDRAuto::DoHDRMain(QString FilePath, Image<unsigned char> *imgs, Image<double> &HDRimg, double* B_shutter, int nImage)
{
	qDebug() << "========== HDR Main ==========";
	CountTime.restart();
	//*--- init weight function -----------ok
	double wt[256];
	int count = 127;
	for (int i = 0; i<256; i++)
		if (i<128)
			wt[i] = (double)i / 16256.0;
		else
		{
			wt[i] = (double)count / 16256.0;
			count--;
		}

	//*--- init z
	srand(time(NULL));
	const int pixelNum = 100;
	int img_row = imgs[0].width;
	int img_col = imgs[0].height;
	unsigned char *_z[3][pixelNum];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < pixelNum; j++)
			_z[i][j] = new unsigned char[nImage];
	unsigned char *__z[3][pixelNum];
	unsigned char **z[3];
	z[0] = __z[0];
	z[1] = __z[1];
	z[2] = __z[2];
	for (int i = 0; i<pixelNum; i++){
		z[0][i] = _z[0][i];
		z[1][i] = _z[1][i];
		z[2][i] = _z[2][i];
	}
	bool Bool_CurrentAreaPoint = false;
	int tx, ty;

	for (int i = 0; i<pixelNum; i++)
	{
		do 
		{
			tx = rand() % img_col;
			ty = rand() % img_row;
			if (tx < 0) tx *= -1;
			if (ty < 0) ty *= -1;
		} while (!Check_In_Area(ty,tx) || !IsUseAblePoint(imgs,nImage,tx,ty));
		for (int j = 0; j < nImage; j++)
		{
			//printf("j=%d\n", j);
			z[0][i][j] = imgs[j].R[tx][ty];
			z[1][i][j] = imgs[j].G[tx][ty];
			z[2][i][j] = imgs[j].B[tx][ty];
		}
	}
	qDebug("init z");

	//*--- Get Response function 
	double *x[3] = {};

	//前半段到256的解給g 後半段的解給LE 
	//double _G[3][256];
	double *G[3];
	G[0] = new double[256];
	G[1] = new double[256];
	G[2] = new double[256];

	//double _LE[3][pixelNum];
	double *LE[3];
	LE[0] = new double[pixelNum];
	LE[1] = new double[pixelNum];
	LE[2] = new double[pixelNum];

	for (int i = 0; i<3; i++)
	{
		//解x ,z ,曝光值陣列B, l_size = 100 ,weight function,pixelNum,i=R G B 三個
		Gsolve(&x[i], z[i], B_shutter, 100, wt, pixelNum, nImage);
		for (int j = 0; j<256 + pixelNum; j++){
			if (j<256)
				G[i][j] = x[i][j];
			else
				LE[i][j - 256] = x[i][j];
		}
	}
	printf("solved\n");

	HDR_image(HDRimg, imgs, nImage, wt, B_shutter, G);

	qDebug("HDR finish");

	//* write HDR
	HDRimg.writeImage((FilePath + "HDR").toLocal8Bit().data(), I_HDR_T);
	qDebug("write HDR image done");

	/*if (DebugMode)
	{
	Image<unsigned char> toneLDR = HDRimg.tomemapping();
	toneLDR.writeImage((FilePath + "tonemappedLDR").toLocal8Bit().data(), I_JPG_T);
	printf("write tonemapped image done\n");

	Image<unsigned char> LDRimg = HDRimg.toByte();
	LDRimg.lsRGBTosRGB_255();
	LDRimg.writeImage((FilePath + "LDR").toLocal8Bit().data(), I_JPG_T);
	qDebug("write LDR image done");
	}*/

	//*--- release memory
	delete[] G[0];
	delete[] G[1];
	delete[] G[2];
	delete[] LE[0];
	delete[] LE[1];
	delete[] LE[2];
	delete[] x[0];
	delete[] x[1];
	delete[] x[2];

	qDebug() << "DoHDRMain => " << CountTime.elapsed() / 1000.0 << " s";
}


void HDRAuto::FillGrayColor(Image<double> *img, int i, int j, int Color)
{
	img->R[i][j] = Color;
	img->G[i][j] = Color;
	img->B[i][j] = Color;
}
void HDRAuto::ImgToGray(Image<double> *Mask, Image<unsigned char> *imgs, int PickNumber)
{
	for (int i = 0; i < Mask->height; i++)
		for (int j = 0; j < Mask->width; j++)
		{
			int Color = ((double)imgs[PickNumber].R[i][j] * 0.2989 +
				(double)imgs[PickNumber].G[i][j] * 0.587 +
				(double)imgs[PickNumber].B[i][j] * 0.114);
			if (Color > 255)
				Color = 255;
			FillGrayColor(Mask, i, j, Color);
		}
}
void HDRAuto::ImageBinarization(Image<double> *Mask, double tempInt)
{
	//////////////////////////////////////////////////////////////////////////
	// Histogram 假設RGB是0~255
	//////////////////////////////////////////////////////////////////////////
	/*int HistogramResult[256];
	for (int i = 0; i < 256; i++)
		HistogramResult[i] = 0;
	for (int i = 0; i < Mask->height; i++)
		for (int j = 0; j < Mask->width; j++)
			if (Check_In_Area(j, i))
				HistogramResult[(int)Mask->R[i][j]] ++;*/
	/*
	//////////////////////////////////////////////////////////////////////////
	//	判斷 是不是在半徑內的點 TotalPixelNumber
	//////////////////////////////////////////////////////////////////////////
	int TotalPixelNumber = 0;// Mask->width * Mask->height;
	for (int i = 0; i < imgs->height; i++)
		for (int j = 0; j < imgs->width; j++)
			if (Check_In_Area(j, i))
				TotalPixelNumber++;
	qDebug() << "TotalPixelNumber =>" << TotalPixelNumber;

	//////////////////////////////////////////////////////////////////////////
	// 做Histogram
	//////////////////////////////////////////////////////////////////////////
	int tempInt = 0;
	for (int i = 0; i < 256; i++)
		if ((double)TotalPixelNumber * threshold >= tempInt)
			tempInt += HistogramResult[i];
		else
		{
			tempInt = i;						//紀錄是哪一個Index
			break;
		}

	for (int i = 0; i < 256; i++)
		qDebug() << "i = "<< i << " = " << HistogramResult[i];*/
	qDebug() << "Threshold Color => " << tempInt;
	for (int i = 0; i < Mask->height; i++)
		for (int j = 0; j < Mask->width; j++)
			if (Mask->R[i][j] < tempInt)
				FillGrayColor(Mask, i, j, 0);
			else
				FillGrayColor(Mask, i, j, 255);
}
bool HDRAuto::IsInQueue(QVector<QVector2D> *Temp, int x, int y)
{
	if (Temp->size() == 0)
		return false;

	int Index = 0;
	while (Temp->size() > Index)
	{
		if ((*Temp)[Index].x() == x && (*Temp)[Index].y() == y)
			return true;
		Index++;
	}
	return false;
}
bool HDRAuto::IsSameColor(Image<double> *img, QVector2D a, QVector2D b)
{
	if (img->R[(int)a.y()][(int)a.x()] == img->R[(int)b.y()][(int)b.x()] &&
		img->G[(int)a.y()][(int)a.x()] == img->G[(int)b.y()][(int)b.x()] &&
		img->B[(int)a.y()][(int)a.x()] == img->B[(int)b.y()][(int)b.x()])
		return true;
	return false;
}
bool HDRAuto::CheckConditionInQueue(QVector<QVector2D> *Temp, int **PointTable, int dx, int dy, Image<double> *img)
{
	if ((*Temp)[0].x() + dx >= 0			&& (*Temp)[0].y() + dy >= 0 &&								// 左上 邊界
		(*Temp)[0].x() + dx < img->width	&& (*Temp)[0].y() + dy < img->height &&						// 右下 邊界
		PointTable[(int)(*Temp)[0].x() + dx][(int)(*Temp)[0].y() + dy] == -1 &&							// 被走過
		IsSameColor(img, (*Temp)[0], QVector2D((*Temp)[0].x() + dx, (*Temp)[0].y() + dy)) &&			// 是否是同一個顏色
		!IsInQueue(Temp, (*Temp)[0].x() + dx, (*Temp)[0].y() + dy))
		return true;
	return false;
}
void HDRAuto::Grouping(Image<double> *img)
{
	//////////////////////////////////////////////////////////////////////////
	//  分群 要先把整個陣列，先把
	//////////////////////////////////////////////////////////////////////////
	// 建一個Table 確定全部都沒有走過
	int **PointCheckTable = new int *[img->width];			// -1 代表沒走過，0代表跟邊界一樣顏色
	for (int i = 0; i < img->width; i++)
	{
		PointCheckTable[i] = new int[img->height];
		for (int j = 0; j < img->height ; j++)
			PointCheckTable[i][j] = -1;
	}

	// 再來用 Queue 做BFS (如果適用Stack，就是教DFS)
	int lastX = 0, lastY = 0, lastID = 0;
	bool TrackFinish = false;
	QVector<QVector<QVector2D> *> GroupingResult;
	QVector<QVector2D> *TempQueue;			// 要存到GroupResult的結果
	QVector<QVector2D> *Temp = new QVector<QVector2D>();				// BFS 的過程
	while (!TrackFinish)
	{
		Temp->push_back(QVector2D(lastX, lastY));
		TempQueue = new QVector<QVector2D>();
		while (Temp->size() != 0)
		{
			PointCheckTable[(int)(*Temp)[0].x()][(int)(*Temp)[0].y()] = lastID;

			if (CheckConditionInQueue(Temp, PointCheckTable, -1, 0, img))	// Left Pixel
			{
				Temp->push_back(QVector2D((*Temp)[0].x() - 1, (*Temp)[0].y()));
				TempQueue->push_back(QVector2D((*Temp)[0].x() - 1, (*Temp)[0].y()));
			}
			if (CheckConditionInQueue(Temp, PointCheckTable, 0, -1, img))	// Top Pixel
			{
				Temp->push_back(QVector2D((*Temp)[0].x(), (*Temp)[0].y() - 1));
				TempQueue->push_back(QVector2D((*Temp)[0].x(), (*Temp)[0].y() - 1));
			}
			if (CheckConditionInQueue(Temp, PointCheckTable, 1, 0, img))	// Right Pixel
			{
				Temp->push_back(QVector2D((*Temp)[0].x() + 1, (*Temp)[0].y()));
				TempQueue->push_back(QVector2D((*Temp)[0].x() + 1, (*Temp)[0].y()));
			}
			if (CheckConditionInQueue(Temp, PointCheckTable, 0, 1, img))	// Down Pixel
			{
				Temp->push_back(QVector2D((*Temp)[0].x(), (*Temp)[0].y() + 1));
				TempQueue->push_back(QVector2D((*Temp)[0].x(), (*Temp)[0].y() + 1));
			}
			Temp->pop_front();
		}
		GroupingResult.push_back(TempQueue);
		//qDebug() << "lastID = " << lastID << lastX << lastY;

		bool StartX = false;
		bool findResult = false;
		for (int j = lastY; j < img->height && !findResult; j++)
		{
			for (int i = (StartX ? 0 : lastX); i < img->width && !findResult; i++)
				if (PointCheckTable[i][j] == -1)
				{
					lastX = i;
					lastY = j;
					findResult = true;
				}
				else if (i == img->width - 1 && j == img->height - 1)
					TrackFinish = true;
				StartX = true;
		}
		lastID++;
	}

	//////////////////////////////////////////////////////////////////////////
	// 分群之後填顏色
	//////////////////////////////////////////////////////////////////////////
	lastID--;	// 把最後家的扣掉
	int *B = new int[lastID - 1];

	QTime time = QTime::currentTime();
	qsrand((uint)time.msec());
	for (int i = 1; i < lastID; i++)
		B[i -1] = qrand() * 101 % 255;

	//////////////////////////////////////////////////////////////////////////
	// 開始填顏色，分群，如果總面積不超過 0.2 就把他填
	//////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < img->height; i++)
		for (int j = 0; j < img->width; j++)
			if (PointCheckTable[j][i] != 0)
				if (GroupingResult[PointCheckTable[j][i]]->size() > Radius * Radius * 0.2)		// 把大的區塊的顏色塗白色
				{
					img->R[i][j] = 255;
					img->G[i][j] = 255;
					img->B[i][j] = 255;
				}
				else																			// 把小的區塊 填上面邊一格的顏色
				{
					img->R[i][j] = img->R[i - 1][j];
					img->G[i][j] = img->G[i - 1][j];
					img->B[i][j] = img->B[i - 1][j];
				}

}
void HDRAuto::FindMaskArea(QString FilePath, Image<double> &HDRimgs, Image<unsigned char> imgs[], int nImage)
{
	qDebug() << "========== Fink Mask Area ==========";
	CountTime.restart();
	int PickNumber = nImage /2;
	qDebug() << "Pick imgs[" << PickNumber << "] to find mask";
	Image<double> *Mask =new Image<double>(imgs[PickNumber].height, imgs[PickNumber].width);

	ImgToGray(Mask, imgs, PickNumber);
	ImageBinarization(Mask, 7);
	Mask = Mask->blur<double>(9);
	ImageBinarization(Mask, 127);
	Grouping(Mask);
	
	if (DebugMode)
		Mask->writeImage((FilePath + "Mask").toLocal8Bit().data(), I_JPG_T);


	qDebug() << "FinkMaskArea => " << CountTime.elapsed() / 1000.0 << " s";
}

void HDRAuto::OpenFileEvent()
{
	//讀 txt 檔
	qDebug() << "========== Open File ==========";
	QString FilePath = QFileDialog::getOpenFileName(this,
		tr("Open Text"),
		"C:/Users/Graphics/Desktop/SourceTree/HDR-Auto/Win32/Release/HDR Pictures/",
		tr("Text Files(*.txt)"));
	QFile file(FilePath);
	qDebug() << "Path: " << FilePath;
	if (file.open(QIODevice::ReadOnly))
	{
		//讀檔案
		QTextStream *ss = new QTextStream(&file);
		int nImage;

		(*ss) >> nImage;
		(*ss) >> Radius;
		Image<unsigned char> *imgs = new Image<unsigned char>[nImage];
		double *B_shutter = new double[nImage];

		QStringList strlist = FilePath.split("/");
		FilePath = "";
		for (int i = 0; i<strlist.length() - 1; i++)
			FilePath += strlist[i] + "/";
		qDebug() << "Location => "<< FilePath;

		// HDR Image
		Image<double> HDRimg;
		if (ReadImageAndCut(FilePath, nImage, ss, imgs, B_shutter))
		{
			DoHDRMain(FilePath, imgs, HDRimg, B_shutter, nImage);
			FindMaskArea(FilePath, HDRimg, imgs, nImage);
			qDebug() << "========== Success ==========";
		}
		else
			qDebug() << "========== Fail ==========";
		file.close();
	}
}