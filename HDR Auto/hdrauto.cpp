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

	if (!QDir(path + "Cut").exists())
		QDir().mkdir(path + "Cut");

	for (int i = 0; i < nImage; i++)
	{
		(*ss) >> r;
		(*ss) >> B_shutter[i];
		QImage temp(path + r);

		if (temp.width() == 2592 && temp.height() == 1728)
		{
			qDebug(qPrintable(r + " Norma Size => 2592 x 1728"));

			temp = temp.copy(567, 153, 1453, 1453);
			temp.save(path + "Cut/" + r);
			imgs[i].readImage((path + "Cut/" + r).toLocal8Bit().data());
			imgs[i].sRGBTolsRGB_255();
			B_shutter[i] = log(1. / B_shutter[i]);
		}
		else if (temp.width() == 1453 && temp.height() == 1453)
		{
			qDebug(qPrintable(r + " Norma Size => 1453 x 1453"));
			imgs[i].readImage((path + r).toLocal8Bit().data());
			imgs[i].sRGBTolsRGB_255();
			B_shutter[i] = log(1. / B_shutter[i]);
		}
		else
		{
			qDebug("ReadImageAndCut error");
			return false;
		}
	}

	//Delete Temp File
	if (!DebugMode)
		if (!QDir(path + "Cut").removeRecursively())
			qDebug() << "Remove Temp Directory error";

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
	//printf("HDR_image\n");

	int img_col = img[0].height;
	int img_row = img[0].width;

	//EnergyMap
	//Image<double> EnergyMap(img_col, img_row);
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
void HDRAuto::DoHDRMain(QString FilePath, Image<unsigned char> *imgs, double* B_shutter, int nImage)
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
	for (int i = 0; i<pixelNum; i++)
	{
		//printf("i=%d, w=%d, h=%d\n", i, img_row, img_col);
		int tx = rand() % img_col;
		int ty = rand() % img_row;
		if (tx < 0) tx *= -1;
		if (ty < 0) ty *= -1;
		//printf("i=%d, tx=%d, ty=%d\n", i, tx, ty);
		for (int j = 0; j< nImage; j++)
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

	//*--- get HDR image
	Image<double> HDRimg;
	HDR_image(HDRimg, imgs, nImage, wt, B_shutter, G);

	qDebug("HDR finish");

	//* write HDR
	HDRimg.writeImage((FilePath + "HDR").toLocal8Bit().data(), I_HDR_T);
	qDebug("write HDR image done");

	if (DebugMode)
	{
		Image<unsigned char> toneLDR = HDRimg.tomemapping();
		toneLDR.writeImage((FilePath + "tonemappedLDR").toLocal8Bit().data(), I_JPG_T);
		printf("write tonemapped image done\n");

		Image<unsigned char> LDRimg = HDRimg.toByte();
		LDRimg.lsRGBTosRGB_255();
		LDRimg.writeImage((FilePath + "LDR").toLocal8Bit().data(), I_JPG_T);
		qDebug("write LDR image done");
	}

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

	qDebug() << "DoHDRMain " << QString::number(CountTime.elapsed() / 1000.0, 'f', 3) << " s";
}
void HDRAuto::OpenFileEvent()
{
	//讀 txt 檔
	qDebug() << "========== Open File ==========";
	QString FilePath = QFileDialog::getOpenFileName(this,
		tr("Open Text"),
		"C:/Users/Graphics/Desktop/HDR Auto/Win32/Release/HDR Pictures/",
		tr("Text Files(*.txt)"));
	QFile file(FilePath);
	qDebug() << "Path: " << FilePath;
	if (file.open(QIODevice::ReadOnly))
	{
		//讀檔案
		QTextStream *ss = new QTextStream(&file);
		int nImage;

		(*ss) >> nImage;
		Image<unsigned char> *imgs = new Image<unsigned char>[nImage];
		double *B_shutter = new double[nImage];

		QStringList strlist = FilePath.split("/");
		FilePath = "";
		for (int i = 0; i<strlist.length() - 1; i++)
			FilePath += strlist[i] + "/";
		qDebug() << "Location => "<< FilePath;

		if (ReadImageAndCut(FilePath, nImage, ss, imgs, B_shutter))
			DoHDRMain(FilePath, imgs, B_shutter, nImage);
		file.close();
		qDebug() << "========== Success ==========";
	}
}