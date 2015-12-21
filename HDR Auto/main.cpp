#include "hdrauto.h"
#include <QtWidgets/QApplication>
int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	HDRAuto w;
	w.show();
	return a.exec();
}
