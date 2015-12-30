#pragma once

#include <png.h>
extern "C"{
#include <jpeglib.h>
}
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "FileHDR.hpp"


template<class MapT>
MapT **newTex(int height, int width){
	MapT **d = new MapT*[height];
	for(int i=0; i<height; i++)
		d[i] = new MapT[width];
	
	return d;
}

template<class MapT>
void delTex(MapT **r, int height){
	for(int i=0; i<height; i++)
		delete[] r[i];
	delete[] r;
}

#define I_PNG_T 0
#define I_JPG_T 1
#define I_HDR_T 2

template<class T>
class Image{
	public:
	char *name;            //file name
	int pixelFormat;       //0: gray, 1:RGB
	int width, height;
	T **R;
	T **G;                 //green or gray
	T **B;
	//T *data;               //data for texture
	
	Image(): name(NULL), pixelFormat(-1), width(0), height(0), R(NULL), G(NULL), B(NULL)
	{}
	
	Image(const Image<T> &r){
		if(r.name != NULL){
			int namel = strlen(r.name);
			this->name = new char[namel];
			for(int i=0; i<namel; i++){
				this->name[i] = r.name[i];
			}
		}
		else  this->name = NULL;
		
		this->init(r.height, r.width, r.pixelFormat);
		
		if(r.pixelFormat){
			for(int i=0; i<height; i++){
				for(int j=0; j<width; j++){
					this->R[i][j] = r.R[i][j];
					this->G[i][j] = r.G[i][j];
					this->B[i][j] = r.B[i][j];
				}
			}
		}
		else{
			for(int i=0; i<height; i++){
				for(int j=0; j<width; j++){
					this->G[i][j] = r.G[i][j];
				}
			}
		}
	}
	
	Image(char *imgPath){
		if(!readImage(imgPath)){
			printf("read file %s fail\n", imgPath);
			clear();
		}
	}
	
	Image(int h, int w, int format =1){
		name = NULL;
		init(h, w, format);
	}
	
	~Image(){
		clear();
	}
	
	//bool readImage(char *imgPath);
	//void clear();
	
	private: 
	//char *getNameFromPath(char *r);
	int dwanJy(const char *a, const char *b) const { //判斷兩字串是否相等
	    while(*b != 0){
	        if(*a != *b){return 0;}
	        a++, b++;
		}
	    return 1;
	}
	//bool ReadPngData( FILE* rFile);
	//bool ReadJpgData( FILE* rFile);
	
	//implememt
	public:
	bool readImage(char *imgPath){
		char *postfix = getNameFromPath(imgPath);
		
		//get file then check
	    FILE *fp = fopen( imgPath, "rb" );
	    if( !fp ) {
	    	printf("can't open file : %s\n'", imgPath);
			return false;            //檔案開啟失敗傳回 -1
		}
		
		//read file
		bool d;
		if(dwanJy(postfix, "png")){
			d = ReadPngData(fp);
		}
		else if(dwanJy(postfix, "jpg") || dwanJy(postfix, "JPG")){
			d = ReadJpgData(fp);
		}
		else{
			printf("unknown file format: %s\n", postfix);
			d = false;
		}
		fclose( fp );
		
		delete[] postfix;
		return d;
	}
	
	bool writeImage(char *fileName, int type){
		char oName[10000] = {};
		FILE *fp = NULL;
		
		switch(type){
		case I_JPG_T:
			sprintf(oName, "%s.jpg", fileName);
			fp = fopen(oName, "wb");
			return WriteJpg(fp);
		case I_HDR_T:
			sprintf(oName, "%s.hdr", fileName);
			fp = fopen(oName, "wb");
			return WriteHDR(fp);
		default:
			printf("unknown file format: %d\n", type);
			return false;
		}
		
	    return false;
	}
	
	void clear(){
		if(name != NULL){
			delete[] name;
			name = NULL;
		}
		if(R != NULL){
			for(int i=0; i<height; i++){
		    	delete[] R[i];
		    }
		    delete[] R;
		    R = NULL;
		}
		if(G != NULL){
			for(int i=0; i<height; i++){
		    	delete[] G[i];
		    }
		    delete[] G;
		    G = NULL;
		}
		if(B != NULL){
			for(int i=0; i<height; i++){
		    	delete[] B[i];
		    }
		    delete[] B;
		    B = NULL;
		}
		width  = 0;
		height = 0;
		pixelFormat = -1;
	}
	
	void init(int h, int w, int format){
		width = w;
		height= h;
		pixelFormat = format;
		
		if(format){
			//RGB
			R = new T*[height];
			G = new T*[height];
			B = new T*[height];
			for(int i=0; i<height; i++){
				R[i] = new T[width];
				G[i] = new T[width];
				B[i] = new T[width];
			}
		}
		else{
			//Gray
			R = NULL;
			G = new T*[height];
			B = NULL;
			for(int i=0; i<height; i++){
				G[i] = new T[width];
			}
		}
	}
	
	void print(){
		printf("image info\n");
		printf("name   : %s\n", name);
		printf("format : %s\n", pixelFormat ? "RGB" : "Gray");
		printf("width  : %d\n", width);
		printf("height : %d\n", height);
	}
	
	// 0-255 -> 0.0-1.0
	void selfToFloat(){
		for(int i=0; i<height; i++){
			for(int j=0; j<width; j++){
				R[i][j] = R[i][j] / (T)255;
				G[i][j] = G[i][j] / (T)255;
				B[i][j] = B[i][j] / (T)255;
			}
		}
	}
	
	// 0.0-1.0 -> 0-255
	void selfToByte(){
		for(int i=0; i<height; i++){
			for(int j=0; j<width; j++){
				R[i][j] = (T)(255.0 * R[i][j]);
				G[i][j] = (T)(255.0 * G[i][j]);
				B[i][j] = (T)(255.0 * B[i][j]);
			}
		}
	}
	
	// 0.0 ~ 1.0
	void sRGBTolsRGB(){
		for(int i=0; i<height; i++){
			for(int j=0; j<width; j++){
				double tr = (double)R[i][j];
				double tg = (double)G[i][j];
				double tb = (double)B[i][j];
				
				if ( tr > 0.04045 ) tr = pow((( tr + 0.055 ) / 1.055 ), 2.4);
				else                tr = tr / 12.92;
				if ( tr > 0.04045 ) tg = pow((( tg + 0.055 ) / 1.055 ), 2.4);
				else                tg = tg / 12.92;
				if ( tr > 0.04045 ) tb = pow((( tb + 0.055 ) / 1.055 ), 2.4);
				else                tb = tb / 12.92;
				
				R[i][j] = (T)tr;
				G[i][j] = (T)tg;
				B[i][j] = (T)tb;
			}
		}
	}
	
	// 0 ~ 255
	void sRGBTolsRGB_255(){
		for(int i=0; i<height; i++){
			for(int j=0; j<width; j++){
				double tr = (double)R[i][j] / 255.0;
				double tg = (double)G[i][j] / 255.0;
				double tb = (double)B[i][j] / 255.0;
				
				if ( tr > 0.04045 ) tr = pow((( tr + 0.055 ) / 1.055 ), 2.4);
				else                tr = tr / 12.92;
				if ( tr > 0.04045 ) tg = pow((( tg + 0.055 ) / 1.055 ), 2.4);
				else                tg = tg / 12.92;
				if ( tr > 0.04045 ) tb = pow((( tb + 0.055 ) / 1.055 ), 2.4);
				else                tb = tb / 12.92;
				
				R[i][j] = (T)(tr * 255.0);
				G[i][j] = (T)(tg * 255.0);
				B[i][j] = (T)(tb * 255.0);
			}
		}
	}
	
	// 0.0 ~ 1.0
	void lsRGBTosRGB(){
		for(int i=0; i<height; i++){
			for(int j=0; j<width; j++){
				double tr = (double)R[i][j];
				double tg = (double)G[i][j];
				double tb = (double)B[i][j];
				
				if ( tr > 0.0031308 ) tr = 1.055 * ( pow(tr, ( 1.0 / 2.4)) ) - 0.055;
				else                  tr = 12.92 * tr;
				if ( tg > 0.0031308 ) tg = 1.055 * ( pow(tg, ( 1.0 / 2.4)) ) - 0.055;
				else                  tg = 12.92 * tg;
				if ( tb > 0.0031308 ) tb = 1.055 * ( pow(tb, ( 1.0 / 2.4)) ) - 0.055;
				else                  tb = 12.92 * tb;
				
				R[i][j] = (T)tr;
				G[i][j] = (T)tg;
				B[i][j] = (T)tb;
			}
		}
	}
	
	// 0 ~ 255
	void lsRGBTosRGB_255(){
		for(int i=0; i<height; i++){
			for(int j=0; j<width; j++){
				double tr = (double)R[i][j] / 255.0;
				double tg = (double)G[i][j] / 255.0;
				double tb = (double)B[i][j] / 255.0;
				
				if ( tr > 0.0031308 ) tr = 1.055 * ( pow(tr, ( 1.0 / 2.4)) ) - 0.055;
				else                  tr = 12.92 * tr;
				if ( tg > 0.0031308 ) tg = 1.055 * ( pow(tg, ( 1.0 / 2.4)) ) - 0.055;
				else                  tg = 12.92 * tg;
				if ( tb > 0.0031308 ) tb = 1.055 * ( pow(tb, ( 1.0 / 2.4)) ) - 0.055;
				else                  tb = 12.92 * tb;
				
				R[i][j] = (T)(tr * 255.0);
				G[i][j] = (T)(tg * 255.0);
				B[i][j] = (T)(tb * 255.0);
			}
		}
	}
	
	void normalize(){
		double max = 0;
		for(int i=0; i<height; i++){
			for(int j=0; j<width; j++){
				if(R[i][j] > max) max = R[i][j];
				if(G[i][j] > max) max = G[i][j];
				if(B[i][j] > max) max = B[i][j];
			}
		}
		
		for(int i=0; i<height; i++){
			for(int j=0; j<width; j++){
				R[i][j] /= max;
				G[i][j] /= max;
				B[i][j] /= max;
			}
		}
	}
	
	// 0.0-1.0 -> 0-255
	Image<unsigned char> toByte(){
		Image<unsigned char> d(height, width, pixelFormat);
		if(pixelFormat){
			for(int i=0; i<height; i++){
				for(int j=0; j<width; j++){
					d.R[i][j] = (unsigned char)(255.0 * (R[i][j] > 1.0 ? 1.0 : R[i][j]));
					d.G[i][j] = (unsigned char)(255.0 * (G[i][j] > 1.0 ? 1.0 : G[i][j]));
					d.B[i][j] = (unsigned char)(255.0 * (B[i][j] > 1.0 ? 1.0 : B[i][j]));
				}
			}
		}
		else{
			for(int i=0; i<height; i++){
				for(int j=0; j<width; j++){
					d.G[i][j] = (unsigned char)(255.0 * (G[i][j] > 1.0 ? 1.0 : G[i][j]));
				}
			}
		}
		return d;
	}
	
	void writeMat(char *fileName){
		FILE *fp = fopen(fileName, "w");
		if( !fp){
			printf("can't open %s\n", fileName);
    		return ;
		}
		
		for(int i=0;i<height; i++){
			for(int j=0; j<height; j++){
				fprintf(fp, "%lg %lg %lg ", (double)R[i][j], (double)G[i][j], (double)B[i][j]);
			}
			fprintf(fp, "\n");
		}
		
		fclose(fp);
	}
	
	template<class RT>
	Image<RT>* blur(int s, Image<RT> *result=NULL){
		int n = s*2+1;
		if(result == NULL) result = new Image<RT>(height, width);
		Image<double> tmp(height, width);
		double **filter = newTex<double>(n, n);
		gaussianFilter(filter, n);
		
		double total = 0;
		for(int i=0; i<n; i++) total += filter[s][i];
		
		//v
		#pragma omp parallel for
		for(int i=0; i<height; i++){
			for(int j=0; j<width; j++){
				double r=0, g=0, b=0;
				for(int _i=0; _i<n; _i++){
					int id = i - s + _i;
					if(id < 0) id *= -1;
					else if(id >= height) id = height - (id - height) -1;
					r += R[id][j] * filter[s][_i];
					g += G[id][j] * filter[s][_i];
					b += B[id][j] * filter[s][_i];
				}
				tmp.R[i][j] = r/total;
				tmp.G[i][j] = g/total;
				tmp.B[i][j] = b/total;
			}
		}
		
		//h
		#pragma omp parallel for
		for(int j=0; j<width; j++){
			for(int i=0; i<height; i++){
				double r=0, g=0, b=0;
				for(int _j=0; _j<n; _j++){
					int id = j - s + _j;
					if(id < 0) id *= -1;
					else if(id >= width) id = width - (id - width) -1;
					r += tmp.R[i][id] * filter[_j][s];
					g += tmp.G[i][id] * filter[_j][s];
					b += tmp.B[i][id] * filter[_j][s];
				}
				result->R[i][j] = r/total;
				result->G[i][j] = g/total;
				result->B[i][j] = b/total;
			}
		}
		return result;
	}
	
	Image<unsigned char> tomemapping(){
		double a = 0.18;
		const int level = 7;
		int phi = 10;
		double xlon = 0.001;
		
		/*normalize
		double _max = 0, _min = 1000000.0, _interval;
		for (int i=0;i<height;i++){
			for (int j=0;j<width;j++){
				if(R[i][j] < _min) _min = R[i][j];
				if(G[i][j] < _min) _min = G[i][j];
				if(B[i][j] < _min) _min = B[i][j];
				if(R[i][j] > _max) _max = R[i][j];
				if(G[i][j] > _max) _max = G[i][j];
				if(B[i][j] > _max) _max = B[i][j];
			}
		}
		_interval = _max - _min;
		for (int i=0;i<height;i++){
			for (int j=0;j<width;j++){
				R[i][j] = (R[i][j] - _min) / _interval;
				G[i][j] = (G[i][j] - _min) / _interval;
				B[i][j] = (B[i][j] - _min) / _interval;
			}
		}//*/
		
		int img_h = height;
		int img_w = width;
		long double avgLw = 0;
		Image<unsigned char> result(img_h, img_w);
		
		double **Lw = newTex<double>(img_h, img_w);
		double **L = newTex<double>(img_h, img_w);
		
		for (int i=0;i<img_h;i++)
		{
			for (int j=0;j<img_w;j++)
			{
				Lw[i][j] = 0.2126 * R[i][j] + 0.7152 * G[i][j] + 0.0722 * B[i][j]; 
				if(Lw[i][j] < 0.0){
					//printf("LE zero!\n");
					Lw[i][j] = 0.0;
				}
				else if(!std::isfinite(Lw[i][j]) || Lw[i][j] > 1.0e3){
					//printf("nan!\n");
					Lw[i][j] = 1.0e3;
				}
				
				avgLw += std::log(Lw[i][j] + 0.0001);
			}
		}
	
		avgLw = std::exp(avgLw / (img_h*img_w));
		printf("avglw= %lg\n", (double)avgLw);
	
		for (int i=0;i<img_h;i++)
		{
			for (int j=0;j<img_w;j++)
			{
				L[i][j]=a*Lw[i][j] / avgLw;
			}
		}
	
		double **blurred[level+1];
		for(int s=0; s<=level; s++){
			blurred[s] = newTex<double>(img_h, img_w);
		}
		
		for(int s=0; s<=level; s++){
			const int n = s*2+1;
			double **tmp = newTex<double>(img_h, img_w);
			double **_filter = new double *[n];
			for (int i = 0; i < n; i++)
				_filter[i] = new double [n];

			double **filter = new double *[n];
			for(int i=0; i<n; i++) filter[i] = _filter[i];
			gaussianFilter(filter, n);
			double total = 0;
			for(int i=0; i<n; i++) total += filter[s][i];
			
			//v
			#pragma omp parallel for
			for(int i=0; i<img_h; i++){
				for(int j=0; j<img_w; j++){
					
					double pxsum = 0;
					for(int _i=0; _i<n; _i++){
						int id = i - s + _i;
						if(id < 0) id *= -1;
						else if(id >= img_h) id = img_h - (id - img_h) -1;
						pxsum += L[id][j] * filter[s][_i];
					}
					tmp[i][j] = pxsum/total;
				}
			}
			
			//h
			#pragma omp parallel for
			for(int j=0; j<img_w; j++){
				for(int i=0; i<img_h; i++){
					
					double pxsum = 0;
					for(int _j=0; _j<n; _j++){
						int id = j - s + _j;
						if(id < 0) id *= -1;
						else if(id >= img_w) id = img_w - (id - img_w) -1;
						pxsum += tmp[i][id] * filter[_j][s];
					}
					blurred[s][i][j] = pxsum/total;
				}
			}
			delTex<double>(tmp , img_h);
		}
		
		double two = 2;
		double tmp = pow(two,phi)*a;
	
		int **sk = newTex<int>(img_h, img_w);
		
		for (int i=0;i<img_h ;i++)
		{
			for (int j=0 ;j<img_w ; j++)
			{
				double dis = 1;
				sk[i][j] = -1;
				for (int s = 0 ;s<level -1 ;s++)
				{
					double s2 = (s+1) * (s+1);
					dis = (blurred[s][i][j]-blurred[s+1][i][j]) / (tmp/(s2)+blurred[s][i][j]);
					
					if(abs(dis)<xlon)
					{
						sk[i][j] = s;
						break;
					}
				}
				if (sk[i][j]<0)
				{
					sk[i][j]=level;
				}
			}
		}
	
		//double **LR = newTex<double>(img_h, img_w);
		for (int i=0;i<img_h;i++)
		{
			for (int j= 0;j<img_w;j++)
			{
				L[i][j] /= (1+blurred[ sk[i][j] ][i][j]);  //要寫的
			}
		}
		
		for (int i =0;i<img_h;i++)
		{
			for (int j=0;j<img_w;j++)
			{
				double r = R[i][j] / Lw[i][j] * L[i][j];
				double g = G[i][j] / Lw[i][j] * L[i][j];
				double b = B[i][j] / Lw[i][j] * L[i][j];
				//*
				if ( r > 0.0031308 ) r = 1.055 * ( pow(r, ( 1.0 / 2.4)) ) - 0.055;
				else                 r = 12.92 * r;
				if ( g > 0.0031308 ) g = 1.055 * ( pow(g, ( 1.0 / 2.4)) ) - 0.055;
				else                 g = 12.92 * g;
				if ( b > 0.0031308 ) b = 1.055 * ( pow(b, ( 1.0 / 2.4)) ) - 0.055;
				else                 b = 12.92 * b;//*/
				
				result.R[i][j] = (unsigned char)(r >= 1.0 ? 255 : r <= 0 ? 0 : r*255);
				result.G[i][j] = (unsigned char)(g >= 1.0 ? 255 : g <= 0 ? 0 : g*255);
				result.B[i][j] = (unsigned char)(b >= 1.0 ? 255 : b <= 0 ? 0 : b*255);
			}
		}
		
		//del
		for(int s=1; s<=level; s++){
			delTex<double>(blurred[s], img_h);
		}
		///↓要寫的
		//imwrite(result_L,'result_G.bmp');
		//figure, imshow(result_L);
		//*/
		
		delTex<double>(Lw, img_h);
		delTex<double>(L , img_h);
	
		return result;
	}
	
	Image<T>& operator= (const Image<T>& r){
		this->clear();
		
		if(r.name != NULL){
			int namel = strlen(r.name);
			this->name = new char[namel];
			for(int i=0; i<namel; i++){
				this->name[i] = r.name[i];
			}
		}
		else  this->name = NULL;
		
		this->init(r.height, r.width, r.pixelFormat);
		
		if(r.pixelFormat){
			for(int i=0; i<height; i++){
				for(int j=0; j<width; j++){
					this->R[i][j] = r.R[i][j];
					this->G[i][j] = r.G[i][j];
					this->B[i][j] = r.B[i][j];
				}
			}
		}
		else{
			for(int i=0; i<height; i++){
				for(int j=0; j<width; j++){
					this->G[i][j] = r.G[i][j];
				}
			}
		}
		return *this;
	}//*/
	
	
	private:
	char *getNameFromPath(const char *r){
		int i=0, idot=-1, ipwd, length;
		while(r[i] != 0) i++;
		length = i;
		i--;
		while(i >=0 && r[i] != '/' && r[i] != '\\'){
			if(idot < 0 && r[i]=='.') idot = i;
			i--;
		}
		ipwd = i;
		
		char *d = new char[length-idot];
		name = new char[idot-ipwd];
		int top=0;
		for(top = 0, i=idot+1; i < length; i++, top++){
			d[top] = r[i];
		}
		d[top] = 0;
		for(top = 0, i=ipwd+1; i < idot; i++, top++){
			name[top] = r[i];
		}
		name[top] = 0;
		return d;
	}
	
	bool ReadPngData( FILE* rFile){
	    png_structp png_ptr  = NULL;
		png_infop   info_ptr = NULL;
	
		//init png read
		try{
		    png_ptr  = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);        //創建 png 讀取結構
		    info_ptr = png_create_info_struct(png_ptr);                               //png 文件信息結構
		    png_init_io(png_ptr, rFile);
		    png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_EXPAND | PNG_TRANSFORM_STRIP_ALPHA, 0); //讀取png，忽略A
		}
		catch(...){
			return false;
		}
	    
	    //read data
	    width = png_get_image_width( png_ptr, info_ptr );         //獲得圖片寬度
	    height = png_get_image_height( png_ptr, info_ptr );       //獲得圖片高度
	    int color_type = png_get_color_type( png_ptr, info_ptr ); //獲得圖片顏色類型
	    int img_size = width * height;                            //計算 pixel 數 
	    png_bytep *row_point = png_get_rows( png_ptr, info_ptr ); //讀取 pixel
	    pixelFormat = (color_type%4 > 0 ? 1 : 0);                 //彩色1，灰階0。忽略 type 3 
	
	    //將讀取到的數字按規定格式讀取 
	    if(pixelFormat){
	    	//RGB
		    R = new T*[height];
		    G = new T*[height];
		    B = new T*[height];
		    for(int i=0; i<height; i++){
		    	R[i] = new T[width];
		    	G[i] = new T[width];
		    	B[i] = new T[width];
		    }
		    
		    for( int x = 0; x < height; x++ ){
		        for( int y = 0; y < width; y++ ){
		            R[x][y] = (T)row_point[x][y*3 + 0];    //R
		            G[x][y] = (T)row_point[x][y*3 + 1];    //G
		            B[x][y] = (T)row_point[x][y*3 + 2];    //B
		        }
		    }
	    }
	    else{
	    	//Gray
	    	R = NULL;
	    	G = new T*[height];
	    	B = NULL;
	    	for(int i=0; i<height; i++){
		    	G[i] = new T[width];
		    }
		    
		    for( int x = 0; x < height; x++ ){
		        for( int y = 0; y < width; y++ ){
		            G[x][y] = (T)row_point[x][y];        //G
		        }
		    }
	    }
	    
	    png_destroy_read_struct(&png_ptr, &info_ptr, 0);
	    
	    return true;
	}
	
	bool ReadJpgData( FILE* rFile){
		struct jpeg_decompress_struct cinfo;
		struct jpeg_error_mgr jerr;
		
		cinfo.err = jpeg_std_error(&jerr);
		jpeg_create_decompress(&cinfo);
		
		//read info
		jpeg_stdio_src(&cinfo, rFile);
		jpeg_read_header(&cinfo, TRUE);
		height = cinfo.image_height;
		width  = cinfo.image_width;
		int pixelBytes = cinfo.num_components;
		pixelFormat = (pixelBytes == 1 ? 0 : 1);
		
		unsigned char *m_bgra = new unsigned char[width*height*pixelBytes];
		memset(m_bgra, 0, width*height*pixelBytes);
		
		jpeg_start_decompress(&cinfo);
		
		JSAMPROW row_pointer[1];
		
		//read data
		while(cinfo.output_scanline < height){
			/*row_pointer[0] = &m_bgr[(cinfo.output_height-cinfo.output_scanline-1)
			 *cinfo.image_width*cinfo.num_components];*/
			row_pointer[0] = &m_bgra[cinfo.output_scanline*width*pixelBytes];
			//printf("%d ", cinfo.output_scanline);
			jpeg_read_scanlines(&cinfo, row_pointer, 1);
		}
		if(pixelFormat){
	    	//RGB
		    R = new T*[height];
		    G = new T*[height];
		    B = new T*[height];
		    for(int i=0; i<height; i++){
		    	R[i] = new T[width];
		    	G[i] = new T[width];
		    	B[i] = new T[width];
		    }
		    
		    for( int x = 0; x < height; x++ ){
		        for( int y = 0; y < width; y++ ){
		            R[x][y] = (T)m_bgra[x*width*3 + y*3 + 0];    //R
		            G[x][y] = (T)m_bgra[x*width*3 + y*3 + 1];    //G
		            B[x][y] = (T)m_bgra[x*width*3 + y*3 + 2];    //B
		        }
		    }
		}
		else{
			//Gray
	    	R = NULL;
	    	G = new T*[height];
	    	B = NULL;
	    	for(int i=0; i<height; i++){
		    	G[i] = new T[width];
		    }
		    
		    for( int x = 0; x < height; x++ ){
		        for( int y = 0; y < width; y++ ){
		            G[x][y] = (T)m_bgra[x*width + y];    //G
		        }
		    }
		}
		
		jpeg_finish_decompress(&cinfo);
		jpeg_destroy_decompress(&cinfo);
		
		delete[] m_bgra;
		return true;
	}
	
	bool WriteJpg( FILE *cFile){
		if(cFile == NULL) return false;
		
		struct jpeg_compress_struct cinfo;
		struct jpeg_error_mgr jerr; 
	    
	    cinfo.err = jpeg_std_error(&jerr); 
	    jpeg_create_compress(&cinfo);
	
		jpeg_stdio_dest(&cinfo, cFile);
	
	    cinfo.image_width = width;
	    cinfo.image_height = height;
	    if(pixelFormat){
	    	cinfo.input_components = 3;
	    	cinfo.in_color_space = JCS_RGB;
	    }
	    else{
	    	cinfo.input_components = 1;
	    	cinfo.in_color_space = JCS_GRAYSCALE;
	    }
		jpeg_set_defaults(&cinfo);
	
	    jpeg_start_compress(&cinfo, FALSE);
		
	    JSAMPROW *row_pointer = new JSAMPROW[height];
	    unsigned char *stride;
	    if(pixelFormat){    //RGB
			stride = (unsigned char *)malloc( height * width * 3);
		    for (int i=0; i<height; i++) {
		    	unsigned char *nrow = &stride[i * width * 3];
		    	row_pointer[i] = &stride[i * width * 3];
				for(int j=0; j<width; j++){
					nrow[j*3  ] = R[i][j];
					nrow[j*3+1] = G[i][j];
					nrow[j*3+2] = B[i][j];
				}
			}
	    }
		else{               //G
			stride = (unsigned char *)malloc( height * width);
		    for (int i=0; i<height; i++) {
		    	unsigned char *nrow = &stride[i * width];
		    	row_pointer[i] = &stride[i * width];
				for(int j=0; j<width; j++){
					nrow[j] = G[i][j];
				}
			}
		}
	
	    int d = jpeg_write_scanlines(&cinfo, row_pointer, height);
	
	    jpeg_finish_compress(&cinfo);
	    jpeg_destroy_compress(&cinfo);
	    
	    fclose(cFile);
	    
	    delete[] row_pointer;
	    delete[] stride;
	
	    return d == height;
	}
	
	bool WriteHDR(FILE *cFile){
		if(cFile == NULL) return false;
		
		// Write trivial hdr header
		fprintf(cFile,"#?RADIANCE\n");
		fprintf(cFile,"FORMAT=32-bit_rle_rgbe\n");
		fprintf(cFile,"\n");
		fprintf(cFile,"-Y %d +X %d\n", height, width);
		
		for (int y = 0; y < height; y++){
			FileHDR::fwritescan<T>(R[y], G[y], B[y], width, cFile);
		}
		
		fclose(cFile);

		return true;
	}
	
	double ComputeGaussian(double n, double theta){
	    return ((1.0 / sqrt(2 * 3.1415926536 * theta)) * exp(-(n * n) / (2 * theta * theta)));
	}
	
	unsigned long long int Pascal_triangle(unsigned long long int n, unsigned long long int r){
		unsigned long long int x = 1; // N階
		unsigned long long int y = 1; // N-R階
		unsigned long long int z = 1; // R階
		unsigned long long int ans = 0; //結果
		
		if ( n==0 || r==n ) return 1;
		else
		{
			for (unsigned long long int i = 1; i<n + 1; i++)
			{
				x =x* i;
			}
			for (unsigned long long int j = 1; j<n - r + 1; j++)
			{
				y=y*j;
			}
			for (unsigned long long int o = 1; o<r + 1; o++)
			{
				z=z*o;
			}
			ans= x / (y*z);
			return ans;
		}
	}
	
	void gaussianFilter (double** filter, int N) {
		//int sumN=0;
		for(int i=0;i<N;i++)
		{
			unsigned long long int v = Pascal_triangle(N - 1, i);
			filter[0][i] = (double)v;
			//sumN+=v;
		}
		for(int j=1;j<N;j++)
		{
			unsigned long long int v = Pascal_triangle(N - 1, j);
			filter[j][0] = (double)v;
			//sumN+=v;
		}
		for(int k=1;k<N;k++)
		{
			for(int l=1;l<N;l++)
			{
				filter[k][l]=filter[k][0]*filter[0][l];
				//sumN += (int)filter[k][l];
			}
		}
		/*
		for(int k=0;k<N;k++)
		{
			for(int l=0;l<N;l++)
			{
				filter[k][l]=filter[k][l] / sumN;
			}
		}//*/
	}
};

