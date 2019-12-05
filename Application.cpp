#define _CRT_SECURE_NO_DEPRECATE
#include "Application.h"
#include <FL/fl_file_chooser.H>
#include "Gui.h"
#include <iostream>
#include <list>
#include <algorithm>
#include "Matrix.h"

extern Gui *gui;
Image curImage; // make inImage a global, need it in other classes
Volume vol;

std::list <Image> ImageHistoryList; 


// these are declared in transFunc.cpp
extern TransferFunction transFunc[4];
extern int maxTransFunc;              
extern int activeTransFunc;       
extern float transFuncColor[4][3]; 
int histogram[3][256];

Matrix rotX(3, 3), rotY(3, 3), rotZ(3, 3);
int volumeSwitch = 0;  //0: xray; 1: mip;
int imageHeight, imageWidth;

// help functions

int getPixel(int i, int j) {
	return (i * curImage.nx + j) * curImage.ncolorChannels;
}

double getDensity(int i, int j, int k) {
	return vol.data[i * (vol.ny * vol.nz) + j * (vol.nz) + k];
}

void loadRotationMatrix() {

	double radx = gui->sliderX->value();
	double rady = gui->sliderY->value();
	double radz = gui->sliderZ->value();
	/*
	printf("radx: %f\n", radx);
	printf("rady: %f\n", rady);
	printf("radz: %f\n", radz);
	*/
	radx = (radx * acos(-1.0)) / 180.0;
	rady = (rady * acos(-1.0)) / 180.0;
	radz = (radz * acos(-1.0)) / 180.0;
	/*
	printf("radx: %f\n", radx);
	printf("rady: %f\n", rady);
	printf("radz: %f\n", radz);
	*/
	//system("pause");
	rotX.data[rotX.getIndex(0, 0)] = 1;
	rotX.data[rotX.getIndex(0, 1)] = 0;
	rotX.data[rotX.getIndex(0, 2)] = 0;

	rotX.data[rotX.getIndex(1, 0)] = 0;
	rotX.data[rotX.getIndex(1, 1)] = cos(radx);
	rotX.data[rotX.getIndex(1, 2)] = -sin(radx);

	rotX.data[rotX.getIndex(2, 0)] = 0;
	rotX.data[rotX.getIndex(2, 1)] = sin(radx);
	rotX.data[rotX.getIndex(2, 2)] = cos(radx);

	rotY.data[rotY.getIndex(0, 0)] = cos(rady);
	rotY.data[rotY.getIndex(0, 1)] = 0;
	rotY.data[rotY.getIndex(0, 2)] = sin(rady);

	rotY.data[rotY.getIndex(1, 0)] = 0;
	rotY.data[rotY.getIndex(1, 1)] = 1;
	rotY.data[rotY.getIndex(1, 2)] = 0;

	rotY.data[rotY.getIndex(2, 0)] = -sin(rady);
	rotY.data[rotY.getIndex(2, 1)] = 0;
	rotY.data[rotY.getIndex(2, 2)] = cos(rady);

	rotZ.data[rotZ.getIndex(0, 0)] = cos(radz);
	rotZ.data[rotZ.getIndex(0, 1)] = -sin(radz);
	rotZ.data[rotZ.getIndex(0, 2)] = 0;

	rotZ.data[rotZ.getIndex(1, 0)] = sin(radz);
	rotZ.data[rotZ.getIndex(1, 1)] = cos(radz);
	rotZ.data[rotZ.getIndex(1, 2)] = 0;

	rotZ.data[rotZ.getIndex(2, 0)] = 0;
	rotZ.data[rotZ.getIndex(2, 1)] = 0;
	rotZ.data[rotZ.getIndex(2, 2)] = 1;
}

void addImageHistory() {
	Image tempImage;
	tempImage.data = new unsigned char[curImage.n * curImage.ncolorChannels];
	tempImage.n = curImage.n;
	tempImage.ncolorChannels = curImage.ncolorChannels;
	tempImage.nx = curImage.nx;
	tempImage.ny = curImage.ny;

	for (int i = 0; i < curImage.n * curImage.ncolorChannels; i++) {
		tempImage.data[i] = curImage.data[i];
	}

	ImageHistoryList.push_back(tempImage);
}

// the constructor method for the Application class

Application::Application()
{
  // initialize the image data structure
  curImage.nx=curImage.ny=curImage.n=curImage.ncolorChannels=0;

  // add more initialization here:

}

// the method that gets executed when the readFile callback is called

void Application::ReadFile()
{
   FILE *fp;
   char imageType[3],str[100];
   int dummy;
   int i,j;

   char *file = fl_file_chooser("Pick a file to READ from", "*.{pgm,ppm,vol}", "");
   if(file == NULL)
		return;

   // Read PGM image file with filename "file"

   // The PGM file format for a GREYLEVEL image is:
   // P5 (2 ASCII characters) <CR>
   // two ASCII numbers for nx ny (number of rows and columns <CR>
   // 255 (ASCII number) <CR>
   // binary pixel data with one byte per pixel

   // The PGM file format for a COLOR image is:
   // P6 (2 ASCII characters) <CR>
   // two ASCII numbers for nx ny (number of rows and columns <CR>
   // 255 (ASCII number) <CR>
   // binary pixel data with three bytes per pixel (one byte for eacg RGB)

    fp=fopen(file,"rb");

    // free memory from old image, if any
    if(curImage.n>0)
    {
	  delete[] curImage.data;
	}

    // read the first ASCII line to find out if we read a color image or
	// a greylevel image

    fgets(str,100,fp);
	sscanf(str,"%s",imageType);

	//OutputDebugStringA(str);

    if(!strncmp(imageType,"P5",2)) // greylevel image 
    {
       curImage.ncolorChannels=1;
       maxTransFunc=1;  // have read a greylevel image: need only one transfer function
       activeTransFunc=0;    // transFunc 0 is active (all the time)
       transFuncColor[0][0]=1.0; transFuncColor[0][1]=0.0; transFuncColor[0][2]=0.0;  // draw transfer function curve in red 
	} 
	else if(!strncmp(imageType,"P6",2)) // color image 
    {
	   curImage.ncolorChannels=3;
       maxTransFunc=1;  // have read a color image: if you want 3 transfer functions (HSV) change this value to 3
       activeTransFunc=0;    // transFunc 0 is active (for now)
       transFuncColor[0][0]=1.0; transFuncColor[0][1]=0.0; transFuncColor[0][2]=0.0;  // draw transfer function curve in red 
       transFuncColor[1][0]=0.0; transFuncColor[1][1]=1.0; transFuncColor[1][2]=0.0;  // draw transfer function curve in red 
       transFuncColor[2][0]=0.0; transFuncColor[2][1]=0.0; transFuncColor[2][2]=1.0;  // draw transfer function curve in red 
	}

	// skip comments embedded in header
    fgets(str,100,fp);  
	while(str[0]=='#')
		fgets(str,100,fp);

    // read image dimensions 
    sscanf(str,"%d %d",&curImage.nx,&curImage.ny);

	// read the next line into dummy variable
    fgets(str,100,fp);  
	 
   	// allocate the memory
	curImage.n=curImage.nx*curImage.ny;

    // read the image data 
	curImage.data = new unsigned char [curImage.n*curImage.ncolorChannels];
	fread(curImage.data,sizeof(unsigned char),curImage.n*curImage.ncolorChannels,fp);
	
  	// unfortunately OpenGL displays the image upside-down
 	// we have to flip it at read time.
    FlipImage(&curImage);

	fclose(fp);

    // call the window drawing routines to display the image curImage
	// and to draw image-related transFuncs

	Caluate_Histogram();
	addImageHistory();
	//---------------------------------------

	gui->DisplayWindow->redraw();
	gui->EditorWindow->redraw();
}


void Application::WriteFile()
{
   FILE *fp;
   char imageType[3],str[100];
   int dummy,i;
   char *file = fl_file_chooser("Specify a filename to WRITE to", "*.{vol,ppm,pgm}", "");
   if(file == NULL)
		return;

   // Write PGM image file with filename "file"

   // The PGM file format for a GREYLEVEL image is:
   // P5 (2 ASCII characters) <CR>
   // two ASCII numbers for nx ny (number of rows and columns <CR>
   // 255 (ASCII number) <CR>
   // binary pixel data with one byte per pixel

   // The PGM file format for a COLOR image is:
   // P6 (2 ASCII characters) <CR>
   // two ASCII numbers for nx ny (number of rows and columns <CR>
   // 255 (ASCII number) <CR>
   // binary pixel data with three bytes per pixel (one byte for each R,G,B)

    fp=fopen(file,"wb");

    // write the first ASCII line with the file type
	if(curImage.ncolorChannels==1)
   	  fprintf(fp,"P5\n"); //greylevel image
    else if(curImage.ncolorChannels==3)
      fprintf(fp,"P6\n");  // color image

    // write image dimensions 
    fprintf(fp,"%d %d\n",curImage.nx,curImage.ny);

	// write '255' to the next line 
    fprintf(fp,"255\n");

	// since we flipped the image upside-down when we read it
 	// we have to write it upside-down so it's stored the right way
    for(i=curImage.ny-1;i>=0;i--)
  	  fwrite(&curImage.data[i*curImage.nx*curImage.ncolorChannels],sizeof(unsigned char),curImage.nx*curImage.ncolorChannels,fp);

	fclose(fp);
}

// flips an image upside down
// you will not have to change anything here

void Application::FlipImage(Image *img)
{
    int i,j,k,rowOffsetSrc,rowOffsetDest,columnOffset;
    unsigned char ctmp;

	for(i=0;i<img->ny/2;i++)
	{
	   rowOffsetSrc=i*img->nx*img->ncolorChannels;
	   rowOffsetDest=(img->ny-1-i)*img->nx*img->ncolorChannels;
       for(j=0;j<img->nx;j++)
	   {
		   columnOffset=j*img->ncolorChannels;
		   for(k=0;k<img->ncolorChannels;k++)
		   {
			   ctmp=img->data[rowOffsetSrc+columnOffset+k];
               img->data[rowOffsetSrc+columnOffset+k]=img->data[rowOffsetDest+columnOffset+k];
               img->data[rowOffsetDest+columnOffset+k]=ctmp;
		   }
	   }
	}
}



// put your application routines here:

void Application::GUI_Update(){
	for (int i = 0; i < ImageHistoryList.back().n * ImageHistoryList.back().ncolorChannels; i++) {
		curImage.data[i] = transFunc[activeTransFunc][ImageHistoryList.back().data[i]];
	}

	Caluate_Histogram();
	addImageHistory();

	gui->DisplayWindow->redraw();
	gui->EditorWindow->redraw();
}
void Application::GUI_GreyEffect() {
	for (int y = 0; y < ImageHistoryList.back().ny; y++) {
		for (int x = 0; x < ImageHistoryList.back().nx; x++) {
			int greyVal = 0;
			for (int c = 0; c < ImageHistoryList.back().ncolorChannels; c++) {
				greyVal = ImageHistoryList.back().data[getPixel(y, x) + c] + greyVal;
			}
			greyVal = greyVal / 3;
			for (int c = 0; c < ImageHistoryList.back().ncolorChannels; c++) {
				curImage.data[getPixel(y, x) + c] = greyVal;
			}
		}
	}

	Caluate_Histogram();
	addImageHistory();

	gui->DisplayWindow->redraw();
	gui->EditorWindow->redraw();
}

void Application::GUI_AverageSmooth(){
	for (int y = 0; y < ImageHistoryList.back().ny; y++) {
		for (int x = 0; x < ImageHistoryList.back().nx; x++) {
			for (int c = 0; c < ImageHistoryList.back().ncolorChannels; c++) { //loop though RGB
				int sum = 0;
				for (int ty = -1; ty < 2; ty++) {
					for (int tx = -1; tx < 2; tx++) {
						// the case on the edge
						int tempy = y + ty;
						int tempx = x + tx;
						if (tempy < 0) tempy = ImageHistoryList.back().ny + ty;
						if (tempx < 0) tempx = ImageHistoryList.back().nx + tx;

						//sum up R || G || B splitly
						sum = ImageHistoryList.back().data[getPixel(tempy, tempx) + c] + sum;
					}
				}
				sum = sum / 9;
				curImage.data[getPixel(y, x) + c] = sum;
			}
		}
	}

	Caluate_Histogram();
	addImageHistory();

	gui->DisplayWindow->redraw();
	gui->EditorWindow->redraw();
}

void Application::GUI_MedianSmooth(){
	for (int y = 0; y < ImageHistoryList.back().ny; y++) {
		for (int x = 0; x < ImageHistoryList.back().nx; x++) {
			for (int c = 0; c < ImageHistoryList.back().ncolorChannels; c++) { //loop though RGB
				int sample[9];
				int counter = 0;
				int median = -1;
				for (int ty = -1; ty < 2; ty++) {
					for (int tx = -1; tx < 2; tx++) {
						// the case on the edge
						int tempy = y + ty;
						int tempx = x + tx;
						if (tempy < 0) tempy = ImageHistoryList.back().ny + ty;
						if (tempx < 0) tempx = ImageHistoryList.back().nx + tx;

						//add up R || G || B splitly to array
						sample[counter] = ImageHistoryList.back().data[getPixel(tempy, tempx) + c];
						counter++;
					}
				}
				int n = sizeof(sample) / sizeof(sample[0]);
				std::sort(sample, sample + n);
				if (n % 2 != 0)
					median = sample[n / 2];
				else
					median = (sample[(n - 1) / 2] + sample[n / 2]) / 2.0;

				curImage.data[getPixel(y, x) + c] = median;
			}
		}
	}

	Caluate_Histogram();
	addImageHistory();

	gui->DisplayWindow->redraw();
	gui->EditorWindow->redraw();
}

void Application::GUI_GaussianSmooth(){
	// init gaussian filliter
	double sum = 0;
	double gaussianFilliter[5][5] = {
		{1, 4, 7, 4, 1},
		{4, 16, 26, 16, 4},
		{7, 26, 41, 26, 7},
		{4, 16, 26, 16, 4},
		{1, 4, 7, 4, 1}
	};

	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			sum += gaussianFilliter[i][j];
		}
	}

	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			gaussianFilliter[i][j] /= sum;
		}
	}

	//apply gaussian filliter
	for (int y = 0; y < ImageHistoryList.back().ny; y++) {
		for (int x = 0; x < ImageHistoryList.back().nx; x++) {
			for (int c = 0; c < ImageHistoryList.back().ncolorChannels; c++) { //loop though RGB
				int sum = 0;
				for (int ty = -2; ty < 3; ty++) {
					for (int tx = -2; tx < 3; tx++) {
						// the case on the edge
						int tempy = y + ty;
						int tempx = x + tx;
						if (tempy < 0) tempy = ImageHistoryList.back().ny + ty;
						if (tempx < 0) tempx = ImageHistoryList.back().nx + tx;

						//sum up R || G || B splitly
						sum = ImageHistoryList.back().data[getPixel(tempy, tempx) + c] * gaussianFilliter[ty +2][tx + 2] + sum;

					}
				}
				curImage.data[getPixel(y, x) + c] = sum;
			}
		}
	}

	Caluate_Histogram();
	addImageHistory();

	gui->DisplayWindow->redraw();
	gui->EditorWindow->redraw();

}

void Application::GUI_EdgeDetect(){
	// init sobel filliter
	Image tempImage;
	tempImage.data = new unsigned char[curImage.n * curImage.ncolorChannels];
	tempImage.n = curImage.n;
	tempImage.ncolorChannels = curImage.ncolorChannels;
	tempImage.nx = curImage.nx;
	tempImage.ny = curImage.ny;

	double Sobelx[3][3] = {
		{-1, 0, 1},
		{-2, 0, 2},
		{-1, 0, 1},
	};
	double Sobely[3][3] = {
		{1, 2, 1},
		{0, 0, 0},
		{-1, 2, -1},
	};

	for (int i = 0; i < ImageHistoryList.back().n * ImageHistoryList.back().ncolorChannels; i++) {
		tempImage.data[i] = ImageHistoryList.back().data[i];
	}

	//----------------------------------------//
	// uncomment for grey scale edge detection//
	//----------------------------------------//
	/*for (int y = 0; y < ImageHistoryList.back().ny; y++) {
		for (int x = 0; x < ImageHistoryList.back().nx; x++) {
			int greyVal = 0;
			for (int c = 0; c < ImageHistoryList.back().ncolorChannels; c++) {
				greyVal = ImageHistoryList.back().data[getPixel(y, x) + c] + greyVal;
			}
			greyVal = greyVal / 3;
			for (int c = 0; c < ImageHistoryList.back().ncolorChannels; c++) {
				tempImage.data[getPixel(y, x) + c] = greyVal;
			}
		}
	}*/
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	//apply sobel filliter
	for (int y = 0; y < tempImage.ny; y++) {
		for (int x = 0; x < tempImage.nx; x++) {
			for (int c = 0; c < tempImage.ncolorChannels; c++) { //loop though RGB
				int sumx = 0;
				int sumy = 0;
				for (int ty = -1; ty < 2; ty++) {
					for (int tx = -1; tx < 2; tx++) {
						// the case on the edge
						int tempy = y + ty;
						int tempx = x + tx;
						if (tempy < 0) tempy = tempImage.ny + ty;
						if (tempx < 0) tempx = tempImage.nx + tx;

						//sum up R || G || B splitly
						sumx = tempImage.data[getPixel(tempy, tempx) + c] * Sobelx[ty + 1][tx + 1] + sumx;
						sumy = tempImage.data[getPixel(tempy, tempx) + c] * Sobely[ty + 1][tx + 1] + sumy;
					}
				}
				if (sumx < 0) sumx = 0;
				if (sumy < 0) sumy = 0;

				curImage.data[getPixel(y, x) + c] = sqrt(pow(min(sumx,255), 2.0) + pow(min(sumy, 255), 2.0));
			}
		}
	}

	Caluate_Histogram();
	addImageHistory();

	gui->DisplayWindow->redraw();
	gui->EditorWindow->redraw();
}

void Application::GUI_UndoButton(){
	if (ImageHistoryList.size() > 1) {
		ImageHistoryList.pop_back();
		for (int i = 0; i < ImageHistoryList.back().n * ImageHistoryList.back().ncolorChannels; i++) {
			curImage.data[i] = ImageHistoryList.back().data[i];
		}

		Caluate_Histogram();

		gui->DisplayWindow->redraw();
		gui->EditorWindow->redraw();
	}
	else {
		std::cout << "you already reach the last layer of your image history" <<std::endl;
	}
}

void Application::GUI_RestoreButton(){
	for (int i = 0; i < ImageHistoryList.front().n * ImageHistoryList.front().ncolorChannels; i++) {
		curImage.data[i] = ImageHistoryList.front().data[i];
	}

	Caluate_Histogram();
	addImageHistory();

	gui->DisplayWindow->redraw();
	gui->EditorWindow->redraw();
}

void Application::Caluate_Histogram() {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 256; j++) {
			histogram[i][j] = 0;
		}
	}
	for (int y = 0; y < curImage.ny; y++) {
		for (int x = 0; x < curImage.nx; x++) {
			for (int c = 0; c < curImage.ncolorChannels; c++) { //loop though RGB
				histogram[c][curImage.data[getPixel(y, x) + c]] += 1;
			}
		}
	}

	/*for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 256; j++) {
			std::cout << histogram[i][j] << " ";
		}
		std::cout << std::endl;
	}*/
}

//------------------------------------------------------------
//--------------------Volume Rendering------------------------
//------------------------------------------------------------

void Application::GUI_ReadVolumeFile() {
	FILE* fp;
	char imageType[20], str[100];

char* file = fl_file_chooser("Pick a file to READ from", "*.{pgm,ppm,vol}", "");
if (file == NULL)
return;

fp = fopen(file, "rb");
fgets(str, 100, fp);
sscanf(str, "%s", imageType);
if (!strncmp(imageType, "P7", 2)) { // volume data
	// skip comments embedded in header
	fgets(str, 100, fp);
	while (str[0] == '#')
		fgets(str, 100, fp);

	// read volume dimensions 
	sscanf(str, "%d %d %d", &vol.nx, &vol.ny, &vol.nz);
	vol.n = vol.nx * vol.ny * vol.nz;

	fgets(str, 100, fp);

	vol.data = (unsigned char*)malloc(vol.n * sizeof(unsigned char));
	fread(vol.data, sizeof(unsigned char), vol.n, fp);
}
fclose(fp);

std::cout << "welcome to the club";
std::cout << vol.data;

GUI_VolumeDraw();
};

//void vIntersectRaywithVolumeBoundingBox(double p[3], double ray[3], double& t_front, double& t_back)
//{
//	double bbx0, bbx1, bby0, bby1, bbz0, bbz1;	// bounding box
//	double dev;
//	double t1, t2, mint, maxt;
//
//	t_front = -1e10;
//	t_back = 1e10;
//	dev = 0.0001;
//	if (useboundingbox) {
//		//		bbx0 = bbx[0] + dev;	bbx1 = bbx[1] - dev;
//		//		bby0 = bby[0] + dev;	bby1 = bby[1] - dev;
//		//		bbz0 = bbz[0] + dev;	bbz1 = bbz[1] - dev;
//		bbx0 = bbx[0] - dev;	bbx1 = bbx[1] + dev;
//		bby0 = bby[0] - dev;	bby1 = bby[1] + dev;
//		bbz0 = bbz[0] - dev;	bbz1 = bbz[1] + dev;
//	}
//	else {
//		bbx0 = 0 + dev;		bbx1 = vol.nx - 1 - dev;
//		bby0 = 0 + dev;		bby1 = vol.ny - 1 - dev;
//		bbz0 = 0 + dev;		bbz1 = vol.nz - 1 - dev;
//	}
//
//	if (fabs(ray[0]) <= 1e-5) {
//		if (p[0]<bbx0 || p[0]>bbx1)
//			return;
//	}
//	if (fabs(ray[1]) <= 1e-5) {
//		if (p[1]<bby0 || p[1]>bby1)
//			return;
//	}
//	if (fabs(ray[2]) <= 1e-5) {
//		if (p[2]<bbz0 || p[2]>bbz1)
//			return;
//	}
//	if (fabs(ray[0]) > 1e-5) {	// intersect with bbx
//		t1 = (bbx0 - p[0]) / ray[0];
//		t2 = (bbx1 - p[0]) / ray[0];
//		mint = min(t1, t2);
//		maxt = max(t1, t2);
//		t_front = max(mint, t_front);
//		t_back = min(maxt, t_back);
//	}
//	if (fabs(ray[1]) > 1e-5) {	// intersect with bby
//		t1 = (bby0 - p[1]) / ray[1];
//		t2 = (bby1 - p[1]) / ray[1];
//		mint = min(t1, t2);
//		maxt = max(t1, t2);
//		t_front = max(mint, t_front);
//		t_back = min(maxt, t_back);
//	}
//	if (fabs(ray[2]) > 1e-5) {	// intersecgt with bbz
//		t1 = (bbz0 - p[2]) / ray[2];
//		t2 = (bbz1 - p[2]) / ray[2];
//		mint = min(t1, t2);
//		maxt = max(t1, t2);
//		t_front = max(mint, t_front);
//		t_back = min(maxt, t_back);
//	}
//	return;
//}

void Application::GUI_xRay() {
	volumeSwitch = 0;
	GUI_VolumeDraw();
};
void Application::GUI_MIP() {
	volumeSwitch = 1;
	GUI_VolumeDraw();
};

void Application::GUI_VolumeDraw() {
	imageHeight = gui->sliderWidth->value();
	imageWidth = gui->sliderHeight->value();

	double s[] = { 0,0,2 };
	double p00[3];
	double ray[] = { 0,0,-1 };
	double u[] = { 1,0,0 };
	double v[] = { 0,1,0 };

	loadRotationMatrix();
	Matrix mu(1, 3), mv(1, 3), mray(1, 3), ms(1, 3);
	for (int i = 0; i < 3; i++) {
		mu.data[mu.getIndex(0, i)] = u[i];
		mv.data[mv.getIndex(0, i)] = v[i];
		mray.data[mray.getIndex(0, i)] = ray[i];
		ms.data[ms.getIndex(0, i)] = s[i];
	}
	mu = ((mu.multiply(rotX)).multiply(rotY)).multiply(rotZ);
	mv = ((mv.multiply(rotY)).multiply(rotZ)).multiply(rotX);
	mray = ((mray.multiply(rotX)).multiply(rotY)).multiply(rotZ);
	ms = ((ms.multiply(rotX)).multiply(rotY)).multiply(rotZ);
	double uv[3];

	for (int i = 0; i < 3; i++) {
		u[i] = mu.data[mu.getIndex(0, i)];
		v[i] = mv.data[mv.getIndex(0, i)];
		ray[i] = mray.data[mray.getIndex(0, i)];
		s[i] = ms.data[ms.getIndex(0, i)];
		uv[i] = -(u[i] + v[i]);
	}

	for (int i = 0; i < 3; i++) {
		p00[i] = s[i] - (((double)(imageWidth / 2.0) / (double)(imageWidth))) * u[i];
	}
	for (int i = 0; i < 3; i++) {
		p00[i] = p00[i] - (((double)(imageHeight / 2.0) / (double)(imageHeight))) * v[i];
	}


	delete[] curImage.data;

	curImage.ncolorChannels = 1;
	curImage.nx = imageHeight;
	curImage.ny = imageWidth;
	curImage.n = curImage.nx * curImage.ny;
	curImage.data = new unsigned char[curImage.ncolorChannels * curImage.n];

	for (int i = 0; i < curImage.nx; i++) {
		for (int j = 0; j < curImage.ny; j++) {
			double pij[3];
			double tempu[3];
			double tempv[3];
			double temp[3];
			for (int k = 0; k < 3; k++) {
				tempu[k] = ((double)(i) / imageHeight) * u[k];
			}

			for (int k = 0; k < 3; k++) {
				tempv[k] = ((double)(j) / imageWidth) * v[k];
			}

			for (int k = 0; k < 3; k++) {
				temp[k] = tempu[k] + tempv[k];;
			}

			for (int k = 0; k < 3; k++) {
				pij[k] = p00[k] + temp[k];
			}


			double tbegin, tend;
			bool f = false;

			double sum = 0;
			double maxval = 0;
			double step = 0.1;

			for (double t = 0; t <= 50; t += step) {
				double sampleloc[3];
				double lx, ly, lz, p;
				double w, v, u;
				double fdens;

				sampleloc[0] = pij[0] + t * ray[0];
				sampleloc[1] = pij[1] + t * ray[1];
				sampleloc[2] = pij[2] + t * ray[2];

				lx = (1.00 / (vol.nx - 1));
				ly = (1.00 / (vol.ny - 1));
				lz = (1.00 / (vol.nz - 1));
				p = ((sampleloc[0] + 0.5) / lx) * lx;
				int i = p * (vol.nx - 1);
				p = ((sampleloc[1] + 0.5) / ly) * ly;
				int j = p * (vol.ny - 1);
				p = ((sampleloc[2] + 0.5) / lz) * lz;
				int k = p * (vol.nz - 1);
				if (i < 0 || j < 0 || k < 0) {
					if (f)	break;
					continue;
				}
				if (i >= vol.nx - 1 || j >= vol.ny - 1 || k >= vol.nz - 1) {
					if (f)	break;
					continue;
				}
				w = sampleloc[0] - (i * lx);
				v = sampleloc[1] - (j * ly);
				u = sampleloc[2] - (k * lz);
				f = true;
				fdens =
					getDensity(i, j, k) * (1 - u) * (1 - v) * (1 - w) +
					getDensity(i + 1, j, k) * u * (1 - v) * (1 - w) +
					getDensity(i + 1, j, k + 1) * u * (1 - v) * (w)+
					getDensity(i, j, k + 1) * (1 - u) * (1 - v) * w +
					getDensity(i, j + 1, k) * (1 - u) * v * (1 - w) +
					getDensity(i + 1, j + 1, k) * v * u * (1 - w) +
					getDensity(i + 1, j + 1, k + 1) * v * w * u +
					getDensity(i, j + 1, k + 1) * (1 - u) * v * w;

				if (volumeSwitch == 0) {
					sum = sum + fdens * step;
				}
				else if (volumeSwitch == 1) {
					maxval = max(maxval, fdens);
				}
			}
			if (volumeSwitch == 0) {
				curImage.data[getPixel(j, i)] = round(sum);
			}
			else if (volumeSwitch == 1) {
				curImage.data[getPixel(j, i)] = round(maxval);
			}
		}
	}
	gui->DisplayWindow->redraw();
}