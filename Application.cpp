#define _CRT_SECURE_NO_DEPRECATE
#include "Application.h"
#include <FL/fl_file_chooser.H>
#include "Gui.h"
#include <iostream>
#include <list>
#include <algorithm>

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

// help functions

int getPixel(int i, int j) {
	return (i * curImage.nx + j) * curImage.ncolorChannels;
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