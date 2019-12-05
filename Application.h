
#ifndef __Application_H
#define __Application_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// put your global typedefs here

typedef struct{
	unsigned char *data;  // an image of bytes
	int *intData;         // an image of integers
    float *fdata;         // an image of floats
    int nx,ny,n;          // image dimensions
    int ncolorChannels;   // number of color channels in the image (1 or 3)
} Image;

typedef struct{
	unsigned char *data;  // a volume of bytes
    int nx,ny,nz,n;       // image dimensions (nz is for a volume, a "3D image")
} Volume;

class Application {
public:
  Application();

  void ReadFile();
  void WriteFile();
  void FlipImage(Image *img);
 
// put your application prototypes here: 

  void GUI_Update();
  void GUI_GreyEffect();
  void GUI_AverageSmooth();
  void GUI_MedianSmooth();
  void GUI_GaussianSmooth();
  void GUI_EdgeDetect();
  void GUI_UndoButton();
  void GUI_RestoreButton();
  void Caluate_Histogram();

  void GUI_ReadVolumeFile();
  void GUI_xRay();
  void GUI_MIP();

  /*void GUI_SliderX();
  void GUI_SliderY();
  void GUI_SliderZ();
  void GUI_SliderWidth();
  void GUI_SliderHeight();*/
  void GUI_VolumeDraw();
};

#endif
