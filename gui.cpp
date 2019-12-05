// generated by Fast Light User Interface Designer (fluid) version 1.0305

#include "gui.h"

void Gui::cb_readFile_i(Fl_Menu_*, void*) {
  app->ReadFile();
}
void Gui::cb_readFile(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_readFile_i(o,v);
}

void Gui::cb_writeFile_i(Fl_Menu_*, void*) {
  app->WriteFile();
}
void Gui::cb_writeFile(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_writeFile_i(o,v);
}

void Gui::cb_updateGUI_i(Fl_Menu_*, void*) {
  app->GUI_Update();
}
void Gui::cb_updateGUI(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_updateGUI_i(o,v);
}

void Gui::cb_greyEffect_i(Fl_Menu_*, void*) {
  app->GUI_GreyEffect();
}
void Gui::cb_greyEffect(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_greyEffect_i(o,v);
}

void Gui::cb_averageSmooth_i(Fl_Menu_*, void*) {
  app->GUI_AverageSmooth();
}
void Gui::cb_averageSmooth(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_averageSmooth_i(o,v);
}

void Gui::cb_medianSmooth_i(Fl_Menu_*, void*) {
  app->GUI_MedianSmooth();
}
void Gui::cb_medianSmooth(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_medianSmooth_i(o,v);
}

void Gui::cb_gaussianSmooth_i(Fl_Menu_*, void*) {
  app->GUI_GaussianSmooth();
}
void Gui::cb_gaussianSmooth(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_gaussianSmooth_i(o,v);
}

void Gui::cb_edgeDetect_i(Fl_Menu_*, void*) {
  app->GUI_EdgeDetect();
}
void Gui::cb_edgeDetect(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_edgeDetect_i(o,v);
}

void Gui::cb_undoButton_i(Fl_Menu_*, void*) {
  app->GUI_UndoButton();
}
void Gui::cb_undoButton(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_undoButton_i(o,v);
}

void Gui::cb_readVolumeFile_i(Fl_Menu_*, void*) {
  app->GUI_ReadVolumeFile();
}
void Gui::cb_readVolumeFile(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_readVolumeFile_i(o,v);
}

void Gui::cb_xRay_i(Fl_Menu_*, void*) {
  app->GUI_xRay();
}
void Gui::cb_xRay(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_xRay_i(o,v);
}

void Gui::cb_MIP_i(Fl_Menu_*, void*) {
  app->GUI_MIP();
}
void Gui::cb_MIP(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_MIP_i(o,v);
}

void Gui::cb_restoreButton_i(Fl_Menu_*, void*) {
  app->GUI_RestoreButton();
}
void Gui::cb_restoreButton(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_restoreButton_i(o,v);
}

void Gui::cb_exitButton_i(Fl_Menu_*, void*) {
  exit(1);
}
void Gui::cb_exitButton(Fl_Menu_* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_exitButton_i(o,v);
}

Fl_Menu_Item Gui::menu_menuBar[] = {
 {"File", 0,  0, 0, 64, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {"Read", 0,  (Fl_Callback*)Gui::cb_readFile, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {"Write", 0,  (Fl_Callback*)Gui::cb_writeFile, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {0,0,0,0,0,0,0,0,0},
 {"Image Processing", 0,  0, 0, 64, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {"Update", 0,  (Fl_Callback*)Gui::cb_updateGUI, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {"grey effect", 0,  (Fl_Callback*)Gui::cb_greyEffect, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {"Average Smooth", 0,  (Fl_Callback*)Gui::cb_averageSmooth, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {"Median Smooth", 0,  (Fl_Callback*)Gui::cb_medianSmooth, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {"Gaussian Smooth", 0,  (Fl_Callback*)Gui::cb_gaussianSmooth, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {"Edge Detect", 0,  (Fl_Callback*)Gui::cb_edgeDetect, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {"Undo", 0,  (Fl_Callback*)Gui::cb_undoButton, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {0,0,0,0,0,0,0,0,0},
 {"Volume Rendering", 0,  0, 0, 64, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {"Read Volume File", 0,  (Fl_Callback*)Gui::cb_readVolumeFile, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {"X-Ray", 0,  (Fl_Callback*)Gui::cb_xRay, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {"MIP", 0,  (Fl_Callback*)Gui::cb_MIP, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
 {0,0,0,0,0,0,0,0,0},
 {"Restore", 0,  (Fl_Callback*)Gui::cb_restoreButton, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 4},
 {"Exit", 0,  (Fl_Callback*)Gui::cb_exitButton, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 1},
 {0,0,0,0,0,0,0,0,0}
};
Fl_Menu_Item* Gui::fileMenu = Gui::menu_menuBar + 0;
Fl_Menu_Item* Gui::readFile = Gui::menu_menuBar + 1;
Fl_Menu_Item* Gui::writeFile = Gui::menu_menuBar + 2;
Fl_Menu_Item* Gui::imageProcess = Gui::menu_menuBar + 4;
Fl_Menu_Item* Gui::updateGUI = Gui::menu_menuBar + 5;
Fl_Menu_Item* Gui::greyEffect = Gui::menu_menuBar + 6;
Fl_Menu_Item* Gui::averageSmooth = Gui::menu_menuBar + 7;
Fl_Menu_Item* Gui::medianSmooth = Gui::menu_menuBar + 8;
Fl_Menu_Item* Gui::gaussianSmooth = Gui::menu_menuBar + 9;
Fl_Menu_Item* Gui::edgeDetect = Gui::menu_menuBar + 10;
Fl_Menu_Item* Gui::undoButton = Gui::menu_menuBar + 11;
Fl_Menu_Item* Gui::volumeRendering = Gui::menu_menuBar + 13;
Fl_Menu_Item* Gui::readVolumeFile = Gui::menu_menuBar + 14;
Fl_Menu_Item* Gui::xRay = Gui::menu_menuBar + 15;
Fl_Menu_Item* Gui::MIP = Gui::menu_menuBar + 16;
Fl_Menu_Item* Gui::restoreButton = Gui::menu_menuBar + 18;
Fl_Menu_Item* Gui::exitButton = Gui::menu_menuBar + 19;

void Gui::cb_sliderX_i(Fl_Value_Slider*, void*) {
  app->GUI_VolumeDraw();
}
void Gui::cb_sliderX(Fl_Value_Slider* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_sliderX_i(o,v);
}

void Gui::cb_sliderY_i(Fl_Value_Slider*, void*) {
  app->GUI_VolumeDraw();
}
void Gui::cb_sliderY(Fl_Value_Slider* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_sliderY_i(o,v);
}

void Gui::cb_sliderZ_i(Fl_Value_Slider*, void*) {
  app->GUI_VolumeDraw();
}
void Gui::cb_sliderZ(Fl_Value_Slider* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_sliderZ_i(o,v);
}

void Gui::cb_sliderWidth_i(Fl_Value_Slider*, void*) {
  app->GUI_VolumeDraw();
}
void Gui::cb_sliderWidth(Fl_Value_Slider* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_sliderWidth_i(o,v);
}

void Gui::cb_sliderHeight_i(Fl_Value_Slider*, void*) {
  app->GUI_VolumeDraw();
}
void Gui::cb_sliderHeight(Fl_Value_Slider* o, void* v) {
  ((Gui*)(o->parent()->user_data()))->cb_sliderHeight_i(o,v);
}

Gui::Gui() {
  { MainWindow = new Fl_Double_Window(985, 555, "Scientific Visualization Projects");
    MainWindow->user_data((void*)(this));
    { menuBar = new Fl_Menu_Bar(0, 0, 985, 25, "menuBar");
      menuBar->menu(menu_menuBar);
    } // Fl_Menu_Bar* menuBar
    { EditorWindow = new CEditorWindow(15, 25, 385, 350, "EditorWindow");
      EditorWindow->box(FL_NO_BOX);
      EditorWindow->color(FL_BACKGROUND_COLOR);
      EditorWindow->selection_color(FL_BACKGROUND_COLOR);
      EditorWindow->labeltype(FL_NORMAL_LABEL);
      EditorWindow->labelfont(0);
      EditorWindow->labelsize(14);
      EditorWindow->labelcolor(FL_FOREGROUND_COLOR);
      EditorWindow->align(Fl_Align(FL_ALIGN_CENTER));
      EditorWindow->when(FL_WHEN_RELEASE);
    } // CEditorWindow* EditorWindow
    { DisplayWindow = new CDisplayWindow(415, 25, 570, 530, "DisplayWindow");
      DisplayWindow->box(FL_NO_BOX);
      DisplayWindow->color(FL_BACKGROUND_COLOR);
      DisplayWindow->selection_color(FL_BACKGROUND_COLOR);
      DisplayWindow->labeltype(FL_NORMAL_LABEL);
      DisplayWindow->labelfont(0);
      DisplayWindow->labelsize(14);
      DisplayWindow->labelcolor(FL_FOREGROUND_COLOR);
      DisplayWindow->align(Fl_Align(FL_ALIGN_CENTER));
      DisplayWindow->when(FL_WHEN_RELEASE);
    } // CDisplayWindow* DisplayWindow
    { sliderX = new Fl_Value_Slider(15, 390, 295, 25, "x");
      sliderX->type(1);
      sliderX->minimum(-180);
      sliderX->maximum(180);
      sliderX->step(1);
      sliderX->callback((Fl_Callback*)cb_sliderX);
    } // Fl_Value_Slider* sliderX
    { sliderY = new Fl_Value_Slider(15, 445, 295, 25, "y");
      sliderY->type(1);
      sliderY->minimum(-180);
      sliderY->maximum(180);
      sliderY->step(1);
      sliderY->callback((Fl_Callback*)cb_sliderY);
    } // Fl_Value_Slider* sliderY
    { sliderZ = new Fl_Value_Slider(15, 500, 295, 25, "z");
      sliderZ->type(1);
      sliderZ->minimum(-180);
      sliderZ->maximum(180);
      sliderZ->step(1);
      sliderZ->callback((Fl_Callback*)cb_sliderZ);
    } // Fl_Value_Slider* sliderZ
    { sliderWidth = new Fl_Value_Slider(330, 380, 25, 145, "width");
      sliderWidth->minimum(256);
      sliderWidth->maximum(512);
      sliderWidth->step(1);
      sliderWidth->value(256);
      sliderWidth->callback((Fl_Callback*)cb_sliderWidth);
    } // Fl_Value_Slider* sliderWidth
    { sliderHeight = new Fl_Value_Slider(375, 380, 25, 145, "height");
      sliderHeight->minimum(256);
      sliderHeight->maximum(512);
      sliderHeight->step(1);
      sliderHeight->value(256);
      sliderHeight->callback((Fl_Callback*)cb_sliderHeight);
    } // Fl_Value_Slider* sliderHeight
    MainWindow->end();
  } // Fl_Double_Window* MainWindow
  app=new Application();
}

void Gui::show() {
  MainWindow->show();
  EditorWindow->show();
  DisplayWindow->show();
}
