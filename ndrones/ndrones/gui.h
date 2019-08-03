#pragma once
#ifndef GUI_H
#define GUI_H
#include <iostream>
#include <vector>
#include <queue>

// FLTK
#define WIN32
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Draw.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_JPEG_Image.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Menu_Bar.H>

// Local includes
#include "base.h"
#include "config.h"

// GLOBALS
static Fl_Input        *G_inptitle = NULL;
static Fl_Input        *G_inpfilter = NULL;
static Fl_Input        *G_inpdirname = NULL;
static Fl_File_Chooser *G_chooser = NULL;


// The main class that display images and drawings
class Canvas : public Fl_Group {
protected:
	// Callback function that handles events through message e, inherited from the group class of FLTK
	int handle(int e);

	// The draw line function that uses standard coord
	void fl_normal_line(float x, float y, float x1, float y1);

	// Draw mouse coords in small black rectangle
	void drawCoords();

	// Draw the axes and the unit square
	void drawAxes();

	// Draw the agents
	void drawAgents();

	// Draw the packages
	void drawPackages();

	// Draw the targets;
	void drawTargets();

	// The main draw function of canvas
	void draw();

	// Convert scenario coords to canvas coord
	void scenToCanvasCoord(float &x, float &y);

	// Convert canvas coords to scenario coord
	void canvasToScenCoord(float &x, float &y);

public:
	// Box to display the background image
	Fl_Box * imageBox;
	static int canvasWidth;
	static int canvasHeight;

	// The canvas to render the ink
	float*** canvas;

	// Axes and unit square information
	int axisLength;
	int origin;
	int squareLength;

	// The scenario to be demo'd
	static Scenario scenario;

	Canvas(int X, int Y, int W, int H, const char *L);
	~Canvas();
	int mouseX, mouseY;
};


class GUI {
public:
	static Canvas *canvas;

	// Buttons
	Fl_Button *exitBu;
	Fl_Button *runStepBu;
	Fl_Button *runAllBu;
	Fl_Button *randomBu;
	Fl_Button *clearBu;
	Fl_Menu_Bar *menuBar;

	Fl_Window *win;

	static std::string canvasFileName;

	// Big redraw signal makes the canvas wipe everything, it's a bit slower so normally the small signal is preferred
	static bool bigRedrawSignal;

	GUI(int winWidth, int winHeight);
	~GUI();

	// Callback function for exit button
	static void exitCallback(Fl_Widget*w, void*data);
	// Callback function for browse button
	static void browseCallback(Fl_Widget*w, void*data);
	// Callback function for save button
	static void saveResultCallback(Fl_Widget*w, void*data);

};

#endif