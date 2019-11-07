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
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Menu_Bar.H>

// Local includes
#include "base.h"
#include "config.h"
#include "scenario.h"

// GLOBALS
static Fl_Input        *G_inptitle = NULL;
static Fl_Input        *G_inpfilter = NULL;
static Fl_Input        *G_inpdirname = NULL;
static Fl_File_Chooser *G_chooser = NULL;


enum AnimationMode{DRONEANI, PACKAGEANI};
enum GridVisualizationMode {UNDEFINED, STPMAP, DRONEUSAGEMAP, DEPOTUSAGEMAP};


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

	// Draw grid lines
	void drawGridLines();

	// Draw grid points
	void drawGridPoints();

	// Draw the agents
	void drawAgents();

	void drawLineAnis(std::vector<LineAnimation> &anis, float timeElapsed);

	// Draw animations in the scenario
	void drawAnimation();

	// Draw the packages
	void drawPackages();

	// Draw the targets;
	void drawTargets();

	// Draw the discrete points
	void drawDiscretePoints();

	// The main draw function of canvas
	void draw();

	// Convert scenario coords to canvas coord
	void scenToCanvasCoord(float &x, float &y);

	// Convert canvas coords to scenario coord
	void canvasToScenCoord(float &x, float &y);

	// The timer function to activate redraw
	static void timerCallback(void *userData);

public:
	// Box to display the background image
	Fl_Box * imageBox;
	static int canvasWidth;
	static int canvasHeight;
	static int canvasX;
	static int canvasY;

	static float marginScale;
	static int margin;

	// The canvas to render the ink
	float*** canvas;

	// Axes and unit square information
	int axisLength;
	int origin;
	int squareLength;

	// The scenario to be demo'd
	static Scenario scenario;
	static float timer;

	// Animation
	AnimationMode aniMode;

	// Grid visualization
	GridVisualizationMode gridVisMode;
	std::vector<Color> gridColor;

	// Compute the colors for grid points based on one of the GridVisualizationMode
	void computeColorMap();

	Canvas(int X, int Y, int W, int H, const char *L);
	~Canvas();
	int mouseX, mouseY;
};


class GUI {
public:
	static Canvas *canvas;

	// Buttons
	Fl_Button *runAllBu;
	Fl_Button *droneAniBu;
	Fl_Button *gridVisBu;
	Fl_Button *scriptBu;
	Fl_Button *clearBu;
	Fl_Menu_Bar *menuBar;
	Fl_Round_Button *droneAniOptBu;
	Fl_Round_Button *packageAniOptBu;
	Fl_Round_Button *shortestPathMapOptBu; // Radio button for shorest path map
	Fl_Round_Button *droneUsageMapOptBu; // Radio button for drone usage map
	Fl_Round_Button *depotUsageMapOptBu; // Radio button for depot usage (useful only if there are more than one depot)
	Fl_Window *win;

	static std::string canvasFileName;

	static bool drawSignal;

	// Big redraw signal makes the canvas wipe everything, it's a bit slower so normally the small signal is preferred
	static bool bigRedrawSignal;

	static bool drawGridSignal;

	GUI(int winWidth, int winHeight);
	~GUI();

	// Callback function for exit button
	static void exitCallback(Fl_Widget*w, void*data);
	// Callback function for browse button
	static void browseCallback(Fl_Widget*w, void*data);
	// Callback function for save button
	static void saveResultCallback(Fl_Widget*w, void*data);

	// Callback for the solver
	static void solverCallback(Fl_Widget*w, void*data);

	// Callback to set draw signal = true
	static void drawSignalCallback(Fl_Widget*w, void*data);
	// Callback to set draw signal for the grid = true
	static void drawGridSignalCallback(Fl_Widget*w, void*data);
	// Callback to run a script
	static void runScriptCallback(Fl_Widget*w, void*data);
	// Callback to set animation mode
	static void aniRadioCallback(Fl_Widget*w, void*data);

	// Callback to set the grid visualization mode
	static void gridVisRadioCallback(Fl_Widget*w, void*data);
};

#endif