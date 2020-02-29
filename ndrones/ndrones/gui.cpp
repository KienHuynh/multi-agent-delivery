#include "gui.h"

int Canvas::canvasWidth = 0;
int Canvas::canvasHeight = 0;
int Canvas::canvasX = 0;
int Canvas::canvasY = 0;
float Canvas::marginScale = 0;
int Canvas::margin = 0;
Scenario Canvas::scenario;

bool GUI::bigRedrawSignal = 0;
bool GUI::drawSignal = false;
bool GUI::drawGridSignal = false;
ZoomMode GUI::zoomMode = ZoomMode::ALL;
std::string GUI::canvasFileName = "";
Canvas* GUI::canvas = NULL;


void GUI::exitCallback(Fl_Widget*w, void*data) {
	exit(0);
}


void GUI::browseCallback(Fl_Widget*w, void*data) {
	if (!G_chooser) {
		G_chooser = new Fl_File_Chooser("", "", Fl_File_Chooser::SINGLE, "");
	}

	G_chooser->directory(NULL);
	G_chooser->filter(NULL);
	G_chooser->type(Fl_File_Chooser::SINGLE);
	G_chooser->label(NULL);
	G_chooser->show();

	// Block until user picks something.
	//     (The other way to do this is to use a callback())
	//
	while (G_chooser->shown()) {
		Fl::wait();
	}

	// Print the results
	if (G_chooser->value() == NULL) {
		fprintf(stderr, "(User hit 'Cancel')\n");
		return;
	}

	fprintf(stderr, "DIRECTORY: '%s'\n", G_chooser->directory());
	fprintf(stderr, "    VALUE: '%s'\n", G_chooser->value());
	fprintf(stderr, "    COUNT: %d files selected\n", G_chooser->count());


	// Flag the signal so that the image box will redraw
	Canvas::scenario.loadFile(G_chooser->value());
	canvas->redraw();
	fprintf(stderr, "--------------------\n");
}


void GUI::saveResultCallback(Fl_Widget*w, void*data) {
	if (!G_chooser) {
		G_chooser = new Fl_File_Chooser("", "", Fl_File_Chooser::CREATE, "");
	}

	G_chooser->directory(NULL);
	G_chooser->filter(NULL);
	G_chooser->type(Fl_File_Chooser::CREATE);
	G_chooser->label(NULL);
	G_chooser->show();

	// Block until user picks something.
	//     (The other way to do this is to use a callback())
	//
	while (G_chooser->shown()) {
		Fl::wait();
	}

	// Print the results
	if (G_chooser->value() == NULL) {
		fprintf(stderr, "(User hit 'Cancel')\n");
		return;
	}

	fprintf(stderr, "DIRECTORY: '%s'\n", G_chooser->directory());
	fprintf(stderr, "    VALUE: '%s'\n", G_chooser->value());
	fprintf(stderr, "    COUNT: %d files selected\n", G_chooser->count());
	fprintf(stderr, "--------------------\n");

	Canvas::scenario.writeSolution(G_chooser->value());
}


void GUI::solverCallback(Fl_Widget *w, void *data) {
	Canvas::scenario.solve();
	Canvas::scenario.aniStart = true;
	canvas->redraw();
}


void GUI::drawSignalCallback(Fl_Widget *w, void *data) {
	drawSignal = true;
	bigRedrawSignal = true;
	if (canvas->aniMode == DRONEANI) Canvas::scenario.createDroneAnimation();
	if (canvas->aniMode == PACKAGEANI) Canvas::scenario.createPackageAnimation();
	auto now = std::chrono::system_clock::now().time_since_epoch();
	int mili = std::chrono::duration_cast<std::chrono::milliseconds>(now).count();
	Canvas::scenario.timer = ((float)mili) / 1000;
}


void GUI::zoomToggleCallback(Fl_Widget *w, void *data) {
	GUI::zoomMode = (ZoomMode) (1 - (int) GUI::zoomMode);
	fl_color(255);
	fl_rectf(0, 0, Canvas::canvasWidth+100, Canvas::canvasHeight+100);
	canvas->redraw();
}


void Canvas::computeColorMap() {
	float s = scenario.points[0].bestTime;
	float l = s;

	for (int i = 1; i < scenario.points.size(); i++) {
		if (s > scenario.points[i].bestTime) s = scenario.points[i].bestTime;
		if (l < scenario.points[i].bestTime) l = scenario.points[i].bestTime;
	}

	// Compute color map
	for (int i = 0; i < scenario.points.size(); i++) {
		Color color(255, 255, 255);
		if (gridVisMode == STPMAP) {
			color = Scenario::agentQueueColorMap(scenario.points[i].agentQueue);	
		}
		if (gridVisMode == STIMEMAP) {
			color = Scenario::bestTimeColorMap(s, l, scenario.points[i].bestTime);
		}
		if (gridVisMode == DEPOTUSAGEMAP) {
			int depot = -1;
			for (int p = 0; p < scenario.packages.size(); p++) {
				if (scenario.points[i].pointQueue.size() == 0) continue;
				if (scenario.points[i].pointQueue[0] == scenario.packages[p].loc) {
					depot = p;
					break;
				}
			}
			
			if (depot > -1)	color = Scenario::depotColorMap(scenario.packages[depot]);
		}
		gridColor.push_back(color);
	}
}

void GUI::drawGridSignalCallback(Fl_Widget *w, void *data) {
	if (canvas->gridVisMode == UNDEFINED) return;
	drawGridSignal = !drawGridSignal;
	canvas->computeColorMap();
};


void GUI::runScriptCallback(Fl_Widget *w, void *data) {

}


void GUI::aniRadioCallback(Fl_Widget*w, void*data) {
	canvas->aniMode = (AnimationMode) (int)data;
}


void GUI::gridVisRadioCallback(Fl_Widget*w, void*data) {
	canvas->gridVisMode = (GridVisualizationMode)(int)data;
}


void Canvas::timerCallback(void *userData) {
	Canvas *o = (Canvas*)userData;
	o->redraw();
	Fl::repeat_timeout(0.001, timerCallback, userData);
}


void Canvas::fl_normal_line(float x, float y, float x1, float y1) {
	fl_line(x, canvasHeight - y, x1, canvasHeight - y1);
}


// Function to handle mouse events such as drag, release, etc.
int Canvas::handle(int e) {
	int ret = Fl_Group::handle(e);
	switch (e) {
	// Mouse events
	case FL_DRAG: {
		ret = 0;
		break;
	}
	case FL_RELEASE: {
		float x = (int)Fl::event_x();
		float y = (int)Fl::event_y();

		// Convert x and y to normalized coord
		y = canvasHeight - y;
		x = (x - origin) / (float)squareLength;
		y = (y - origin) / (float)squareLength;

		if (!(x < 0 || y < 0 || x > 1 || y > 1)) {
			
		}

		ret = 1;
		redraw();
		break;
	}
	case FL_PUSH: {
		// Has to return 1 so to make it clear that this widget "wants" to be clicked.
		// Otherwise FL_DRAG and such will be ignored
		ret = 1;
		break;
	}
	case FL_ENTER:
	{
		mouseX = (int)Fl::event_x();
		mouseY = (int)Fl::event_y();
		ret = 1;                // FL_ENTER: must return(1) to receive FL_MOVE
		break;
	}

	case FL_MOVE:
	{
		mouseX = (int)Fl::event_x();
		mouseY = (int)Fl::event_y();
		// FL_MOVE: mouse movement causes 'user damage' and redraw..
		damage(FL_DAMAGE_USER1);
		ret = 1;
		break;
	}

	// Keyboard events
	case FL_FOCUS:
	case FL_UNFOCUS:
		ret = 1;                // enables receiving keyboard events
		break;
	case FL_SHORTCUT:           // incase widget that isn't ours has focus
	case FL_KEYDOWN:
	{
		int key = Fl::event_key();
		if (key == 'z') {
			//TODO: Implement undo
			this->redraw();
			//this->draw();
			ret = 1;
		}
		else {
			ret = 0;
		}
		if (key == 'y') {
			//TODO: Implement redo
			ret = 0;
		}

		break;
	}

	}
	
	return(ret);
}


void Canvas::scenToCanvasCoord(float &x, float &y) {
	float xDiff, yDiff;
	float minX, minY, maxX, maxY;
	if (GUI::zoomMode == ZoomMode::ALL) {
		minX = scenario.minX;
		maxX = scenario.maxX;
		minY = scenario.minY;
		maxY = scenario.maxY;

		xDiff = maxX - minX;
		yDiff = maxY - minY;
	}
	if (GUI::zoomMode == ZoomMode::ST) {
		minX = scenario.packages[0].loc.x;
		maxX = scenario.packages[0].loc.x;
		minY = scenario.packages[0].loc.y;
		maxY = scenario.packages[0].loc.y;
		for (auto package : scenario.packages) {
			if (minX > package.loc.x) minX = package.loc.x;
			if (maxX < package.loc.x) maxX = package.loc.x;
			if (minY > package.loc.y) minY = package.loc.y;
			if (maxY < package.loc.y) maxY = package.loc.y;
		}
		for (auto target : scenario.targets) {
			if (minX > target.loc.x) minX = target.loc.x;
			if (maxX < target.loc.x) maxX = target.loc.x;
			if (minY > target.loc.y) minY = target.loc.y;
			if (maxY < target.loc.y) maxY = target.loc.y;
		}
		xDiff = maxX - minX;
		yDiff = maxY - minY;

		if (xDiff > yDiff) {
			yDiff = xDiff;
			minY = (maxY + minY - xDiff) / 2;
		}
		else {
			xDiff = yDiff;
			minX = (maxX + minX - yDiff) / 2;
		}
	}
	

	float diff = xDiff > yDiff ? xDiff : yDiff;
	x = (x + abs(minX))*canvasWidth;
	x /= diff;
	x = x*marginScale + margin + canvasX;
	// (canvasHeight - ((y + abs)*canvasHeight / ((float) (scenario.maxY - scenario.minY))))*marginScale + margin + canvasY;
	y = ((y + abs(minY))*canvasHeight / diff);
	y = canvasHeight - y;
	y = y*marginScale + margin + canvasY;
}


void Canvas::canvasToScenCoord(float &x, float &y) {
	x = ((x - canvasX - margin)/marginScale)*((float)(scenario.maxX - scenario.minX))/((float)canvasWidth);
	y = y - margin - canvasY;
	y = canvasHeight - y;
	y = y * ((float)(scenario.maxY - scenario.minY)) / ((float)canvasHeight);
}


void Canvas::drawPie360(int x, int y, int dx, int dy, char str, int label) {
	fl_pie(x - 5, y - 5, dx, dy, 0, 360);

	fl_color(0, 0, 0);
		
	//if (GUI::zoomMode == ZoomMode::ALL) {
		char s[80];
		if (label < 0) return;
		sprintf_s(s, "%d", label);
		fl_font(FL_HELVETICA, 18);
		fl_draw(s, x + 5, y + 5);
	//}
		
	/*if (GUI::zoomMode == ZoomMode::ST) {
		char s[80];
		if (label < 0) sprintf_s(s, "%c %d, %d", str, x, y);
		else sprintf_s(s, "%c %d, %d, %d", str, label, x, y);
		
		fl_font(FL_HELVETICA, 18);
		fl_draw(s, x + 5, y + 5);
	}*/
}


void Canvas::drawCoords() {
	// Coordinates as a string
	char s[80];
	float x = (float)Fl::event_x();
	float y = (float)Fl::event_y();
	canvasToScenCoord(x, y);
	sprintf_s(s, "x=%.3f y=%.3f", x, y);
	// Black rect
	fl_color(FL_BLACK);
	fl_rectf(Canvas::canvasWidth - 180, this->y(), 180, 25);
	// White text
	fl_color(FL_WHITE);
	fl_font(FL_HELVETICA, 18);
	fl_draw(s, Canvas::canvasWidth - 180 + 7, this->y() + 20);
}


void Canvas::drawAxes() {
	// Draw the square
	char dashStyle[3] = { 10,5,0 };
	fl_color(fl_rgb_color(100, 100, 255));
	fl_line_style(FL_DASH, 2, dashStyle);
	fl_normal_line(origin, origin, origin, origin + squareLength);
	fl_normal_line(origin, origin, origin + squareLength, origin);
	fl_normal_line(origin + squareLength, origin, origin + squareLength, origin + squareLength);
	fl_normal_line(origin, origin + squareLength, origin + squareLength, origin + squareLength);

	// Draw the axes
	fl_color(fl_rgb_color(50, 50, 200));
	fl_line_style(FL_SOLID, 2);

	// Horizontal axis
	fl_normal_line(origin, origin, origin + axisLength, origin);
	// Arrow tip
	fl_normal_line(origin + axisLength, origin, origin + axisLength - 7, origin - 5);
	fl_normal_line(origin + axisLength, origin, origin + axisLength - 7, origin + 5);

	// Vertical axis
	fl_normal_line(origin, origin, origin, origin + axisLength);
	// Arrow tip
	fl_normal_line(origin, origin + axisLength, origin - 5, origin + axisLength - 7);
	fl_normal_line(origin, origin + axisLength, origin + 5, origin + axisLength - 7);

}


void Canvas::drawAgents() {
	for (int i = 0; i < scenario.agents.size(); i++) {
		float x = scenario.agents[i].loc.x;
		float y = scenario.agents[i].loc.y;
		scenToCanvasCoord(x, y);
		float v = scenario.agents[i].v;

		// Assign the color based on the velocity and the max/min speed in the scenario
		int color = 255*(v - scenario.minSpeed) / (scenario.maxSpeed - scenario.minSpeed);
		fl_color(fl_rgb_color(25, 25, color));

		drawPie360(x, y, 10, 10, 'a', scenario.agents[i].v);
	}
}


void Canvas::drawPackages() {
	if (scenario.packageInputMode == SINGLE_POINT) {
		for (int i = 0; i < scenario.packages.size(); i++) {
			float x = scenario.packages[i].loc.x;
			float y = scenario.packages[i].loc.y;
			scenToCanvasCoord(x, y);
			fl_color(200, 25, 25);

			drawPie360(x, y, 10, 10, 't', scenario.packages[i].ID);
		}
	}
	
	if (scenario.packageInputMode == POLY) {
		int n = scenario.packages.size();
		for (int i = 0; i < n; i++) {
			float x0 = scenario.packages[i].loc.x;
			float y0 = scenario.packages[i].loc.y;
			float x1 = scenario.packages[(i + 1) % n].loc.x;
			float y1 = scenario.packages[(i + 1) % n].loc.y;
			scenToCanvasCoord(x0, y0);
			scenToCanvasCoord(x1, y1);
			fl_color(200, 25, 25);
			//fl_normal_line(x0, canvasHeight - y0, x1, canvasHeight - y1);
			fl_line(x0, y0, x1, y1);

			// Draw target ID
			fl_color(0, 0, 0);
			if (scenario.packages[i].ID >= 0 && i == 0) {
				char s[80];
				sprintf_s(s, "%d", scenario.packages[i].ID);
				fl_font(FL_HELVETICA, 18);
				fl_draw(s, x0 + 5, y0 + 5);
			}
		}
	}
}


void Canvas::drawTargets() {
	if (scenario.targetInputMode == SINGLE_POINT) {
		for (int i = 0; i < scenario.targets.size(); i++) {
			float x = scenario.targets[i].loc.x;
			float y = scenario.targets[i].loc.y;
			scenToCanvasCoord(x, y);
			fl_color(25, 200, 25);
			fl_pie(((int)x) - 5, ((int)y) - 5, 10, 10, 0, 360);

			// Draw target ID
			fl_color(0, 0, 0);
			if (scenario.targets[i].ID >= 0) {
				char s[80];
				sprintf_s(s, "%d", scenario.targets[i].ID);
				fl_font(FL_HELVETICA, 18);
				fl_draw(s, x + 5, y + 5);
			}
		}
	}

	if (scenario.targetInputMode == POLY) {
		int n = scenario.targets.size();
		for (int i = 0; i < n; i++) {
			float x0 = scenario.targets[i].loc.x;
			float y0 = scenario.targets[i].loc.y;
			float x1 = scenario.targets[(i + 1) % n].loc.x;
			float y1 = scenario.targets[(i + 1) % n].loc.y;
			scenToCanvasCoord(x0, y0);
			scenToCanvasCoord(x1, y1);
			fl_color(25, 200, 25);
			//fl_normal_line(x0, canvasHeight - y0, x1, canvasHeight - y1);
			fl_line(x0, y0, x1, y1);

			// Draw target ID
			fl_color(0, 0, 0);
			if (scenario.targets[i].ID >= 0 && i == 0) {
				char s[80];
				sprintf_s(s, "%d", scenario.targets[i].ID);
				fl_font(FL_HELVETICA, 18);
				fl_draw(s, x0 + 5, y0 + 5);
			}
		}
	}
}


void Canvas::drawObstacles() {
	for (auto o : scenario.obstacles) {
		for (auto t : o.triIdx) {
			float x0 = o.points[t[0]].x;
			float y0 = o.points[t[0]].y;
			float x1 = o.points[t[1]].x;
			float y1 = o.points[t[1]].y;
			float x2 = o.points[t[2]].x;
			float y2 = o.points[t[2]].y;
			scenToCanvasCoord(x0, y0);
			scenToCanvasCoord(x1, y1);
			scenToCanvasCoord(x2, y2);
			
			fl_color(250, 200, 100);
			fl_polygon(x0, y0, x1, y1, x2, y2);
		}
		
		for (int i = 0; i < o.hullPtIdx.size(); i++) {
			int j0 = o.hullPtIdx[i];
			int j1 = o.hullPtIdx[(i + 1) % o.hullPtIdx.size()];
			float x0 = o.points[j0].x;
			float y0 = o.points[j0].y;
			float x1 = o.points[j1].x;
			float y1 = o.points[j1].y;
			scenToCanvasCoord(x0, y0);
			scenToCanvasCoord(x1, y1);
			
			fl_color(200, 100, 250);
			fl_line_style(FL_SOLID, 2);
			fl_line(x0, y0, x1, y1);
		}
	}
}


void Canvas::drawDiscretePoints() {
	for (int i = 0; i < scenario.points.size(); i++) {
		float x = scenario.points[i].p.x;
		float y = scenario.points[i].p.y;
		scenToCanvasCoord(x, y);
		fl_color(200, 200, 200);
		fl_pie(((int)x) - 5, ((int)y) - 5, 10, 10, 0, 360);
	}
}


void Canvas::drawLineAnis(std::vector<LineAnimation> &anis, float timeElapsed) {
	for (int i = 0; i < anis.size(); i++) {
		LineAnimation ani = anis[i];

		// Special treatment for the first loop
		Point2D newP;
		Point2D prevP;
		bool draw = false;
		if (timeElapsed >= ani.startTime && timeElapsed <= ani.endTime && anis[i].active == false) {
			if (anis[i].prevTimer < 0) {
				anis[i].prevTimer = timeElapsed;
				continue;
			}

			anis[i].active = true;
			newP = ((ani.end - ani.start) * (timeElapsed - ani.startTime)
				/ ani.duration) + ani.start;
			prevP = ani.start;
			draw = true;
		}

		else if (timeElapsed >= ani.startTime && timeElapsed <= ani.endTime && anis[i].active == true) {
			newP = ((ani.end - ani.start) * (timeElapsed - ani.startTime)
				/ ani.duration) + ani.start;
			prevP = ((ani.end - ani.start) * (ani.prevTimer - ani.startTime)
				/ ani.duration) + ani.start;
			draw = true;
		}

		// Last loop
		if (timeElapsed > ani.endTime && ani.active == true) {
			anis[i].active = false;

			newP = ani.end;
			prevP = ((ani.end - ani.start) * (ani.prevTimer - ani.startTime)
				/ ani.duration) + ani.start;
			draw = true;
		}

		// Time passed without the ani activated because the speed of the ani was too fast for the comp
		if (timeElapsed > ani.endTime && ani.active == false) {
			newP = ani.end;
			prevP = ani.start;
			draw = true;
		}

		if (draw) {
			fl_color(fl_rgb_color(ani.color[0], ani.color[1], ani.color[2]));
			fl_line_style(FL_SOLID, 4);
			scenToCanvasCoord(newP.x, newP.y);
			scenToCanvasCoord(prevP.x, prevP.y);
			// newP.y = canvasHeight - newP.y;
			// prevP.y = canvasHeight - prevP.y;
			fl_line(newP.x, newP.y, prevP.x, prevP.y);
		}
	}
}


// TODO: Should be in a separate animation class
void Canvas::drawAnimation() {
	if (!scenario.aniStart || !GUI::drawSignal) return;
	auto now = std::chrono::system_clock::now().time_since_epoch();
	int mili = std::chrono::duration_cast<std::chrono::milliseconds>(now).count();
	float timeElapsed = ((float)mili / 1000) - scenario.timer;
	
	if (aniMode == DRONEANI) drawLineAnis(scenario.droneAnis, timeElapsed);
	if (aniMode == PACKAGEANI) drawLineAnis(scenario.packageAnis, timeElapsed);
}


void Canvas::drawGridLines() {

}


void Canvas::drawGridPoints() {
	if (!GUI::drawGridSignal) return;
	for (int i = 0; i < scenario.points.size(); i++) {
		if (scenario.points[i].isDesignatedPoint) continue;
		
		fl_color(fl_rgb_color(gridColor[i].r, gridColor[i].g, gridColor[i].b));
		float x = scenario.points[i].p.x;
		float y = scenario.points[i].p.y;
		scenToCanvasCoord(x, y);
		fl_pie(((int)x) - 2.5, ((int)y) - 2.5, 5, 5, 0, 360);
	}
}


void Canvas::draw() {

	if (GUI::bigRedrawSignal) {
		// Not sure how do you clear drawn objects in FLTK so I just draw a white rectangle instead
		fl_color(255);
		fl_rectf(0, 0, Canvas::canvasWidth + 100, Canvas::canvasHeight + 100);
		GUI::bigRedrawSignal = 0;
	}

	// Let group draw itself
	Fl_Group::draw();
	{
		static int redraws = 0;
		fl_color(FL_BLACK);
		fl_font(FL_COURIER, 80);
	}

	drawGridPoints();
	drawObstacles();
	drawAgents();
	drawPackages();
	drawTargets();
	drawAnimation();
	
}


Canvas::Canvas(int X, int Y, int W, int H, const char *L = 0) : Fl_Group(X, Y, W, H, L) {
	canvas = new float**[H];

	// Create a white background
	imageBox = new Fl_Box(X, Y, W, H, "");
	imageBox->color(255);
	color(255);

	axisLength = ((W < H ? W : H)) * 0.9;
	origin = ((W < H ? W : H)) * 0.05;
	squareLength = ((W < H ? W : H)) * 0.8;

	Fl::add_timeout(0.001, timerCallback, (void*)this);

	aniMode = DRONEANI;
	gridVisMode = UNDEFINED;
}


Canvas::~Canvas() {
	delete imageBox;
}


GUI::GUI(int winWidth, int winHeight) {
	// Init the main window
	win = new Fl_Window(winWidth, winHeight);
	win->color(255);

	// Init the top menu bar
	int menuBarWidth = 1080;
	int menuBarHeight = 25;

	// Init the canvas for drawing / displaying points
	Canvas::canvasWidth = winHeight - menuBarHeight - 50;
	Canvas::canvasHeight = winHeight - menuBarHeight - 50;
	Canvas::canvasX = 25;
	Canvas::canvasY = menuBarHeight + 25;
	Canvas::marginScale = 0.96;
	Canvas::margin = 7;

	// Init the width/height unit of the buttons
	int yButtonUnit = winHeight / 16;
	int xButtonUnit = winWidth / 8;

	// Init the menu bar and create related menu options
	Fl_Menu_Bar *menu = new Fl_Menu_Bar(0, 0, winWidth, menuBarHeight);
	menu->add("File/Open", FL_CTRL + 'o', browseCallback);
	menu->add("File/Save", FL_CTRL + 's', saveResultCallback);

	// Assign callbacks to corresponding buttons
	runAllBu = new Fl_Button(Canvas::canvasWidth + xButtonUnit, menuBarHeight + yButtonUnit, 160, yButtonUnit/2, "Solve");
	runAllBu->callback(solverCallback);

	droneAniBu = new Fl_Button(Canvas::canvasWidth + xButtonUnit, menuBarHeight + yButtonUnit * 2, 160, yButtonUnit / 2, "Animate");
	droneAniBu->callback(drawSignalCallback);

	// Group of radio buttons which allow user to choose animation mode [droneani]
	Fl_Group* aniOptionGroup = new Fl_Group(Canvas::canvasWidth + xButtonUnit, menuBarHeight + yButtonUnit * 2.8, 180, 50); 
	aniOptionGroup->box(FL_UP_FRAME);

	droneAniOptBu = new Fl_Round_Button(Canvas::canvasWidth + xButtonUnit + 5, menuBarHeight + yButtonUnit * 3, 10, 10, "Animate drones path");
	droneAniOptBu->type(102);
	droneAniOptBu->down_box(FL_ROUND_DOWN_BOX);
	droneAniOptBu->callback(aniRadioCallback, (void*)DRONEANI);
	
	packageAniOptBu = new Fl_Round_Button(Canvas::canvasWidth + xButtonUnit + 5, menuBarHeight + yButtonUnit * 3.5, 10, 10, "Animate package path");
	packageAniOptBu->type(102);
	packageAniOptBu->down_box(FL_ROUND_DOWN_BOX);
	packageAniOptBu->callback(aniRadioCallback, (void*)PACKAGEANI);

	aniOptionGroup->end();

	gridVisBu = new Fl_Button(Canvas::canvasWidth + xButtonUnit, menuBarHeight + yButtonUnit * 4.5, 160, yButtonUnit / 2, "Visualiza grid");
	gridVisBu->callback(drawGridSignalCallback);

	// Group of radio buttons which allow user to choose grid visualization mod [gridvis]
	Fl_Group* gridVisOptionGroup = new Fl_Group(Canvas::canvasWidth + xButtonUnit, menuBarHeight + yButtonUnit * 5.3, 180, 75);
	gridVisOptionGroup->box(FL_UP_FRAME);

	shortestPathMapOptBu = new Fl_Round_Button(Canvas::canvasWidth + xButtonUnit + 5, menuBarHeight + yButtonUnit * 5.5, 10, 10, "Shortest path map");
	shortestPathMapOptBu->type(102);
	shortestPathMapOptBu->down_box(FL_ROUND_DOWN_BOX);
	shortestPathMapOptBu->callback(gridVisRadioCallback, (void*)STPMAP);

	droneUsageMapOptBu = new Fl_Round_Button(Canvas::canvasWidth + xButtonUnit + 5, menuBarHeight + yButtonUnit * 6, 10, 10, "Depot usage map");
	droneUsageMapOptBu->type(102);
	droneUsageMapOptBu->down_box(FL_ROUND_DOWN_BOX);
	droneUsageMapOptBu->callback(gridVisRadioCallback, (void*)DEPOTUSAGEMAP);

	depotUsageMapOptBu = new Fl_Round_Button(Canvas::canvasWidth + xButtonUnit + 5, menuBarHeight + yButtonUnit * 6.5, 10, 10, "Shortest time map");
	depotUsageMapOptBu->type(102);
	depotUsageMapOptBu->down_box(FL_ROUND_DOWN_BOX);
	depotUsageMapOptBu->callback(gridVisRadioCallback, (void*)STIMEMAP);
	gridVisOptionGroup->end();

	// TODO: replace this with actual mouse zooming function
	zoomBu = new Fl_Button(Canvas::canvasWidth + xButtonUnit, menuBarHeight + yButtonUnit * 7.5, 160, yButtonUnit / 2, "Zoom in/out");
	zoomBu->callback(zoomToggleCallback);

	// Create  the actual canvas
	canvas = new Canvas(Canvas::canvasX, Canvas::canvasY, Canvas::canvasWidth, Canvas::canvasHeight, 0);

	// Init the color palette
	Palette::createPalette();

	win->resizable(win);
	win->show();
}


GUI::~GUI() {
	delete canvas;
	delete runAllBu;
	delete clearBu;
	delete menuBar;
	delete win;
}