#include "base.h"


// Color
Color::Color() {
	r = 255;
	g = 255;
	b = 255;
}


Color::Color(unsigned char _r, unsigned char _g, unsigned char _b) {
	r = _r;
	g = _g;
	b = _b;
}


Color Color::HSV2RGB(float h, float s, float v) {
	float c = s * v;
	float h_ = h / 60.0;
	h_ = h_ < 2 ? h_ : h_ - (float) 2*(((int)h_)/2);
	float x = c * (1.0 - abs(h_ - 1.0));
	float m = v - c;
	float r, g, b;
	if (h < 60) {
		r = c;
		g = x;
		b = 0;
	}
	else if (h < 120) {
		r = x;
		g = c;
		b = 0;
	}
	else if (h < 180) {
		r = 0;
		g = c;
		b = x;
	}
	else if (h < 240) {
		r = 0;
		g = x;
		b = c;
	}
	else if (h < 300) {
		r = x;
		g = 0;
		b = c;
	}
	else if (h < 360) {
		r = c;
		g = 0;
		b = x;
	}
	r = (r + m) * 255;
	g = (g + m) * 255;
	b = (b + m) * 255;

	return Color((unsigned char) r, (unsigned char) g, (unsigned char) b);
}


// Palette
Color Palette::palette[cfg::numColor];
Color Palette::huePalette[cfg::numColor];

void Palette::createPalette() {
	// Create the random palette
	for (int i = 0; i < cfg::numColor; i++) {
		float h = ((float)i) * 360.0 / (float)cfg::numColor;
		float s = ((i + 1) % 2 == 0) ? 1 : 0.7;
		float v = (i %  2 == 0) ? 0.9 : 0.4;
		Color c = Color::HSV2RGB((int)h, s, v);
		palette[i] = c;
	}

	// Shuffle them up
	for (int i = 0; i < cfg::numColor; i++) {
		srand(i);
		int j = rand() % cfg::numColor;
		Color c = palette[i];
		palette[i] = palette[j];
		palette[j] = c;
	}

	// Create the hue palette
	for (int i = 0; i < cfg::numColor; i++) {
		float h = ((float)i) * 360.0 / (float)cfg::numColor;
		float s = 0.8;
		float v = 0.7;
		Color c = Color::HSV2RGB((int)h, s, v);
		huePalette[i] = c;
	}
}


//TODO: Separate basic classes from scenario
Point2D::Point2D() {
	x = -1;
	y = -1;
}


Point2D::Point2D(float _x, float _y) {
	x = _x;
	y = _y;
}


float Point2D::l2Distance(Point2D a, Point2D b) {
	float tmp = (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y);
	return sqrt(tmp);
}


float Point2D::abs(Point2D a) {
	return sqrt(a.x*a.x + a.y*a.y);
}


Point2D Point2D::operator+(Point2D const &p) {
	Point2D tmp;
	tmp.x = x + p.x;
	tmp.y = y + p.y;
	return tmp;
}


Point2D Point2D::operator-(Point2D const &p) {
	Point2D tmp;
	tmp.x = x - p.x;
	tmp.y = y - p.y;
	return tmp;
}


Point2D Point2D::operator/(float const f) {
	Point2D tmp;
	tmp.x = x / f;
	tmp.y = y / f;
	return tmp;
}


Point2D Point2D::operator*(float const f) {
	Point2D tmp;
	tmp.x = x * f;
	tmp.y = y * f;
	return tmp;
}


bool Point2D::operator == (Point2D &obj) {
	if (x == obj.x && y == obj.y) return true;
	else return false;
}


bool Point2D::operator == (Point2D const &obj) {
	if (x == obj.x && y == obj.y) return true;
	else return false;
}


bool PointState::bestTimeLT(PointState a, PointState b) {
	return a.bestTime < b.bestTime;
}


PointState::PointState(Point2D _p) {
	p = _p;
	bestTime = -1;
	isDesignatedPoint = false;
}


DesignatedPoint::DesignatedPoint(Point2D _p) {
	loc = _p;
	currentLoc = _p;
	ID = -1;
}


DesignatedPoint::DesignatedPoint(Point2D _p, int _id) {
	loc = _p;
	currentLoc = _p;
	ID = _id;
}


DesignatedPoint::DesignatedPoint() {

}


Agent::Agent(int _ID, int _x, int _y, float _v) {
	ID = _ID;
	loc0.x = _x;
	loc0.y = _y;
	loc = loc0;
	v = _v;

	delay = 0;
	finTime = 0;

	maxFuel = 0;
	fuel = 0;

	isAvail = true;
	orderOfEx = 0;
	prevMatching = -1;
	currentMatching = -1;
	nextMatching = -1;
}


Agent::Agent(int _ID, Point2D _loc0, float _v) {
	ID = _ID;
	loc0 = _loc0;
	loc = loc0;
	v = _v;

	delay = 0;
	finTime = 0;

	maxFuel = 0;
	fuel = 0;

	isAvail = true;
	orderOfEx = 0;
	prevMatching = -1;
	currentMatching = -1;
	nextMatching = -1;
}


Agent::Agent(int _ID, int _x, int _y, float _v, float _delay) {
	ID = _ID;
	
	loc0.x = _x;
	loc0.y = _y;
	loc = loc0;
	v = _v;

	delay = _delay;
	finTime = 0;

	maxFuel = 0;
	fuel = 0;

	isAvail = true;
	orderOfEx = 0;
	prevMatching = -1;
	currentMatching = -1;
	nextMatching = -1;
}


Agent::Agent(int _ID, Point2D _loc0, float _v, float _delay) {
	ID = _ID;
	
	loc0 = _loc0;
	loc = loc0;
	v = _v;

	delay = _delay;
	finTime = 0;

	maxFuel = 0;
	fuel = 0;

	isAvail = true;
	orderOfEx = 0;
	prevMatching = -1;
	currentMatching = -1;
	nextMatching = -1;
}


float Agent::timing(Point2D p) {
	float d = Point2D::l2Distance(p, loc);
	return d / v;
}


float Agent::timing(Point2D i, Point2D j) {
	float d = Point2D::l2Distance(i, j);
	return d / v;
}


bool Agent::operator<(Agent const &obj) {
	return v < obj.v;
}


bool Agent::operator==(Agent const &obj) {
	if (loc0 == obj.loc0 &&
		loc == obj.loc &&
		v == obj.v &&
		delay == obj.delay &&
		maxFuel == obj.maxFuel &&
		fuel == obj.fuel &&
		isAvail == obj.isAvail) return true;

	return false;
}