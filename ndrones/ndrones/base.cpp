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


float Point2D::l2norm(Point2D a) {
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


LinkedPoint2D::LinkedPoint2D(Point2D p) {
	x = p.x;
	y = p.y;
	prev = NULL;
	next = NULL;
}


LinkedPoint2D::LinkedPoint2D(float _x, float _y) {
	x = _x;
	y = _y;
	prev = NULL;
	next = NULL;
}


void LinkedPoint2D::operator=(Point2D const &p) {
	x = p.x;
	y = p.y;
	prev = NULL;
	next = NULL;
}


SimplePolygon::SimplePolygon(std::vector<LinkedPoint2D> _points) {
	points = _points;

	points[0].prev = points.size() - 1;
	points[0].next = 1;
	for (int i = 1; i < points.size(); i++) {
		points[i].prev = i - 1;
		points[i].next = (i + 1) % points.size();
	}
}


SimplePolygon::SimplePolygon(std::vector<Point2D> _points) {
	std::vector<LinkedPoint2D> test(_points.begin(), _points.end());
	points = test;
	points[0].prev = points.size() - 1;
	points[0].next = 1;
	for (int i = 1; i < points.size(); i++) {
		points[i].prev = i - 1;
		points[i].next = (i + 1) % points.size();
	}
}


std::vector<Point2D> SimplePolygon::findCVHull() {
	return std::vector<Point2D>();
}


float SimplePolygon::Area2(int a, int b, int c) {
	return (points[b].x - points[a].x)*(points[c].y - points[a].y) -
		(points[c].x - points[a].x)*(points[b].y - points[a].y);
}


bool SimplePolygon::left(int a, int b, int c) {
	return Area2(a, b, c) > 0;
}


bool SimplePolygon::leftOn(int a, int b, int c) {
	return Area2(a, b, c) >= 0;
}


bool SimplePolygon::colinear(int a, int b, int c) {
	return Area2(a, b, c) == 0;
}


bool SimplePolygon::inCone(int a, int b) {
	int n = points.size();
	int a0 = points[a].prev;
	int a1 = points[a].next;
	
	if (leftOn(a, a1, a0)) {
		return left(a, b, a0) && left(b, a, a1);
	}

	return !(leftOn(a, b, a1) && leftOn(b,a,a0));
}
	

bool SimplePolygon::between(int a, int b, int c) {
	int ba, ca;
	if (!colinear(a, b, c)) return false;

	if (points[a].x != points[b].x)
		return
		((points[a].x <= points[c].x) && (points[c].x <= points[b].x)) ||
		((points[a].x >= points[c].x) && (points[c].x >= points[b].x));
	else
		return
		((points[a].y <= points[c].y) && (points[c].y <= points[b].y)) ||
		((points[a].y >= points[c].y) && (points[c].y >= points[b].y));
}


bool SimplePolygon::intersectProp(int a, int b, int c, int d) {
	if (colinear(a, b, c) ||
		colinear(a, b, d) ||
		colinear(c, d, a) ||
		colinear(c, b, d)) return false;

	return 
		(!(left(a, b, c)) ^ !(left(a, b, d))) &&
		(!(left(c, d, a)) ^ !(left(c, d, b)));
}


bool SimplePolygon::intersect(int a, int b, int c, int d) {
	if (intersectProp(a, b, c, d)) {
		return true;
	}
	else if (
		between(a, b, c) || between(a, b, d) ||
		between(c, d, a) || between(c, d, b)
		)
		return true;

	return false;
}


bool SimplePolygon::diagonalie(int a, int b) {
	int c = 0, c1;
	do {
		c1 = points[c].next;
		if ((c != a) && (c1 != a) &&
			(c != b) && (c1 != b) &&
			intersect(a, b, c, c1)) return false;
		c = points[c].next;
	} while (c != 0);
	return true;
}


bool SimplePolygon::diagonal(int a, int b) {
	return inCone(a, b) && inCone(b, a) && diagonalie(a, b);
}


void SimplePolygon::earInit() {
	int v0, v1, v2;
	v1 = 0;
	do {
		v2 = points[v1].next;
		v0 = points[v1].prev;
		points[v1].isEar = diagonal(v0, v2);
		v1 = points[v1].next;
	} while (v1 != 0);
}


void SimplePolygon::triangulate() {
	int v0, v1, v2, v3, v4;
	std::vector<LinkedPoint2D> pointsCopy = points;
	int n = pointsCopy.size();
	int start = 0;
	while (n > 3) {
		v2 = start;
		do {
			if (pointsCopy[v2].isEar) {
				// v1, v2, v3 is a triangle


				v3 = pointsCopy[v2].next;
				v4 = pointsCopy[v3].next;
				v1 = pointsCopy[v2].prev;
				v0 = pointsCopy[v1].prev;
				
				pointsCopy[v1].isEar = diagonal(v0, v3);
				pointsCopy[v3].isEar = diagonal(v1, v4);
				

				pointsCopy[v1].next = v3;
				pointsCopy[v3].prev = v1;
				
				start = v3;
				n--;
				break;
			}
			v2 = pointsCopy[v2].next;
		} while (v2 != start);
	}
}