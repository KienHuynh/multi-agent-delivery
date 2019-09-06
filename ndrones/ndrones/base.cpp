#include "base.h"


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


PointState::PointState(Point2D _p) {
	p = _p;
	bestTime = -1;
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