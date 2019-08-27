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


Agent::Agent(int _x, int _y, float _v) {
	loc0.x = _x;
	loc0.y = _y;
	loc = loc0;
	v = _v;

	delay = 0;

	maxFuel = 0;
	fuel = 0;

	isAvail = true;
}


Agent::Agent(Point2D _loc0, float _v) {
	loc0 = _loc0;
	loc = loc0;
	v = _v;

	delay = 0;

	maxFuel = 0;
	fuel = 0;

	isAvail = true;
}


Agent::Agent(int _x, int _y, float _v, float _delay) {
	loc0.x = _x;
	loc0.y = _y;
	loc = loc0;
	v = _v;

	delay = _delay;

	maxFuel = 0;
	fuel = 0;

	isAvail = true;
}


Agent::Agent(Point2D _loc0, float _v, float _delay) {
	loc0 = _loc0;
	loc = loc0;
	v = _v;

	delay = _delay;

	maxFuel = 0;
	fuel = 0;

	isAvail = true;
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