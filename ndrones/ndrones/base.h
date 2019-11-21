/*
 * base.h
 * This file contains basic classes and entities needed for the main solvers
*/

#pragma once
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <iostream>
#include <chrono>

#include "config.h"
#include "util.h"

constexpr auto PI = 3.14156;

enum Solver {ECLD_2D_DYNAMIC};


// The main for Agent
class Agent;


// Class to manage color, value range is 0-255
class Color {
public:
	unsigned char r;
	unsigned char g;
	unsigned char b;
	Color();
	Color(unsigned char, unsigned char, unsigned char);
	//static Color HSV2RGB(float h, float s, float v);
	static Color HSV2RGB(float h, float s, float v);
};


class Palette {
public:
	static Color palette[cfg::numColor];
	static Color huePalette[cfg::numColor];
	static void createPalette();
};


// Basic 2D point class
class Point2D {
public:
	float x;
	float y;

	Point2D();
	Point2D(float, float);

	void copyFrom(Point2D);

	// Compute L2 distance
	static float l2Distance(Point2D, Point2D);
	static float abs(Point2D);

	Point2D operator + (Point2D const &obj);
	Point2D operator - (Point2D const &obj);
	Point2D operator / (float const);
	Point2D operator * (float const);
	bool operator == (Point2D const &obj);
	bool operator == (Point2D &obj);
};


// This class is used to store point state such as best time / agent queue, which is needed to store the solution
class PointState {
public:
	Point2D p;

	// Is true if this is a DesignatedPoint (package, target, etc.). Assigned when loading file.
	bool isDesignatedPoint;

	std::vector<Agent> agentQueue;
	// Queue of points that were previously visited by the above agents
	// Example:
	// Agent a1->a2
	// PQueue p1->p2
	// This means that a1 visited the package at p1
	// a1 then moved to p2 to meet a2
	// a2 then moved to *this* point to meet some other drones
	std::vector<Point2D> pointQueue;
	float bestTime;

	PointState(Point2D _p);

	// Compare if the bestTime (so far) of point a is less than point b
	static bool bestTimeLT(PointState, PointState);
};


class Agent {
public:
	int ID;
	int orderOfEx;
	int prevMatching;
	int currentMatching;
	int nextMatching;

	// Initial location
	Point2D loc0;
	// Current location
	Point2D loc;

	// Velocity
	float v;

	// Initial delay time
	float delay;

	// Finishing time
	float finTime;

	// Max amount of fuel
	float maxFuel;
	// Current fuel
	float fuel;

	// Occupied by a task or not
	bool isAvail;

	Agent(int _ID, int _x, int _y, float _v);
	Agent(int _ID, Point2D _loc0, float _v);
	Agent(int _ID, int _x, int _y, float _v, float _delay);
	Agent(int _ID, Point2D, float _v, float _delay);

	// Compute the time to go to point p based on current loc and velocity
	float timing(Point2D);
	// Compute the time to go from point i to j
	float timing(Point2D, Point2D);

	bool operator < (Agent const &obj);
	bool operator == (Agent const &obj);
};


// Specific points that the agents need to travel from/to
class DesignatedPoint {
public:
	// Staring location
	Point2D loc;
	Point2D currentLoc;
	// ID of this point, used for package-target matching
	int ID;
	// The id reference of this point in the main grid (i.e. the member std::vector<PointState> points in scenario)
	int gridRef;

	DesignatedPoint();
	DesignatedPoint(Point2D);
	DesignatedPoint(Point2D, int);
};


// Basic line animation
// This works by specifying a start time and end time of an animation, along with a start point and end point
// Each time the clock tick, the drawing function will interpolate a line segment between the start and end
// There can be previous animation to an animation, which makes the current animation waits until the previous one is done
class LineAnimation {
public:
	Point2D start;
	Point2D end;
	int color[3];
	float startTime;
	float endTime;
	float duration;

	// Store the previous timer, used to calculate the current line segment (span from prevPoint to currentPoint)
	float prevTimer;
	// This is a list of previous animations this animation has to wait
	std::vector<int> prevAni;
	// Specify if the current animation is active or not
	bool active;

	LineAnimation();
	void setColor(int, int, int);
};