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


class Agent;
class LinkedPoint2D;

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
	static float l2norm(Point2D);

	Point2D operator + (Point2D const &obj);
	Point2D operator - (Point2D const &obj);
	Point2D operator / (float const);
	Point2D operator * (float const);
	void operator = (LinkedPoint2D const &obj);
	bool operator == (Point2D const &obj);
	bool operator == (Point2D &obj);

	// Compute 2 * the area of the triangle a, b, c
	// @param[in] Point2D a
	// @param[in] Point2D b
	// @param[in] Point2D c
	// @param[out] float
	static float Area2(Point2D a, Point2D b, Point2D c);

	// Return true of c is to the left of (a,b)
	// @param[in] Point2D a
	// @param[in] Point2D b
	// @param[in] Point2D c
	// @param[out] bool
	static bool left(Point2D a, Point2D b, Point2D c);

	// Return true of c is to the left of or on (a,b)
	// The "On" here is basically infeasible to compute because all coordinates are float
	// But this is still a different function from left( ) just in case
	// @param[in] Point2D a
	// @param[in] Point2D b
	// @param[in] Point2D c
	// @param[out] bool
	static bool leftOn(Point2D a, Point2D b, Point2D c);

	// Check if three points are colinear
	// @param[in] Point2D a
	// @param[in] Point2D b
	// @param[in] Point2D c
	// @param[out] bool
	static bool colinear(Point2D a, Point2D b, Point2D c);

	// Check if c lies on (a, b) and is between a and b
	// @param[in] Point2D a
	// @param[in] Point2D b
	// @param[in] Point2D c
	// @param[out] bool
	static bool between(Point2D a, Point2D b, Point2D c);

	// Check if two segments 
	// @param[in] Point2D a
	// @param[in] Point2D b
	// @param[in] Point2D c
	// @param[in] Point2D d
	// @param[out] bool
	static bool intersectProp(Point2D a, Point2D b, Point2D c, Point2D d);

	// Check if (a, b) intersects with (c, d)
	// @param[in] Point2D a
	// @param[in] Point2D b
	// @param[in] Point2D c
	// @param[in] Point2D d
	// @param[out] bool
	static bool intersect(Point2D a, Point2D b, Point2D c, Point2D d);
};


// This class is used to store point state such as best time / agent queue, which is needed to store the solution
class PointState {
public:
	Point2D p;

	// Is true if this is a DesignatedPoint (package, target, etc.). Assigned when loading file.
	bool isDesignatedPoint;

	// Is on obstacle boundary
	bool isOb;

	// For shortest path calculation
	bool visited;
	float dist;
	int prev;

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
	// @param[in] Point2D a
	// @param[out] float
	float timing(Point2D);

	// Compute the time to go from point i to j
	// @param[in] Point2D a
	// @param[in] Point2D b
	// @param[out] float
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


class LinkedPoint2D : public Point2D {
public:
	LinkedPoint2D(Point2D p);
	LinkedPoint2D(float x, float y);
	int prev, next;
	bool isEar;
	bool hullPoint;
	void operator = (Point2D const &obj);
};


enum PointLocation {
	INPOLY,
	INPOCKET,
	OUTOFHULL
};


// Simple polygon class
class SimplePolygon {
public:
	SimplePolygon(std::vector<LinkedPoint2D>);
	SimplePolygon(std::vector<Point2D>);
	
	std::vector<LinkedPoint2D> points;
	float perimeter;
	// Hull point index, storing indices of hull points
	std::vector<int> hullPtIdx;
	// Triangle index, storing indices of all triangles after triangulation
	std::vector<std::vector<int>> triIdx;
	// Pocket triangle index, storing indices of all triangles after triangulation of each pocket
	std::vector<std::vector<int>> pocketTriIdx;

	// Using Melkman's algorithm
	void findCVHull();

	// Find the pockets
	void findPocketTriangles();

	// Check if point is inside polygon
	// @param[in] Point2D p
	// @param[out] bool
	bool contain(Point2D p);

	// Check if a line segment intersect the polygon
	// @param[in] Point2D a
	// @param[in] Point2D b
	bool segIntersect(Point2D a, Point2D b);

	// Return true if (a, b) lies inside the polygon
	// @param[in] Point2D a
	// @param[in] Point2D b
	// @param[out] bool
	bool inCone(int a, int b);

	// Check if (a, b) cuts any edge of the polygon
	// Return true if (a, b) lies inside the polygon
	// @param[in] int a, index of the point a in this.points
	// @param[in] int b, index of the point b in this.points
	// @param[out] bool
	bool diagonalie(int a, int b);

	// Return true if (a, b) is a diagonal
	// Return true if (a, b) lies inside the polygon
	// @param[in] int a, index of the point a in this.points
	// @param[in] int b, index of the point b in this.points
	// @param[out] bool
	bool diagonal(int a, int b);

	// Assign isEar = true for all vertices that are ear
	void earInit();

	// Triangulation for drawing (ear-clipping)
	void triangulate();
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