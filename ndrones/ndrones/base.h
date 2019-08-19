#pragma once
#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <iostream>
#include <chrono>

#include "config.h"

enum Solver {ECLD_2D_DYNAMIC};

class Agent;

class Point2D {
public:
	float x;
	float y;

	Point2D();
	Point2D(float, float);

	void copy(Point2D);
	// Compute L2 distance
	static float l2Distance(Point2D, Point2D);
	static float abs(Point2D);

	Point2D operator + (Point2D const &obj);
	Point2D operator - (Point2D const &obj);
	Point2D operator / (float const);
	Point2D operator * (float const);
};


class PointState {
public:
	Point2D p;
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
};


class Agent {
public:
	// Initial location
	Point2D loc0;
	// Current location
	Point2D loc;

	// Velocity
	float v;

	// Max amount of fuel
	float maxFuel;
	// Current fuel
	float fuel;

	Agent(int _x, int _y, float _v);
	Agent(Point2D _loc0, float _v);

	// Compute the time to go to point p based on current loc and velocity
	float timing(Point2D);
	// Compute the time to go from point i to j
	float timing(Point2D, Point2D);

	bool operator < (Agent const &obj);
};


class DesignatedPoint {
public:
	// Staring location
	Point2D loc;
	Point2D currentLoc;
	// ID of this point, used for package-target matching
	int ID;
	// The id reference of this point in the main grid (i.e. the member std::vector<PointState> points in scenario)
	int gridRef;

	DesignatedPoint(Point2D);
	DesignatedPoint(Point2D, int);
};


class LineAnimation {
public:
	Point2D start;
	Point2D end;
	int color[3];
	float startTime;
	float endTime;
	float duration;
	float prevTimer;
	bool active;

	LineAnimation();
	void setColor(int, int, int);
};


class Scenario {
public:
	// TODO: use bit flag
	// n drones deliver 1 package from point (s) to point (s)
	// n drones deliver m packages from point (s) to point (s), all designated points have the same roles
	// n drones deliver m packages from point (s) to point (s), each point is unique
	int problemType;
	int solverType;

	std::vector<Agent> agents;
	std::vector<std::vector<LineAnimation>> anis;

	float maxSpeed, minSpeed;
	int minX, maxX, minY, maxY;

	// Animation related members
	float timer;
	bool aniStart;

	// The grid
	std::vector<PointState> points;
	// Store the packages
	std::vector<DesignatedPoint> packages;
	// Store the id of the point of the package
	std::vector<DesignatedPoint> targets;

	Scenario();

	// Load scenario from file
	void loadFile(const char*);

	// Check if a point is a package
	// TODO: use better data structure?
	bool isPackage(int);
	
	// Euclidean 2D (problem) type 0, dynamic solver for 1 drone and 1 package
	// Approximation of opt.
	void ecld2DType0DynamicSolve11();
	// Euclidean 2D (problem) type 0, dynamic solver for n drones and m packages
	// Approximation of opt.
	void ecld2DType0DynamicSolveNM();
	// Euclidean 2D (problem) type 1, dynamic solver for n drones and m packages
	// Heuristic
	void ecld2DType1DynamicSolveNM();

	// Create an animation based on the solution
	void createAnimation();
	void solve();

private:
	// This is used to load target points or package points
	// Load points such as packages and targets
	// If a point already exists in the grid above, re-use it
	// Otherwise, add a new point to the grid
	// @param myfile ifstream object with the input file already loaded
	// @param problemType problem type, see inputDescription.txt
	// @param dPoints the array of DesignatedPoint to store the package/target locations
	// @param nDPoint number of designated points
	// @param points the point grid
	void loadDesignatedPoint(std::ifstream &myfile, int problemType, std::vector<DesignatedPoint> &dPoints, 
		int nDPoint, std::vector<PointState> &points);
	// Load package regions and target regions as polygons
	// These polygons will be discretized by their edges
	// TODO: Allow this function to add polygons with different IDs
	// @param myfile ifstream object with the input file already loaded
	// @param problemType problem type, see inputDescription.txt
	// @param dPoints the array of DesignatedPoint to store the package/target locations
	// @param nDPoint number of designated points
	// @param points the point grid
	void loadDesignatedPolygon(std::ifstream &myfile, int problemType, std::vector<DesignatedPoint> &dPoints, 
		int nDPoint, std::vector<PointState> &points);
};