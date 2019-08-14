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
	Point2D loc;
	Point2D currentLoc;
	DesignatedPoint(Point2D);
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

	std::vector<PointState> points;
	std::vector<DesignatedPoint> packages;
	// Store the id of the point of the package
	std::vector<int> packageIdx;
	std::vector<DesignatedPoint> targets;
	std::vector<int> targetIdx;

	Scenario();

	// Load scenario from file
	void loadFile(const char*);

	// Check if a point is a package
	// TODO: use better data structure?
	bool isPackage(int);
	
	// Euclidean 2D dynamic solver for 1-drone-1-package
	void ecld2DDynamicSolve11();
	void ecld2DDynamicSolveNN();

	// Create an animation based on the solution
	void createAnimation();
	void solve();

private:
	// This is used to load target points or package points
	void loadDesignatedPoint(std::ifstream &myfile, std::vector<DesignatedPoint> &dPoints, int nDPoint, std::vector<int> &idx,
		std::vector<PointState> &points);
	void loadDesignatedPolygon(std::ifstream &myfile, std::vector<DesignatedPoint> &dPoints, int nDPoint, std::vector<int> &idx,
		std::vector<PointState> &points);
};