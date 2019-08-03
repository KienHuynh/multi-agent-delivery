#pragma once
#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <iostream>

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

	Point2D operator - (Point2D const &obj);
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


class Package {
public:
	Point2D loc0;
	Point2D loc;

	Package(Point2D);
};


class Target {
public:
	Point2D loc;
	
	Target(Point2D);
};


class Scenario {
public:
	std::vector<Agent> agents;
	float maxSpeed, minSpeed;
	int minX, maxX, minY, maxY;

	std::vector<PointState> points;
	std::vector<Package> packages;
	// Store the id of the point of the package
	std::vector<int> packageIdx;
	std::vector<Target> targets;
	std::vector<int> targetIdx;

	// Load scenario from file
	void loadFile(const char*);

	// Euclidean 2D dynamic solver for 1-drone-1-package
	void ecld2DDynamicSolve11();
	void solve(Solver);
};

