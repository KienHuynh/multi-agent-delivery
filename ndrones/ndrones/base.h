#pragma once
#include <vector>
#include <fstream>

enum Solver {ECLD_2D_DYNAMIC};

class Agent;

class Point2D {
public:
	float x;
	float y;

	Point2D();
	Point2D(float, float);

	void copy(Point2D);
};


class PointState {
public:
	Point2D p;
	std::vector<Agent> agentQueue;
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
	std::vector<Point2D> points;
	std::vector<Package> packages;
	std::vector<Target> targets;

	void loadFile(const char*);

	void ecld2DDynamicSolve();
	void solve(Solver);
};

