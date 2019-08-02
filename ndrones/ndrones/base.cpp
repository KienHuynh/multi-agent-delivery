#include "base.h"


Point2D::Point2D() {
	x = -1;
	y = -1;
}


Point2D::Point2D(float _x, float _y) {
	x = _x;
	y = _y;
}


PointState::PointState(Point2D _p) {
	p = _p;
	bestTime = -1;
}


Package::Package(Point2D _loc0) {
	loc0 = _loc0;
	loc = loc0;
}


Target::Target(Point2D _loc) {
	loc = _loc;
}


Agent::Agent(int _x, int _y, float _v) {
	loc0.x = _x;
	loc0.y = _y;
	loc = loc0;
	v = _v;

	maxFuel = 0;
	fuel = 0;
}


Agent::Agent(Point2D _loc0, float _v) {
	loc0 = _loc0;
	loc = loc0;
	v = _v;

	maxFuel = 0;
	fuel = 0;
}


void Scenario::loadFile(const char* fname) {
	std::ifstream myfile;
	myfile.open(fname, std::ios::in);

	int minX = 0, minY = 0, maxX = 0, maxY = 0;
	int nPackage = 0, nTarget = 0, nAgent = 0;

	myfile >> minX >> maxX >> minY >> maxY;
	for (int i = minX; i < maxX; i++) {
		for (int j = minY; j < maxY; j++) {
			points.push_back(Point2D(i, j));
		}
	}

	myfile >> nPackage;
	for (int i = 0; i < nPackage; i++) {
		int x, y;
		myfile >> x >> y;
		packages.push_back(Package(Point2D(x, y)));
	}

	myfile >> nTarget;
	for (int i = 0; i < nTarget; i++) {
		int x, y;
		myfile >> x >> y;
		targets.push_back(Target(Point2D(x, y)));
	}

	myfile >> nAgent;
	maxSpeed = -1;
	minSpeed = -1;
	for (int i = 0; i < nAgent; i++) {
		int x, y;
		float v;
		myfile >> x >> y >> v;
		if (maxSpeed == -1) {
			maxSpeed = v;
			minSpeed = v;
		}
		if (v > maxSpeed) {
			maxSpeed = v;
		}
		if (v < minSpeed) {
			minSpeed = v;
		}
		//Point2D p(x, y);
		//Agent a(p, v);
		Agent a(x, y, v);
		agents.push_back(a);
	}
}


void Scenario::ecld2DDynamicSolve() {

}


void Scenario::solve(Solver solver) {
	if (targets.size() == 1) {
		if (packages.size() == 1) {
			if (solver == Solver::ECLD_2D_DYNAMIC) {
				ecld2DDynamicSolve();
			}
		}
	}
}