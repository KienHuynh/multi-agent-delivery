#include "base.h"


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


Point2D Point2D::operator-(Point2D const &p) {
	Point2D tmp;
	tmp.x = x - p.x;
	tmp.y = y - p.y;
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


void Scenario::loadFile(const char* fname) {
	std::ifstream myfile;
	myfile.open(fname, std::ios::in);

	int nPackage = 0, nTarget = 0, nAgent = 0;
	int stepX, stepY;
	myfile >> minX >> maxX >> stepX >> minY >> maxY >> stepY;
	// TODO: Make this more efficient
	for (int i = minX; i < maxX; i += stepX) {
		for (int j = minY; j < maxY; j += stepY) {
			points.push_back(PointState(Point2D(i, j)));
		}
	}

	myfile >> nPackage;
	for (int i = 0; i < nPackage; i++) {
		int x, y;
		myfile >> x >> y;
		packages.push_back(Package(Point2D(x, y)));
		bool in = false;
		for (int p = 0; p < points.size(); p++) {
			Point2D tmp = points[p].p;
			if (tmp.x == packages[packages.size() - 1].loc0.x &&
				tmp.y == packages[packages.size() - 1].loc0.y) {
				packageIdx.push_back(p);
				in = true;
			}
		}
		if (in == false) {
			points.push_back(PointState(packages[packages.size()-1].loc0));
			packageIdx.push_back(points.size() - 1);
		}
	}

	myfile >> nTarget;
	for (int i = 0; i < nTarget; i++) {
		int x, y;
		myfile >> x >> y;
		targets.push_back(Target(Point2D(x, y)));

		bool in = false;
		for (int p = 0; p < points.size(); p++) {
			Point2D tmp = points[p].p;
			if (tmp.x == targets[targets.size() - 1].loc.x &&
				tmp.y == targets[targets.size() - 1].loc.y) {
				targetIdx.push_back(p);
				in = true;
			}
		}
		if (in == false) {
			points.push_back(PointState(targets[targets.size() - 1].loc));
			targetIdx.push_back(points.size() - 1);
		}
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
	// Sort the agents by speed
	std::sort(agents.begin(), agents.end());
}


float max(float x, float y) {
	if (x > y) return x;
	else return y;
}


void Scenario::ecld2DDynamicSolve11() {
	// Find the drone that can get to the package in the shortest time
	float bestTime = agents[0].timing(packages[0].loc0);
	int bestAgent = 0;
	for (int k = 1; k < agents.size(); k++) {
		float timeTmp = agents[k].timing(packages[0].loc0);
		if (timeTmp < bestTime) {
			bestTime = timeTmp;
			bestAgent = k;
		}
	}

	// Assign best time, agent queue and point queue to all points in the scenario
	for (int i = 0; i < points.size(); i++) {
		points[i].bestTime = agents[bestAgent].timing(points[i].p);
		points[i].agentQueue.push_back(agents[bestAgent]);
		points[i].pointQueue.push_back(packages[0].loc0);
	}

	// No point in handing the package to a slower drone
	for (int k = bestAgent+1; k < agents.size(); k++) {
		std::vector<PointState> pointsCopy = points;
		for (int i = 0; i < points.size(); i++) {
			// Find the best handoff point j so that it can travel to i faster
			int best_j = -1;
			if (i % 10 == 0) std::cout << i << std::endl;
			for (int j = 0; j < pointsCopy.size(); j++) {
				float timeTo_j = agents[k].timing(pointsCopy[j].p);
				// Factor waiting time
				timeTo_j = max(timeTo_j, pointsCopy[j].bestTime);

				float time_jTo_i = agents[k].timing(pointsCopy[j].p, pointsCopy[i].p);
				if (timeTo_j + time_jTo_i < pointsCopy[i].bestTime) {
					
					points[i].bestTime = timeTo_j + time_jTo_i;
					best_j = j;
				}
			}
			if (best_j > -1) {
				points[i].agentQueue = pointsCopy[best_j].agentQueue;
				points[i].agentQueue.push_back(agents[k]);
				points[i].pointQueue = pointsCopy[best_j].pointQueue;
				points[i].pointQueue.push_back(pointsCopy[best_j].p);
			}
		}
	}
}


void Scenario::solve(Solver solver) {
	if (targets.size() == 1) {
		if (packages.size() == 1) {
			if (solver == Solver::ECLD_2D_DYNAMIC) {
				ecld2DDynamicSolve11();
				int a = 1;
				a = a + 1;
			}
		}
	}
}