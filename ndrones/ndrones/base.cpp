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


void Scenario::loadDesignatedPoint(std::ifstream &myfile, int problemType, std::vector<DesignatedPoint> &dPoints, 
	int nDPoint, std::vector<PointState> &points) {

	for (int i = 0; i < nDPoint; i++) {
		int x, y, id;
		if (problemType == 0) {
			myfile >> x >> y;
			dPoints.push_back(DesignatedPoint(Point2D(x, y)));
		}
		else if (problemType == 1) {
			myfile >> x >> y >> id;
			dPoints.push_back(DesignatedPoint(Point2D(x, y), id));
		}
		
		bool in = false;
		for (int p = 0; p < points.size(); p++) {
			Point2D tmp = points[p].p;
			// If the new designated point already exists in the grid
			// Do not add a new point into the grid
			if (tmp.x == dPoints[dPoints.size() - 1].loc.x &&
				tmp.y == dPoints[dPoints.size() - 1].loc.y) {
				dPoints[dPoints.size() - 1].gridRef = p;
				in = true;
				break;
			}
		}

		// Add a new point into the grid
		if (in == false) {
			points.push_back(PointState(dPoints[dPoints.size() - 1].loc));
			dPoints[dPoints.size() - 1].gridRef = points.size() - 1;
		}
	}
}


void Scenario::loadDesignatedPolygon(std::ifstream &myfile, int problemType, std::vector<DesignatedPoint> &dPoints, 
	int nDPoint, std::vector<PointState> &points) {
	std::vector<Point2D> vList;
	for (int i = 0; i < nDPoint; i++) {
		int x, y;
		myfile >> x >> y;
		vList.push_back(Point2D(x, y));
	}

	for (int i=0; i < vList.size(); i++) {
		Point2D start = vList[i % nDPoint];
		Point2D end = vList[(i+1) % nDPoint];
		// Interpolation
		for (int s = 0; s < cfg::polySamplingRate; s++) {
			float rate = ((float)s) / ((float)cfg::polySamplingRate);
			Point2D t1 = end - start;
			Point2D t2 = t1 * rate;
			Point2D middle = (end - start)*rate + start;

			dPoints.push_back(DesignatedPoint(middle));
			bool in = false;
			for (int p = 0; p < points.size(); p++) {
				Point2D tmp = points[p].p;
				if (tmp.x == dPoints[dPoints.size() - 1].loc.x &&
					tmp.y == dPoints[dPoints.size() - 1].loc.y) {
					dPoints[dPoints.size() - 1].gridRef = p;
					in = true;
				}
			}
			if (in == false) {
				points.push_back(PointState(dPoints[dPoints.size() - 1].loc));
				dPoints[dPoints.size() - 1].gridRef = points.size() - 1;
			}
		}
	}
}


void Scenario::loadFile(const char* fname) {

	std::ifstream myfile;
	myfile.open(fname, std::ios::in);
	
	myfile >> problemType;
	myfile >> solverType;

	// Specify whether the inputs will be in regional (polygon) format or discrete points format
	int packageInputMode = 0;
	int targetInputMode = 0;

	int nPackage = 0, nTarget = 0, nAgent = 0, nPVertex = 0;
	int stepX, stepY;
	myfile >> minX >> maxX >> stepX >> minY >> maxY >> stepY;
	// TODO: Make this more efficient
	// Add points to the data pool
	for (int i = minX; i < maxX; i += stepX) {
		for (int j = minY; j < maxY; j += stepY) {
			points.push_back(PointState(Point2D(i, j)));
		}
	}

	myfile >> packageInputMode;
	myfile >> nPackage;
	if (packageInputMode == 0) {
		loadDesignatedPoint(myfile, problemType, packages, nPackage, points);
	}
	else if (packageInputMode == 1) {
		loadDesignatedPolygon(myfile, problemType, packages, nPackage, points);
	}

	myfile >> targetInputMode;
	myfile >> nTarget;
	if (targetInputMode == 0) {
		loadDesignatedPoint(myfile, problemType, targets, nTarget, points);
	}
	else if (targetInputMode == 1) {
		loadDesignatedPolygon(myfile, problemType, targets, nTarget, points);
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


bool Scenario::isPackage(int id) {
	for (int ip = 0; ip < packages.size(); ip++) {
		if (id == packages[ip].gridRef) {
			return true;
		}
	}
	return false;
}


void Scenario::ecld2DType0DynamicSolveNM() {
	// First run with the slowest drone
	// Compute the best time it takes for this drone to get to a package and fly to any other point on the grid
	for (int i = 0; i < points.size(); i++) {
		if (isPackage(i)) continue;

		float bestTime = -1;
		int bestPackageID = packages[0].gridRef;
		for (int ip = 0; ip < packages.size(); ip++) {
			PointState ps_ip = points[packages[ip].gridRef];
			float time = agents[0].timing(ps_ip.p) + agents[0].timing(ps_ip.p, points[i].p);

			if (bestTime == -1 or bestTime > time) {
				bestTime = time;
				bestPackageID = packages[ip].gridRef;
			}
		}

		points[i].agentQueue.push_back(agents[0]);
		points[i].bestTime = bestTime;
		points[i].pointQueue.push_back(points[bestPackageID].p);
	}

	for (int k = 1; k < agents.size(); k++) {
		std::vector<PointState> pointsCopy = points;
		for (int i = 0; i < points.size(); i++) {
			if (isPackage(i)) continue;
			
			int best_j = -1;
			if (i % 10 == 0) std::cout << i << std::endl;
			// Find the best handoff point j so that it can travel to i faster
			for (int j = 0; j < pointsCopy.size(); j++) {

				float timeTo_j = agents[k].timing(pointsCopy[j].p);
				// Factor in waiting time
				timeTo_j = max(timeTo_j, pointsCopy[j].bestTime);

				float time_jTo_i = agents[k].timing(pointsCopy[j].p, pointsCopy[i].p);
				if (timeTo_j + time_jTo_i < points[i].bestTime) {
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


void Scenario::ecld2DType0DynamicSolve11() {
	// Find the drone that can get to the package in the shortest time
	float bestTime = agents[0].timing(packages[0].currentLoc);
	int bestAgent = 0;
	for (int k = 1; k < agents.size(); k++) {
		float timeTmp = agents[k].timing(packages[0].currentLoc);
		if (timeTmp < bestTime) {
			bestTime = timeTmp;
			bestAgent = k;
		}
	}

	// Assign best time, agent queue and point queue to all points in the scenario
	for (int i = 0; i < points.size(); i++) {
		points[i].bestTime = bestTime + agents[bestAgent].timing(packages[0].loc, points[i].p);
		points[i].agentQueue.push_back(agents[bestAgent]);
		points[i].pointQueue.push_back(packages[0].loc);
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
				// Factor in waiting time
				timeTo_j = max(timeTo_j, pointsCopy[j].bestTime);
				 
				float time_jTo_i = agents[k].timing(pointsCopy[j].p, pointsCopy[i].p);
				if (timeTo_j + time_jTo_i < points[i].bestTime) {
					if (i == 401) {
						std::cout << j << std::endl;
					}
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


void Scenario::ecld2DType1DynamicSolveNM() {
	std::vector<int> activeID;
	for (int p = 0; p < packages.size(); p++) {
		for (int t = 0; t < targets.size(); t++) {
			if (packages[p].ID == targets[t].ID) {
				activeID.push_back(targets[t].ID);
				break;
			}
		}
	}

	while (true) {

	}

	for (int _id = 0; _id < activeID.size(); _id++) {

	}
}


void LineAnimation::setColor(int c0, int c1, int c2) {
	color[0] = c0;
	color[1] = c1;
	color[2] = c2;
}


LineAnimation::LineAnimation() {
	prevTimer = -1;
}


void Scenario::createAnimation() {
	if (problemType == 0) {
		int bestTargetID = targets[0].gridRef;
		float bestTime = points[targets[0].gridRef].bestTime;
		for (int i = 1; i < targets.size(); i++) {
			if (bestTime > points[targets[i].gridRef].bestTime) {
				bestTargetID = targets[i].gridRef;
				bestTime = points[targets[i].gridRef].bestTime;
			}
		}

		std::vector<Agent> agentQ = points[bestTargetID].agentQueue;
		std::vector<Point2D> pointQ = points[bestTargetID].pointQueue;
		for (int i = 0; i < agentQ.size(); i++) {
			int color = 255 * (agentQ[i].v - minSpeed) / (maxSpeed - minSpeed);
			LineAnimation tmpAni0, tmpAni1;
			std::vector<LineAnimation> tmpAni;
			if (i < agentQ.size() - 1) {
				tmpAni0.setColor(25, 25, color);
				tmpAni1.setColor(25, 25, color);
				
				tmpAni0.start = agentQ[i].loc0;
				tmpAni0.end = pointQ[i];
				tmpAni0.startTime = 0;
				tmpAni0.endTime = tmpAni0.startTime + 
					agentQ[i].timing(tmpAni0.start, tmpAni0.end);
				tmpAni0.duration = tmpAni0.endTime - tmpAni0.startTime;
				

				tmpAni1.start = tmpAni0.end;
				tmpAni1.end = pointQ[i + 1];
				tmpAni1.startTime = tmpAni0.endTime;
				tmpAni1.endTime = tmpAni1.startTime + 
					agentQ[i].timing(tmpAni1.start, tmpAni1.end);
				tmpAni1.duration = tmpAni1.endTime - tmpAni1.startTime;
			}
			// Special treatment for last agent
			// TODO: trim this down later, redundant code
			else {
				tmpAni0.setColor(25, 25, color);
				tmpAni1.setColor(25, 25, color);

				tmpAni0.start = agentQ[i].loc0;
				tmpAni0.end = pointQ[i];
				tmpAni0.startTime = 0;
				tmpAni0.endTime = tmpAni0.startTime +
					agentQ[i].timing(tmpAni0.start, tmpAni0.end);
				tmpAni0.duration = tmpAni0.endTime - tmpAni0.startTime;

				tmpAni1.start = tmpAni0.end;
				tmpAni1.end = points[bestTargetID].p;
				tmpAni1.startTime = tmpAni0.endTime;
				tmpAni1.endTime = tmpAni1.startTime +
					agentQ[i].timing(tmpAni1.start, tmpAni1.end);
				tmpAni1.duration = tmpAni1.endTime - tmpAni1.startTime;
			}
			tmpAni.push_back(tmpAni0);
			tmpAni.push_back(tmpAni1);
			anis.push_back(tmpAni);
		}
		// Do another run to add waiting time
		for (int i = 1; i < anis.size(); i++) {
			float prevAgentAniEndTime = anis[i - 1][anis[i - 1].size() - 1].endTime;
			if (anis[i][0].endTime < prevAgentAniEndTime) {
				float diff = prevAgentAniEndTime - anis[i][0].endTime;
				
				// Only add wait time after handoff
				for (int j = 1; j < anis[i].size(); j++) {
					anis[i][j].startTime += diff;
					anis[i][j].endTime += diff;
				}

			}
		}
	}
}


Scenario::Scenario() {
	timer = 0;
	aniStart = 0;
}


void Scenario::solve() {
	if (problemType == 0) {
		if (solverType == 0) {
			//ecld2DType0DynamicSolve11();
			ecld2DType0DynamicSolveNM();
			int a = 1;
			a = a + 1;
		}
	}
	if (problemType == 1) {
		if (solverType == 1) {

		}
	}
}