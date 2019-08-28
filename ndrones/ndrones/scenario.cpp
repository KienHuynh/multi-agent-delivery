#include "scenario.h"


void Scenario::loadDesignatedPoint(
	std::ifstream &myfile,
	int problemType,
	std::vector<DesignatedPoint> &dPoints,
	int nDPoint,
	std::vector<PointState> &points) {

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


void Scenario::loadDesignatedPolygon(
	std::ifstream &myfile,
	int problemType,
	std::vector<DesignatedPoint> &dPoints,
	int nDPoint,
	std::vector<PointState> &points) {

	std::vector<Point2D> vList;
	for (int i = 0; i < nDPoint; i++) {
		int x, y;
		myfile >> x >> y;
		vList.push_back(Point2D(x, y));
	}

	for (int i = 0; i < vList.size(); i++) {
		Point2D start = vList[i % nDPoint];
		Point2D end = vList[(i + 1) % nDPoint];
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


bool Scenario::isPackage(int id, std::vector<DesignatedPoint> _packages) {
	for (int ip = 0; ip < _packages.size(); ip++) {
		if (id == _packages[ip].gridRef) {
			return true;
		}
	}
	return false;
}


void Scenario::ecld2DType0DynamicNMCommon(
	std::vector<Agent> _agents,
	std::vector<DesignatedPoint> _packages,
	std::vector<PointState> &_points) {
	// First run with the slowest drone
	// Compute the best time it takes for this drone to get to a package and fly to any other point on the grid
	for (int i = 0; i < _points.size(); i++) {
		if (isPackage(i, _packages)) continue;

		float bestTime = -1;
		int bestPackageID = -1;
		for (int ip = 0; ip < _packages.size(); ip++) {
			PointState ps_ip = _points[_packages[ip].gridRef];
			float time = _agents[0].timing(ps_ip.p) + _agents[0].timing(ps_ip.p, _points[i].p);

			if (bestTime == -1 || bestTime > time) {
				bestTime = time;
				bestPackageID = _packages[ip].gridRef;
			}
		}

		_points[i].agentQueue.push_back(_agents[0]);
		_points[i].bestTime = bestTime;
		_points[i].pointQueue.push_back(_points[bestPackageID].p);
	}

	for (int k = 1; k < _agents.size(); k++) {
		std::vector<PointState> pointsCopy = _points;
		for (int i = 0; i < _points.size(); i++) {
			if (isPackage(i, _packages)) continue;

			int best_j = -1;
			if (i % 100 == 0) std::cout << i << std::endl;
			// Find the best handoff point j so that it can travel to i in the shortest time
			for (int j = 0; j < pointsCopy.size(); j++) {

				float timeTo_j = _agents[k].timing(pointsCopy[j].p);
				// Factor in waiting time
				timeTo_j = max(timeTo_j, pointsCopy[j].bestTime);

				float time_jTo_i = _agents[k].timing(pointsCopy[j].p, pointsCopy[i].p);
				if (timeTo_j + time_jTo_i < _points[i].bestTime) {
					_points[i].bestTime = timeTo_j + time_jTo_i;
					best_j = j;
				}
			}
			if (best_j > -1) {
				_points[i].agentQueue = pointsCopy[best_j].agentQueue;
				_points[i].agentQueue.push_back(_agents[k]);
				_points[i].pointQueue = pointsCopy[best_j].pointQueue;
				_points[i].pointQueue.push_back(pointsCopy[best_j].p);
			}
		}
	}
}


void Scenario::ecld2DType0DynamicNM() {
	ecld2DType0DynamicNMCommon(agents, packages, points);
}


void Scenario::ecld2DType0Dynamic11() {
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
	for (int k = bestAgent + 1; k < agents.size(); k++) {
		std::vector<PointState> pointsCopy = points;
		for (int i = 0; i < points.size(); i++) {
			// Find the best handoff point j so that it can travel to i in the shortest time
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



float Scenario::findBestTimeFromTargets(std::vector<PointState> _points, std::vector<DesignatedPoint> _targets) {
	float bestTime = _points[_targets[0].gridRef].bestTime;
	for (int t = 1; t < _targets.size(); t++) {
		if (bestTime > _points[_targets[t].gridRef].bestTime) {
			bestTime = _points[_targets[t].gridRef].bestTime;
		}
	}
	return bestTime;
}


DesignatedPoint Scenario::findBestTarget(std::vector<PointState> _points, std::vector<DesignatedPoint> _targets) {
	float bestTime = _points[_targets[0].gridRef].bestTime;
	int best_t = 0;
	for (int t = 1; t < _targets.size(); t++) {
		if (bestTime > _points[_targets[t].gridRef].bestTime) {
			bestTime = _points[_targets[t].gridRef].bestTime;
			best_t = t;
		}
	}
	return targets[best_t];
}


// TODO: optimize this
float Scenario::maxValWithoutK(float *arr, int size, int k) {
	if (size == 1) return -1;
	float m = 0;
	float *copy = new float[size];
	
	for (int i = 0; i < size; i++) copy[i] = arr[i];
	std::sort(copy, copy + size);
	
	if (arr[k] == copy[size - 1]) m = copy[size - 2];
	else m = copy[size - 1];
	
	delete copy;

	return m;
}


// TODO: Implement target / package cluster with unique ID for easier management
void Scenario::ecld2DType1DynamicNM() {

	// An ID is active if there exists a package-target matching for it
	// Assume there are A active ids
	std::vector<int> activeID;
	for (int p = 0; p < packages.size(); p++) {
		for (int t = 0; t < targets.size(); t++) {
			if (packages[p].ID == targets[t].ID) {
				if (activeID.size() > 0 &&
					std::find(activeID.begin(), activeID.end(), packages[p].ID) != activeID.end()) continue;
				activeID.push_back(targets[t].ID);
				break;
			}
		}
	}

	// This stores the agent queue for each matching
	bestAgentQueues = new std::vector<Agent>[activeID.size()];
	bestPointQueues = new std::vector<Point2D>[activeID.size()];
	bestTargets = new DesignatedPoint[activeID.size()];

	std::vector<PointState> pointsCopy;
	std::vector<Agent> agentsCopy, agentsCopyWithout_k;

	// Split agents into groups, each belong to a matching
	// i.e. the conflict solving stage
	for (int k = 0; k < agents.size(); k++) {

		// Store the best delivery time for each matching, assume that this matching uses all available drones
		float* bestTimes = new float[activeID.size()];
		// Store the best delivery time for each matching without the "k" agent
		float* bestTimesWithout_k = new float[activeID.size()];

		// Compute the "opt" solution for each matching with the same pool of agents
		// There might be agent(s) appearing in more than one matching
		for (int i = 0; i < activeID.size(); i++) {
			int matchID = activeID[i];
			agentsCopy = agents;
			agentsCopyWithout_k = agents;

			// Remove all drones (up to k-1) that don't belong to this particular matching
			for (int j = k-1; j >= 0; j--) {
				if ((std::find(bestAgentQueues[i].begin(), bestAgentQueues[i].end(), agents[j])) == bestAgentQueues[i].end()) {
					agentsCopy.erase(agentsCopy.begin() + j);
					agentsCopyWithout_k.erase(agentsCopyWithout_k.begin() + j);
				}
			}

			// Remove the k-th drone for agentsCopyWithout_k
			agentsCopyWithout_k.erase(
				std::find(agentsCopyWithout_k.begin(), agentsCopyWithout_k.end(), agents[k])
			);

			std::vector<DesignatedPoint> packagesOfID;
			std::vector<DesignatedPoint> targetsOfID;
			std::copy_if(packages.begin(), packages.end(), std::back_inserter(packagesOfID),
				[&](DesignatedPoint p) {return p.ID == matchID; });
			std::copy_if(targets.begin(), targets.end(), std::back_inserter(targetsOfID),
				[&](DesignatedPoint p) {return p.ID == matchID; });

			// Solve for all available agents
			pointsCopy = points;
			ecld2DType0DynamicNMCommon(agentsCopy, packagesOfID, pointsCopy);
			bestTimes[i] = findBestTimeFromTargets(pointsCopy, targetsOfID);

			// Solve for all available agents but the k-th one
			if (agentsCopyWithout_k.size() > 0) {
				pointsCopy = points;
				ecld2DType0DynamicNMCommon(agentsCopyWithout_k, packagesOfID, pointsCopy);
				bestTimesWithout_k[i] = findBestTimeFromTargets(pointsCopy, targetsOfID);
			}
			else {
				// Only happens in the last loop
				bestTimesWithout_k[i] = INFINITY;
			}
		}
		
		// Compare the results. Assign the k-th agent to a matching so that the overall time after assignment would increase with the smallest amount
		float overallTime = *std::max_element(bestTimes, bestTimes + activeID.size());
		float newOverallTime = -1;
		int assignment = -1;

		for (int i = 0; i < activeID.size(); i++) {
			// Check if k-th drone is used
			if (bestTimes[i] != bestTimesWithout_k[i]) {
				float tmp = maxValWithoutK(bestTimesWithout_k, activeID.size(), i);
				if (assignment == -1 || newOverallTime > tmp) {
					assignment = i;
					newOverallTime = tmp;
				}
			}
		}

		if (assignment != -1) bestAgentQueues[assignment].push_back(agents[k]);
	}
	
	// Re-compute the solution with the above assignments
	for (int i = 0; i < activeID.size(); i++) {
		int matchID = activeID[i];
		pointsCopy = points;
		
		std::vector<DesignatedPoint> packagesOfID;
		std::vector<DesignatedPoint> targetsOfID;
		std::copy_if(packages.begin(), packages.end(), std::back_inserter(packagesOfID),
			[&](DesignatedPoint p) {return p.ID == matchID; });
		std::copy_if(targets.begin(), targets.end(), std::back_inserter(targetsOfID),
			[&](DesignatedPoint p) {return p.ID == matchID; });

		ecld2DType0DynamicNMCommon(bestAgentQueues[i], packagesOfID, pointsCopy);
		bestTargets[i] = findBestTarget(pointsCopy, targetsOfID);
		float bestTime = findBestTimeFromTargets(pointsCopy, targetsOfID);
		bestPointQueues[i] = pointsCopy[bestTargets[i].gridRef].pointQueue;
		// Shouldn't be different
		bestAgentQueues[i] = pointsCopy[bestTargets[i].gridRef].agentQueue;
	}

	// TODO: Insert unused/available agent into the list as a new one according to its speed & delay time
	// Index: k + 1
	int a = 1;
	a++;
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
	// TODO: Update this to use the new solution variables
	if (problemType == 0) {
		int bestTargetID = targets[0].gridRef;
		float bestTime = points[targets[0].gridRef].bestTime;
		for (const DesignatedPoint target : targets) {
			if (bestTime > points[target.gridRef].bestTime) {
				bestTargetID = target.gridRef;
				bestTime = points[target.gridRef].bestTime;
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

	if (problemType == 1) {

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
			ecld2DType0DynamicNM();
			int a = 1;
			a = a + 1;
		}
	}
	if (problemType == 1) {
		if (solverType == 1) {
			ecld2DType1DynamicNM();
		}
	}
}