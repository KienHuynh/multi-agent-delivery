#include "scenario.h"


void Scenario::loadFile(const char* fileName) {
	ScenarioIO::loadFile(fileName, *this);
}


void Scenario::updateMinMaxXY(float x, float y) {
	if (x < minX) minX = x;
	if (x > maxX) maxX = x;
	if (y < minY) minY = y;
	if (y > maxY) maxY = y;
}


void Scenario::writeSolution(const char* fileName) {
	ScenarioIO::writeSolution(fileName, *this);
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


Color Scenario::agentQueueColorMap(std::vector<Agent> _agents) {
	int seed = 0;
	for (int i = 0; i < _agents.size(); i++) {
		seed += _agents[i].ID*_agents[i].loc.x + _agents[i].loc.y;
	}
	srand(seed);
	int id = (rand() + 1111) % cfg::numColor;
	return Palette::palette[id];
}


Color Scenario::agentColorMap(Agent agent) {
	srand(agent.ID);
	return Palette::palette[rand() % cfg::numColor];
}


Color Scenario::bestTimeColorMap(float s, float l, float t) {
	int paletteIndex = (int)(t * cfg::numColor / (l - s));
	return Palette::huePalette[paletteIndex];
}


Color Scenario::depotColorMap(DesignatedPoint dp) {
	srand(dp.ID*(int)dp.loc.x + (int)dp.loc.y);
	int x = rand() % cfg::numColor;
	/*std::cout << x << " " << (int)Palette::palette[x].r << " " <<
		(int) Palette::palette[x].g << " " <<
		(int)Palette::palette[x].b << " " <<std::endl;*/
	return Palette::palette[x];
}


bool Scenario::containPoint(Point2D p) {
	for (int i = 0; i < points.size(); i++) {
		if (points[i].p == p) return true;
	}
	return false;
}


std::vector<Point2D> Scenario::indices2Point2D(std::vector<int> indices) {
	std::vector<Point2D> v;
	for (auto i : indices) {
		v.push_back(points[i].p);
	}
	return v;
}


float Scenario::timing(Agent a, int i) {
	ShortestPath stp = geodesicL2Distance(a.gridRef, i);
	return (stp.length / a.v);
}


float Scenario::timing(Agent a, int i, int j) {
	ShortestPath stp = geodesicL2Distance(i, j);
	return (stp.length / a.v);
}


float Scenario::timing(Agent agent, Point2D a, Point2D b) {
	int ai, bi;
	for (int k = 0; k < points.size(); k++) {
		if (points[k].p == a) ai = k;
		if (points[k].p == b) bi = k;
	}
	return timing(agent, ai, bi);
}


bool Scenario::isDesignatedPoint(PointState ps) {
	return ps.isDesignatedPoint;
}


Scenario::ShortestPath Scenario::geodesicL2Distance(int a, int b) {
	ShortestPath stp;
	std::map<std::vector<int>, ShortestPath>::iterator it =
		stpMap.find(std::vector<int>({ a, b }));

	// Straight line
	if (it == stpMap.end()) {
		stp.path.push_back(b);
		stp.length = Point2D::l2Distance(points[a].p, points[b].p);
	}
	else {
		stp = it->second;
	}

	return stp;
}


Scenario::ShortestPath Scenario::geodesicL2Distance(Point2D a, Point2D b) {
	int ai, bi;
	for (int k = 0; k < points.size(); k++) {
		if (points[k].p == a) ai = k;
		if (points[k].p == b) bi = k;
	}
	return geodesicL2Distance(ai, bi);
}


void Scenario::constructSTPMap() {
	// If there is no obstacle, we can use L2 distance for everything
	// and there is no need for a graph
	if (obs.size() == 0) return;

	// Construct the graph
	// NOTE: at this point, all points inside the obstacles should have been removed
	// The SimplePolygon::segIntersect function still check for the cases where
	// a point lies inside the polygon just in case
	EdgeList* edgeListArr = new EdgeList[points.size()];
	for (int i = 0; i < points.size(); i++) {
		EdgeList edgeList_i;
		for (int j = i + 1; j < points.size(); j++) {
			EdgeList edgeList_j;
			if (!points[i].isOb && !points[j].isOb) continue;
			if (i == j) continue;

			Point2D pi = points[i].p;
			Point2D pj = points[j].p;
			bool isEdge = true;
			// |V|^2 * K 
			// TODO: modify segIntersect so that it takes logK
			for (auto ob : obs) {
				if (i == 88 && j == 91) {
					int a = 0;
					a++;
				}

				// If (i,j) intersects with any obstacle, it's not a graph edge
				if (ob.segIntersect(pi, pj)) {
					isEdge = false;
					break;
				}

				// If (i,j) is not intersecting, but is actually an obstacle's diagonal
				// it's not an edge either
				else if (ob.diagonal(pi, pj)) {
					isEdge = false;
					break;
				}
			}

			if (isEdge) {
				edgeList_i.e.push_back(j);
				edgeList_i.w.push_back(Point2D::l2Distance(pi, pj));
				
				edgeList_j.e.push_back(i);
				edgeList_j.w.push_back(Point2D::l2Distance(pi, pj));
				edgeListArr[j].e.insert(edgeListArr[j].e.end(), edgeList_j.e.begin(), edgeList_j.e.end());
				edgeListArr[j].w.insert(edgeListArr[j].w.end(), edgeList_j.w.begin(), edgeList_j.w.end());
			}
		}
		edgeListArr[i].e.insert(edgeListArr[i].e.end(), edgeList_i.e.begin(), edgeList_i.e.end());
		edgeListArr[i].w.insert(edgeListArr[i].w.end(), edgeList_i.w.begin(), edgeList_i.w.end());
	}

	// Dijkstra's algorithm
	std::vector<int> q0;
	for (int i = 0; i < points.size(); i++) {
		q0.push_back(i);
		edgeList.push_back(edgeListArr[i]);
	}

	for (int i = 0; i < points.size(); i++) {
		if (i % 10 == 0) std::cout << i << std::endl;
		std::vector<PointState> pointsCopy = points;
		std::vector<int> q = q0;
		pointsCopy[i].dist = 0;

		//q.push_back(i);
		while (q.size() > 0) {
			int qu = 0;
			float minDist = pointsCopy[q[0]].dist;
			for (int j = 0; j < q.size(); j++) {
				if (pointsCopy[q[j]].dist < minDist) {
					minDist = pointsCopy[q[j]].dist;
					qu = j;
				}
			}

			int u = q[qu];
			q.erase(q.begin() + qu);

			for (int j = 0; j < edgeList[u].e.size(); j++) {
				int v = edgeList[u].e[j];
				float alt = pointsCopy[u].dist + edgeList[u].w[j];
				if (alt < pointsCopy[v].dist) {
					pointsCopy[v].dist = alt;
					pointsCopy[v].prev = u;
				}
			}
		}

		// Paste shortest paths to the map
		for (int j = 0; j < points.size(); j++) {
			if (j == i) continue;

			Point2D pi = points[i].p;
			Point2D pj = points[j].p;
			bool isEdge = true;
			// |V|^2 * K 
			for (auto ob : obs) {
				if (ob.segIntersect(pi, pj)) {
					isEdge = false;
					break;
				}
				else if (ob.diagonal(pi, pj)) {
					isEdge = false;
					break;
				}
			}
			if (!isEdge) {
				std::vector<int> key({ i,j });
				ShortestPath stp;
				stp.length = 0;
				stp.path.push_back(j);
				int k = j;
				while (pointsCopy[k].prev != i) {
					stp.path.insert(stp.path.begin(), pointsCopy[k].prev);
					k = pointsCopy[k].prev;
				}
				stp.length = pointsCopy[j].dist;

				stpMap.insert(std::pair<std::vector<int>, ShortestPath>(key, stp));
			}
		}
	}

	delete[] edgeListArr;
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
			int ps_ipGridRef = _packages[ip].gridRef;
			float time = timing(_agents[0], ps_ipGridRef) +
				timing(_agents[0], ps_ipGridRef, i) + 
				_agents[0].delay;
			// float time = _agents[0].timing(ps_ip.p) + _agents[0].timing(ps_ip.p, _points[i].p) + _agents[0].delay;

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
		std::cout << k << std::endl;
		std::vector<PointState> pointsCopy = _points;

		// Sort the pointsCopy by their best time
		std::vector<size_t> indices = argSort(pointsCopy);
		
		for (int i = 0; i < _points.size(); i++) {
			if (isPackage(i, _packages)) continue;

			// Find index of point i in the sorted list
			// int sorted_i = findVectorIndexWithBestTime(pointsCopy, _points[i]);

			int best_j = -1;
			// Find the best handoff point j so that it can travel to i in the shortest time
			for (auto j: indices) {
				if (i == j) continue;

				if (pointsCopy[j].bestTime > _points[i].bestTime) break;
				float timeTo_j = timing(_agents[k], j) + _agents[k].delay;

				// Factor in waiting time
				timeTo_j = max(timeTo_j, pointsCopy[j].bestTime);
				float time_jTo_i = timing(_agents[k], j, i);
				
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
	activeID.push_back(packages[0].ID);
	ecld2DType0DynamicNMCommon(agents, packages, points);
	bestAgentQueues = new std::vector<Agent>;
	bestPointQueues = new std::vector<Point2D>;
	bestTargets = new DesignatedPoint;
	bestTimes = new float;

	DesignatedPoint bestTarget = findBestTarget(points, targets);

	bestAgentQueues[0] = points[bestTarget.gridRef].agentQueue;
	bestPointQueues[0] = points[bestTarget.gridRef].pointQueue;
	bestPointQueues[0].push_back(bestTarget.loc);
	bestTargets[0] = bestTarget;
	bestTimes[0] = points[bestTarget.gridRef].bestTime;
	overallTime = bestTimes[0];

	makespan = points[bestTarget.gridRef].bestTime;
	std::cout << makespan << std::endl;
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
	return _targets[best_t];
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


void Scenario::updateReusedAgent(std::vector<Agent> & _agents, std::vector<Point2D> _pointQueue) {
	float prevTime = 0;
	for (int i = 0; i < _agents.size(); i++) {
		float timeToPoint_i = max(prevTime, _agents[i].timing(_pointQueue[i]) + _agents[i].delay);
		float time_iToHandOff = _agents[i].timing(_pointQueue[i], _pointQueue[i + 1]);

		prevTime = timeToPoint_i + time_iToHandOff;
		if (i < _agents.size() - 1) {
			prevTime = max(prevTime, _agents[i + 1].delay + _agents[i + 1].timing(_pointQueue[i + 1]));
		}

		_agents[i].delay = prevTime;
		_agents[i].loc0 = _pointQueue[i + 1];
		_agents[i].loc = _pointQueue[i + 1];
	}
}


bool Scenario::containAgentID(std::vector<Agent> _agents, Agent a) {
	for (const auto & _agent : _agents) {
		if (_agent.ID == a.ID) return true;
	}
	return false;
}


bool Scenario::containAgent(std::vector<Agent> _agents, Agent a) {
	for (const auto & _agent : _agents) {
		if (_agent.ID == a.ID &&
			_agent.delay == a.delay &&
			a.loc == _agent.loc) return true;
	}
	return false;
}


bool Scenario::containAgentAfterOrder(std::vector<Agent> _agents, Agent a) {
	for (const auto & _agent : _agents) {
		if (_agent.ID == a.ID &&
			a.orderOfEx < _agent.orderOfEx) return true;
	}
	return false;
}


// Return false if there is no conflict to resolve, telling the caller that it should just let the assignments stay the same
bool Scenario::conflictResolve(
	std::vector<PointState> _points,
	std::vector<Agent> _agents,
	std::vector<DesignatedPoint>* _packagesOfID,
	std::vector<DesignatedPoint>* _targetsOfID,
	std::vector<int> _activeID,
	int** _agentAssignment,
	float* _bestTimes,
	int &bestMatchingInd,
	int &bestAgentInd) {

	//// TODO: improve this block of codes
	//bool stop = false;
	//for (int k = 0; k < agents.size() && !stop; k++) {
	//	int _bestAgentInd = -1;
	//	int _bestMatchingInd = -1;
	//	int assignmentCount = 0;
	//	for (int a = 0; a < activeID.size() && !stop; a++) {
	//		if (_agentAssignment[a][k] == 1) {
	//			_bestAgentInd = k;
	//			_bestMatchingInd = a;
	//			assignmentCount += 1;
	//		}
	//		if (assignmentCount == 2) {
	//			stop = true;
	//			break;
	//		}
	//	}
	//	if (assignmentCount == 1) {
	//		bestMatchingInd = _bestMatchingInd;
	//		bestAgentInd = _bestAgentInd;
	//		return false;
	//	}
	//}

	float** bestTimesWithout_k = new float*[_activeID.size()];
	// This stores the overall time (i.e. max of all delivery time of all matchings) when 
	// agent k is assigned to matching i (other matchings DO NOT HAVE this agent)
	float** overallTime_k = new float*[_activeID.size()];
	for (int a = 0; a < activeID.size(); a++) {
		bestTimesWithout_k[a] = new float[_agents.size()];
		overallTime_k[a] = new float[_agents.size()];
		for (int k = 0; k < _agents.size(); k++) {
			bestTimesWithout_k[a][k] = -1;
			overallTime_k[a][k] = -1;
		}
	}

	bool isConflict = false;
	for (int k = 0; k < _agents.size(); k++) {
		std::vector<int> matchingInd;
		for (int a = 0; a < _activeID.size(); a++) {
			if (_agentAssignment[a][k] == 1) matchingInd.push_back(a);
		}

		// Compute the solution without the i-th agent for each matching found above
		if (matchingInd.size() <= 1) continue;
		else {
			for (const auto a : matchingInd) {
				std::vector<Agent> agentsCopy = _agents;
				std::vector<PointState> pointsCopy = _points;

				// Removed shared agents
				removeSharedAgents(bestAgentQueues, a, agentsCopy);

				if (agentsCopy.size() < 1) continue;

				for (int i = 0; i < bestAgentQueues[a].size(); i++) {
					int vectorInd = findVectorIndexWithID(agentsCopy, bestAgentQueues[a][i]);
					if (vectorInd < 0) continue;
					agentsCopy[vectorInd] = bestAgentQueues[a][i];
				}

				int vectorIndex = findVectorIndexWithID(agentsCopy, _agents[k]);
				if (vectorIndex < 0) continue;

				agentsCopy.erase(agentsCopy.begin() + vectorIndex);

				if (agentsCopy.size() == 0) {
					bestTimesWithout_k[a][k] = INFINITY;
				}
				else {
					ecld2DType0DynamicNMCommon(agentsCopy, _packagesOfID[a], pointsCopy);
					bestTimesWithout_k[a][k] = findBestTimeFromTargets(pointsCopy, _targetsOfID[a]);
				}
			}

			for (const auto a1 : matchingInd) {
				float maxVal = -1;
				for (const auto a2 : matchingInd) {
					if (a1 == a2) continue;
					if (maxVal < bestTimesWithout_k[a2][k]) {
						maxVal = bestTimesWithout_k[a2][k];
					}
				}
				overallTime_k[a1][k] = maxVal;
			}

			isConflict = true;
		}
	}

	// Find the best matching and agent index, i.e. the index with lowest overall time
	float minVal = -1;
	for (int a = 0; a < activeID.size(); a++) {
		for (int k = 0; k < agents.size(); k++) {
			if (overallTime_k[a][k] == -1) continue;
			if (minVal == -1 || minVal > overallTime_k[a][k]) {
				bestMatchingInd = a;
				bestAgentInd = k;
			}
		}
	}

	return isConflict;
}


int Scenario::findVectorIndexWithID(std::vector<Agent> _agents, Agent a) {
	for (int k = 0; k < _agents.size(); k++) {
		if (_agents[k].ID == a.ID) return k;
	}
	return -1;
}


int Scenario::findVectorIndexFull(std::vector<Agent> _agents, Agent a) {
	for (int k = 0; k < _agents.size(); k++) {
		if (_agents[k].ID == a.ID &&
			_agents[k].delay == a.delay &&
			_agents[k].orderOfEx == a.orderOfEx) return k;
	}
	return -1;
}


int Scenario::findVectorIndexWithBestTime(std::vector<PointState> _points, PointState p) {
	int min_i = 0;
	int max_i = points.size();
	int i = -1;
	while (min_i != max_i) {
		i = (min_i + max_i) / 2;
		if (_points[i].bestTime == p.bestTime) return i;
		if (_points[i].bestTime < p.bestTime) {
			min_i = i;
		}
		else {
			max_i = i;
		}
		if (min_i == max_i) return i;
	}
	return -1;
}


void Scenario::removeSharedAgents(std::vector<Agent>* queues, int id, std::vector<Agent> &agents) {
	for (int a = 0; a < activeID.size(); a++) {
		if (a == id) continue;
		for (int i = 0; i < queues[a].size(); i++) {
			int removalIndex = findVectorIndexFull(agents, queues[a][i]);
			if (removalIndex < 0) continue;
			agents.erase(agents.begin() + removalIndex);
		}
	}
}


bool Scenario::compareOrderOfEx(Agent a, Agent b) {
	return a.orderOfEx < b.orderOfEx;
}


void Scenario::removeGapAgents(std::vector<Agent> &_agents, int _id, int _orderOfEx) {
	for (int k = _agents.size() - 1; k >= 0; k--) {
		if (_agents[k].orderOfEx > _orderOfEx && _agents[k].ID == _id) _agents.erase(_agents.begin() + k);
	}
}


void Scenario::bagAgentsByOrder(std::vector<Agent>* _agentQueues, std::vector<Agent>*& bag) {
	// Inserting all agents into the orderOfExecutionQueue
	delete[] bag;
	bag = new std::vector<Agent>[agents.size()];
	for (int a = 0; a < activeID.size(); a++) {
		for (int k = 0; k < bestAgentQueues[a].size(); k++) {
			int vectorInd = findVectorIndexWithID(agents, bestAgentQueues[a][k]);
			bag[vectorInd].push_back(bestAgentQueues[a][k]);
		}
	}

	// Sort the orderOfExecutionQueue
	for (int k = 0; k < agents.size(); k++) {
		std::sort(bag[k].begin(), bag[k].end(),
			compareOrderOfEx);
	}
}


bool Scenario::equalAgentQueue(std::vector<Agent>* qA, std::vector<Agent>* qB) {
	for (int a = 0; a < activeID.size(); a++) {
		if (qA[a].size() == qB[a].size()) {
			for (int k = 0; k < qA[a].size(); k++) {
				if (qA[a][k].ID != qB[a][k].ID ||
					qA[a][k].orderOfEx != qB[a][k].orderOfEx ||
					qA[a][k].delay != qB[a][k].delay) return false;
			}
		}
		else return false;
	}
	return true;
}


bool Scenario::missingQueue(std::vector<Agent>* qs) {
	for (int a = 0; a < activeID.size(); a++) {
		if (qs[a].size() == 0) return true;
	}
	return false;
}


void Scenario::ecld2DType1DynamicNM() {

	// An ID is active if there exists a package-target matching for it
	// Assume there are A active ids
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
	std::vector<Agent>* prevBestAgentQueues = new std::vector<Agent>[activeID.size()];
	bestAgentQueues = new std::vector<Agent>[activeID.size()];
	bestPointQueues = new std::vector<Point2D>[activeID.size()];

	bestTargets = new DesignatedPoint[activeID.size()];
	bestTimes = new float[activeID.size()];

	float* prevBestTimes = new float[activeID.size()];
	for (int i = 0; i < activeID.size(); i++) {
		prevBestTimes[i] = INFINITY;
	}

	// Storing packages and targets of corresponding ID
	std::vector<DesignatedPoint>* packagesOfID = new std::vector<DesignatedPoint>[activeID.size()];
	std::vector<DesignatedPoint>* targetsOfID = new std::vector<DesignatedPoint>[activeID.size()];
	for (int a = 0; a < activeID.size(); a++) {
		std::copy_if(packages.begin(), packages.end(), std::back_inserter(packagesOfID[a]),
			[&](DesignatedPoint p) {return p.ID == activeID[a]; });
		std::copy_if(targets.begin(), targets.end(), std::back_inserter(targetsOfID[a]),
			[&](DesignatedPoint p) {return p.ID == activeID[a]; });
	}

	std::vector<PointState> pointsCopy;
	std::vector<Agent> newAgents, agentsCopy, agentsCopyWithout_k;

	newAgents = agents;

	/*for (int k = 0; k < agents.size(); k++) {
		orderOfExecutionQueue[k].push_back(agents[k]);
	}*/

	// Stopping flag, is true if there is no change
	bool stop = false;
	bool mainLoopStop = false;
	int loop = 0;
	while (!mainLoopStop) {
		std::vector<Agent>* orderOfExecutionQueue;
		orderOfExecutionQueue = new std::vector<Agent>[agents.size()];

		std::cout << loop << std::endl;
		loop++;

		stop = false;

		int innerLoop = 0;
		while (!stop) {
			std::cout << "Inner loop: " << innerLoop << std::endl;
			innerLoop++;
			float* bestTimes = new float[activeID.size()];
			// This table store 0/1 values, indicating if an agent is used in a matching or not
			// Used for conflict resolving
			int** agentAssignment = new int*[activeID.size()];
			for (int a = 0; a < activeID.size(); a++) {
				agentAssignment[a] = new int[agents.size()];
				for (int k = 0; k < agents.size(); k++) {
					agentAssignment[a][k] = 0;
				}
			}

			for (int a = 0; a < activeID.size(); a++) {
				int matchID = activeID[a];
				agentsCopy = newAgents;

				// Remove all agents that is already belong to another matching
				removeSharedAgents(bestAgentQueues, a, agentsCopy);

				// Inserting its old agent from bestAgentQueues[a]
				// agentsCopy.insert(agentsCopy.end(), bestAgentQueues[a].begin(), bestAgentQueues[a].end());
				// std::sort(agentsCopy.begin(), agentsCopy.end());

				for (int i = 0; i < bestAgentQueues[a].size(); i++) {
					int vectorInd = findVectorIndexWithID(agentsCopy, bestAgentQueues[a][i]);
					if (vectorInd < 0) {
						agentsCopy.push_back(bestAgentQueues[a][i]);
					}
					else {
						agentsCopy[vectorInd] = bestAgentQueues[a][i];
					}
				}
				std::sort(agentsCopy.begin(), agentsCopy.end());

				if (agentsCopy.size() == 0) continue;

				pointsCopy = points;
				ecld2DType0DynamicNMCommon(agentsCopy, packagesOfID[a], pointsCopy);
				bestTimes[a] = findBestTimeFromTargets(pointsCopy, targetsOfID[a]);
				DesignatedPoint bestTarget = findBestTarget(pointsCopy, targetsOfID[a]);
				std::vector<Agent> agentQueue_a = pointsCopy[bestTarget.gridRef].agentQueue;

				for (int i = 0; i < agentQueue_a.size(); i++) {
					int vectorInd = findVectorIndexWithID(newAgents, agentQueue_a[i]);
					// If this matching already has this agent (assigned from a previous loop)
					// There is no point in adding it now
					if (containAgentID(bestAgentQueues[a], agentQueue_a[i])) continue;
					agentAssignment[a][vectorInd] = 1;
				}
			}

			int bestMatchingInd, bestAgentInd;
			bool isConflict = conflictResolve(
				points,
				newAgents,
				packagesOfID,
				targetsOfID,
				activeID,
				agentAssignment,
				bestTimes,
				bestMatchingInd,
				bestAgentInd);

			if (isConflict) {
				bestAgentQueues[bestMatchingInd].push_back(newAgents[bestAgentInd]);
			}
			else stop = true;

			for (int a = 0; a < activeID.size(); a++) {
				delete[] agentAssignment[a];
			}
			delete[] agentAssignment;
			delete bestTimes;
		}

		// Compute the solution for each matching using agents from bestAgentQueues
		// and remaining agents that weren't assigned to anything (they are either generating
		// no conflict or they have no use)
		std::vector<Agent> tmpAgents = agents;
		for (int a = 0; a < activeID.size(); a++) {
			pointsCopy = points;
			agentsCopy = newAgents;

			// Remove all agents that is already belong to another matching
			removeSharedAgents(bestAgentQueues, a, agentsCopy);

			// Inserting its old agent from bestAgentQueues[a]
			// agentsCopy.insert(agentsCopy.end(), bestAgentQueues[a].begin(), bestAgentQueues[a].end());
			// std::sort(agentsCopy.begin(), agentsCopy.end());

			for (int i = 0; i < bestAgentQueues[a].size(); i++) {
				int vectorInd = findVectorIndexWithID(agentsCopy, bestAgentQueues[a][i]);
				if (vectorInd < 0) {
					agentsCopy.push_back(bestAgentQueues[a][i]);
				}
				else {
					agentsCopy[vectorInd] = bestAgentQueues[a][i];
				}
			}
			std::sort(agentsCopy.begin(), agentsCopy.end());

			if (agentsCopy.size() == 0) continue;

			ecld2DType0DynamicNMCommon(agentsCopy, packagesOfID[a], pointsCopy);
			DesignatedPoint bestTarget = findBestTarget(pointsCopy, targetsOfID[a]);
			bestAgentQueues[a] = pointsCopy[bestTarget.gridRef].agentQueue;

			std::vector<Point2D> pointQueueTmp = pointsCopy[bestTarget.gridRef].pointQueue;
			pointQueueTmp.push_back(pointsCopy[bestTarget.gridRef].p);
			std::vector<Agent> agentQueueTmp = bestAgentQueues[a];
			updateReusedAgent(agentQueueTmp, pointQueueTmp);
			for (int k = 0; k < bestAgentQueues[a].size(); k++) {
				bestAgentQueues[a][k].finTime = agentQueueTmp[k].delay;
			}

			// Update newAgents list
			// Multiple matchings can share one agent
			// The agent state will be updated using the information from the matching that uses it last
			for (int k = 0; k < agentQueueTmp.size(); k++) {
				int vectorIndex = findVectorIndexWithID(tmpAgents, agentQueueTmp[k]);
				if (tmpAgents[vectorIndex].delay < agentQueueTmp[k].delay) {
					agentQueueTmp[k].orderOfEx++;
					tmpAgents[vectorIndex] = agentQueueTmp[k];
				}
			}
		}

		bagAgentsByOrder(bestAgentQueues, orderOfExecutionQueue);

		for (int k = 0; k < tmpAgents.size(); k++) {
			for (int i = 1; i < orderOfExecutionQueue[k].size(); i++) {
				// Checking for discrepancy in time between two executions of an agent
				// This would only happen if:
				// Agent A is executed for the (i-1)-th time for matching 1, with fin time t1
				// Agent A is then used for another matching for matching 2, order of execution is i-th, using delay time t1
				// Matching 1 gets assigned with a new drone, this affects fin time t1, maing it t1'
				// Now there is a difference between finTime of A at (i-1)-th and delayTime of A at i-th.
				if (orderOfExecutionQueue[k][i].delay != orderOfExecutionQueue[k][i - 1].finTime) {
					Agent tmp = orderOfExecutionQueue[k][i - 1];

					// Remove all instances of this agent k with orderOfEx i-th onward
					for (int a = 0; a < activeID.size(); a++) {
						if (containAgentAfterOrder(bestAgentQueues[a], tmp)) {
							// Completely clear this queue
							bestAgentQueues[a].clear();
						}
					}
					break;
				}
			}
		}

		bagAgentsByOrder(bestAgentQueues, orderOfExecutionQueue);

		// Detect a gap in order of execution
		for (int k = 0; k < tmpAgents.size(); k++) {
			int gapIndex = -1;
			for (int i = 0; i < orderOfExecutionQueue[k].size(); i++) {
				if (orderOfExecutionQueue[k][i].orderOfEx - gapIndex > 1) {
					break;
				}
				gapIndex++;
			}
			// Remove all agents in the "best" queue with index > gapIndex
			for (int a = 0; a < activeID.size(); a++) {
				removeGapAgents(bestAgentQueues[a], tmpAgents[k].ID, gapIndex);
			}

			// Put the agent back into tmpAgents as new
			if (tmpAgents[k].orderOfEx - gapIndex > 1) {
				if (gapIndex < 0) tmpAgents[k] = agents[k];
				else tmpAgents[k] = orderOfExecutionQueue[k][gapIndex];
			}
		}

		for (int a = 0; a < activeID.size(); a++) {
			if (bestAgentQueues[a].size() == 0) continue;
			std::cout << "Matching " << a << ": ";
			for (int k = 0; k < bestAgentQueues[a].size(); k++) {
				std::cout << bestAgentQueues[a][k].ID << "-" << bestAgentQueues[a][k].orderOfEx << " ";
			}
			std::cout << std::endl;
		}

		if (loop == 1) {
			for (int a = 0; a < activeID.size(); a++) {
				prevBestAgentQueues[a] = bestAgentQueues[a];
			}
		}
		else {
			if (equalAgentQueue(prevBestAgentQueues, bestAgentQueues) &&
				!missingQueue(prevBestAgentQueues)) {
				mainLoopStop = true;
				stop = true;
				break;
			}
			for (int a = 0; a < activeID.size(); a++) {
				prevBestAgentQueues[a] = bestAgentQueues[a];
			}
		}

		newAgents = tmpAgents;

		delete[] orderOfExecutionQueue;
	}

	// Re-compute the solution with the above assignments
	std::vector<DesignatedPoint> allTargets;
	overallTime = -1;
	for (int a = 0; a < activeID.size(); a++) {
		int matchID = activeID[a];
		pointsCopy = points;

		ecld2DType0DynamicNMCommon(bestAgentQueues[a], packagesOfID[a], pointsCopy);
		bestTargets[a] = findBestTarget(pointsCopy, targetsOfID[a]);
		bestTimes[a] = findBestTimeFromTargets(pointsCopy, targetsOfID[a]);
		if (overallTime == -1 || overallTime < bestTimes[a]) overallTime = bestTimes[a];
		bestPointQueues[a] = pointsCopy[bestTargets[a].gridRef].pointQueue;
		bestPointQueues[a].push_back(bestTargets[a].loc);
		// Shouldn't be different
		bestAgentQueues[a] = pointsCopy[bestTargets[a].gridRef].agentQueue;

		// Debug purpose
		allTargets.push_back(findBestTarget(pointsCopy, targetsOfID[a]));
	}
	std::cout << "Overall time: " << overallTime << std::endl;

	// TODO: replace dynamic array with a better data struct
	delete[] packagesOfID;
	delete[] targetsOfID;
}


void LineAnimation::setColor(int c0, int c1, int c2) {
	color[0] = c0;
	color[1] = c1;
	color[2] = c2;
}


LineAnimation::LineAnimation() {
	prevTimer = -1;
}


void Scenario::createDroneAnimation() {
	if ((problemType & (TWODIM | EUCLID | DISCRETE)) == (TWODIM | EUCLID | DISCRETE)) {
		for (int a = 0; a < activeID.size(); a++) {
			int matchID = activeID[a];
			int aniSize = droneAnis.size();
			float bestTime = points[bestTargets[a].gridRef].bestTime;
			float totalDelay = 0;
			
			for (int i = 0; i < bestAgentQueues[a].size(); i++) {
				int color = 255 * (bestAgentQueues[a][i].v - minSpeed) / (maxSpeed - minSpeed);

				// This is used to store every animation
				// There is no time-wise relation between the elements
				std::vector<LineAnimation> tmpAni;

				// This animation is for the i-th drone to fly from its loc0 to the meeting point
				// There shouldn't be any delay
				// First, get the shorest path from loc0 to the meeting point
				ShortestPath stp = geodesicL2Distance(bestAgentQueues[a][i].loc0,
					bestPointQueues[a][i]);
				std::vector<Point2D> stpPoint = indices2Point2D(stp.path);
				stpPoint.insert(stpPoint.begin(), bestAgentQueues[a][i].loc0);
				
				float prevTime = bestAgentQueues[a][i].delay;
				for (int j = 0; j < stpPoint.size()-1; j++) {
					LineAnimation tmpAni0;
					tmpAni0.setColor(25, 25, color);
					tmpAni0.start = stpPoint[j];
					tmpAni0.end = stpPoint[j+1];
					tmpAni0.startTime = prevTime;
					tmpAni0.endTime = tmpAni0.startTime +
						timing(bestAgentQueues[a][i], tmpAni0.start, tmpAni0.end);
					tmpAni0.duration = tmpAni0.endTime - tmpAni0.startTime;
					tmpAni.push_back(tmpAni0);

					prevTime = tmpAni0.endTime;
				}

				// This animation is for the i-th drone to fly from the previous meeting point
				// to a new one
				// There could be a delay here as this drone might have to wait for a
				// slower one to arrive
				stp = geodesicL2Distance(tmpAni[tmpAni.size()-1].end,
					bestPointQueues[a][i+1]);
				stpPoint = indices2Point2D(stp.path);
				stpPoint.insert(stpPoint.begin(), tmpAni[tmpAni.size() - 1].end);
				
				prevTime = totalDelay > prevTime ?
					totalDelay : prevTime;

				for (int j = 0; j < stpPoint.size()-1; j++) {
					LineAnimation tmpAni0;
					tmpAni0.setColor(25, 25, color);
					tmpAni0.start = stpPoint[j];
					tmpAni0.end = stpPoint[j + 1];
					tmpAni0.startTime = prevTime;
					tmpAni0.endTime = tmpAni0.startTime +
						timing(bestAgentQueues[a][i], tmpAni0.start, tmpAni0.end);
					tmpAni0.duration = tmpAni0.endTime - tmpAni0.startTime;
					tmpAni.push_back(tmpAni0);

					prevTime = tmpAni0.endTime;
				}
				totalDelay = prevTime;
				droneAnis.insert(droneAnis.end(), tmpAni.begin(), tmpAni.end());
				// These are only specific to this problem only because we know that each agent only generates 2 animations
				/*tmpAni1.prevAni.push_back(i * 2 + aniSize);
				if (i > 0) tmpAni1.prevAni.push_back((i - 1) * 2 + 1 + aniSize);
				droneAnis.push_back(tmpAni0);
				droneAnis.push_back(tmpAni1);*/
			}
			// Do another run to add waiting time
			/*for (int i = 1; i < droneAnis.size(); i++) {
				if (droneAnis[i].prevAni.size() == 0) continue;
				float maxPrevTime = droneAnis[droneAnis[i].prevAni[0]].endTime;

				for (int j = 1; j < droneAnis[i].prevAni.size(); j++) {
					if (maxPrevTime < droneAnis[droneAnis[i].prevAni[j]].endTime) maxPrevTime = droneAnis[droneAnis[i].prevAni[j]].endTime;
				}

				float diff = maxPrevTime - droneAnis[i].startTime;
				if (diff < 0) continue;
				droneAnis[i].startTime += diff;
				droneAnis[i].endTime += diff;
			}*/
		}
	}
}


void Scenario::createPackageAnimation() {
	if ((problemType & (TWODIM | EUCLID | DISCRETE)) == (TWODIM | EUCLID | DISCRETE)) {
		for (int a = 0; a < activeID.size(); a++) {
			int matchID = activeID[a];
			int aniSize = packageAnis.size();

			for (int i = 0; i < bestAgentQueues[a].size(); i++) {
				LineAnimation tmpAni0, tmpAni1;
				std::vector<LineAnimation> tmpAni;

				tmpAni1.setColor(255, 0, 0);

				tmpAni0.start = bestAgentQueues[a][i].loc0;
				tmpAni0.end = bestPointQueues[a][i];
				tmpAni0.startTime = bestAgentQueues[a][i].delay;
				tmpAni0.endTime = tmpAni0.startTime +
					bestAgentQueues[a][i].timing(tmpAni0.start, tmpAni0.end);
				tmpAni0.duration = tmpAni0.endTime - tmpAni0.startTime;

				tmpAni1.start = tmpAni0.end;
				tmpAni1.end = bestPointQueues[a][i + 1];
				tmpAni1.startTime = tmpAni0.endTime;
				tmpAni1.endTime = tmpAni1.startTime +
					bestAgentQueues[a][i].timing(tmpAni1.start, tmpAni1.end);
				tmpAni1.duration = tmpAni1.endTime - tmpAni1.startTime;

				// These are only specific to this problem only because we know that each agent only generates 2 animations
				// tmpAni1.prevAni.push_back(i * 2 + aniSize);
				if (i > 0) tmpAni1.prevAni.push_back(i - 1 + aniSize);
				// packageAnis.push_back(tmpAni0);
				packageAnis.push_back(tmpAni1);
			}
			for (int i = 1; i < packageAnis.size(); i++) {
				if (packageAnis[i].prevAni.size() == 0) continue;
				float maxPrevTime = packageAnis[packageAnis[i].prevAni[0]].endTime;

				for (int j = 1; j < packageAnis[i].prevAni.size(); j++) {
					if (maxPrevTime < packageAnis[packageAnis[i].prevAni[j]].endTime) maxPrevTime = packageAnis[packageAnis[i].prevAni[j]].endTime;
				}

				float diff = maxPrevTime - packageAnis[i].startTime;
				if (diff < 0) continue;
				packageAnis[i].startTime += diff;
				packageAnis[i].endTime += diff;
			}
		}
	}
}


Scenario::Scenario() {
	timer = 0;
	aniStart = 0;
	makespan = 0;
}


void Scenario::solve() {
	if ((problemType & (TWODIM | EUCLID | DISCRETE)) == (TWODIM | EUCLID | DISCRETE)) {
		if ((problemType & (SINGLE_ID)) == SINGLE_ID) ecld2DType0DynamicNM();
		else ecld2DType1DynamicNM();
	}
}