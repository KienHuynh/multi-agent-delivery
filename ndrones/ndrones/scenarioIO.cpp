#include "scenarioIO.h"


std::istream& operator >> (std::istream &input, DesignatedPointInputMode& f) {
	std::string tmp;
	input >> tmp;
	if (stringEqual(tmp, "single_point")) f = DesignatedPointInputMode::SINGLE_POINT;
	if (stringEqual(tmp, "poly")) f = DesignatedPointInputMode::POLY;

	return input;
}


std::istream& operator >> (std::istream &input, ProblemType& p) {
	std::string tmp;
	input >> tmp;
	if (stringEqual(tmp, "onedim")) p = ONEDIM;
	if (stringEqual(tmp, "twodim")) p = TWODIM;
	if (stringEqual(tmp, "euclid"))	p = EUCLID;
	if (stringEqual(tmp, "graph")) p = GRAPH;
	if (stringEqual(tmp, "discrete")) p = DISCRETE;
	if (stringEqual(tmp, "single_target")) p = SINGLE_TARGET;
	if (stringEqual(tmp, "single_id")) p = SINGLE_ID;

	return input;
}


void ScenarioIO::loadDesignatedPoint(
	std::ifstream &myfile,
	std::vector<DesignatedPoint> &dPoints,
	Scenario &scenario,
	int nDPoint) {

	for (int i = 0; i < nDPoint; i++) {
		int x, y, id;
		if ((scenario.problemType & SINGLE_ID) != 0) {
			myfile >> x >> y;
			dPoints.push_back(DesignatedPoint(Point2D(x, y), -1));
		}
		else if ((scenario.problemType & SINGLE_ID) == 0) {
			myfile >> x >> y >> id;
			dPoints.push_back(DesignatedPoint(Point2D(x, y), id));
		}

		bool in = false;
		for (int p = 0; p < scenario.points.size(); p++) {
			Point2D tmp = scenario.points[p].p;
			// If the new designated point already exists in the grid
			// Do not add a new point into the grid
			if (tmp.x == dPoints[dPoints.size() - 1].loc.x &&
				tmp.y == dPoints[dPoints.size() - 1].loc.y) {
				dPoints[dPoints.size() - 1].gridRef = p;
				scenario.points[p].isDesignatedPoint = true;
				in = true;
				break;
			}
		}

		// Add a new point into the grid
		if (in == false) {
			scenario.points.push_back(PointState(dPoints[dPoints.size() - 1].loc));
			scenario.points[scenario.points.size() - 1].isDesignatedPoint = true;
			dPoints[dPoints.size() - 1].gridRef = scenario.points.size() - 1;
		}
	}
}


void ScenarioIO::loadDesignatedPolygon(
	std::ifstream &myfile,
	std::vector<DesignatedPoint> &dPoints,
	Scenario &scenario,
	int nDPoint) {

	std::vector<DesignatedPoint> vList;
	for (int i = 0; i < nDPoint; i++) {
		int x, y, id = -1;
		if ((scenario.problemType & (SINGLE_ID)) == SINGLE_ID) myfile >> x >> y;
		else myfile >> x >> y >> id;
		vList.push_back(DesignatedPoint(Point2D(x, y), id));
	}

	for (int i = 0; i < vList.size(); i++) {
		Point2D start = vList[i % nDPoint].loc;
		Point2D end = vList[(i + 1) % nDPoint].loc;
		int id = vList[i % nDPoint].ID;
		// Interpolation
		for (int s = 0; s < cfg::polySamplingRate; s++) {
			float rate = ((float)s) / ((float)cfg::polySamplingRate);
			Point2D t1 = end - start;
			Point2D t2 = t1 * rate;
			Point2D middle = (end - start)*rate + start;

			dPoints.push_back(DesignatedPoint(middle, id));

			// Check if this designated point already belong to the grid
			// If not, add it to the point list
			bool in = false;
			for (int p = 0; p < scenario.points.size(); p++) {
				Point2D tmp = scenario.points[p].p;
				if (tmp.x == dPoints[dPoints.size() - 1].loc.x &&
					tmp.y == dPoints[dPoints.size() - 1].loc.y) {
					dPoints[dPoints.size() - 1].gridRef = p;
					scenario.points[p].isDesignatedPoint = true;
					in = true;
				}
			}
			if (in == false) {
				scenario.points.push_back(PointState(dPoints[dPoints.size() - 1].loc));
				scenario.points[scenario.points.size() - 1].isDesignatedPoint = true;
				dPoints[dPoints.size() - 1].gridRef = scenario.points.size() - 1;
			}
		}
	}
}


void ScenarioIO::loadFile(const char* fname, Scenario &scenario) {
	scenario.outputFileName = std::string(fname);
	scenario.outputFileName = scenario.outputFileName.substr(scenario.outputFileName.find_last_of("/\\") + 1);
	std::string::size_type const p(scenario.outputFileName.find_last_of('.'));
	scenario.outputFileName = scenario.outputFileName.substr(0, p);

	std::ifstream myfile;
	myfile.open(fname, std::ios::in);

	// Reading problem type
	scenario.problemType = ProblemType::UNKNOWN;
	std::string str;
	getline(myfile, str);
	std::istringstream ss(str);

	ProblemType tmp;
	while (ss >> tmp) {
		scenario.problemType = (ProblemType)(scenario.problemType | tmp);
	}

	int nPackage = 0, nTarget = 0, nAgent = 0, nPVertex = 0;
	float stepX, stepY;
	myfile >> scenario.minX >> scenario.maxX >> stepX >> scenario.minY >> scenario.maxY >> stepY;

	// TODO: Make this more efficient
	// Add points to the grid
	if (cfg::stepX > 0) stepX = cfg::stepX;
	if (cfg::stepY > 0) stepY = cfg::stepY;

	for (float i = scenario.minX; i < scenario.maxX; i += stepX) {
		for (float j = scenario.minY; j < scenario.maxY; j += stepY) {
			scenario.points.push_back(PointState(Point2D(i, j)));
		}
	}

	myfile >> scenario.packageInputMode;
	myfile >> nPackage;
	if (scenario.packageInputMode == SINGLE_POINT) {
		loadDesignatedPoint(myfile, scenario.packages, scenario, nPackage);
	}
	else if (scenario.packageInputMode == POLY) {
		loadDesignatedPolygon(myfile, scenario.packages, scenario, nPackage);
	}

	myfile >> scenario.targetInputMode;
	myfile >> nTarget;
	if (scenario.targetInputMode == SINGLE_POINT) {
		loadDesignatedPoint(myfile, scenario.targets, scenario, nTarget);
	}
	else if (scenario.targetInputMode == POLY) {
		loadDesignatedPolygon(myfile, scenario.targets, scenario, nTarget);
	}

	myfile >> nAgent;
	scenario.maxSpeed = -1;
	scenario.minSpeed = -1;
	for (int i = 0; i < nAgent; i++) {
		int x, y;
		float v;
		myfile >> x >> y >> v;
		if (scenario.maxSpeed == -1) {
			scenario.maxSpeed = v;
			scenario.minSpeed = v;
		}
		if (v > scenario.maxSpeed) {
			scenario.maxSpeed = v;
		}
		if (v < scenario.minSpeed) {
			scenario.minSpeed = v;
		}

		Agent a(i, x, y, v);
		scenario.agents.push_back(a);
	}
	// Sort the agents by speed
	std::sort(scenario.agents.begin(), scenario.agents.end());
}


void ScenarioIO::writeSolution(const char *outputFile, Scenario scenario) {
	std::ofstream myfile;
	myfile.open(outputFile);

	myfile << scenario.overallTime << std::endl;
	for (int a = 0; a < scenario.activeID.size(); a++) {
		int pairID = scenario.activeID[a];
		myfile << pairID << " " << scenario.bestTimes[a] << std::endl;
		myfile << scenario.bestAgentQueues[a].size() << std::endl;
		for (int k = 0; k < scenario.bestAgentQueues[a].size(); k++) {
			std::cout << "";
			myfile << scenario.bestAgentQueues[a][k].loc0.x << " " << 
				scenario.bestAgentQueues[a][k].loc0.y << " " <<
				scenario.bestAgentQueues[a][k].v << std::endl;
		}
		myfile << (scenario.bestPointQueues[a].size() + 1) << std::endl;
		for (int p = 0; p < scenario.bestPointQueues[a].size(); p++) {
			myfile << scenario.bestPointQueues[a][p].x << " " << scenario.bestPointQueues[a][p].y << std::endl;
		}
		myfile << scenario.bestTargets[a].loc.x << " " << scenario.bestTargets[a].loc.y << std::endl;
	}

	myfile.close();
}