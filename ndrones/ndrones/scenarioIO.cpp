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


std::istream& operator >> (std::istream &input, SamplingMethod& p) {
	std::string tmp;
	input >> tmp;
	if (stringEqual(tmp, "apgrid")) p = APGRID;
	if (stringEqual(tmp, "grid")) p = GRID;
	if (stringEqual(tmp, "circular")) p = CIRCULAR;

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
	SamplingMethod sm = UNSPECIFIED;
	while (ss >> tmp) {
		scenario.problemType = (ProblemType)(scenario.problemType | tmp);
	}

	int nPackage = 0, nTarget = 0, nAgent = 0, nPVertex = 0;

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

	// TODO: Make this more efficient
	// Add points to the grid
	myfile >> sm;
	if (sm == SamplingMethod::APGRID) {
		generateAPGrid(myfile, scenario);
	}
	if (sm == SamplingMethod::GRID) {
		generateGrid(myfile, scenario);
	}
	if (sm == SamplingMethod::CIRCULAR) {
		generateCircularGrid(myfile, scenario);
	}
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

// TODO: Don't generate points outside of the convex hull (?) for both grid method
void ScenarioIO::generateAPGrid(std::ifstream &myFile, Scenario& scenario) {
	float stepX, stepY;
	myFile >> scenario.minX >> scenario.maxX >> stepX >> scenario.minY >> scenario.maxY >> stepY;
	if (cfg::stepX > 0) stepX = cfg::stepX;
	if (cfg::stepY > 0) stepY = cfg::stepY;

	for (float i = scenario.minX; i < scenario.maxX; i += stepX) {
		for (float j = scenario.minY; j < scenario.maxY; j += stepY) {
			scenario.points.push_back(PointState(Point2D(i, j)));
		}
	}
}


void ScenarioIO::generateGrid(std::ifstream &myFile, Scenario& scenario) {
	int nH, nW;
	myFile >> nH >> nW;
	
	float stDistance = Point2D::l2Distance(scenario.packages[0].loc, scenario.targets[0].loc);
	float longestDis = 0;
	for (auto agent : scenario.agents) {
		for (auto target : scenario.targets) {
			if (longestDis < Point2D::l2Distance(agent.loc, target.loc))
				longestDis = Point2D::l2Distance(agent.loc, target.loc);
		}
		for (auto package : scenario.packages) {
			if (longestDis < Point2D::l2Distance(agent.loc, package.loc))
				longestDis = Point2D::l2Distance(agent.loc, package.loc);
		}
	}
	for (auto package : scenario.packages) {
		for (auto target : scenario.targets) {
			if (longestDis < Point2D::l2Distance(package.loc, target.loc))
				longestDis = Point2D::l2Distance(package.loc, target.loc);
		}
	}

	Point2D s = scenario.packages[0].loc;
	Point2D t = scenario.targets[0].loc;
	Point2D v = (t - s) / Point2D::l2norm(t - s);
	float theta = atan(v.y / v.x); // For grid rotation
	float y0 = -longestDis;
	float x0 = 0;
	scenario.minX = s.x;
	scenario.minY = s.y;
	scenario.maxX = s.x;
	scenario.maxY = s.y;
	for (int i = 0; i < nH; i++) {
		for (int j = 0; j < nW; j++) {
			float x1 = x0 + ((float)j / (float)(nW-1)) * stDistance;
			float y1 = y0 + ((float)i / (float)nH) * 2 * longestDis;
			float x, y;
			// Rotate & Translate
			x = cos(theta)*x1 - sin(theta)*y1 + s.x;
			y = sin(theta)*x1 + cos(theta)*y1 + s.y;
			scenario.points.push_back(Point2D(x, y));

			if (x < scenario.minX) scenario.minX = x;
			if (x > scenario.maxX) scenario.maxX = x;
			if (y < scenario.minY) scenario.minY = y;
			if (y > scenario.maxY) scenario.maxY = y;
		}
		std::cout << std::endl;
	}
}


void ScenarioIO::generateCircularGrid(std::ifstream &myFile, Scenario& scenario) {
	float nTheta;
	float nR;
	myFile >> nTheta >> nR;
	if (cfg::nTheta > 0) nTheta = cfg::nTheta;
	if (cfg::nR > 0) nR = cfg::nR;

	// Compare pairwise distance of every relevant pairs
	// TODO: if there might be more pairs, use the O(nlogn) algorithm
	float longestDis = 0;
	for (auto agent : scenario.agents) {
		for (auto target : scenario.targets) {
			if (longestDis < Point2D::l2Distance(agent.loc, target.loc))
				longestDis = Point2D::l2Distance(agent.loc, target.loc);
		}
		for (auto package : scenario.packages) {
			if (longestDis < Point2D::l2Distance(agent.loc, package.loc))
				longestDis = Point2D::l2Distance(agent.loc, package.loc);
		}
	}
	for (auto package : scenario.packages) {
		for (auto target : scenario.targets) {
			if (longestDis < Point2D::l2Distance(package.loc, target.loc))
				longestDis = Point2D::l2Distance(package.loc, target.loc);
		}
	}
	float r = longestDis / nR;

	std::vector<Point2D> centerList;
	for (auto agent : scenario.agents) {
		centerList.push_back(agent.loc);
	}
	for (auto package : scenario.packages) {
		centerList.push_back(package.loc);
	}
	for (auto target : scenario.targets) {
		centerList.push_back(target.loc);
	}

	scenario.minX = centerList[0].x;
	scenario.minY = centerList[0].y;
	scenario.maxX = centerList[0].x;
	scenario.maxY = centerList[0].y;
	
	float dTheta = 360.0 / (float)nTheta;
	for (auto center : centerList) {
		for (float theta_i = 0; theta_i < 360; theta_i += dTheta) {
			for (float r_i = r; r_i < (longestDis + r); r_i += r) {
				float x = cos(theta_i * PI / 180.0)*r_i + center.x;
				float y = sin(theta_i * PI / 180.0)*r_i + center.y;
				scenario.points.push_back(Point2D(x, y));

				if (x < scenario.minX) scenario.minX = x;
				if (x > scenario.maxX) scenario.maxX = x;
				if (y < scenario.minY) scenario.minY = y;
				if (y > scenario.maxY) scenario.maxY = y;
			}
		}
	}
}