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
	if (stringEqual(tmp, "loggrid")) p = LOGGRID;

	return input;
}


std::istream& operator >> (std::istream &input, ObstacleType& t) {
	std::string tmp;
	input >> tmp;
	if (stringEqual(tmp, "type1")) t = TYPE1;
	if (stringEqual(tmp, "type2")) t = TYPE2;
	if (stringEqual(tmp, "type3")) t = TYPE3;
	if (stringEqual(tmp, "type4")) t = TYPE4;
	if (stringEqual(tmp, "type5")) t = TYPE5;

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

		scenario.updateMinMaxXY(scenario.points[scenario.points.size() - 1].p.x,
			scenario.points[scenario.points.size() - 1].p.y);
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
			scenario.updateMinMaxXY(scenario.points[scenario.points.size() - 1].p.x,
				scenario.points[scenario.points.size() - 1].p.y);
		}
	}
}


void ScenarioIO::loadAgent(std::ifstream &myfile,
	Scenario &scenario,
	int nAgents) {
	std::string str;
	getline(myfile, str);
	std::istringstream ss(str);

	scenario.maxSpeed = -1;
	scenario.minSpeed = -1;

	for (int i = 0; i < nAgents; i++) {
		float x, y;
		float v;
		ObstacleType obsType = TYPE0;
		int combinedObsType = TYPE0;

		std::string str;
		getline(myfile, str);
		std::istringstream ss(str);
		
		ss >> x >> y >> v;
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

		while (ss >> obsType) {
			combinedObsType = obsType | combinedObsType;
		}

		Agent a(i, x, y, v, combinedObsType);

		Point2D p(x, y);
		if (!scenario.containPoint(p)) {
			scenario.points.push_back(p);
			a.gridRef = scenario.points.size() - 1;
		}
		scenario.agents.push_back(a);
		scenario.updateMinMaxXY(x, y);
	}
}


void ScenarioIO::loadObstacle(std::ifstream &myfile,
	Scenario &scenario,
	int nObs) {

	std::vector<SimplePolygon> obstacles;
	for (int i = 0; i < nObs; i++) {
		std::vector<Point2D> points;
		int nPoint = 0;
		ObstacleType t;
		myfile >> nPoint >> t;

		for (int j = 0; j < nPoint; j++) {
			Point2D p;
			myfile >> p.x >> p.y;
			points.push_back(p);
		}

		SimplePolygon polygon(points);
		// Triangulate it (for drawing and coloring, mostly)
		polygon.triangulate();
		// Find convex hull
		polygon.findCVHull();
		polygon.type = t;

		obstacles.push_back(polygon);
	}
	scenario.obs = obstacles;

	// Eliminate points inside the obstacles
	for (int oi = 0; oi < obstacles.size(); oi++) {
		SimplePolygon o = obstacles[oi];
		int count = 0;
		for (int i = scenario.points.size()-1; i >= 0; i--)
		{
			if (o.contain(scenario.points[i].p)) {
				scenario.points.erase(scenario.points.begin() + i);
				count++;
			}
		}

		// Add points on the boundary of each obstacle as potential handoff points
		float nNewPoint = sqrt((float)count);
		for (int i = 0; i < o.points.size(); i++) {
			Point2D p0(o.points[i].x, o.points[i].y);
			Point2D p1;
			p1 = o.points[(i + 1) % o.points.size()];

			if (!scenario.containPoint(p0)) {
				PointState ps(p0);
				ps.isOb = true;
				ps.obIndex = oi;
				scenario.points.push_back(ps);
			}
			
			//float length = Point2D::l2Distance(p0, p1);
			//int nNewPoint_i = (int) (length * nNewPoint / o.perimeter);
			//for (int j = 0; j < nNewPoint_i; j++) {
			//	float ratio = ((float)(j + 1)) / ((float)(nNewPoint_i + 2));
			//	float newX = p0.x + ratio * (p1.x - p0.x);
			//	float newY = p0.y + ratio * (p1.y - p0.y);
			//	Point2D newP(newX, newY);

			//	/*float eps = o.perimeter / 100;
			//	for (int r = 0; r < 4; r++) {
			//		Point2D newPtmp = newP;
			//		if (r == 0) {
			//			newPtmp.x += eps;
			//			newPtmp.y += eps;
			//		}
			//		if (r == 1) {
			//			newPtmp.x -= eps;
			//			newPtmp.y += eps;
			//		}
			//		if (r == 2) {
			//			newPtmp.x += eps;
			//			newPtmp.y -= eps;
			//		}
			//		if (r == 3) {
			//			newPtmp.x -= eps;
			//			newPtmp.y -= eps;
			//		}
			//		if (!o.contain(newPtmp)) {
			//			newP = newPtmp;
			//			break;
			//		}
			//	}*/

			//	if (!scenario.containPoint(newP)) {
			//		PointState ps(Point2D(newX, newY));
			//		ps.isOb = true;
			//		// ps.obIndex = oi;
			//		scenario.points.push_back(ps);
			//	}
			//}
		}
	}

	// Construct the graph and pre-compute pairwise distance between all necessary points
	scenario.constructSTPMap();
}


void ScenarioIO::loadFile(const char* fname, Scenario &scenario) {
	scenario.outputFileName = std::string(fname);
	scenario.outputFileName = scenario.outputFileName.substr(scenario.outputFileName.find_last_of("/\\") + 1);
	std::string::size_type const p(scenario.outputFileName.find_last_of('.'));
	scenario.outputFileName = scenario.outputFileName.substr(0, p);

	int nPackage = 0, 
		nTarget = 0, 
		nAgent = 0, 
		nPVertex = 0,
		nObs = 0;

	std::ifstream myfile;
	myfile.open(fname, std::ios::in);

	// Reading problem type
	scenario.problemType = ProblemType::UNKNOWN;
	std::string str;
	getline(myfile, str);
	std::istringstream ss(str);

	ProblemType pt;
	SamplingMethod sm = SamplingMethod::UNSPECIFIED;
	while (ss >> pt) {
		scenario.problemType = (ProblemType)(scenario.problemType | pt);
	}


	// Read package input mode
	myfile >> scenario.packageInputMode;
	myfile >> nPackage;
	// Read packages' info
	if (scenario.packageInputMode == SINGLE_POINT) {
		loadDesignatedPoint(myfile, scenario.packages, scenario, nPackage);
	}
	else if (scenario.packageInputMode == POLY) {
		loadDesignatedPolygon(myfile, scenario.packages, scenario, nPackage);
	}

	// Read target input mode
	myfile >> scenario.targetInputMode;
	myfile >> nTarget;
	// Read targets' info
	if (scenario.targetInputMode == SINGLE_POINT) {
		loadDesignatedPoint(myfile, scenario.targets, scenario, nTarget);
	}
	else if (scenario.targetInputMode == POLY) {
		loadDesignatedPolygon(myfile, scenario.targets, scenario, nTarget);
	}

	// Read agents' info
	myfile >> nAgent;
	loadAgent(myfile, scenario, nAgent);

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
	if (sm == SamplingMethod::LOGGRID) {
		generateLogGrid(myfile, scenario);
	}

	std::string tmp;
	myfile >> tmp;
	myfile >> nObs;
	loadObstacle(myfile, scenario, nObs);

	

	myfile.close();
}


void ScenarioIO::writeSolution(const char *outputFile, Scenario scenario) {
	std::ofstream myfile;
	myfile.open(outputFile);

	// TODO
	/*myfile << scenario.overallTime << std::endl;
	myfile << scenario.points.size() << std::endl;
	for (int a = 0; a < scenario.activeID.size(); a++) {
		int pairID = scenario.activeID[a];
		myfile << pairID << " " << scenario.bestTimes[a] << std::endl;
		myfile << scenario.bestAgentQueues[a].size() << std::endl;
		for (int k = 0; k < scenario.bestAgentQueues[a].size(); k++) {
			myfile << scenario.bestAgentQueues[a][k].loc0.x << " " << 
				scenario.bestAgentQueues[a][k].loc0.y << " " <<
				scenario.bestAgentQueues[a][k].v << std::endl;
		}
		myfile << (scenario.bestPointQueues[a].size() + 1) << std::endl;
		for (int p = 0; p < scenario.bestPointQueues[a].size(); p++) {
			myfile << scenario.bestPointQueues[a][p].x << " " << scenario.bestPointQueues[a][p].y << std::endl;
		}
		myfile << scenario.bestTargets[a].loc.x << " " << scenario.bestTargets[a].loc.y << std::endl;
	}*/

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
	float maxY;
	myFile >> nH >> nW >> maxY;
	float longestDis = 0;
	float stDistance = Point2D::l2Distance(scenario.packages[0].loc, scenario.targets[0].loc);
	if (maxY <= 0) {
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
	}
	else {
		longestDis = maxY;
	}
	

	Point2D s = scenario.packages[0].loc;
	Point2D t = scenario.targets[0].loc;
	Point2D v = (t - s) / Point2D::l2norm(t - s);
	float theta = atan(v.y / v.x); // For grid rotation
	float y0 = -longestDis;
	float x0 = 0;
	scenario.updateMinMaxXY(s.x, s.y);

	for (int i = 0; i < nH; i++) {
		for (int j = 0; j < nW; j++) {
			float x1 = x0 + ((float)j / (float)(nW-1)) * stDistance;
			float y1 = y0 + ((float)i / (float)nH) * 2 * longestDis;
			float x, y;
			// Rotate & Translate
			x = cos(theta)*x1 - sin(theta)*y1 + s.x;
			y = sin(theta)*x1 + cos(theta)*y1 + s.y;
			scenario.points.push_back(Point2D(x, y));
			//std::cout << x << " " << y << std::endl;
			scenario.updateMinMaxXY(x, y);
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

	scenario.updateMinMaxXY(centerList[0].x, centerList[0].y);
	
	float dTheta = 360.0 / (float)nTheta;
	for (auto center : centerList) {
		for (float theta_i = 0; theta_i < 360; theta_i += dTheta) {
			for (float r_i = r; r_i < (longestDis + r); r_i += r) {
				float x = cos(theta_i * PI / 180.0)*r_i + center.x;
				float y = sin(theta_i * PI / 180.0)*r_i + center.y;
				scenario.points.push_back(Point2D(x, y));
				scenario.updateMinMaxXY(x, y);
			}
		}
	}
}


void ScenarioIO::generateLogGrid(std::ifstream &myFile, Scenario& scenario) {
	int nH, nW;
	float f, maxY;
	myFile >> nH >> nW >> f >> maxY;
	float longestDis = 0;
	float stDistance = Point2D::l2Distance(scenario.packages[0].loc, scenario.targets[0].loc);
	if (maxY <= 0) {
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
	}
	else {
		longestDis = maxY;
	}

	Point2D s = scenario.packages[0].loc;
	Point2D t = scenario.targets[0].loc;
	Point2D v = (t - s) / Point2D::l2norm(t - s);
	float theta = atan(v.y / v.x); // For grid rotation
	float y0 = 0;
	float x0 = 0;
	scenario.updateMinMaxXY(s.x, s.y);
	
	for (int j = 0; j < nW; j++) {
		bool stop = false;
		for (int i = -nH/2; i <= nH /2; i++) {
			float x1 = x0 + ((float)j / (float)(nW - 1)) * stDistance;
			float y1 = 0;
			if (i != 0) y1 = y0 + (i/abs(i)) * ((1-pow(f,abs(i)-1)) / (1-f)) / (float)nH * longestDis;
			
			if (abs(y1) > longestDis) continue;
			
			float x, y;
			// Rotate & Translate
			x = cos(theta)*x1 - sin(theta)*y1 + s.x;
			y = sin(theta)*x1 + cos(theta)*y1 + s.y;
			scenario.points.push_back(Point2D(x, y));
			scenario.updateMinMaxXY(x, y);
		}
	}
}