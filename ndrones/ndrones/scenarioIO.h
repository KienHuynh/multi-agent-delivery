#pragma once
#include "scenario.h"


class Scenario;
enum ProblemType;
enum DesignatedPointInputMode;

// APGRID: axis-parallel grid
// GRID: s-t aligned grid
// RECTUNION: union of rectangles
// CIRCULAR: circular sampling method
// LOGGRID: s-t aligned log grid
enum SamplingMethod {
	UNSPECIFIED = 0,
	APGRID = 1,
	GRID = 2,
	RECTUNION = 3,
	CIRCULAR = 4,
	LOGGRID = 5
};

// Some input operators for ProblemType and DesignatedPointInputMode to make life easier
std::istream& operator >> (std::istream &input, ProblemType& p);
std::istream& operator >> (std::istream &input, DesignatedPointInputMode& f);


class ScenarioIO {
public:
	// Load scenario from file.
	// @param[in] const char* fileName.
	// @param[out] Sceanrio scenario.
	static void loadFile(const char*, Scenario &scenario);

	// Write solution to a file.
	// @param[in] const char* fileName.
	// @param[out] Sceanrio scenario.
	static void writeSolution(const char* outputFile, Scenario scenario);

	// This is used to load target points or package points.
	// Load points such as packages and targets.
	// If a point already exists in the grid above, re-use it.
	// Otherwise, add a new point to the grid.
	// @param[out] myfile ifstream object with the input file already loaded.
	// @param[out] std::vector<DesignatedPoint> dPpoints the point list to store the DesignatedPoint (s).
	// @param[in] Sceanrio scenario.
	// @param[in] int nDPoint number of designated points to be loaded (also specified in the input file).
	static void loadDesignatedPoint(
		std::ifstream &myfile,
		std::vector<DesignatedPoint> &dPpoints,
		Scenario &scenario,
		int nDPoint);

	// Load package regions and target regions as polygons.
	// These polygons will be discretized by their edges.
	// @param[out] myfile ifstream object with the input file already loaded.
	// @param[out] std::vector<DesignatedPoint> dPpoints the point list to store the DesignatedPoint (s).
	// @param[in] Scenario scenario.
	// @param[in] int nDPoint number of designated points to be loaded (also specified in the input file).
	static void loadDesignatedPolygon(
		std::ifstream &myfile,
		std::vector<DesignatedPoint> &dPpoints,
		Scenario &scenario,
		int nDPoint);

	//
	static void loadObstacle(
		std::ifstream &myfile,
		std::vector<SimplePolygon> &obstacles,
		Scenario &scenario,
		int nObs
	);


private:
	// Generate points using grid scheme. This grid is axis parallel.
	// @param[out] std::ifstream &myfile.
	// @param[out] Scenario scenario, storing the points.
	static void generateAPGrid(std::ifstream &myfile, Scenario&);

	// Generate points using a better grid scheme.
	// This grid is limited to the union of package-target strips.
	// @param[out] std::ifstream &myfile.
	// @param[out] Scenario scenario, storing the points.
	static void generateGrid(std::ifstream &myfile, Scenario&);

	// Generate points using circular grid scheme
	// @param[out] std::ifstream &myfile.
	// @param[out] Scenario scenario, storing the points
	static void generateCircularGrid(std::ifstream &myfile, Scenario&);

	// Generate points using log grid scheme
	// @param[out] std::ifstream &myfile.
	// @param[out] Scenario scenario, storing the points
	static void generateLogGrid(std::ifstream &myfile, Scenario&);
};