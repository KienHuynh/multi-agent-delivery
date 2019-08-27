#pragma once
#include "base.h"

class Scenario {
public:
	// TODO: use bit flag
	// n drones deliver 1 package from point (s) to point (s)
	// n drones deliver m packages from point (s) to point (s), all designated points have the same roles
	// n drones deliver m packages from point (s) to point (s), each point is unique
	int problemType;
	int solverType;

	std::vector<Agent> agents;
	std::vector<std::vector<LineAnimation>> anis;

	float maxSpeed, minSpeed;
	int minX, maxX, minY, maxY;

	// Animation related members.
	float timer;
	bool aniStart;

	// The grid.
	std::vector<PointState> points;
	// Store the packages.
	std::vector<DesignatedPoint> packages;
	// Store the id of the point of the package.
	std::vector<DesignatedPoint> targets;

	Scenario();

	// Load scenario from file.
	void loadFile(const char*);

	// Check if a point is a package.
	// TODO: use better data structure?
	bool isPackage(int);

	// Euclidean 2D (problem) type 0, dynamic solver for 1 drone and 1 package.
	// Approximation of opt..=
	void ecld2DType0Dynamic11();

	// Euclidean 2D (problem) type 0, dynamic solver for n drones and m packages.
	// Approximation of opt.
	void ecld2DType0DynamicNM();

	// The actual body of ecld2DType0DynamicNM.
	// It is used so that it can be shared with more complicated functions.
	// @param[in] _agents the list of agents to be used.
	// @param[in] _packages the list of packages to be used. All of them share the same role.
	// @param[out] _points the grid (sample) points which will also be used to store solutions. If you don't want your data points to be overwritten, use of copy of your data points before calling this.
	void ecld2DType0DynamicNMCommon(
		std::vector<Agent> _agents,
		std::vector<DesignatedPoint> _packages,
		std::vector<PointState> &_points);

	// Euclidean 2D (problem) type 0, dynamic solver for n drones and m packages.
	// Approximation of opt.
	// @param[in] matchID the ID of the package-target matching we want to compute.
	void ecld2DType0DynamicNM(int matchID);

	// Euclidean 2D (problem) type 0, dynamic solver for n drones and m packages.
	// Approximation of opt.
	void ecld2DType1DynamicNM();

	// Create an animation based on the solution.
	void createAnimation();
	void solve();

private:
	// This is used to load target points or package points.
	// Load points such as packages and targets.
	// If a point already exists in the grid above, re-use it.
	// Otherwise, add a new point to the grid.
	// @param[in] myfile ifstream object with the input file already loaded.
	// @param[in] problemType problem type, see inputDescription.txt.
	// @param[out] dPoints the array of DesignatedPoint to store the package/target locations.
	// @param[out] nDPoint number of designated points.
	// @param[in] points the point grid.
	void loadDesignatedPoint(
		std::ifstream &myfile,
		int problemType,
		std::vector<DesignatedPoint> &dPoints,
		int nDPoint,
		std::vector<PointState> &points);

	// Load package regions and target regions as polygons.
	// These polygons will be discretized by their edges.
	// TODO: Allow this function to add polygons with different IDs.
	// @param[in] myfile ifstream object with the input file already loaded.
	// @param[in] problemType problem type, see inputDescription.txt.
	// @param[out] dPoints the array of DesignatedPoint to store the package/target locations.
	// @param[out] nDPoint number of designated points.
	// @param[in] points the point grid.
	void loadDesignatedPolygon(
		std::ifstream &myfile,
		int problemType,
		std::vector<DesignatedPoint> &dPoints,
		int nDPoint,
		std::vector<PointState> &points);
};