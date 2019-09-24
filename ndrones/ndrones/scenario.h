/*
 * scenario.h
 * This file contains the main scenario class
*/
#pragma once

#include <array>

#include "base.h"

enum DesignatedPointInputMode {
	SINGLE_POINT, POLY
};


// TODO: consider the following case
// There are multiple targets and IDs, some of the targets might have shared IDs. 
// We must deliver the package (does not matter where they come from) to EACH of these target.
enum ProblemType {
	UNKNOWN = 0,
	ONEDIM = 1,
	TWODIM = 2,
	EUCLID = 4,
	GRAPH = 8,
	DISCRETE = 16,
	// CONTINUOUS: 5th bit = 0
	SINGLE_TARGET = 32,
	//MULTI_TARGET: 6th bit = 0
	SINGLE_ID = 64
};


std::istream& operator >> (std::istream &input, ProblemType& p);


std::istream& operator >> (std::istream &input, DesignatedPointInputMode& f);

// The core class of the project
// Store the scenario: data points, designated points (packages, targets), agents
// Include the solvers to the problems and creation of animations
class Scenario {
public:
	// TODO: use bit flag
	// n drones deliver 1 package from point (s) to point (s)
	// n drones deliver m packages from point (s) to point (s), all designated points have the same roles
	// n drones deliver m packages from point (s) to point (s), each point is unique
	ProblemType problemType;
	int solverType;

	std::vector<Agent> agents;
	std::vector<LineAnimation> droneAnis;
	std::vector<LineAnimation> packageAnis;

	float maxSpeed, minSpeed;
	int minX, maxX, minY, maxY;
	float makespan;

	// Animation related members.
	float timer;
	bool aniStart;

	// The grid.
	std::vector<PointState> points;
	// Store the packages.
	std::vector<DesignatedPoint> packages;
	// Store the targets.
	std::vector<DesignatedPoint> targets;
	// Store the list of matching id (s)
	std::vector<int> activeID;

	// Specify whether the inputs will be in regional (polygon) format or discrete points format
	DesignatedPointInputMode packageInputMode;
	DesignatedPointInputMode targetInputMode;

	// Store the package polygon
	std::vector<Point2D> targetPoly;
	// Store the target polygon
	std::vector<Point2D> packagePoly;

	std::string outputFileName;

	// Variables storing the solution
	std::vector<Agent>* bestAgentQueues;
	std::vector<Point2D>* bestPointQueues;
	DesignatedPoint* bestTargets;
	float* bestTimes;
	float overallTime;

	Scenario();

	// Load scenario from file.
	void loadFile(const char*);

	// Write solution to a file
	void writeSolution(const char* outputFile);

	// Check if a point is one of the packages
	// TODO: use better data structure?
	// @param[in] int index of a grid point.
	// @param[in] _packages the list of packages.
	// @return bool
	bool isPackage(int, std::vector<DesignatedPoint>);

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
	// This will only use packages and targets of matchID to compute.
	// Approximation of opt.
	// @param[in] matchID the ID of the package-target matching we want to compute.
	// @param[in] _agents the list of agents to be used.
	// @param[out] _points the grid (sample) points which will also be used to store solutions. If you don't want your data points to be overwritten, use of copy of your data points before calling this.
	void ecld2DType0DynamicNM(
		int matchID,
		std::vector<Agent> _agents,
		std::vector<PointState> &_points);

	// Euclidean 2D (problem) type 0, dynamic solver for n drones and m packages.
	// Approximation of opt.
	void ecld2DType1DynamicNM();

	// Create an animation based on the solution for the drones.
	void createDroneAnimation();
	// Create an animation for the flight paths of the packages.
	void createPackageAnimation();
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
		ProblemType problemType,
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
	// @param[out] points the point grid.
	// @param[out] poly the polygon to be stored for drawing.
	void loadDesignatedPolygon(
		std::ifstream &myfile,
		ProblemType problemType,
		std::vector<DesignatedPoint> &dPoints,
		int nDPoint,
		std::vector<PointState> &points,
		std::vector<Point2D> &poly);

	// There might be multiple targets of the same matching / role (i.e. our agents can deliver one package to any of them
	// Find the best delivery among all targets
	// @param[in] _points the point grid
	// @param[in] _targets the list of targets
	// @return bestTime float
	float findBestTimeFromTargets(std::vector<PointState> _points, std::vector<DesignatedPoint> _targets);

	// There might be multiple targets of the same matching / role (i.e. our agents can deliver one package to any of them
	// Find the target with the best time
	// @param[in] _points the point grid
	// @param[in] _targets the list of targets
	// @return bestTarget DesignatedPoint
	DesignatedPoint findBestTarget(std::vector<PointState> _points, std::vector<DesignatedPoint> _targets);

	// Utility function, find the max value of an array without considering the k-th element
	// @param[in] arr the pointer to the array
	// @param[in] size it's size
	// @param[in] k the index to be avoided
	// @return float
	float maxValWithoutK(float *arr, int size, int k);

	void updateReusedAgent(std::vector<Agent> & _agents, std::vector<Point2D> _pointQueue);

	// Check if a list of agents contain another agent using their ID
	bool containAgentID(std::vector<Agent> _agents, Agent a);

	// Check if a list of agents contain another agent (every variable has to be the same)
	bool containAgent(std::vector<Agent> _agents, Agent a);

	// Resolving conflict between all package-target matchings that share some agents.
	// @param[in] _points the point grid.
	// @param[in] _agents the agent list.
	// @param[in] _packageOfID array of list of packages, separated based on their IDs.
	// @param[in] _targetsOfID array of list of targets, separated based on their IDs.
	// @param[in] _ativeID list of active ID.
	// @param[in] _agentAssignment matching-agent assignment table, 0 at [i][j] means agent j was assigned to the i-th matching, 1 if it is assigned.
	// @param[in] _bestTimes list of best time computed for each matching when all agents are available for it
	// @param[out] bestMatchingInd the best matching that will be used after all conflicts are resolved.
	// @param[out] bestAgentInd the best agent to be assigned to the above matching
	// @return isConflict true if there was a conflict, false if there was not.
	bool conflictResolve(
		std::vector<PointState> _points,
		std::vector<Agent> _agents,
		std::vector<DesignatedPoint>* _packagesOfID,
		std::vector<DesignatedPoint>* _targetsOfID,
		std::vector<int> _activeID,
		int** _agentAssignment,
		float* _bestTimes,
		int &bestMatchingInd,
		int &bestAgentInd);

	// Given a vector of agents _agents, find the location of Agent a in _agents, using its ID
	int findVectorIndexWithID(std::vector<Agent> _agents, Agent a);
	// Given a vector of agents _agents, find the location of Agent a in _agents, using its ID and delayTime
	int findVectorIndexFull(std::vector<Agent> _agents, Agent a);

	// Remove any agents in _agents if it also exists in queues[i]
	// However, if i == id, then no one will be removed
	void removeSharedAgents(std::vector<Agent>* queues, int id, std::vector<Agent> &_agents);

	// Return true if a.orderOfEx > b.orderOfEx
	static bool compareOrderOfEx(Agent a, Agent b);

	// Remove all agents in the list with orderOfEx greater than _orderOfEx
	void removeGapAgents(std::vector<Agent> &_agents, int _id, int _orderOfEx);

	// Check if this list contains agent a with bigger orer of execution
	bool containAgentAfterOrder(std::vector<Agent> _agents, Agent a);

	// Bag agents into some queues and sort them by order of executions
	void bagAgentsByOrder(std::vector<Agent>* _agentQueues, std::vector<Agent>*& bag);

	// Check if two queues of agents are equal
	bool equalAgentQueue(std::vector<Agent>* qA, std::vector<Agent>* qB);

	// Check if there is an empty queue here
	bool missingQueue(std::vector<Agent>* qs);
};