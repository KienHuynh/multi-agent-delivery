#pragma once
namespace cfg {
	const float timeStep = 0.1;

	// Step size for sampling the points on polygon edges
	const int polySamplingRate = 3;

	// Step size for sampling the grid points
	// If these two are > 0, they will override the values in the input file
	// For debug purpose only
	const float stepX = -1;
	const float stepY = -1;
	// Circular grid parameters
	const float theta = -1;
	const float nR = -1;

	// Number of colors for the palette
	const int numColor = 2048;
}