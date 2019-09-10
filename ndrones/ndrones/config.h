#pragma once
namespace cfg {
	const float timeStep = 0.1;

	// Step size for sampling the points on polygon edges
	const int polySamplingRate = 3;

	// Step size for sampling the grid points
	// If these two are > 0, they will override the values in the input file
	const float stepX = -1;
	const float stepY = -1;
}