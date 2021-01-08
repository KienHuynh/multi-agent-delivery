#include <iostream>
#include <cxxopts.hpp>

// Local includes
#include "base.h"
#include "config.h"
#include "scenario.h"

int main(int argc, char **argv) {
	cxxopts::Options options(argv[0], "Options");
	options
		.positional_help("[optional args]")
		.show_positional_help();

	options.add_options()
		("i, input", "Path to input file", cxxopts::value<std::string>())
		("o, output", "Path to output file", cxxopts::value<std::string>())
		;
	auto cmd = options.parse(argc, argv);

	std::string inputFile = cmd["i"].as<std::string>();
	std::string outputFile = cmd["o"].as<std::string>();

	Scenario scenario;
	scenario.loadFile(inputFile.c_str());
	scenario.solve();
	scenario.writeSolution(outputFile.c_str());

}
