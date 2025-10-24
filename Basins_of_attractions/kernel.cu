#include "basinsHOST.h"
#include "systems.cuh"

#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <conio.h>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>

const std::string BASINS_OUTPUT_PATH = "C:/Users/user/Documents";

void runBasinsPerformanceTests() {
	std::ofstream resultsFile(std::string(BASINS_OUTPUT_PATH) + "/basins_performance_results.csv");
	if (!resultsFile) {
		std::cerr << "Error: Failed to open file for writing results!" << std::endl;
		return;
	}

	resultsFile << "Parameter,CT,Resolution,Library,ExecutionTime_ms" << std::endl;

	std::vector<int> resolutions = { 100,200,400,600,800,1000 };  // Resolution tests
	std::vector<int> modelingTimes = { 500, 1000, 1500, 2000, 2500, 3000 };  // Simulation time tests

	numb params[5]{ 0.5, 0.1665, 1.4,  15.552, 2 };
	numb init[3]{ 0, 0, 0, };
	numb ranges[4]{ -6, 6, -6, 6 };
	int indicesOfMutVars[2]{ 0, 1 };
	const int custom_block_size = 32;

	std::cout << "\n===== Test 1: Influence of resolution on execution time =====\n";

	for (int modelingTime : modelingTimes) {

		for (int resolution : resolutions) {
			std::cout << "\nTesting with resolution = " << resolution << " ModelingTime = " << modelingTime << std::endl;

			std::cout << "  Running Basins::basinsOfAttraction_2..." << std::endl;
			auto start1 = std::chrono::high_resolution_clock::now();
			long long duration1 = 0;

			try {
				int time[3];
				Basins::basinsOfAttraction_2(
					500,                // System simulation time
					resolution,         // Diagram resolution
					0.01,               // Integration step
					sizeof(init) / sizeof(numb),   // Number of initial conditions
					init,               // Array of initial conditions
					ranges,
					indicesOfMutVars,
					1,                  // Equation index for the diagram
					100000000,          // Maximum value
					modelingTime,       // Time to simulate
					params,             // Parameters
					sizeof(params) / sizeof(numb),  // Number of parameters
					1,                  // Multiplier
					0.05,               // Epsilon for DBSCAN
					custom_block_size,
					std::string(BASINS_OUTPUT_PATH) + "/basins_res_test_" + std::to_string(resolution) + ".csv",
					time
				);

				auto end1 = std::chrono::high_resolution_clock::now();
				duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();

				resultsFile << "Resolution, " << modelingTime << "," << resolution << ",Basins," << duration1 << std::endl;
				std::cout << "    Execution time: " << duration1 << " ms" << std::endl;
			}
			catch (const std::exception& e) {
				std::cerr << "    Error: " << e.what() << std::endl;
				resultsFile << "Resolution with CT = " << modelingTime << "," << resolution << ",Basins,ERROR" << std::endl;
			}

			std::cout << "  Running old_library::basinsOfAttraction_2..." << std::endl;
			auto start2 = std::chrono::high_resolution_clock::now();

			resultsFile.flush();
		}

	}

	resultsFile.close();
	std::cout << "\nPerformance tests completed. Results saved to performance_results.csv" << std::endl;
}

int main()
{
	size_t startTime = std::clock();
	numb h = (numb)0.01;

#ifdef USE_SYSTEM_FOR_BASINS
	numb params[5]{ 0.5, 0.1665, 1.4,  15.552, 2 };
	numb init[3]{ 0, 0, 0, };
	numb ranges[4]{ -6, 6, -6, 6 };
	int indicesOfMutVars[2]{ 0, 1 };
	const int custom_block_size = 256;
	//runBasinsPerformanceTests();
	 {
	  std::cout << "Start basins" << std::endl;
	  auto start = std::chrono::high_resolution_clock::now();
	  int time[3];
	  Basins::basinsOfAttraction_2(
	  	700,       // CT
	  	300,       // Resolution
	  	h,         // time step
	  	sizeof(init) / sizeof(numb),   // amount of init conditions
	  	init,         // init conditions
	  	ranges,			// parameters range
	  	indicesOfMutVars, // indices of butual variables
	 	1,          // Index of the equation to use for plotting the diagram
	 	100000000,  // Maximum value (by absolute value); above this the system is considered "diverged"
	 	500,       // Time that will be simulated before computing the diagram
	 	params,     // Parameters
	 	sizeof(params) / sizeof(numb),  // Number of parameters
	 	1,          // Multiplier that reduces time and computation load (only every 'preScaller' point will be computed)
	 	0.05,       // Epsilon for the DBSCAN algorithm
		custom_block_size,
	  	std::string(BASINS_OUTPUT_PATH) + "/bas.csv",
		time
	  );
	  auto end = std::chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	  std::cout << "Time taken: " << duration << " milliseconds" << std::endl;
	  std::cout << "Time taken: for system " << time[0] << " ms --- for dbscan " << time[1] << " ms --- at all " << time[2] << " ms" << std::endl;

	  }

	//runPerformanceTests();
#endif

#ifdef USE_SYSTEM_FOR_BASINS_2
	numb params[5]{ 0.5, 0.1, 1.4,  15.552, 2 };
	numb init[3]{ 0, 0, 0, };

	{
		std::cout << "Start basins" << std::endl;
		auto start = std::chrono::high_resolution_clock::now();
		runBasinsPerformanceTests();

		//Basins::basinsOfAttraction_2(
		//	200,        // Simulation time of the system
		//	100,        // Diagram resolution
		//	0.01,       // Integration step
		//	sizeof(init) / sizeof(numb),   // Number of initial conditions (equations in the system)
		//	init,       // Array of initial conditions
		//	new numb[4] { -200, 200, -60, 60 },
		//	new int[2] { 0, 1 },
		//	1,          // Index of the equation to use for plotting the diagram
		//	100000000,  // Maximum value (by absolute value); above this the system is considered "diverged"
		//	1000,       // Time that will be simulated before computing the diagram
		//	params,     // Parameters
		//	sizeof(params) / sizeof(numb),  // Number of parameters
		//	1,          // Multiplier that reduces time and computation load (only every 'preScaller' point will be computed)
		//	0.05,       // Epsilon for the DBSCAN algorithm
		//	std::string(BASINS_OUTPUT_PATH) + "/bas_2.csv"
		//);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "Time taken: " << duration << " milliseconds" << std::endl;
	}
#endif


	std::cout << "Time taken: " << (std::clock() - startTime) / (numb)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

	return 0;
}