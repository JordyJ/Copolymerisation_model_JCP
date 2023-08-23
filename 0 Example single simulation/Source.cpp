#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include <math.h>
#include <thread>
#include "PolymerSystemIO.h"
// Include windows.h if compiling and running on windows
//#include <windows.h>
using namespace std;
/* 
Developed by Jordan Juritz
Imperial College London
October 2021 */

int main(int argc, char* argv[]) {
	
	cout << argv[0] << endl; // Read in input file containing parameters

	PolymerSystemIO PolySystem; // Create a polymer system object

	if (argc > 1) {
		PolySystem.inputfilename = argv[1]; 
	}

	PolySystem.setup(); // Initialise the simulation
	
	// Optional below: 
	// Can be used to generate a non-empty initial condition. 
	// The following generates a dense brush of polymers
	//for (size_t i = 0; i < PolySystem.length_limit; i++)
	//{
	//	PolySystem.generate_initial_condition(0, i, 1, i); 
	//}
	
	
	//PolySystem.generate_initial_condition(0, 9, 8, 0);
	
	//PolySystem.generate_initial_condition(0, PolySystem.length_limit-1, PolySystem.length_limit -1, 0);
	
	double trip_time = 0;
	
	// Print the initial state of the sim
	cout << "\n" << PolySystem.time << "%\n";//
	PolySystem.print_state_vertically();
	
	// Simulation termination conditions
	//while (PolySystem.time >= 0 && PolySystem.time < PolySystem.time_limit) // Use this line to terminate at a time time_limit
	while (PolySystem.time >= 0 && PolySystem.pols_written < PolySystem.max_pols) // Use this to terminate when a fixed number of polymers max_pols have been written 
	{
		// Simulate one step		
 		PolySystem.grow();
	}

	cout << "\n\n Finished at simulation time T: " << PolySystem.time << endl;

	PolySystem.outputfile.close();

}