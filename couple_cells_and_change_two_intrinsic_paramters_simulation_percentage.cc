#include <iostream>
#include <fstream>
#include <iomanip> 
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <string>
#include <thread>
#include <chrono>
#include <future>
#include <functional>
#include "couple_cells_simulation_infinite_intrinsic_changes.h"

using namespace std;
#define PI 3.14159265


int test(){
	int lol = 0;
	return lol;
}

int main(int argc, char** argv){


int first_cell;
int second_cell;
float g_from_cell_1_to_cell_2;
float g_from_cell_2_to_cell_1;
float V_syn_reverse_cell_1 = -78;
float V_syn_reverse_cell_2 = -78;
float synapse_V_half = 45;
float synapse_V_slope = -2;
float synapse_tau = 100;
int change_max_cond_of_cell_0[1];
int change_max_cond_of_cell_1[1];
char Filename_V [999];
float change_max_cond_to_0[1];
float change_max_cond_to_1[1];
float g_value_range[2];
int short_simulation_time  = 505000;   // in ms 
double g_vec[8][2001];
char stim_file [999];
FILE *stimfile;
sprintf (stim_file, "conductances_new_set_2000.dat");
stimfile = fopen(stim_file,"r");
read_stimulus(stimfile,g_vec);
fclose(stimfile);



float min_period_threshold = 100;
float max_period_threshold = 800;
float min_phase_difference_threshold = 0.47;
float max_phase_difference_threshold = 0.53;
float min_burst_exclusion_threshold = 0.95;

/* parameters that are relevant for running this script multiple times */
// The process number we give this specific script
int process_number = atoi(argv[1]);
// The number of pairs the script should scan through
const int number_of_pairs = 200; // 500;
// The pair from the working file this script should start with  
unsigned int start_with_pair_number = process_number * number_of_pairs;

// The name of the file we load the shuffled pairs from
string load_file = "./real_chosen_sub_10.dat";

int cell_vec[2][2000];
float g_value_vec[2][2000];

// First we randomly pick two cells from our 2000 cells (unique combination)
	
	// Here instead of regenerating the random pairs for each run, we now pick 
	// from a file with pre-generated combinations

	
	ifstream input_file(load_file);
	if (input_file){

		string line;
		string temp;
		int counter;
    
		for (int i = 0; getline(input_file, line); ++i){
						
						stringstream os;
						os << line;
						counter = 0;
						while (getline(os, temp,' ')){

							if (counter == 0){
						    	cell_vec[0][i] = stoi(temp);
						    	// cout << " " << first_cell_arr[counter_for_extracted_pairs] << endl; 
			    			}
			    			if (counter == 1){
						    	cell_vec[1][i] = stoi(temp);
			    			}
			    			if (counter == 2){
					    		g_value_vec[0][i] = stof(temp);
					    	}
					    	if (counter == 3){
					    		g_value_vec[1][i] = stof(temp);

					    	}
					    	counter ++;
					    }

		}







    
    input_file.close();
	}
	else{
		cout << "could not find file " << endl;
	}





float V_start_1 = -38;
float V_start_2 = -42;


int number_of_Threads = 31;
char Filename_sp [number_of_Threads][999];

float ranges[7][2] = {{800,1200},{0,6},{0,12},{20,130},{20,140},{90,120},{0.0,2}}; // defines the ranges in which we originally set the conductances
int runnning_processes = 0;

thread Threads[number_of_Threads];
std::packaged_task<int()> The_Tasks[number_of_Threads];
std::future<int> outputs[number_of_Threads];


float change_it_to_1;
float change_it_to_2;

int first_cell_arr[number_of_Threads];
int second_cell_arr[number_of_Threads];
float change_it_to_arr[number_of_Threads][2];
int what_to_change_arr[number_of_Threads][2];
float g_from_cell_2_to_cell_1_arr[number_of_Threads];
float g_from_cell_1_to_cell_2_arr[number_of_Threads];
int got_output;
int changing_parameter_maximum = 31; 
int angular_maximum = 15;
float calculated_changing_paramter;

int what_to_change_1 = atoi(argv[2]);
int what_to_change_2 = atoi(argv[4]);



// First loop - all circuits that are investigated

for (int cell_pair = process_number * number_of_pairs; cell_pair < process_number * number_of_pairs + number_of_pairs; ++cell_pair)
{	
	first_cell = cell_vec[0][cell_pair];
	second_cell = cell_vec[1][cell_pair];
	cout << " " << first_cell << " " << second_cell << " "  << endl;

	g_from_cell_1_to_cell_2 = g_value_vec[0][cell_pair];
	g_from_cell_2_to_cell_1 = g_value_vec[1][cell_pair];
	
		// second loop - degrees from 0 to 180
		for (int angular_parameter = 0; angular_parameter < angular_maximum ; ++angular_parameter)
		{	

			if ( what_to_change_1 < what_to_change_2 && (angular_parameter == 0 || angular_parameter == 14) )
			{
				continue;
			}
			if (!(angular_parameter == 0 || angular_parameter == 7 || angular_parameter == 14))
			{
				continue;
			}	
			cout << cell_pair << ": " << first_cell << " " << second_cell << " "  << angular_parameter << endl; 

			// third loop - the change in the parameter
			for (int changing_parameter = 0; changing_parameter < changing_parameter_maximum ; ++changing_parameter)
			{	

				change_it_to_1 = ceil(changing_parameter * 1000 * 2 * g_vec[1+what_to_change_1][first_cell] /float(changing_parameter_maximum-1))/1000.;


				calculated_changing_paramter = (changing_parameter- (changing_parameter_maximum-1)/2 )/float((changing_parameter_maximum-1)/2) * sin(-PI/2. + PI * angular_parameter/float(angular_maximum-1));
				

				change_it_to_2 = ceil(calculated_changing_paramter * 1000 * g_vec[1+what_to_change_2][first_cell] + 1000*g_vec[1+what_to_change_2][first_cell])/1000.;


				
				

				
				int dummy_process_number = cell_pair/25;
				sprintf (Filename_sp[runnning_processes], "./%s_%s/set_%d/%d_%d_sp_%s_%g_%s_%g_g10_%g_g01_%g.dat",argv[3],argv[5],dummy_process_number,first_cell,second_cell,argv[3],change_it_to_1,argv[5],change_it_to_2,g_from_cell_2_to_cell_1,g_from_cell_1_to_cell_2);

				// cout << Filename_sp[runnning_processes];

				first_cell_arr[runnning_processes] = first_cell ;
				second_cell_arr[runnning_processes] = second_cell;
				change_it_to_arr[runnning_processes][0] = change_it_to_1;
				change_it_to_arr[runnning_processes][1] = change_it_to_2 ;
				what_to_change_arr[runnning_processes][0] = what_to_change_1;
				what_to_change_arr[runnning_processes][1] = what_to_change_2;
				g_from_cell_2_to_cell_1_arr[runnning_processes] = g_from_cell_2_to_cell_1;
				g_from_cell_1_to_cell_2_arr[runnning_processes] = g_from_cell_1_to_cell_2;
				// New way to try
				The_Tasks[runnning_processes] = packaged_task<int()> (bind(simulation, 
						  first_cell_arr[runnning_processes] , second_cell_arr[runnning_processes], V_start_1, V_start_2,
				          g_from_cell_1_to_cell_2_arr[runnning_processes],g_from_cell_2_to_cell_1_arr[runnning_processes], 
				          V_syn_reverse_cell_1, V_syn_reverse_cell_2,
				          synapse_V_half, synapse_V_slope, synapse_tau,
				          synapse_V_half, synapse_V_slope, synapse_tau,   
				          Filename_V, Filename_sp[runnning_processes], true, true, false, false,
				          // change_max_cond_of_cell_0, change_max_cond_to, 
				          what_to_change_arr[runnning_processes], change_it_to_arr[runnning_processes], 2,
				          change_max_cond_of_cell_1, change_max_cond_to_1, 0 ,
				          short_simulation_time,
				          min_period_threshold, max_period_threshold,
						  min_phase_difference_threshold, max_phase_difference_threshold,
						  min_burst_exclusion_threshold));
				outputs[runnning_processes] = The_Tasks[runnning_processes].get_future();
				Threads[runnning_processes] = thread(move(The_Tasks[runnning_processes]));
				

			


				runnning_processes ++;
				this_thread::sleep_for(chrono::milliseconds(2));



				if (runnning_processes >= number_of_Threads){
					for (int i = 0; i < number_of_Threads; ++i)
					{
						Threads[i].join();
						got_output = outputs[i].get();
						
					}
					runnning_processes = 0;

				}


			// dritter Loop über die unterschiedlichen Parameter die geändert werden können 	
			}

		// Zweiter Loop über alle g_werte die verändert werden können [V_half, tau_syn, g_A, g_H, g_CaS, g_CaT]
		}


// Erster Loop über alle Zellen die untersucht werden sollen	
}




return 0;
}

