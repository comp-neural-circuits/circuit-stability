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
#include <random>

#include<vector>
using namespace std;

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
char Filename_V [999];
float change_to = 0;
float g_value_range[2];
int short_simulation_time  = 505000;   // in ms 
double g_vec[8][2001];
char stim_file [999];
FILE *stimfile;
sprintf (stim_file, "./conductances_dataset_2000.dat");
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
const int number_of_pairs = 50; // 500;
// The pair from the working file this script should start with  
unsigned int start_with_pair_number = process_number * number_of_pairs;

// The name of the file we load the shuffled pairs from
string load_file = "./selected_circuits_8757.dat";

int cell_vec[2][9000];
float g_value_vec[2][9000];

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

const int number_of_Threads = 31;
char Filename_sp [number_of_Threads][999];

float ranges[7][2] = {{800,1200},{0,6},{0,12},{20,130},{20,140},{90,120},{0.0,2}}; // defines the ranges in which we originally set the conductances
int runnning_processes = 0;

thread Threads[number_of_Threads];
std::packaged_task<int()> The_Tasks[number_of_Threads];
std::future<int> outputs[number_of_Threads];


float change_it_to;

int first_cell_arr[number_of_Threads];
int second_cell_arr[number_of_Threads];
int change_max_cond_of_cell_0[1] = {-1};
int change_max_cond_of_cell_1[1] = {-1};
float change_max_cond_to_0[1] = {0};
float change_max_cond_to_1[1] = {0};
float change_it_to_arr[number_of_Threads][1];
int what_to_change_arr[number_of_Threads][1];
float g_from_cell_2_to_cell_1_arr[number_of_Threads];
float g_from_cell_1_to_cell_2_arr[number_of_Threads];
int got_output;
int changing_parameter_maximum = 31; 



// First Loop - all circuits that are investigated in this run

for (int cell_pair = start_with_pair_number; cell_pair < start_with_pair_number + number_of_pairs; ++cell_pair)
{	

	first_cell = cell_vec[0][cell_pair];
	second_cell = cell_vec[1][cell_pair];
	cout << cell_pair <<  " " << first_cell << " " << second_cell << " "  << endl;

	g_from_cell_1_to_cell_2 = g_value_vec[0][cell_pair];
	g_from_cell_2_to_cell_1 =  g_value_vec[1][cell_pair];


	//swap the cells - to ensure that every cell can be the perturbed cell 
	//int rand_int =  rand() % 2 + 1;
	//if (rand_int == 2){
	//	first_cell = cell_vec[1][cell_pair];
	//	second_cell = cell_vec[0][cell_pair];
	//	g_from_cell_1_to_cell_2 = g_value_vec[1][cell_pair];
	//	g_from_cell_2_to_cell_1 =  g_value_vec[0][cell_pair];
	//}
	
	

	// second loop pver all possible perturbations
	for (int what_to_change = 0; what_to_change < 7; ++what_to_change)
	{	
		if (!(what_to_change == 1 or what_to_change == 6)){
			cout << "continued " << what_to_change << endl;
			continue;
		}
		// Third loop over the paramter=range
		for (int changing_parameter = 0; changing_parameter < changing_parameter_maximum ; ++changing_parameter)
		{	

			change_it_to = ceil(changing_parameter * 1000 * 2 * g_vec[1+what_to_change][first_cell] /float(changing_parameter_maximum-1))/1000.;

			



				sprintf (Filename_V, "not_relevant.dat");
				switch (what_to_change)
				{
					case 0:
					sprintf (Filename_sp[runnning_processes], "./mes_%d/Na/%d_%d_sp_Na_%g_g10_%g_g01_%g.dat", process_number + 40 ,first_cell,second_cell,change_it_to,g_from_cell_2_to_cell_1,g_from_cell_1_to_cell_2);
					break;
					case 1:
					sprintf (Filename_sp[runnning_processes], "./mes_%d/CaT/%d_%d_sp_CaT_%g_g10_%g_g01_%g.dat", process_number,first_cell,second_cell,change_it_to,g_from_cell_2_to_cell_1,g_from_cell_1_to_cell_2);
					break;
					case 2:
					sprintf (Filename_sp[runnning_processes], "./mes_%d/CaS/%d_%d_sp_CaS_%g_g10_%g_g01_%g.dat", process_number +40,first_cell,second_cell,change_it_to,g_from_cell_2_to_cell_1,g_from_cell_1_to_cell_2);
					break;
					case 3:
					sprintf (Filename_sp[runnning_processes], "./mes_%d/A/%d_%d_sp_A_%g_g10_%g_g01_%g.dat", process_number +40,first_cell,second_cell,change_it_to,g_from_cell_2_to_cell_1,g_from_cell_1_to_cell_2);
					break;
					case 4:
					sprintf (Filename_sp[runnning_processes], "./mes_%d/KCa/%d_%d_sp_KCa_%g_g10_%g_g01_%g.dat", process_number +40 ,first_cell,second_cell,change_it_to,g_from_cell_2_to_cell_1,g_from_cell_1_to_cell_2);
					break;
					case 5:
					sprintf (Filename_sp[runnning_processes], "./mes_%d/Kd/%d_%d_sp_Kd_%g_g10_%g_g01_%g.dat", process_number +40,first_cell,second_cell,change_it_to,g_from_cell_2_to_cell_1,g_from_cell_1_to_cell_2);
					break;
					case 6:
					sprintf (Filename_sp[runnning_processes], "./mes_%d/H/%d_%d_sp_H_%g_g10_%g_g01_%g.dat", process_number,first_cell,second_cell,change_it_to,g_from_cell_2_to_cell_1,g_from_cell_1_to_cell_2);
					break;
				}



				cout << Filename_sp[runnning_processes] << endl;

				first_cell_arr[runnning_processes] = first_cell ;
				second_cell_arr[runnning_processes] = second_cell;
				change_it_to_arr[runnning_processes][0] = change_it_to;
				what_to_change_arr[runnning_processes][0] = what_to_change;
				g_from_cell_2_to_cell_1_arr[runnning_processes] = g_from_cell_2_to_cell_1;
				g_from_cell_1_to_cell_2_arr[runnning_processes] = g_from_cell_1_to_cell_2;

				// New way to try
				The_Tasks[runnning_processes] = packaged_task<int()> (bind(simulation, first_cell_arr[runnning_processes] , second_cell_arr[runnning_processes], V_start_1, V_start_2,
				          g_from_cell_1_to_cell_2_arr[runnning_processes],g_from_cell_2_to_cell_1_arr[runnning_processes], 
				          V_syn_reverse_cell_1, V_syn_reverse_cell_2,
				          synapse_V_half, synapse_V_slope, synapse_tau,
				          synapse_V_half, synapse_V_slope, synapse_tau,   
				          Filename_V, Filename_sp[runnning_processes], false, true, false, false,
				          what_to_change_arr[runnning_processes], change_it_to_arr[runnning_processes], 1,
				          change_max_cond_of_cell_1, change_max_cond_to_1, 0,
				          short_simulation_time,
				          min_period_threshold, max_period_threshold,
						  min_phase_difference_threshold, max_phase_difference_threshold,
						  min_burst_exclusion_threshold));
				outputs[runnning_processes] = The_Tasks[runnning_processes].get_future();
				Threads[runnning_processes] = thread(move(The_Tasks[runnning_processes]));


				runnning_processes ++;
				this_thread::sleep_for(chrono::milliseconds(200));

				if (runnning_processes >= number_of_Threads){
					for (int i = 0; i < number_of_Threads; ++i)
					{
						Threads[i].join();
						got_output = outputs[i].get();
					}
					runnning_processes = 0;

				}

		}
		// Third loop over the paramter=range
	}
	// second loop pver all possible perturbations
	}


// First Loop - all circuits that are investigated in this run

return 0;
}

