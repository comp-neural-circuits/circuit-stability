#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip> 
#include <time.h>
#include <sstream>
#include <thread>
#include <chrono>
#include <future>
#include <string>
//#include <boost/filesystem.hpp>
#include "couple_cells_simulation_infinite_intrinsic_changes.h"
// #include <gsl/gsl_randist.h>
// #include <gsl/gsl_rng.h>
#include <random>
#include <algorithm>
#include <iterator>
#include <functional>
#include <vector>





using namespace std;

int main(int argc, char** argv){

	/* DEFINING PARAMETERS */

		// Defines the number of threads the program should use
		const int number_of_Threads = 31;
		

		/* parameters that are relevant for running this script multiple times */
			// The process number we give this specific script
			int process_number = atoi(argv[1]);
			// The number of pairs the script should scan through
			int number_of_pairs = 500; // 500;
			// The pair from the working file this script should start with  
			unsigned int start_with_pair_number = process_number * number_of_pairs;

		// The name of the file we load the shuffled pairs from
		string stim_file = "./all_comb_of_the_750_cell_pairs.txt";
		

		/* prameters relevant for the simulation */
		

			// The simulation time for the fast scan
			int short_simulation_time  = 45000;   // in ms 
			// the simulaion time for the consollidation scan
			int long_simulation_time  = 405000; // 1205000; // in ms
			                                     
			float V_syn_reverse_cell_0 = -78;
			float V_syn_reverse_cell_1 = -78;
			float synapse_V_half = 45;
			float synapse_V_slope = -2;
			float synapse_tau = 100;
			float V_start_0 = -38;
			float V_start_1 = -42;

		/* parameters for spike train evaluation */
			float min_period_threshold = 100;
			float max_period_threshold = 800;
			float min_phase_difference_threshold = 0.47;
			float max_phase_difference_threshold = 0.53;
			float min_burst_exclusion_threshold = 0.95;


		/* How to scan through the g_plane */

			// This value gives the minimal amount of circuits that should be 
			// found before we stop scanning the g_plane for a specific pair
			int valid_circuit_results_cut_off = 2;
			// This is the actual value of the maximal g_syn
			// because we always add 1 to avoid 0-connections
			const int max_g_syn_integer = 40; // 32
			// step-length of the scan through the g-plane
			const float step_length = 0.025;


			bool make_the_long_run = false;


		/* Generate output files */

			bool supress_V_file_output = true;
			bool supress_Sp_file_output = true;
			bool supress_V_file_output_long_run = true;
			bool supress_Sp_file_output_long_run = false;

	/* DEFINING PARAMETERS */



char Filename_sp [number_of_Threads][999];
char Filename_V [number_of_Threads][999];
char bookkeeping_file[999];


/*For random number generation:*/
  // gsl_rng * r;
  // r = gsl_rng_alloc (gsl_rng_default);
  // gsl_rng_set(r, time(NULL)); //for random seed: time(NULL)
  // /*done*/

// INTRODUCING THE BOOKKEEPING FILES
// 

  	char bookkeeping_file_0[99];
	sprintf (bookkeeping_file_0,"./set_%d/%d_doc_silent.dat",process_number, process_number);
	FILE *BookFile_circuits_silent;
	BookFile_circuits_silent = fopen(bookkeeping_file_0,"w");
	fprintf(BookFile_circuits_silent,"#This file contains all the g_values for the silent circuits \n");
	fclose(BookFile_circuits_silent);

	char bookkeeping_file_1[99];
	sprintf (bookkeeping_file_1,"./set_%d/%d_doc_non_tested_circuits.dat",process_number, process_number);
	FILE *BookFile_circuits_not_tested;
	BookFile_circuits_not_tested = fopen(bookkeeping_file_1,"w");
	fprintf(BookFile_circuits_not_tested,"#This file contains all the g_values for the circuits that were not tested, because the were discarded already because of the g_value border. The first two values of each line indicate the cell pair, all the following pairs give the integer for g_values. \n");
	fclose(BookFile_circuits_not_tested);

	char bookkeeping_file_2[99];
	sprintf (bookkeeping_file_2,"./set_%d/%d_doc_non_satisfactory_circuits.dat",process_number, process_number);
	FILE *BookFile_circuits_evaluated_as_not_fitting;
	BookFile_circuits_evaluated_as_not_fitting = fopen(bookkeeping_file_2,"w");	
	fprintf(BookFile_circuits_evaluated_as_not_fitting,"#This file contains all the g_values for the circuits that were evaluated to NOT satisfy the conditions and thus beeing discarded. \n");
	fclose(BookFile_circuits_evaluated_as_not_fitting);
	char bookkeeping_file_3[99];
	sprintf (bookkeeping_file_3,"./set_%d/%d_doc_satisfactory_circuits.dat",process_number, process_number);
	FILE *BookFile_circuits_evaluated_as_fitting;
	BookFile_circuits_evaluated_as_fitting = fopen(bookkeeping_file_3,"w");	
	fprintf(BookFile_circuits_evaluated_as_fitting,"#This file contains all the g_values for the circuits that were evaluated to satisfy the conditions and thus beeing possibly selected. \n");
	fclose(BookFile_circuits_evaluated_as_fitting);
	char bookkeeping_file_4[99];
	sprintf (bookkeeping_file_4,"./set_%d/%d_doc_chosen_circuits.dat",process_number, process_number);
	FILE *BookFile_chosen_circuits;
	BookFile_chosen_circuits = fopen(bookkeeping_file_4,"w");	
	fprintf(BookFile_chosen_circuits,"#This file contains all the g_values for the selected circuits. \n");
	fclose(BookFile_chosen_circuits);

// INTRODUCING THE BOOKKEEPING FILES


float ranges[7][2] = {{800,1200},{0,6},{0,12},{20,130},{20,140},{90,120},{0.0,2}}; // defines the ranges in which we originally set the conductances

thread Threads[number_of_Threads];
std::packaged_task<int()> The_Tasks[number_of_Threads];
std::future<int> outputs[number_of_Threads];



int change_max_cond_of_cell_0[1] = {-1};
int change_max_cond_of_cell_1[1] = {-1};
float change_max_cond_to_0[1] = {0};
float change_max_cond_to_1[1] = {0};
int what_to_change_arr[number_of_Threads];
int first_cell_arr[number_of_Threads];
int second_cell_arr[number_of_Threads];
float g_from_cell_1_to_cell_0_arr[number_of_Threads];
float g_from_cell_0_to_cell_1_arr[number_of_Threads];
int actual_g_values_of_this_Thread[number_of_Threads][2];

int got_output;

int first_cell;
int second_cell;
int g_01;
int g_10;
float g_from_cell_0_to_cell_1;
float g_from_cell_1_to_cell_0;
vector<float> working_gsyn_01_vec;
vector<float> working_gsyn_10_vec;
/* The following three variables are used to restrict the searched space of possibilities
it works like this:
when we skim through the space and and say we chose (g_01,g_10) = (2,6) 
Now the simulation over a short period gives that cell 0 does not spike at all
This suggests, that the cell will not spike for (2, X > 6) as well, so our new
g_border_for_10[2] = 6 
*/
float g_border_for_01[max_g_syn_integer];
float g_border_for_10[max_g_syn_integer];



// counts for every pair how many different g_value combinations are in the buffer
int number_of_g_pairs_that_leads_to_desired_output;
bool cancel_fast_scan;



int runnning_processes = 0;
int batch_counter = 0;

auto time_at_simulation_start = std::chrono::system_clock::now();

int cell_vec[2][number_of_pairs];

	
	// First we randomly pick two cells from our 2000 cells (unique combination)
	// Here instead of regenerating the random pairs for each run, we now pick 
	// from a file with pre-generated combinations

	
	ifstream input_file(stim_file);
    for(int i = 0; i < start_with_pair_number; ++i)
    {	
    	int dummy;
    	input_file >> dummy >> dummy;
    }
    if (start_with_pair_number + number_of_pairs > 281625 - 1){
    	number_of_pairs = 281625 - 1 - start_with_pair_number;
    }
    for (int i = 0; i < number_of_pairs; ++i)
    {
    	input_file >> cell_vec[0][i] >> cell_vec[1][i];
    }
    
    input_file.close();


	for (int cell_pair_iterator = 0; cell_pair_iterator < number_of_pairs; ++cell_pair_iterator)
	// begin cell_pair_loopp
	{	
	

		
		first_cell = cell_vec[0][cell_pair_iterator];
		second_cell= cell_vec[1][cell_pair_iterator];
		auto intermediate_time = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = intermediate_time - time_at_simulation_start;
	

		number_of_g_pairs_that_leads_to_desired_output = 0;

		BookFile_circuits_not_tested = fopen(bookkeeping_file_1,"a");	
		fprintf(BookFile_circuits_not_tested,"\n%d %d", first_cell, second_cell);
		fclose(BookFile_circuits_not_tested);

		BookFile_circuits_evaluated_as_not_fitting = fopen(bookkeeping_file_2,"a");	
		fprintf(BookFile_circuits_evaluated_as_not_fitting,"\n%d %d", first_cell, second_cell);
		fclose(BookFile_circuits_evaluated_as_not_fitting);

		BookFile_circuits_silent = fopen(bookkeeping_file_0,"a");	
		fprintf(BookFile_circuits_silent,"\n%d %d", first_cell, second_cell);
		fclose(BookFile_circuits_silent);
		// Then we pick randomly from the possible ranges for the synaptic connections


		int n(0);
		random_device rd;
		mt19937 g(rd());
		vector<int> gsyn_vec;
		gsyn_vec.reserve(max_g_syn_integer*max_g_syn_integer);
		generate_n(back_inserter(gsyn_vec), max_g_syn_integer*max_g_syn_integer, [n]()mutable { return n++; });
		shuffle(gsyn_vec.begin(), gsyn_vec.end(), g);
		// Now we do modulo devision to get the g_value_01, and integer division to get the g_value_10


		// Initialize the g_borders before each run through the g_syn_plane
		for (int ii = 0; ii < max_g_syn_integer; ++ii)
		{
			g_border_for_01[ii] = max_g_syn_integer + 1;
			g_border_for_10[ii] = max_g_syn_integer + 1;
		}

		for (int g_syn_iterator = 0; g_syn_iterator < max_g_syn_integer*max_g_syn_integer; ++g_syn_iterator)
		// begin g_syn_loop
		{	


			// We always add 1 so we don't get results with one connection beeing 0
			g_01 = 1 + int(gsyn_vec[g_syn_iterator]%max_g_syn_integer);
			g_10 = 1 + int(gsyn_vec[g_syn_iterator]/max_g_syn_integer);
			
			

			if( g_01 > g_border_for_01[g_10] || g_10 > g_border_for_10[g_01])
			{
				BookFile_circuits_not_tested = fopen(bookkeeping_file_1,"a");
				fprintf(BookFile_circuits_not_tested," %d %d", g_01,g_10 );
				fclose(BookFile_circuits_not_tested);
				if (g_syn_iterator != (max_g_syn_integer*max_g_syn_integer-1)){
					continue;
				}
			}

			
			int rand_int =  rand() % 100;//gsl_rng_uniform_int(r,100);
			g_from_cell_0_to_cell_1 = step_length * (g_01 + (rand_int - 50)/100.);
			rand_int = rand() % 100;
			g_from_cell_1_to_cell_0 = step_length * (g_10 + (rand_int - 50)/100.);
			actual_g_values_of_this_Thread[runnning_processes][0] = g_01;
			actual_g_values_of_this_Thread[runnning_processes][1] = g_10;



			
			sprintf (Filename_V[runnning_processes], "./set_%d/%d_%d_V_g01_%g_g10_%g.dat",process_number,first_cell,second_cell,g_from_cell_0_to_cell_1,g_from_cell_1_to_cell_0);
			sprintf (Filename_sp[runnning_processes],"./set_%d/%d_%d_sp_g01_%g_g10_%g.dat",process_number,first_cell,second_cell,g_from_cell_0_to_cell_1,g_from_cell_1_to_cell_0);
		
			


			first_cell_arr[runnning_processes] = first_cell ;
			second_cell_arr[runnning_processes] = second_cell;
			g_from_cell_1_to_cell_0_arr[runnning_processes] = g_from_cell_1_to_cell_0;
			g_from_cell_0_to_cell_1_arr[runnning_processes] = g_from_cell_0_to_cell_1;

			// New way to try
			The_Tasks[runnning_processes] = packaged_task<int()> (bind(simulation, 
					  first_cell_arr[runnning_processes] , second_cell_arr[runnning_processes], V_start_0, V_start_1,
			          g_from_cell_0_to_cell_1_arr[runnning_processes],g_from_cell_1_to_cell_0_arr[runnning_processes], 
			          V_syn_reverse_cell_0, V_syn_reverse_cell_1,
			          synapse_V_half, synapse_V_slope, synapse_tau,
					  synapse_V_half, synapse_V_slope, synapse_tau,   
			          Filename_V[runnning_processes], Filename_sp[runnning_processes], true , supress_V_file_output, supress_Sp_file_output, true,
			          change_max_cond_of_cell_0, change_max_cond_to_0, 0,
			          change_max_cond_of_cell_1, change_max_cond_to_1, 0,
			          short_simulation_time,
			          min_period_threshold, max_period_threshold,
					  min_phase_difference_threshold, max_phase_difference_threshold,
					  min_burst_exclusion_threshold));
			outputs[runnning_processes] = The_Tasks[runnning_processes].get_future();
			Threads[runnning_processes] = thread(move(The_Tasks[runnning_processes]));


			runnning_processes ++;
			this_thread::sleep_for(chrono::milliseconds(2));
			

			if (( g_syn_iterator == (max_g_syn_integer*max_g_syn_integer-1)))
			{	
				//This loop ensures, that we stay do not run simulations for two 
				//cell pairs at the same time
				while( runnning_processes < number_of_Threads){
					The_Tasks[runnning_processes] = packaged_task<int()> (do_nothing);
					outputs[runnning_processes] = The_Tasks[runnning_processes].get_future();
					Threads[runnning_processes] = thread(move(The_Tasks[runnning_processes]));
					runnning_processes ++;
				}
			}

			if (runnning_processes >= number_of_Threads)
			// Begin start to evaluate the results from one simulation
			{
				for (int i = 0; i < number_of_Threads; ++i)
				{
					Threads[i].join();
					got_output = outputs[i].get();
					
					
					if (got_output == -1){
					// This means we received a return from the "do_nothing" function
						continue;
					}

					// In case that the return is 1 or 2, one of the cells does not spike at all
					// this means we neglect this try and also block all other possibilities with
					// g_values higher than the ones that already produced no output
					// (we specifically only block the connection that was to strong)
					if( got_output == 1)
					{
						if (g_border_for_10[actual_g_values_of_this_Thread[i][0]] > actual_g_values_of_this_Thread[i][1])
						{
							g_border_for_10[actual_g_values_of_this_Thread[i][0]] = actual_g_values_of_this_Thread[i][1];
						}

						BookFile_circuits_silent = fopen(bookkeeping_file_0,"a");	
						fprintf(BookFile_circuits_silent," %f %f", g_from_cell_0_to_cell_1_arr[i],g_from_cell_1_to_cell_0_arr[i]);
						fclose(BookFile_circuits_silent);
						
					}
					if( got_output == 2)
					{
						if (g_border_for_01[actual_g_values_of_this_Thread[i][1]] > actual_g_values_of_this_Thread[i][0])
						{
							g_border_for_01[actual_g_values_of_this_Thread[i][1]] = actual_g_values_of_this_Thread[i][0];
						}

						BookFile_circuits_silent = fopen(bookkeeping_file_0,"a");	
						fprintf(BookFile_circuits_silent," %f %f", g_from_cell_0_to_cell_1_arr[i],g_from_cell_1_to_cell_0_arr[i]);
						fclose(BookFile_circuits_silent);
						
					}
					if(got_output == 0)
					{	
						working_gsyn_01_vec.push_back(g_from_cell_0_to_cell_1_arr[i]);
						working_gsyn_10_vec.push_back(g_from_cell_1_to_cell_0_arr[i]);
						number_of_g_pairs_that_leads_to_desired_output ++;


						
					}
					if (got_output == 3){

							BookFile_circuits_evaluated_as_not_fitting = fopen(bookkeeping_file_2,"a");	
							fprintf(BookFile_circuits_evaluated_as_not_fitting," %f %f", g_from_cell_0_to_cell_1_arr[i],g_from_cell_1_to_cell_0_arr[i]);
							fclose(BookFile_circuits_evaluated_as_not_fitting);
					}
				}

				runnning_processes = 0;

				if (number_of_g_pairs_that_leads_to_desired_output > valid_circuit_results_cut_off)
				{
					cancel_fast_scan = true;
				}
				else{

					if (( g_syn_iterator == (max_g_syn_integer*max_g_syn_integer-1)))
					{ 	
						if (number_of_g_pairs_that_leads_to_desired_output == 0){
							// In case we are through the whole g_plane 
							cancel_fast_scan = false;
						}
						else{
							cancel_fast_scan = true;
						}
					}
					else{
						cancel_fast_scan = false;
					}
				}

				if (cancel_fast_scan)
				{

				batch_counter ++;
				
				BookFile_circuits_evaluated_as_fitting = fopen(bookkeeping_file_3,"a");	
				fprintf(BookFile_circuits_evaluated_as_fitting,"\n%d %d", 0, number_of_g_pairs_that_leads_to_desired_output);
				fprintf(BookFile_circuits_evaluated_as_fitting," %d %d", first_cell, second_cell);
				for (int ii = 0; ii < number_of_g_pairs_that_leads_to_desired_output; ++ii)
				{
					fprintf(BookFile_circuits_evaluated_as_fitting," %f %f", working_gsyn_01_vec[ii],working_gsyn_10_vec[ii]);
				}
				fclose(BookFile_circuits_evaluated_as_fitting);


				// In this case we just don't need to scan through the g_plane anymore
				// but can continue with the next cell_pair instead
				gsyn_vec.clear();

				working_gsyn_01_vec.clear();
				working_gsyn_10_vec.clear();




				break;
				

				}
				






			// End start to evaluate the results from one simulation
			}


		// This is important to guarantee a clear restart
		gsyn_vec.clear();

		// end g_syn_loop			
		}

	
	if ( (batch_counter == number_of_Threads || cell_pair_iterator == (number_of_pairs-1)) && make_the_long_run)
	// Begin - if - start the longer simulation
	{

		int batch_counter_dummy = batch_counter;

		ifstream input_satisfaction(bookkeeping_file_3);
		int used_already;
		int max_value;
		int counter;
		int counter_for_extracted_pairs = 0;
		string line;
		string temp;

		for (int i = 0; getline(input_satisfaction, line); ++i)
		{
			if (i > 1)
			{
				
				if (counter_for_extracted_pairs == batch_counter){
					break;
				}  
				stringstream os;
				os << line;
				counter = 0;
				used_already = 0;
				while (getline(os, temp,' ')){
					if (counter == 0){
	    				used_already = stoi(temp);
					}
					if (counter == 1){
	    				max_value = stoi(temp);
	    				if (max_value < 0){
	    				// This is either the case when max_value is set to -1 for "its done and approved"
	    				// or when its set to -2 for "we are through all the options and did not find any"
	    					break;
	    				}
					}
					if (counter == 2){
				    	first_cell_arr[counter_for_extracted_pairs] = stoi(temp);
				    	// cout << " " << first_cell_arr[counter_for_extracted_pairs] << endl; 
	    			}
	    			if (counter == 3){
				    	second_cell_arr[counter_for_extracted_pairs] = stoi(temp);
	    			}
	    			if (counter == 2*used_already + 4){
			    		g_from_cell_0_to_cell_1_arr[counter_for_extracted_pairs] = stof(temp);
			    	}
			    	if (counter == 2*used_already + 5){
			    		g_from_cell_1_to_cell_0_arr[counter_for_extracted_pairs] = stof(temp);
			    		counter_for_extracted_pairs ++;
			    		break;

			    	}
			    	counter += 1;
				}

			}
		}
	    input_satisfaction.close();

	    

	    for (int ii = 0; ii < batch_counter; ++ii)
	    {	

	    	sprintf (Filename_V[ii], "./set_%d/%d_%d_V_g01_%g_g10_%g.dat",process_number,first_cell_arr[ii],second_cell_arr[ii],g_from_cell_0_to_cell_1_arr[ii],g_from_cell_1_to_cell_0_arr[ii]);
			sprintf (Filename_sp[ii],"./set_%d/%d_%d_sp_g01_%g_g10_%g.dat",process_number,first_cell_arr[ii],second_cell_arr[ii],g_from_cell_0_to_cell_1_arr[ii],g_from_cell_1_to_cell_0_arr[ii]);
			cout << Filename_sp[ii] << endl;
			The_Tasks[ii] = packaged_task<int()> (bind(simulation, first_cell_arr[ii] , second_cell_arr[ii], V_start_0, V_start_1,
			          g_from_cell_0_to_cell_1_arr[ii],g_from_cell_1_to_cell_0_arr[ii], 
			          V_syn_reverse_cell_0, V_syn_reverse_cell_1,
			          synapse_V_half, synapse_V_slope, synapse_tau,
					  synapse_V_half, synapse_V_slope, synapse_tau,   
			          Filename_V[ii], Filename_sp[ii], true , supress_V_file_output_long_run, supress_Sp_file_output_long_run, true,
			          change_max_cond_of_cell_0, change_max_cond_to_0, 0,
			          change_max_cond_of_cell_1, change_max_cond_to_1, 0,
			          long_simulation_time,
			          min_period_threshold, max_period_threshold,
					  min_phase_difference_threshold, max_phase_difference_threshold,
					  min_burst_exclusion_threshold ));

			outputs[ii] = The_Tasks[ii].get_future();
			Threads[ii] = thread(move(The_Tasks[ii]));


			this_thread::sleep_for(chrono::milliseconds(2));
	    }


	    cout << " wait to finish " <<  endl;
	    bool bookkeeping_bool[number_of_Threads];
	    for (int i = 0; i < batch_counter; ++i)
		{	
			Threads[i].join();
			
			got_output = outputs[i].get(); 
	
			if(got_output == 0)
			{	
				BookFile_circuits_evaluated_as_fitting = fopen(bookkeeping_file_4,"a");	
				fprintf(BookFile_circuits_evaluated_as_fitting,"\n%d %d %f %f", first_cell_arr[i],second_cell_arr[i],g_from_cell_0_to_cell_1_arr[i],g_from_cell_1_to_cell_0_arr[i]);
				fclose(BookFile_circuits_evaluated_as_fitting);
				cout << i << endl;
				bookkeeping_bool[i] = true;
			}
			else
			{	
				cout << " ok " << endl;
				bookkeeping_bool[i] = false;

	
			}

			cout << "beeeing at the end" << endl;
			
		}

		ifstream input_satisfaction_2(bookkeeping_file_3);
    
		counter_for_extracted_pairs = 0;

		ostringstream whole_file;

		for (int i = 0; getline(input_satisfaction_2, line); ++i)
		{
			
			if (i > 1){
			
			stringstream os;
			os << line;
			counter = 0;
			int indicator = 0;
			while (getline(os, temp,' ')){
				int dummy = stoi(temp);
				if (counter == 0){
					indicator = dummy;
				}
				else if (counter == 1){
					if ( dummy >= 0)
					{
						if (bookkeeping_bool[counter_for_extracted_pairs]){
							dummy = -1;

							// this actually means we clean the record out of the memory
							batch_counter_dummy --;
						} 
						else {
							indicator ++;
							if (dummy == indicator){
							// Here I want to make it visible that we went to the end
							// but its no treally necessary because we check only for
							// beeing bigger than the indicator up  when reading out
								dummy = -2;

							// this actually means we clean the record out of the memory
							// here we have to do it as well because we ran out of options
							batch_counter_dummy --;
							}
						}
						counter_for_extracted_pairs ++;
					}
					whole_file << "\n" << indicator << " " << dummy; 
				}
				else{
					whole_file << " " << temp;
				}
				
				    	
		    	counter += 1;
				
			
				}
			}
			if (i == 0)
				{
				whole_file << line << "\n";
				}
		}


	    input_satisfaction_2.close();

	    ofstream Satisfactory_file;
	    Satisfactory_file.open(bookkeeping_file_3);

	    Satisfactory_file << whole_file.str();
		Satisfactory_file.close();
		

		batch_counter = batch_counter_dummy;



	// End - if - start the longer simulation
	}


	// end cell_pair_loopp		
}


return 0;
}