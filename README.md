# Half-center-oscillator

The file 'couple_cells_and_change_synaptic_strength' can be used to scan the space of synaptic connections

The file 'couple_cells_and_change_intrinsic_paramters' can be used to run single perturbations

The file 'couple_cells_and_change_two_intrinsic_paramters_simulation_percentage' can be used to run the double perturbations


 you can compile the files with 

 g++ -c couple_cells_and_change_synaptic_strength.cc -pthread -std=c++11 && g++ couple_cells_and_change_synaptic_strength.o -o scan_space -pthread -std=c++11

 g++ -c couple_cells_and_change_intrinsic_paramters.cc -pthread -std=c++11 && g++ couple_cells_and_change_intrinsic_paramters.o -o neuron_simulation -pthread -std=c++11
 
 g++ -c couple_cells_and_change_two_intrinsic_paramters_simulation_percentage.cc -pthread -std=c++11 && g++ couple_cells_and_change_two_intrinsic_paramters_simulation_percentage.o -o neuron_simulation_double -pthread -std=c++11

 
and rund them via 

./scan_space x , wiht x defining the 'process number', i.e., wich fraction of cells to choose ( 1 = first 500, 2 = 500 - 1000, ...)

It requires a folder structure that seperates between the processes. For each process 'x' the folder needs to be named './set_x'


./neuron_simulation x, with the x being the process number (1 = first 50, 2 = 50 - 100, ...)
This also requires a certain folder structure at the same location of the executed file:
./mes_x with x being the process number and subfolders for all the possible changes (Na, CaT, CaS, A, KCa, Kd, H)

./neuron_simulation_double x a b f g   , with x as the process number (1 = first 200, ..) b and g indicate the name of the perturbation ('Na','CaT','CaS','A','KCa','Kd','H') and a and f indicate the index corresponding of these perturbations in the list ('Na' = 0, 'CaT' = 1, ...)

For these, the folder structure follows the following logic. 
./b_g/set_x wit b and g as defined above, the names of the perturbations and x being the process number. 

The stored files in general report the spike times (therefore they start with 'sp'). These can be evaluated using the python file 'calculate_bExclu_DCratio_PhaseDiff_CorrCoeff' for the single perturbations and 'calculate_2D_bExclu_DCratio_PhaseDiff_CorrCoeff' for the double perturbation 

That evaluates all circuit instances with respecto to the following properties:
'mean spikes per burst (0)','mean spikes per burst (1)', mean period (0)', 'mean period (1)', mean duty cycle (0)', 'mean duty cycle (1)', burst exclusion metric', mean duty cycle ratio', mean phase difference',  'mean phase difference std', correlation coefficient', mean period error (0)', 'mean period error (1)'

it needs a folder that is called './00_dictionaries'


The python script 'calculate_stability.py' can load the generated dictionaries to calculate the stability values.



