# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.widgets
import mpl_toolkits.axes_grid1
import pickle
# import cPickle
import argparse
import seaborn as sns
from pandas import DataFrame
import random
import scipy.integrate as integrate



__author__ = 'See'
 
parser = argparse.ArgumentParser(description='This script needs precalculated dictonaries with pages for the different amounts of perturbations with H, Na, or others.')
parser.add_argument('-b','--burst_detect_mode', help='Specify the burst detection mode. \n 0 means that each burst is count seperately \n 1 means, that it counts as one burst, as long as one cell bursts \n Examples: \n ||| ... ||| ||| ||| ... ... 0: 7 bursts    1: 4 bursts ',required=False, default=2)
args = parser.parse_args()

#plt.rc('text', usetex=True) # Erlaubt die Beschriftung in Latex r'[...]' , bzw r"[...]" 
#plt.rc('font', family='serif')

# fig_of_two,ax_of_two = plt.subplots(1,2, figsize = (25.,15.))
# GLOBAL_cell_0 = 0
# GLOBAL_cell_1 = 0



def load_obj(name_ ):
	with open(name_ + '.pkl', 'rb') as f:
		return pickle.load(f)

def save_obj(obj, name_):
	with open(name_ + '.pkl', 'wb') as f:
		pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)



def measure_of_stability(complete_dictonary_, cell_pair_key_, reference_value_):
	''' This function differentiates between the left and the right side of
	the initial value (i.e. did we decrease or increase the maximum conductance)
	'''
	
	def localization_func(distance, max_x_value):
		X = 1
		# X = 0.25+2*np.exp(-((8*distance)**2)/float(max_x_value)**2)
		return X

	dic_keys_ =  ['mean spikes per burst (0)','mean spikes per burst (1)',
                    'mean period (0)', 'mean period (1)',
                    'mean duty cycle (0)', 'mean duty cycle (1)',
                    'burst exclusion metric',
                    'mean duty cycle ratio',
                    'mean phase difference', 'mean phase difference std',
                    'correlation coefficient',
                    'mean period error (0)', 'mean period error (1)',
                    'at least one cell only spikes','the active cell']


	reference_phase_difference = complete_dictonary_["props"][cell_pair_key_][1]
	this_series_of_phase_differences = complete_dictonary_[dic_keys_[8]][cell_pair_key_]
	this_series_of_burst_exclusion = complete_dictonary_[dic_keys_[6]][cell_pair_key_]
	X_values_of_this_series = complete_dictonary_["X"][cell_pair_key_]

	phase_burst = [this_series_of_phase_differences ,this_series_of_burst_exclusion]


	# First for the left side
	difference_norm = []
	difference_list_phase_difference = []
	
	set_to_nan_later = False
	if np.isnan(reference_phase_difference) and reference_value_ > X_values_of_this_series[int((len(X_values_of_this_series)-1)/2)]:
		reference_phase_difference = 0.5
		set_to_nan_later = True

	for jj in range(0, int((len(X_values_of_this_series)-1)/2)):

		ii = int((len(X_values_of_this_series)-1)/2 - (jj + 1))

		phase = phase_burst[0][ii]
		burst = phase_burst[1][ii]

		dummy_distance = abs(reference_value_ - X_values_of_this_series[ii])
		localization_factor = localization_func(dummy_distance, X_values_of_this_series[-1])

		difference_norm.append(localization_factor)



		if 0.95 < burst < 1.01  and not np.isnan(phase)  and 0 <= phase <= 1:

			# Now this is to look at the change in the phase difference
			dummy_difference = reference_phase_difference-phase
			if dummy_difference > 0:
				normalization = reference_phase_difference
			else: 
				normalization = 1 - reference_phase_difference
			
			difference_list_phase_difference.append(localization_factor*abs(dummy_difference/normalization))

		else:

			difference_list_phase_difference.append(localization_factor)

	return_value = 1 - (np.sum(difference_list_phase_difference))/np.sum(difference_norm), difference_norm, difference_list_phase_difference
	if np.isnan(return_value[0]):
		return_value = (0, difference_norm, difference_list_phase_difference)
	
	left_side_return = return_value

	# For the right side


	if np.isnan(reference_phase_difference) and reference_value_ < X_values_of_this_series[int((len(X_values_of_this_series)-1)/2)]:
		reference_phase_difference = 0.5
		
	if set_to_nan_later:
		reference_phase_difference = float('nan')




	difference_norm = []
	difference_list_phase_difference = []
	for ii in range(int((len(X_values_of_this_series)-1)/2) + 1, len(X_values_of_this_series)):
		phase = phase_burst[0][ii]
		burst = phase_burst[1][ii]

		dummy_distance = abs(reference_value_ - X_values_of_this_series[ii])

		localization_factor = localization_func(dummy_distance, X_values_of_this_series[-1])
		difference_norm.append(localization_factor)
		# print ii, phase, burst
		if 0.95 < burst < 1.01 and not np.isnan(phase) and 0 <= phase <= 1:

			# Now this is to look at the change in the phase difference
			dummy_difference = reference_phase_difference-phase
			if dummy_difference > 0:
				normalization = reference_phase_difference
			else: 
				normalization = 1 - reference_phase_difference
			
			difference_list_phase_difference.append(localization_factor*abs(dummy_difference/normalization))

		else:

			difference_list_phase_difference.append(localization_factor)

		# print difference_list_phase_difference


	return_value = 1 - (np.sum(difference_list_phase_difference))/np.sum(difference_norm), difference_norm, difference_list_phase_difference

	if np.isnan(return_value[0]):
		return_value = (0, difference_norm, difference_list_phase_difference)
		print (" got one"
			)
	right_side_return = return_value






	return ( left_side_return[0], left_side_return[1],left_side_return[2], right_side_return[0],right_side_return[1],right_side_return[2])




def main():





	g_values = [i.strip().split() for i in open("./conductances_dataset_2000.dat").readlines()]
	g_values_float = [[float(ii) for ii in jj[1:]] for jj in g_values]

	
	value_ranges = [[800,1200],[0,6],[0,12],[20,130],[20,140],[90,120],[0.0,2]]


	dic_keys = ['mean spikes per burst (0)','mean spikes per burst (1)',
                    'mean period (0)', 'mean period (1)',
                    'mean duty cycle (0)', 'mean duty cycle (1)',
                    'burst exclusion metric',
                    'mean duty cycle ratio',
                    'mean phase difference', 'mean phase difference std',
                    'correlation coefficient',
                    'mean period error (0)', 'mean period error (1)',
                    'at least one cell only spikes','the active cell']

	path = './'

	intrinsic = ["Na","CaT","CaS","A","KCa","Kd","H"] 
	circuitrinsic = ["V_th","tau_syn","g_01","g_10"]
	change_it = intrinsic[1:2] + intrinsic[-1:] 


	change_it = change_it[:]
	
	

	MEASUREs = {}
	MEASUREs["cell pair"] = []
	MEASUREs["stable for all"] = []
	MEASUREs["Fano Factor"] = []
	MEASUREs["g_syn 01"] = []
	MEASUREs["g_syn 10"] = []
	


	for keys in dic_keys:
		MEASUREs["circuit_values_of_starting_point {}".format(keys)] = []


	for element in change_it:
		MEASUREs["decrease " + element] = []
		MEASUREs["increase " + element] = []
		MEASUREs["g_value cell 0 " + element] = []
		MEASUREs["g_value cell 1 " + element] = []
		MEASUREs["the active cell " + element] = []
		MEASUREs["path " + element] = []
		for aa in [-3,-2,-1,1,2,3]:
			MEASUREs["phase_burst {} {}".format(element,aa)] = []

		MEASUREs[element] = []


	intrinsic_left = ["decrease " + element for element in intrinsic]
	intrinsic_right = ["increase " + element for element in intrinsic]




	for process_number in range(4):#160):
		print (process_number)
		for change_ii, change_this in enumerate(change_it):


			try:
				# check, whether this dictonary has already been calculated

				object_name_ = path + "00_calculations/mes_{}_detected_values_{}_bursts_as_circuit".format(process_number,change_this)
				bb_dic = load_obj(object_name_)
				
					
			except IOError: 

				print ("did not find a dictonary for ", object_name_)

			""" 
				####################################################################################
				This can be used to show the distribution of g-values in the chosen set of circuits
				####################################################################################
			"""

			# fig, axes = plt.subplots(figsize = [9,9])
			# dudes = [bb_dic["props"][key][0] for key in bb_dic["props"]]
			# dudes_0 = [element[1] for element in dudes]
			# dudes_1 = [element[0] for element in dudes]
			# axes.scatter(dudes_0,dudes_1)
			# axes.set_xlabel("g_10 / 0.025")
			# axes.set_ylabel("g_01 / 0.025")
			# X = np.linspace(0,1,100)
			# # axes.plot(X,X)
			# choice = 0.125
			# axes.plot(X,X-choice)
			# axes.plot(X,X+choice)
			# plt.show()
			# return 

			""" 
				####################################################################################
				This can be used to show the distribution of g-values in the chosen set of circuits
				####################################################################################
			"""

			MEASUREs['decrease ' + change_this] = MEASUREs['decrease ' + change_this] + [42 for ii in range(len(MEASUREs["cell pair"])-len(MEASUREs['decrease ' + change_this]))]
			MEASUREs['increase ' + change_this] = MEASUREs['increase ' + change_this] + [42 for ii in range(len(MEASUREs["cell pair"])-len(MEASUREs['increase ' + change_this]))]
			MEASUREs["the active cell " + change_this] = MEASUREs["the active cell " + change_this] + [42 for ii in range(len(MEASUREs["cell pair"])-len(MEASUREs["the active cell " + change_this]))]
			MEASUREs["path " + change_this] = MEASUREs["path " + change_this] + [42 for ii in range(len(MEASUREs["cell pair"])-len(MEASUREs["path " + change_this]))]
			for aa in [-3,-2,-1,1,2,3]:
				MEASUREs["phase_burst {} {}".format(change_this,aa)] = MEASUREs["phase_burst {} {}".format(change_this,aa)] + [42 for ii in range(len(MEASUREs["cell pair"])-len(MEASUREs["phase_burst {} {}".format(change_this,aa)]))]

			MEASUREs["g_value cell 0 " + change_this] = MEASUREs["g_value cell 0 " + change_this] + [42 for ii in range(len(MEASUREs["cell pair"])-len(MEASUREs["g_value cell 0 " + change_this]))]
			MEASUREs["g_value cell 1 " + change_this] = MEASUREs["g_value cell 1 " + change_this] + [42 for ii in range(len(MEASUREs["cell pair"])-len(MEASUREs["g_value cell 1 " + change_this]))]
			MEASUREs[change_this] = MEASUREs[change_this] + [42 for ii in range(len(MEASUREs["cell pair"])-len(MEASUREs[change_this]))]
			
			for ii,(key_2) in enumerate(bb_dic["props"]):
				
							
				X = bb_dic["X"][key_2]
							
				g_values_ = bb_dic["props"][key_2][0]
				circuitrinsic_original_values = [45,100,g_values_[0],g_values_[1]]
				original_intrinsic_value_cell_0 = g_values_float[key_2[0]][intrinsic.index(change_this)]
				original_intrinsic_value_cell_1 = g_values_float[key_2[1]][intrinsic.index(change_this)]
							
				original_phase_diff_value = bb_dic["props"][key_2][1]

				# if key_2 != (1642,1845):
				# 	continue
				print 
				stability_measure_left, norm_it_left, list_it_left, stability_measure_right, norm_it_right, list_it_right = measure_of_stability(bb_dic,key_2,original_intrinsic_value_cell_0)
				# print change_this, stability_measure_left, stability_measure_right

				# This is needed for every list that is the same for all the different possible changes
				list_id = -1
				if key_2 in MEASUREs["cell pair"]:

					list_id = MEASUREs["cell pair"].index(key_2)
				else:
					MEASUREs["cell pair"].append(key_2)

				if list_id == -1:

					MEASUREs['decrease ' + change_this].append(0)
					MEASUREs['increase ' + change_this].append(0)
					MEASUREs["the active cell " + change_this].append(0)
					MEASUREs["path " + change_this].append(0)
					MEASUREs[change_this].append(0)

			# 	continue

					MEASUREs["stable for all"].append(True)
					# print bb_dic[dic_keys[9]][key_2][15],bb_dic[dic_keys[8]][key_2][15]
					MEASUREs['Fano Factor'].append(bb_dic[dic_keys[9]][key_2][15]**2/bb_dic[dic_keys[8]][key_2][15])


					for keys in dic_keys:
						MEASUREs["circuit_values_of_starting_point {}".format(keys)].append(bb_dic[keys][key_2][15])

					MEASUREs["g_value cell 0 " + change_this].append(0)
					MEASUREs["g_value cell 1 " + change_this].append(0)

					MEASUREs["g_syn 01"].append(0)
					MEASUREs["g_syn 10"].append(0)

					for aa in [-3,-2,-1,1,2,3]:
						MEASUREs["phase_burst {} {}".format(change_this,aa)].append(0)

				
				MEASUREs['decrease ' + change_this][list_id] = stability_measure_left
				MEASUREs['increase ' + change_this][list_id] = stability_measure_right
				MEASUREs["the active cell " + change_this][list_id] = bb_dic[dic_keys[14]][key_2]
				try:
					MEASUREs["path " + change_this][list_id] = bb_dic['path'][key_2]
	
				except KeyError as e:
					MEASUREs["path " + change_this][list_id] = None
					

				for aa in [-3,-2,-1,1,2,3]:
					MEASUREs["phase_burst {} {}".format(change_this,aa)][list_id] = (bb_dic[dic_keys[8]][key_2][15+aa],bb_dic[dic_keys[6]][key_2][15+aa])


				MEASUREs[change_this][list_id] = (stability_measure_left + stability_measure_right)/2.


				if MEASUREs["stable for all"][list_id] and MEASUREs[change_this][list_id] > 0.85:
					MEASUREs["stable for all"][list_id] = True
				else:
					MEASUREs["stable for all"][list_id] = False

				MEASUREs["g_value cell 0 " + change_this][list_id] = original_intrinsic_value_cell_0
				MEASUREs["g_value cell 1 " + change_this][list_id] = original_intrinsic_value_cell_1
				
				MEASUREs["g_syn 01"][list_id] = g_values_[0]
				MEASUREs["g_syn 10"][list_id] = g_values_[1]


				





	save_obj(MEASUREs,"./00_calculations/MEASUREs")





	return 0 
  


# check, whether this module runs as a program itself or 
# acts just as an imported module
if __name__ == '__main__':
    # in that case - call the main func
    main()














