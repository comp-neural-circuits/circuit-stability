# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, Data, RealData
from pylab import *
import bExclu_DCratio_PhaseDiff_CorrCoeff as get_measures
from mpl_toolkits.axes_grid1 import make_axes_locatable
import copy
import pickle
import argparse
import platform
from matplotlib.patches import Circle, Arc

CONST_PI = 3.14159265

__author__ = 'See'
 
parser = argparse.ArgumentParser(description='This script calls ( and thus requires! ) the script "bExclu_DCratio_PhaseDiff_CorrCoeff.py".')
parser.add_argument('-b','--burst_detect_mode', help='Specify the burst detection mode. \n 0 means that each burst is count seperately \n 1 means, that it counts as one burst, as long as one cell bursts \n Examples: \n ||| ... ||| ||| ||| ... ... 0: 7 bursts    1: 4 bursts ',required=False, default=2)
parser.add_argument('-v','--argument_for_calc', help='just give me a value',required=False, default=0)
parser.add_argument('-c','--change_two_conductances', help='specify the two conductances you want to be changed',required=True)
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




def calc_everything(cells_, change_1_, change_to_1_,change_2_, change_to_2_, g10_, g01_, set_number):
	# returns a dictionary
	
	annoying_change_to_list = [(0,0),(0,0.001),(0.001,0),(0,-0.001),(-0.001,0),(0.001,0.001),(-0.001,-0.001),(0.001,-0.001),(-0.001,0.001),(0,0)]
	for ee in annoying_change_to_list:
		# print np.nan_to_num(dummy)
		new_change_to_2 = change_to_2_ + ee[0]
		new_change_to_1 = change_to_1_ + ee[1]
		path = './' + "{}_{}".format(change_1_,change_2_) +"/set_{}".format(set_number) + "/{0}_{1}_sp_{2}_{3:g}_{4}_{5:g}_g10_{6}_g01_{7}.dat".format(cells_[0],cells_[1],change_1_,new_change_to_1,change_2_,new_change_to_2,g10_,g01_)
		found_one = False
		if os.path.isfile(path):
			dummy = get_measures.main(burst_detect_mode_ = 2,  path_ = path)
			found_one = True
			break

	if not found_one:
		dummy = get_measures.main(burst_detect_mode_ = 2,  path_ = path)

	return dummy





def main():

	filename = "conductances_dataset_2000.dat"
	readout = [i.strip().split() for i in open(filename).readlines()]
	g_values = [[float(element) for ii,element in enumerate(bigger_element) if ii != 0] for bigger_element in readout]
	
	# Here we decide on our refrence value
	# change_it_1 will be<

	change_1 = args.change_two_conductances[0]
	change_2 = args.change_two_conductances[1]
	change_it = (int(change_1),int(change_2))


	cells_per_batch = 25
	batches = 1
	number_of_angular_perturbations = 15

	poss_perturbations = ["Na","CaT","CaS","A","KCa","Kd","H"]

	process_number = int(args.argument_for_calc)

	print 'startet run from ii = ', process_number*batches, ' to ', process_number*batches+batches
	for ii in range(process_number*batches,process_number*batches+batches):
		print "ii: ", ii
		# The reference value comes first
		name_for_these_pairs = "change_{}_and_{}__process_{}_with_{}_pairs_per_process".format(poss_perturbations[change_it[0]],poss_perturbations[change_it[1]],ii,cells_per_batch)

		
		path_2 = "./00_dictionaries/" 
		object_name_ = path_2 + name_for_these_pairs
		burst_detect_mode_ = int(args.burst_detect_mode)
		
		title = ['mean spikes per burst (0)','mean spikes per burst (1)',
				'mean period (0)', 'mean period (1)',
				'mean duty cycle (0)', 'mean duty cycle (1)',
				'burst exclusion metric',
				'mean duty cycle ratio',
				'mean phase difference',  'mean phase difference std',
				'correlation coefficient',
				'mean period error (0)', 'mean period error (1)']

		try:
			# check, whether the specific pair has already been calculated

			aa_dic = load_obj(object_name_)

		except IOError: 


			filename = "selected_circuits_8757.dat"
			readout = [i.strip().split() for i in open(filename).readlines()]
			corresponding_g_values = [[float(element) for ff,element in enumerate(bigger_element) if ff > 1] for gg,bigger_element in enumerate(readout) if (ii * cells_per_batch) <= gg < (ii * cells_per_batch)  + cells_per_batch]
			poss_pairs = [[int(element) for ff,element in enumerate(bigger_element) if ff < 2] for gg,bigger_element in enumerate(readout) if (ii * cells_per_batch) <= gg < (ii * cells_per_batch)  + cells_per_batch]


			aa_dic = {} 
			number_of_perturbations = 31
			number_of_angular_perturbations = 15
			for kk, cells in enumerate(poss_pairs):
				print ii, kk, cells
				cell_key = (cells[0],cells[1])

				
				aa_dic[cell_key] = {} 
				aa_dic[cell_key]["Y"] = {}
				aa_dic[cell_key]["X"] = {}
				aa_dic[cell_key]["props"] = {}
				aa_dic[cell_key]["reference_X"] = {}
				aa_dic[cell_key]["reference_conductance"] = change_1

				# it it is not yet calculated, do this now
				g01 = corresponding_g_values[kk][0]
				g10 = corresponding_g_values[kk][1]
				

				result = calc_everything(cells, poss_perturbations[change_it[0]], round(g_values[cells[0]][change_it[0]],3) ,poss_perturbations[change_it[1]] , round(g_values[cells[0]][change_it[1]],3) , g10, g01, ii)
				for bb,key_ in enumerate(title):
					aa_dic[cell_key][key_] = {} 
					aa_dic[cell_key][key_]["reference"] = result[bb]
			

				for nn in range(2):
					# This loop goes over the two cells beeing the reference cell (different to the reference value, attention)
				
				
					
					for alpha in range(nn, number_of_angular_perturbations - nn):
						# Here the nn neatly deals witht the fact that we double calculations for the two diagonal lines (at the beginngin and the end of the second run, we calculate them 
						# already with the first cell as a reference)
						# Now we want to split up the calculation regarding the new alpha value - Phi
						
						# We want to skip some angles, because they were chosen to close to each other:

						exclude_those =[
						[1,2,12,13,15],
						[1,2,12,13,15]
						]

						if alpha in exclude_those[nn]:
								continue
								# print " exclude "


						
						
						Angle = -CONST_PI/2. + CONST_PI * alpha/float(number_of_angular_perturbations-1)
						# This value goes between -PI/2 and PI/2


						half = (number_of_perturbations-1)/2
						X_first_half = [ceil(change_to * 1000 * 2 * g_values[cells[0]][change_it[nn]]/float(number_of_perturbations-1))/1000. for change_to in range(0,half)]
						X_second_half = [ceil(change_to * 1000 * 2 * g_values[cells[0]][change_it[nn]]/float(number_of_perturbations-1))/1000. for change_to in range(half+1, number_of_perturbations)]
						X = [X_first_half, X_second_half]
						short = number_of_perturbations
						nn_2 = (nn+1)%2
						Y_first_half = [ ceil( (    (cc - (short -1)/2)/float((short-1)/2) * sin(Angle)) * 1000 * g_values[cells[0]][change_it[nn_2]] + 1000 * g_values[cells[0]][change_it[nn_2]])/1000. for cc in range(0,half) ]
						Y_second_half = [ ceil( (    (cc - (short -1)/2)/float((short-1)/2) * sin(Angle)) * 1000 * g_values[cells[0]][change_it[nn_2]] + 1000 * g_values[cells[0]][change_it[nn_2]])/1000. for cc in range(half+1,short) ]

						Y = [Y_first_half,Y_second_half]

						new_keys = [[315,135,314.3,134.3,312,132,308,128,302,122,293.5,113.5,282.5,102.5,270,90,257.5,77.5,246.5,66.5,238,58,232,52,228,48,225.7,45.7,225,45],[135.7,315.7,138,318,142,322,148.1,328.1,156.5,336.5,167.5,347.5,180,0,193,12.5,203.5,23.5,212,32,218,38,222,42,224.3,44.3]]
				
						for mm in range(2):
							key_3 = float('NaN')
							if nn == 0 and mm == 0 and Angle >= 0:
								key_3 = round(270 + 180 * arctan( ((Y[mm][0-mm]/g_values[cells[0]][change_it[nn_2]])-1) )/CONST_PI, 2)
							if nn == 0 and mm == 0 and Angle < 0:
								key_3 = round(270 - 180 * arctan( -((Y[mm][0-mm]/g_values[cells[0]][change_it[nn_2]])-1) )/CONST_PI, 2)
							if nn == 0 and mm == 1 and Angle >= 0:
								key_3 = round(90 + 180 * arctan( -((Y[mm][0-mm]/g_values[cells[0]][change_it[nn_2]])-1) )/CONST_PI, 2)
							if nn == 0 and mm == 1 and Angle < 0:
								key_3 = round(90 - 180 * arctan(((Y[mm][0-mm]/g_values[cells[0]][change_it[nn_2]])-1) )/CONST_PI, 2)
							if nn == 1 and mm == 0 and Angle >= 0:
								key_3 = round(180 - 180 * arctan( ((Y[mm][0-mm]/g_values[cells[0]][change_it[nn_2]])-1) )/CONST_PI, 2)
							if nn == 1 and mm == 0 and Angle < 0:
								key_3 = round(180 + 180 * arctan( -((Y[mm][0-mm]/g_values[cells[0]][change_it[nn_2]])-1) )/CONST_PI, 2)
							if nn == 1 and mm == 1 and Angle >= 0:
								key_3 = round(0 + 180 * arctan( ((Y[mm][0-mm]/g_values[cells[0]][change_it[nn_2]])-1) )/CONST_PI, 2)
							if nn == 1 and mm == 1 and Angle < 0:
								key_3 = round(360 - 180 * arctan( -((Y[mm][0-mm]/g_values[cells[0]][change_it[nn_2]])-1) )/CONST_PI, 2)
								
							# print key_3, new_keys[nn][(alpha-1*nn)*2+mm] 
							key_3 = new_keys[nn][(alpha-1*nn)*2+mm]

							# print cells[0], cells[1], change_it[nn], change_it[nn_2], key_3, new_keys[nn][(alpha-1*nn)*2+mm]
							
							# print mm, Y[mm][0-mm], (X[mm][0-mm]/g_values[cells[0]][change_it[nn]])-1, change_it[nn],(Y[mm][0-mm]/g_values[cells[0]][change_it[nn_2]])-1,change_it[nn_2], key_3, new_keys[nn][(alpha-1*nn)*2+mm]


							aa_dic[cell_key]["X"][key_3] = X[mm]
							aa_dic[cell_key]["Y"][key_3] = Y[mm]

							aa_dic[cell_key]["reference_X"][key_3] = g_values[cells[0]][change_it[nn]]

							for xx, yy in zip(X[mm],Y[mm]):	
								# print poss_perturbations[change_it[nn]], xx ,poss_perturbations[change_it[nn_2]] , yy	

								result = calc_everything(cells, poss_perturbations[change_it[nn]], xx ,poss_perturbations[change_it[nn_2]] , yy , g10, g01, ii)

								for ll,tt in enumerate(title):
									try:
										aa_dic[cell_key][tt][key_3].append(result[ll])
									except KeyError:
										aa_dic[cell_key][tt][key_3] = [result[ll]]


					# if calculated here, we save the result in a dictonary so we can 
					# reopen it later without further need of calculation 
					
					 	
							
		save_obj(aa_dic, object_name_)

			# each dictonary saved here has the following properties. 
			# the kex "X" gives the x-values of the changed parameter (Na, H, ... depending on thename of the dictionary)
			# the key "props" gives for each cell pair (saved as keys) the corresponding g_values and the original 



		
	return 0 
	  


# check, whether this module runs as a program itself or 
# acts just as an imported module
if __name__ == '__main__':
    # in that case - call the main func
    main()














