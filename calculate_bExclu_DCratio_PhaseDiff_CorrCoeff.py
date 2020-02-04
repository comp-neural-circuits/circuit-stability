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
from matplotlib.patches import Circle


__author__ = 'See'
 
parser = argparse.ArgumentParser(description='This script calls ( and thus requires! ) the script "bExclu_DCratio_PhaseDiff_CorrCoeff.py".')
parser.add_argument('-b','--burst_detect_mode', help='Specify the burst detection mode. \n 0 means that each burst is count seperately \n 1 means, that it counts as one burst, as long as one cell bursts \n Examples: \n ||| ... ||| ||| ||| ... ... 0: 7 bursts    1: 4 bursts ',required=False, default=2)
parser.add_argument('-v','--argument_for_calc', help='just give me a value',required=False, default=0)
args = parser.parse_args()

process_number = int(args.argument_for_calc)
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


def main():

    
    filename = "conductances_dataset_2000.dat"
    readout = [i.strip().split() for i in open(filename).readlines()]
    g_values = [[float(element) for ii,element in enumerate(bigger_element) if ii != 0] for bigger_element in readout]


    set_start = 0
    batches = 4



    for ii in range(process_number*batches,process_number*batches+batches):

        print (ii )
        filename = "selected_circuits_8757.dat"
        readout = [i.strip().split() for i in open(filename).readlines()]
        poss_pairs = [[int(element) for ff,element in enumerate(bigger_element) if ff < 2] for gg,bigger_element in enumerate(readout) if (ii * 12) <= gg < (ii * 12)  + 12]
        corresponding_g_values = [[float(element) for ff,element in enumerate(bigger_element) if ff > 1] for gg,bigger_element in enumerate(readout) if (ii * 12) <= gg < (ii * 12)  + 12]



        NUMBER_OF_PAIR = ii#+(40*set_start)
        print (NUMBER_OF_PAIR)
        name_for_these_pairs = "mes_{}".format(NUMBER_OF_PAIR)
        poss_perturbations = ["Na","CaT","CaS",
                              "A","KCa","Kd","H"
                              ]
        ranges = [[800,1200],[0,6],[0,12],[20,130],[20,140],[90,120],[0.0,2]]

        for ll, pert in enumerate(poss_perturbations):

            print (ll)
            path_2 = "./" 
            object_name_ = path_2 + "/00_calculations/{}_detected_values_{}".format(name_for_these_pairs,pert)
            burst_detect_mode_ = int(args.burst_detect_mode)
            if burst_detect_mode_ == 1:
                object_name_ = path_2 + "/00_calculations/{}_detected_values_{}_bursts_as_one".format(name_for_these_pairs,pert)
            if burst_detect_mode_ == 2:
                object_name_ = path_2 + "/00_calculations/{}_detected_values_{}_bursts_as_circuit".format(name_for_these_pairs,pert)
            
            title = ['mean spikes per burst (0)','mean spikes per burst (1)',
                    'mean period (0)', 'mean period (1)',
                    'mean duty cycle (0)', 'mean duty cycle (1)',
                    'burst exclusion metric',
                    'mean duty cycle ratio',
                    'mean phase difference', 'mean phase difference std',
                    'correlation coefficient',
                    'mean period error (0)', 'mean period error (1)',
                    'at least one cell only spikes','the active cell']
            aa = [[0 for kk in range(len(poss_pairs))] for ii in range(len(title))]
 
            try:
                # check, whether the specific pair has already been calculated

                aa_dic = load_obj(object_name_)
                for ii_, key_ in enumerate(title):
                    aa[ii_] = aa_dic[key_]
            except IOError: 

                print (pert)
                aa_dic = {} 
                aa_dic["X"] = {}
                aa_dic["props"] = {}
                aa_dic["path"] = {}
                number_of_perturbations = 31
                for kk, cells in enumerate(poss_pairs):

                    cell_1 = cells[0]
                    cell_2 = cells[1]
                    change_to_1 = [ceil(mm * 1000 * 2 * g_values[cell_1][ll]/float(number_of_perturbations-1))/1000. for mm in range(number_of_perturbations)]
                    change_to = change_to_1

                   


                    aa_dic["X"][(cell_1,cell_2)] = change_to
                    # print dummy_dic[title[8]][corresponding_g_values[kk][0]][corresponding_g_values[kk][1]]

                    for cc in range(number_of_perturbations):
                        # it it is not yet calculated, do this now
                        g01 = corresponding_g_values[kk][0]
                        g10 = corresponding_g_values[kk][1]
                        if g01 == 0.00:
                            g01 = 0
                        if g10 == 0.00: 
                            g10 = 0
                        
                        Found_File = False
                        annoying_change_to_list = [0,0.001,-0.001,0.002,-0.002]
                        for ee in annoying_change_to_list:
                            new_change_to = change_to[cc] + ee
                            path = './' + name_for_these_pairs + '/' + pert + "/{0}_{1}_sp_{2}_{3:g}_g10_{4}_g01_{5}.dat".format(cell_1,cell_2,pert,new_change_to,g10,g01)
                            if os.path.isfile(path):
                                dummy = get_measures.main(burst_detect_mode_ = burst_detect_mode_,  path_ = path)
                                Found_File = True
                                voltage_channel_path = './' + name_for_these_pairs + '_v/' + pert + "/{0}_{1}_channels_and_voltage_traces_{2}_{3:g}_g10_{4}_g01_{5}.dat".format(cell_1,cell_2,pert,new_change_to,g10,g01)
                                if os.path.isfile(voltage_channel_path):
                                    try:
                                         aa_dic['path'][(cell_1,cell_2)].append(voltage_channel_path)
                                    except KeyError as e:
                                         aa_dic['path'][(cell_1,cell_2)] = [voltage_channel_path]

                                   

                                break
                            else:
                                dummy = get_measures.main(burst_detect_mode_ = burst_detect_mode_,  path_ = path)







                        if not Found_File:
                            print ('not found ',path)
                            try:
                                 aa_dic['path'][(cell_1,cell_2)].append(None)
                            except KeyError as e:
                                 aa_dic['path'][(cell_1,cell_2)] = [None]
                        if cc == 15:
                            aa_dic["props"][(cell_1,cell_2)] = (corresponding_g_values[kk],dummy[8],dummy[2])
                 
                        for bb in range(len(title)):
                            aa[bb][kk] = dummy[bb]

                        # if calculated here, we save the result in a dictonary so we can 
                        # reopen it later without further need of calculation 
                        for ii_, key_ in enumerate(title):
                            try:
                                aa_dic[key_]
                            except KeyError:
                                aa_dic[key_] = {}
                            key_2 =(cell_1,cell_2)
                            try:
                                aa_dic[key_][key_2].append(aa[ii_][kk])
                            except KeyError:
                                aa_dic[key_][key_2] = [aa[ii_][kk]]

    						
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














