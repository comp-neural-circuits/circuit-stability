import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from math import*

def getThreshold(spike_times_):
	''' Input is just a list of spike times '''
	### output is the guessed threshold and 
	### also which method was used ( 1,2,3 )

	ISIs = np.diff(spike_times_)
	sorted_ISI = np.sort(ISIs[:])

	dum_value = [sorted_ISI[0]]


	sp_thr = (np.percentile(ISIs,90) + np.min(ISIs))/2.;
	method = 2
	
	if (abs(sp_thr-sorted_ISI[-1])<10 or abs(sp_thr-sorted_ISI[0])<10) and method == 2:
		method = 3
		sp_thr = 0.99*min(ISIs)


	return (sp_thr, method)


def main(burst_detect_mode_ = 0, path_ = './'):

	at_least_one_cell_only_spikes = False
	The_active_cell = 2
	# 2 means both are active
	filename = path_
	# print(filename)
	if os.path.isfile(filename): 
		readout = [i.strip().split() for i in open(filename).readlines()]
		# readout = [i.strip().split() for i in open(filepath+"cell1_probe_spks.dat".format(cell)).readlines()]
		spike_times_1 = [float(readout[ii][1]) for ii in range(len(readout)) if int(readout[ii][0]) == 0]
		spike_times_2 = [float(readout[ii][1]) for ii in range(len(readout)) if int(readout[ii][0]) == 1]

		cell_1_spikes = True
		cell_2_spikes = True
		if spike_times_1 == []:
			cell_1_spikes = False
			The_active_cell = 1
		if spike_times_2 == []:
			cell_2_spikes = False
			if The_active_cell == 1:
				The_active_cell = -1
			else:
				The_active_cell = 0

		if len(spike_times_1) == 1 or len(spike_times_2) == 1:
			print("############################")
			print("got one here")
			print("############################")
			cell_1_spikes = cell_2_spikes = False
	else:
		cell_1_spikes = False
		cell_2_spikes = False



		#########################################################
		####### 			BURST DETECTION				#########
		#########################################################

		### cell 1
	if burst_detect_mode_ == 0 or burst_detect_mode_ == 2:


		#### First each cell individually ####
		######################################

		### cell 1 ###
		##############

		if cell_1_spikes:




			(sp_thr_1,method_1) =  getThreshold(spike_times_1)

			num_of_spikes_per_burst_1 = [0 for ii in range(30)]
			ii = 0
			num_sp_1 = [1]
			first_spike_1 = [spike_times_1[0]]
			last_spike_1 = []
			n_bursts_1 = 0


			while len(spike_times_1) > ii+1:
				# for each burst
				if (spike_times_1[ii+1] - spike_times_1[ii] < sp_thr_1): # if spikes are close enough in time
					num_sp_1[-1] += 1 	# add that spike to the actual burst
				else:
					last_spike_1.append(spike_times_1[ii])
					n_bursts_1 += 1
					#num_of_spikes_per_burst_1[num_sp_1[-1]] += 1
					num_sp_1.append(1)
					first_spike_1.append(spike_times_1[ii+1])
				ii += 1;

			last_spike_1.append(spike_times_1[ii])

			first_spike_1 = first_spike_1[1:]
			last_spike_1 = last_spike_1[1:]
			num_sp_1 = num_sp_1[1:]


		### cell 2 ###
		##############
		if cell_2_spikes:

			(sp_thr_2,method_2) =  getThreshold(spike_times_2)

			num_of_spikes_per_burst_2 = [0 for ii in range(30)]
			ii = 0
			num_sp_2 = [1]
			first_spike_2 = [spike_times_2[0]]
			last_spike_2 = []
			n_bursts_2 = 0

			while len(spike_times_2) > ii+1:
				# for each burst
				if (spike_times_2[ii+1] - spike_times_2[ii] < sp_thr_2): # if spikes are close enough in time
					num_sp_2[-1] += 1 	# add that spike to the actual burst
				else:
					last_spike_2.append(spike_times_2[ii])
					n_bursts_2 += 1
					#num_of_spikes_per_burst_2[num_sp_2[-1]] += 1
					num_sp_2.append(1)
					first_spike_2.append(spike_times_2[ii+1])
				ii += 1;

			last_spike_2.append(spike_times_2[ii])

			# print(last_spike_2[0:3])
			# print(first_spike_2[0:3])

			first_spike_2 = first_spike_2[1:]
			last_spike_2 = last_spike_2[1:]
			num_sp_2 = num_sp_2[1:]


		#### Now looking at the circuit properties ####
		###############################################

		# Here we want to capture the alternating effects of the circuit
		# such that we want to call everything a burst, that happens between two bursts of the other cell
		# In case there is an overlap of two bursts, this is counted as a new burst as well

		if cell_1_spikes and cell_2_spikes and burst_detect_mode_ == 2:

			# We alredy got the detected bursts from the code above
			bursts_1 = [(aa,bb,cc) for aa,bb,cc in zip(first_spike_1,last_spike_1,num_sp_1)]
			bursts_2 = [(aa,bb,cc) for aa,bb,cc in zip(first_spike_2,last_spike_2,num_sp_2)]



			### First Loop: Run through the detected bursts of cell one and check, if there is a burst of cell 2 between them
			for omega in range(2):
				if omega == 0:
					bursts_a = bursts_1[:]
					bursts_b = bursts_2[:]
				if omega == 1:
					bursts_a = bursts_2[:]
					bursts_b = bursts_1[:]


				first_spike = []
				last_spike = []
				num_sp = []
				last_kk = 0
				found_next_burst = False


				ii = 0
				while ii < len(bursts_a):
					b1_beginning = bursts_a[ii][0]
					b1_ending = bursts_a[ii][1]
					b1_num_sp = bursts_a[ii][2]

					first_spike.append(b1_beginning)
					last_spike.append(b1_ending)
					num_sp.append(b1_num_sp)

					for jj in range(ii + 1,len(bursts_a)):
						next_b1_beginning = bursts_a[jj][0]
						next_b1_ending = bursts_a[jj][1]
						next_b1_num_sp = bursts_a[jj][2]

						for kk in range(last_kk, len(bursts_b)):
							b2_beginning = bursts_b[kk][0]
							b2_ending = bursts_b[kk][1]

							if b1_beginning <= b2_beginning <= next_b1_beginning:
								last_kk = kk
								found_next_burst = True
								break
							if b1_beginning < next_b1_beginning < b2_beginning:
								last_kk = kk
								found_next_burst = False
								break

						if found_next_burst:
							ii = jj - 1
							break
						else:
							if jj == (len(bursts_a) - 1):
								last_spike[-1] = next_b1_ending
								num_sp[-1] += next_b1_num_sp
								ii = (len(bursts_a) - 1)
								break
							else:
								last_spike[-1] = next_b1_ending
								num_sp[-1] += next_b1_num_sp
								continue

					ii += 1

				if omega == 0:
					first_spike_1 = first_spike[:]
					last_spike_1 = last_spike[:]
					num_sp_1 = num_sp[:]

				if omega == 1:
					first_spike_2 = first_spike[:]
					last_spike_2 = last_spike[:]
					num_sp_2 = num_sp[:]



	
	if burst_detect_mode_ == 1:

		if cell_1_spikes and not cell_2_spikes:
			cell_1_spikes = False
		if cell_2_spikes and not cell_1_spikes:
			cell_2_spikes = False

		if cell_2_spikes and cell_1_spikes:
			num_sp_1 = []
			first_spike_1 = []
			last_spike_1 = []
			n_bursts_1 = 0

			num_sp_2 = []
			first_spike_2 = []
			last_spike_2 = []
			n_bursts_2 = 0

			all_spikes = [(int(element[0]),float(element[1])) for element in readout]

			actual_spiking_cell = all_spikes[0][0]

			if actual_spiking_cell == 0:
				last_spike_1.append(all_spikes[0][1])
				first_spike_1.append(all_spikes[0][1])
				num_sp_1.append(0)
			if actual_spiking_cell == 1:
				last_spike_2.append(all_spikes[0][1])
				first_spike_2.append(all_spikes[0][1])
				num_sp_2.append(0)

			for element in all_spikes:
				if element[0] != actual_spiking_cell:
					actual_spiking_cell = element[0]
					if element[0] == 0:
						n_bursts_1 += 1
						last_spike_1.append(element[1])
						first_spike_1.append(element[1])
						num_sp_1.append(0)
					if element[0] == 1:
						n_bursts_2 += 1
						last_spike_2.append(element[1])
						first_spike_2.append(element[1])
						num_sp_2.append(0)

				if element[0] == 0:
					last_spike_1[-1] = element[1]
					num_sp_1[-1] += 1
				if element[0] == 1:

					last_spike_2[-1] = element[1]
					num_sp_2[-1] += 1



			first_spike_1 = first_spike_1[1:]
			last_spike_1 = last_spike_1[1:]
			num_sp_1 = num_sp_1[1:]

			first_spike_2 = first_spike_2[1:]
			last_spike_2 = last_spike_2[1:]
			num_sp_2 = num_sp_2[1:]

			


						

	# ___________________BURST DETECTION____________________#
	#########################################################




	#########################################################
	####### 		SINGLE CELL PROPERTIES 			#########
	#########################################################
	if cell_1_spikes:
		mean_numSPB_1 = np.mean(num_sp_1)
		mean_period_1 = np.mean(np.diff(first_spike_1))
		mean_period_1_err = np.std(np.diff(first_spike_1))
		mean_burst_duration_1 = np.mean([e - s for e,s in zip(last_spike_1,first_spike_1) ])
		mean_DC_1 = mean_burst_duration_1/mean_period_1
	else:
		mean_numSPB_1 = float('NaN')
		mean_period_1=  float('NaN')
		mean_period_1_err = float('NaN')
		mean_DC_1 = float('NaN')
	if cell_2_spikes:
		mean_numSPB_2 = np.mean(num_sp_2)
		mean_period_2 = np.mean(np.diff(first_spike_2))
		mean_period_2_err = np.std(np.diff(first_spike_2))
		mean_burst_duration_2 = np.mean([e - s for e,s in zip(last_spike_2,first_spike_2) ])
		mean_DC_2 = mean_burst_duration_2/mean_period_2
	else:
		mean_numSPB_2 = float('NaN')
		mean_period_2 = float('NaN')
		mean_period_2_err = float('NaN')
		mean_DC_2 = float('NaN')

	# ______________SINGLE CELL PROPERTIES__________________#
	#########################################################






	if cell_1_spikes and cell_2_spikes and (burst_detect_mode_ == 0 or burst_detect_mode_ == 2):






		#########################################################
		####### 		BURST EXCLUSION METRIC 			#########
		#########################################################


		last_ii = 0
		Onet = 0

		# first case:
		# 	s1#####e1 
		# 		       s2######e2

		#second case:
		#				s1####e1
		# 	s2######e2

		# the '-1 ' here is just needed for the calculation of th eduty cycle ratio
		

		for jj in range(len(first_spike_1)-1):
			s1 = first_spike_1[jj]
			e1 = last_spike_1[jj]
			if s1 == e1:
				continue
			for ii in range(last_ii, len(first_spike_2)-1):
				if s1 == e1:
					continue
				s2 = first_spike_2[ii]
				e2 = last_spike_2[ii]
				if s2 > e1:
				# here the first loop needs to catch the second one
				# so no 'last_ii'-update because we still want to use
				# the actual burst 2 for the next comparison  
					break
				else:
					last_ii = ii
					if e2 < s1:
						# now the second loop needs to catch the first, so we just 
						# keep it running
						continue
					else:
						# now e2 >= s1 and e1 >= s2
						# possible cases:
						''' 
						1: 
						 		s1########e1
						 	s2########e2
						2:
								s1#######e1
									s2########e2
						3: 		
								s1########e1
								   s2###e2
						4:
								s1########e1
							s2################e2
						'''
					burst_1 = abs(e1-s1)
					burst_2 = abs(e2-s2)

					### cases 3 and 4
					if (s1 > s2 and e1 < e2) or (s1 <= s2 and e1 > e2):
						Onet += min(burst_1,burst_2)
					else:
					### cases 1 and 2 
						if s1 > s2:
							Onet += burst_1 - (e1-e2)
						else:
							Onet += burst_2 - (e2-e1)

						
					
					

					
					if e2 < e1:
						# here it could be, that there is a second burst of line 2, 
						# which also intersects with burst 1 
						# so we just keep circulating in loop 2
						continue
					if e2 > e1:
						# now here burst 2 ends after burst 1, so 
						# maybe there is another

						break


		overall_1 = sum([last_spike_1[ii]-first_spike_1[ii] for ii in range(len(first_spike_1))])
		# complete time when cell 1 is spiking
		overall_2 =  sum([last_spike_2[ii]-first_spike_2[ii] for ii in range(len(first_spike_2))])
		# complete time when cell 2 is spiking
		T_trial = first_spike_1[-1]-first_spike_1[0]
		# print (overall_1,overall_2, )


		if (overall_1+overall_2>T_trial):
		    # overlap time if random
		    Orand = min(overall_1,overall_2) - 0.5*(T_trial-max(overall_1,overall_2))
		    # minimum possible overlap time
		    Omin = overall_1 + overall_2 - T_trial
		else:
		    Orand = min(overall_1,overall_2)**2/(2*(T_trial-max(overall_1,overall_2)))
		    Omin = 0
		
		Omax = min(overall_1,overall_2) #max possible overlap time

		if Orand == 0:
			# This means, that there are only spikes, no burst and thus no overlap possible
			chi = float('NaN')
			at_least_one_cell_only_spikes = True
		else:
			if (Onet>Orand):
			    chi = (Orand-Onet)/(Omax-Orand)

			else:
			    # print(g01,g10,overall_1,overall_2)
			    chi = (Orand-Onet)/(Orand-Omin)
		
		
		# return chi

		# ________________BURST EXCLUSION METRIC________________#
		#########################################################




		#########################################################
		####### 		DUTY CYCLE RATIO 	 			#########
		#########################################################
		if mean_burst_duration_2 == 0:
			duty_cycle_ratio = float('NaN')
		else:
			duty_cycle_ratio = mean_DC_1/mean_DC_2

		# ________________DUTY CYCLE RATIO _____________________#
		#########################################################




		#########################################################
		####### 		MEAN PHASE DIFFERENCE			#########
		#########################################################



		phase_difference = []
		period = []
		last_ii = 0
		jj = 0


		for jj in range(len(first_spike_1)-1):

			s1 = first_spike_1[jj]
			period.append(first_spike_1[jj + 1] - s1)

			for ii in range(last_ii,len(first_spike_2)):

				s2 = first_spike_2[ii]
				if (s2 < s1 ):
					continue
				else:
					last_ii = ii
					phase_difference.append((s2 - s1)/(first_spike_1[jj+1] - s1))
					break

						
		phase_difference_mean = np.mean(phase_difference)
		phase_difference_std = np.std(phase_difference)
	
		

		# ________________MEAN PHASE DIFFERENCE ________________#
		#########################################################


		

		#########################################################
		####### 		CORRELATION COEFFICIENT			#########
		#########################################################

		a = np.histogram([spike_times_1[ii] - 5000 for ii in range(len(spike_times_1))] ,np.arange(0,1000000,10))
		b = np.histogram([spike_times_2[ii] - 5000 for ii in range(len(spike_times_2))] ,np.arange(0,1000000,10))
		corr_coeff = np.corrcoef(a[0],b[0])[0][1]


		# _____________ CORRELATION COEFFICIENT ________________#
		#########################################################



	else:
		chi = float('NaN')
		duty_cycle_ratio = float('NaN')
		phase_difference_mean = float('NaN')
		phase_difference_std = float('NaN')
		corr_coeff = float('NaN')





	return [mean_numSPB_1,mean_numSPB_2,mean_period_1,mean_period_2,mean_DC_1, mean_DC_2,chi,duty_cycle_ratio,phase_difference_mean,phase_difference_std,corr_coeff,
			mean_period_1_err, mean_period_2_err, at_least_one_cell_only_spikes, The_active_cell]




	"""
	phase_difference = [ [] for ii in range(7)]
	#########   detecting the phase difference to the next burst of spike 2 ###


	for kk in range(7):
		skip = kk

		for hh in range(skip+1):
			jj = hh

			last_ii = 0
			while jj < len(last_spike_1):
				e1 = last_spike_1[jj]
				for ii in range(last_ii, len(first_spike_2)-1):
					e2 = first_spike_2[ii]
					if e2 < e1:
						continue
					else:
						last_ii = ii 
						phase_difference[skip].append(e2-e1)
						# print(e1,first_spike_2[ii], first_spike_2[ii-1],first_spike_2[ii]-e1)
						break
				jj += (skip + 1)
					

	phase_difference_std = [ np.std(phase_difference[ii]) for ii in range(len(phase_difference))]

	(mu, sigma) = norm.fit(phase_difference[0])
	min_sigma = [sigma,0]
	for ii in range(7):
		(mu, sigma) = norm.fit(phase_difference[ii])
		if sigma < min_sigma[0]:
			min_sigma[0] = sigma
			min_sigma[1] = ii


	return min_sigma


	"""



	
	"""
	x,bins,p = plt.hist(phase_difference, bins=30)
	max_height = 0
	for item in p:
		height = item.get_height()/sum(x)
		item.set_height(height)
		if height > max_height:
			max_height = height




	plt.ylim((0,max_height*1.1))
	plt.tight_layout()
	plt.savefig("./try_{}_{}_histograms/{}_{}.jpg".format(cell_1,cell_2,g10,g01))
	# plt.show()
	plt.close('all')

	sum_height = 0
	for item in p:
		height = item.get_height()
		if height > 0.15:
			sum_height += height
	"""


		




# check, whether this module runs as a program itself or 
# acts just as an imported module
if __name__ == '__main__':
    # in that case - call the main func
    main()




