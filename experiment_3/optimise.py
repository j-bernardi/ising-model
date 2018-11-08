import sys, random, time, copy
import numpy as np 
from matplotlib import pyplot as plt

sys.path.insert(0, '..')
import main, exp3


def get_n_smallest(n, dct) :
	"""Returns key with nth smallest mean value in number dict of lists, mean and stdev"""

	#dct : temp : [Cs]
	#find temp with lowest mean Cs
	#return temp of the nth min C, the nth min C, the stdev of the nth min C

	lst = []
	for key in dct :
		lst.append((key, exp3.get_mean(dct[key]), exp3.get_stdev(dct[key])))

	lst.sort(key=lambda x: x[1])

	return lst[n-1][0], lst[n-1][1], lst[n-1][2]

def simulate_dct(T_dct) :
	"""progress for each temperature in the dict
	append the C result to the dictionary"""

	these_ens = []
	print("\nProgressing for first element in dict- return to start")
	lattice.progress_n(n_eq_N)
	temp_count = 0
	for temp in T_dct :

		lattice.new_run(temp)
		temp_count +=1
		print("\nProgressing for T_eq, temp", temp, "-", temp_count, "/", len(T_dct))
		lattice.progress_n(n_eq_T)

		print("Progressing for data")
		lines = lattice.progress_n(n_prog)

		for l in lines :
			splt = l.split(",")
			l_en, l_mag = float(splt[1].strip()), float(splt[2].strip())
			
			these_ens.append(l_en/N ** 2)
	
		var_e = exp3.get_var(these_ens)
	
		C = var_e/(main.k_b*temp**2)
	
		T_dict[temp].append(C)

	return T_dict

def simulate_key(T_dct, temp) :
	"""Simulates only 1 temperature in the temp dict"""

	these_ens = []

	lattice.new_run(temp)

	print("Progressing for T_eq")
	lattice.progress_n(n_eq_T)

	lines = lattice.progress_n(n_prog)

	for l in lines :
		splt = l.split(",")
		l_en, l_mag = float(splt[1].strip()), float(splt[2].strip())
			
		these_ens.append(l_en/N ** 2)
	
	var_e = exp3.get_var(these_ens)
	
	C = var_e/(main.k_b*temp**2)
	
	T_dict[temp].append(C)

	return T_dict

#########
tests = 3

T_range = [2.2, 3.5]
#Reset the simulation parameters
t_mid = exp3.get_mean(T_range)
	
#(Temp: [values])
T_dict = {T_range[0] : [], t_mid : [], T_range[1] : []}

typ = 'r'

#field
h = 0 

n_eq_N = 600

n_eq_T = 200

n_prog = 200

#number of graphs to display
graphs = 0

N = 30
#################

#simulate n times until smallest + stdev < 2nd smallest - stdev

#Initialise
lattice = main.lattice(N, T_range[0], field = h, exch_energy = main.J_e, mag_moment = main.mu_e, typ=typ)

acc_range = 0.1
first = True

#While the T range is too large, chop the range
while abs(T_range[0] - T_range[1] ) > acc_range :

	print("T_dict so far")
	print(T_dict)

	#Simulate for each temperature twice and get a list - initial run only
	if first : 
		
		print("\n***INITIALISING ***\n")
		print("Simulating for all temps twice- initial range")

		for i in range(2) :
		
			#run initial simulations

			print("\nInitial run, temp", i+1, "/2")
	
			simulate_dct(T_dict)
		first = False

	#else bring the new sliced temperature up to speed
	else :
		if len(T_dict[t_mid]) != len(T_dict[max_T]) :
			print("\nBringing new temp up to speed, T =", t_mid)
			
			for i in range(len(T_dict[max_T])) :
			
				print("\nRun", i+1, "/", len(T_dict[max_T]))
			
				T_dict = simulate_key(T_dict, t_mid)

	print("\n*** EXAMINING RANGE", T_range, "***\n")

	print("Current dict:", T_dict)
	
	#Now compare the bottom 2 values to remove the smallest one.
	min_T, min_C, min_C_std = get_n_smallest(1, T_dict)
	next_T, next_C, next_C_std = get_n_smallest(2, T_dict)
	max_T, max_C, max_C_std = get_n_smallest(3, T_dict)

	print("\nmin, next, max T found:", min_T, next_T, max_T)
	print("From Cs", min_C, next_C, max_C, '\n')

	#If not overlapping, remove the bottom one and reslice range, else simulate whole thing more
	if min_T != t_mid and min_C + min_C_std < next_C - next_C_std and min_C + min_C_std < max_C - max_C_std :

		#Success
		print("Removing bottom pt, ", min_T, "- reducing range")
		#remove the min_
		del T_dict[min_T]

		#Slice the new range
		T_range = [next_T, max_T]
		print("New T_range:", T_range)
		
		t_mid = exp3.get_mean(T_range)
		print("New mid temp:", t_mid)

		#Update the dictionary to hold the next 3 positions to examine
		T_dict[t_mid] = []
		
		#Return to the while loop- check the new range

	else :
		print("Simulating again- still overlapping.")
		print("min+std:", min_C + min_C_std, "mid-std", next_C - next_C_std)
		
		#Not yet- run another simulation
		T_dict = simulate_dct(T_dict)
			
		#Update values
		min_T, min_C, min_C_std = get_n_smallest(1, T_dict)
		next_T, next_C, next_C_std = get_n_smallest(2, T_dict)
		max_T, max_C, max_C_std = get_n_smallest(3, T_dict)

print("Final dictionary:")
print(T_dict)
print()
print("Mid range:", t_mid, "+/-", acc_range/2)

			
