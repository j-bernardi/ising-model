import sys, random, time, copy
import numpy as np 
from matplotlib import pyplot as plt

sys.path.insert(0, '../..')
sys.path.insert(0, '..')
import main

#Make a script that prints all plot results to a text file
#if same temp done, take one with the lower stdev
#Store N

#Make script that plots from the text file

#Use it to fill in the peak around T_C quantatively

#maybe use polyfit with high N, display it and find the peak? 

def get_mean(lst) :

	return sum(lst)/len(lst)

def get_stdev(lst) :

	mean = get_mean(lst)
	var = 0
	for l in lst :
		var += (mean - l) ** 2

	return np.sqrt(var/len(lst))

def get_var(lst) :
	""""""
	mean = get_mean(lst)
	var = 0
	for l in lst :
		var += (mean - l) ** 2

	return var

def append_result_to_file(filename, results) :
	"""Append these simulation results to the file, to save them"""

	with open(filename, 'a') as f :

		for line in results :
			f.write(line)

	print("Appended results to file", filename, "in cwd.")

def get_results_from_file(filename) :
	"""Get all the lines from the file, ready for plotting"""

	lines = []
	with open(filename, 'r') as f :

		for ln in f :
			lines.append((ln.strip().split(",")))
	
	lines = np.array(lines)

	NS = lines[:,0]
	ns = lines[:,1]
	Ts = lines[:,2]
	Cs = lines[:,3]
	Xs = lines[:,4]

	print("Read all result from", filename, "in cwd.")

	return NS, ns, Ts, Cs, Xs
			
def run(result_file, temp_list, h, Num, n_eq_N, n_eq_T, n_prog, graphs) :

	if temp_list[0] > 2.27 :
		typ = 'r'
	else :
		typ = 'a+'
	
	lattice = main.lattice(Num, temp_list[0], field = h, exch_energy = main.J_e, mag_moment = main.mu_e, typ=typ)

	print("\nProgressing for initial eq")
	lattice.progress_n(n_eq_N)

	temp_track = 0

	Xs, Cs, result_lines = [], [], []

	for temp in temp_list :

		these_ens, these_mags = [], []

		temp_track += 1
		print("\nTemp =", temp, "-", temp_track, "/", len(temp_list))

		#Init lattice
		lattice.new_run(temp, h)

		if temp != temp_list[0] :

			print("Progressing for T_eq")
			lattice.progress_n(n_eq_T)

		print("Progressing for heat capacity data")
		lines = lattice.progress_n(n_prog, graphs=graphs)

		for l in lines :
			splt = l.split(",")
			l_en, l_mag = float(splt[1].strip()), float(splt[2].strip())
			
			these_ens.append(l_en/Num ** 2)
			these_mags.append(l_mag)

		var_e = get_var(these_ens)
		var_m = get_var(these_mags)

		C = var_e/(main.k_b*temp**2)
		X = main.mu * var_m/(main.k_b*temp)

		Cs.append(C) 
		Xs.append(X)

		result_lines.append(str(Num) +"," + str(n_prog)+"," + str(temp) +"," + str(C) + "," +str(X) + "\n")

		append_result_to_file(result_file, [str(Num) +"," + str(n_prog)+"," + str(temp) +"," + str(C) + "," +str(X) + "\n"])
		print("Appended this run to", result_file)

	return result_lines

if __name__ == "__main__" :

	################################################	
	############# runtime parameters ###############

	sim = 0#Number times to simulate for each pt

	result_file = "h_0.5_results.txt"

	#size of lattice NxN
	Ns = [50]#[40, 60, 70, 80]
	
	#Temperatures to sketch out curve
	temp_list = [0.5, 1.0, 1.5, 1.8, 2, 2.2,  2.4, 2.6, 3.0, 3.5, 4.5]
	temp_list +=  [6, 8, 10, 12, 13, 14, 15, 16, 17, 18, 20, 25, 30, 35, 40]
	temp_list = [x for x in np.linspace(2, 5, num=10)]

	#field
	h = 0.5

	n_eq_N = 300

	n_eq_T = 400

	n_prog = 300

	#number of graphs to display
	graphs = 0
	################################################
	for N in Ns :
		
		print("\n*** RUN FOR N =", N, "***\n")
		
		simulate = copy.copy(sim)
		
		#Run the simulations
		while simulate != 0 :
			print()
			print("***NEW RUN : runs remaining:", simulate, "***")
			
			result_lines = run(result_file, temp_list, h, N, n_eq_N, n_eq_T, n_prog, graphs)	
	
			#append_result_to_file(result_file, result_lines)
	
			simulate -= 1
	
		all_Ns, all_ns, all_temps, all_Cs, all_Xs = get_results_from_file(result_file)
	
		print("Collecting results...")
		
		#Collect C and X data to dictionaries for each temp
		C_dict, X_dict = {}, {}
	
		for i in range(len(all_temps)) :

			if int(all_Ns[i]) == N :
	
				if float(all_temps[i]) in C_dict :
					C_dict[float(all_temps[i])].append(float(all_Cs[i]))
					X_dict[float(all_temps[i])].append(float(all_Xs[i]))
				else :
					C_dict[float(all_temps[i])] = [float(all_Cs[i])]
					X_dict[float(all_temps[i])] = [float(all_Xs[i])]
	
		#Find statistical errors on these
		C_temps, C_plots, C_stdevs = [], [], []
		for C_temp in C_dict :
			C_temps.append(C_temp)
			C_plots.append(get_mean(C_dict[C_temp]))
			C_stdevs.append(get_stdev(C_dict[C_temp]))
	
		X_temps, X_plots, X_stdevs = [], [], []
		for X_temp in X_dict :
			X_temps.append(X_temp)
			X_plots.append(get_mean(X_dict[X_temp]))
			X_stdevs.append(get_stdev(X_dict[X_temp]))
	
		#print(C_dict)
		print("cstevs\n", C_stdevs)
		print("xstevs\n", X_stdevs)
	
		print("Plotting results...")
		
		plt.rc('xtick', labelsize=22) 
		plt.rc('ytick', labelsize=22)
		plt.figure()
		plt.errorbar(C_temps, C_plots, fmt='+', yerr=C_stdevs)
		plt.title("Heat Capacity, N =" + str(N) + ", B = " +str(h), size =22)
		plt.xlabel("Temp", size =22)
		plt.ylabel("Heat Capacity", size =22)
		plt.rc('xtick', labelsize=22) 
		plt.rc('ytick', labelsize=22)
		 
		
		plt.figure()
		plt.errorbar(X_temps, X_plots,fmt='+', yerr=X_stdevs)
		plt.title("Magnetic Susceptibility, N = " + str(N) + ", B = " +str(h), size =22)
		plt.ylabel("Susceptibility", size =22)
		plt.xlabel("Temp", size =22)
		plt.rc('xtick', labelsize=20) 
		plt.rc('ytick', labelsize=20) 
	
	plt.rc('xtick', labelsize=22) 
	plt.rc('ytick', labelsize=22)

	plt.show()


