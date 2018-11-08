import sys
import numpy as np 
from matplotlib import pyplot as plt

sys.path.insert(0, '..')
import main

"""A script to perform runs at each temperature
Tracks the evolution of net energy and spin per site in a graph
Useful for visualisation"""

def get_mean(lst) :

	return sum(lst)/len(lst)

def get_stdev(lst) :

	mean = get_mean(lst)
	var = 0
	for l in lst :
		var += (mean - l) ** 2

	return np.sqrt(var/len(lst))

def add_lines_to_plot(rows, energy = True) :
	"""Plots energy (true) or spin from rows output"""

	if energy :
		tit = "energy" 
	else :
		tit = "spin"

	plt.title("Evolution of net lattice " + tit + " with time.")
	plt.xlabel("Time steps taken")

	for i in range(len(rows)):
		
		splt = rows[i].split(", ")

		for el in range(len(splt)) :
			splt[el] = float(splt[el])

		rows[i] = splt

	rows = np.array(rows)

	steps = rows[:,0]
	energies = rows[:,1]
	log_energies = np.log(abs(energies))
	spins = rows[:,2]
	log_spins = np.log(abs(spins))

	if energy :
		plt.ylabel("Lattice Energy / Joules/(exchange_energy*h_bar^2)")
		plt.plot(steps, energies)
		#plt.figure()
		#plt.plot(steps[1:], log_energies[1:])
	else :
		plt.ylabel("Net integer spin")
		plt.plot(steps, spins)

if __name__ == "__main__" :

	#size of lattice NxN

	#For bad plots
	Ns = [1,2]#[50, 150]

	#type of lattice (r andom, a+ a- ligned)
	typ = 'a+' 
	#field
	h = 0

	#temp_list = [2, 1]
	temp_list = [x for x in np.linspace(0.5, 1.5, num=3)]
	temp_list += [x for x in np.linspace(1.8, 3, num = 8)]
	temp_list += [x for x in np.linspace(3.5, 5, num = 4)]

	#steps - whole lattice progressions - use data after j (e.g. j steps to eq)
	n_eq_N = 300
	n_eq_T = 200
	n = 150
	
	#graphing
	#number of graphs to display
	graphs = 0
	
	#init plots
	en_plot = plt.figure()
	plt.xlabel("Temp", size = 22)
	plt.ylabel("Mean Energy", size = 22)
	plt.title("Mean energy on temperature for various N", size = 22)
	
	mag_plot = plt.figure()
	plt.xlabel("Temp", size = 22)
	plt.ylabel("Spin per site", size = 22)
	plt.title("Mean magnetisation on temperature for various N", size = 22)

	e_leg, m_leg = [], []

	#iterate through temps, Temp in J_e/k_b (set = 1)
	for N in Ns :
		temp_count = 0
		#initial lattice
		lat = main.lattice(N, temp_list[0], field=h, typ=typ)

		print("*** Run for N =", N, "***")
		
		#For this N
		ens, mags = [], []
		stdev_ens, stdev_mags = [], []

		lat.progress_n(n_eq_N, graphs=graphs)

		#Run for each temp and average the last points
		for temp in temp_list :

			these_ens, these_mags = [], []
			#Preserve spin state, regen random numbers and reset temp 
			lat.new_run(temp)
	
			print()
			temp_count += 1 
			print("Running for temp =", ("%.2f" % temp), ":", temp_count, "/", len(temp_list))
			print()

			lat.progress_n(n_eq_T, graphs=graphs)
			
			#Progress the lattice n times 
			lines = lat.progress_n(n, graphs=graphs)
	
			for l in lines :
				splt = l.split(",")
				l_en, l_mag = float(splt[1].strip()), float(splt[2].strip())
				
				these_ens.append(l_en/N ** 2)
				these_mags.append(abs(l_mag))
	
			ens.append(get_mean(these_ens))
			stdev_ens.append(get_stdev(these_ens))
			mags.append(get_mean(these_mags))
			stdev_mags.append(get_stdev(these_mags))
	
		#Write this temp to file
		with open("N_" + str(N) + "_n_" +str(n) + ".txt", 'w') as f :
			
			for i in range(len(temp_list)) :
	
				f.write(str(temp_list[i]) + ", " + str(ens[i]) + ", " + str(mags[i]))

		#plot for this N
		plt.rc('xtick', labelsize=22) 
		plt.rc('ytick', labelsize=22) 
		plt.figure(en_plot.number)
		plt.errorbar(temp_list, ens, yerr=stdev_ens)
		e_leg.append("N_" +str(N) +"_n_" +str(n))
		plt.rc('xtick', labelsize=22) 
		plt.rc('ytick', labelsize=22) 
		
	
		plt.figure(mag_plot.number)
		plt.errorbar(temp_list, mags, yerr=stdev_mags)
		m_leg.append("N_" +str(N) +"_n_" +str(n))
		plt.rc('xtick', labelsize=22) 
		plt.rc('ytick', labelsize=22) 
	
	plt.rc('xtick', labelsize=22) 
	plt.rc('ytick', labelsize=22)
	
	plt.figure(en_plot.number)
	plt.legend(e_leg)
	plt.rc('xtick', labelsize=22) 
	plt.rc('ytick', labelsize=22) 

	plt.figure(mag_plot.number)
	plt.legend(m_leg)
	plt.rc('xtick', labelsize=22) 
	plt.rc('ytick', labelsize=22) 
	
	plt.show()
