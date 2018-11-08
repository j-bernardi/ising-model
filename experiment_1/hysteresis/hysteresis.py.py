import sys, copy
import numpy as np 
from matplotlib import pyplot as plt

sys.path.insert(0, '../..')
sys.path.insert(0, '..')
import main

"""A script to perform runs at each temperature
Tracks the evolution of net energy and spin per site in a graph
Useful for visualisation
Runs forward and backward to check the same"""

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
	N = 50

	#type of lattice (r andom, a+ a- ligned)
	typ = 'r' 
	#field
	hs_init = np.linspace(0, -0.2, num=3)
	hs = [x for x in np.linspace(-0.2, 0.2, num = 11)]

	rev = copy.copy(hs)
	rev.reverse()

	temp = 2

	#steps - whole lattice progressions - use data after j (e.g. j steps to eq)
	n = 200
	n_T_eq = 400
	
	#graphing
	#number of graphs to display
	graphs = 0

	#initial lattice
	lat = main.lattice(N, temp, field=hs[0], typ=typ)

	#Progress from 0 to max

	#iterate through temps, Temp in J_e/k_b (set = 1)
	init_count = 0
	for h in hs_init :

		init_count += 1
		print("h =", h, "-", init_count, "/", len(hs_init))
		
		these_ens, these_mags = [], []
		#Preserve spin state, regen random numbers and reset temp 
		lat.new_run(temp, h)
		lat.progress_n(n_T_eq)

		"""
		#Progress the lattice n times 
		lines = lat.progress_n(n, graphs=graphs)
		
		for l in lines :
			splt = l.split(",")
			l_en, l_mag = float(splt[1].strip()), float(splt[2].strip())
			
			these_ens.append(l_en/(N**2))
			these_mags.append(l_mag)

		ens.append(get_mean(these_ens))
		stdev_ens.append(get_stdev(these_ens))
		mags.append(get_mean(these_mags))
		stdev_mags.append(get_stdev(these_mags))
		"""
	
	#Plot of all the energies
	ens, mags = [], []
	stdev_ens, stdev_mags = [], []
	e_leg, m_leg = ["forward"], ["forward"]

	temp_count = 0
	for h in hs :

		temp_count += 1
		print("h =", h, "-", temp_count, "/", 2*len(hs))
		
		these_ens, these_mags = [], []
		#Preserve spin state, regen random numbers and reset temp 
		lat.new_run(temp, h)
		lat.progress_n(n_T_eq)

		#Progress the lattice n times 
		lines = lat.progress_n(n, graphs=graphs)

		for l in lines :
			splt = l.split(",")
			l_en, l_mag = float(splt[1].strip()), float(splt[2].strip())
			
			these_ens.append(l_en/(N**2))
			these_mags.append(l_mag)

		ens.append(get_mean(these_ens))
		stdev_ens.append(get_stdev(these_ens))
		mags.append(get_mean(these_mags))
		stdev_mags.append(get_stdev(these_mags))
	
	rens, rmags = [], []
	rstdev_ens, rstdev_mags = [], []
	e_leg.append("reverse")
	m_leg.append("reverse")
	
	for h in rev :

		temp_count += 1
		print("h =", h, "-", temp_count, "/", 2*len(rev))
		
		these_ens, these_mags = [], []
		#Preserve spin state, regen random numbers and reset temp 
		lat.new_run(temp, h)
		lat.progress_n(n_T_eq)

		#Progress the lattice n times 
		lines = lat.progress_n(n, graphs=graphs)

		for l in lines :
			splt = l.split(",")
			l_en, l_mag = float(splt[1].strip()), float(splt[2].strip())
			
			these_ens.append(l_en/(N**2))
			these_mags.append(l_mag)

		rens.append(get_mean(these_ens))
		rstdev_ens.append(get_stdev(these_ens))
		rmags.append(get_mean(these_mags))
		rstdev_mags.append(get_stdev(these_mags))


	plt.figure()
	plt.errorbar(hs, ens, yerr=stdev_ens)
	plt.errorbar(rev, rens, yerr=rstdev_ens)
	plt.title("Energy per Site, N = " + str(N) + " T = " + str(temp), size = 22)
	plt.xlabel("Applied field H", size = 22)
	plt.ylabel("Avg energy", size = 22)
	plt.legend(e_leg)
	

	plt.figure()
	plt.errorbar(hs, mags, yerr=stdev_mags)
	plt.errorbar(rev, rmags, yerr=rstdev_mags)
	plt.title("Spin per Site, N = " + str(N) + " T = " + str(temp), size = 22)
	plt.xlabel("Applied field H", size = 22)
	plt.ylabel("Avg spin", size = 22)

	plt.legend(m_leg)

	plt.show()
