import sys
import numpy as np 
from matplotlib import pyplot as plt

sys.path.insert(0, '../..')
import main

"""A script to perform runs at each temperature
Tracks the evolution of net energy and spin per site in a graph
Useful for visualisation"""

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
	N = 100 
	
	#type of lattice (r andom, a+ a- ligned)
	typ = 'a+' 
	#field
	h = 0

	T_i, T_f = 0.1, 5
	#temp_list = [2, 1]
	temp_list = [x for x in np.linspace(T_i, T_f, num=10)]

	#steps - whole lattice progressions
	n = 100
	
	#graphing
	#number of graphs to display
	graphs = 6

	#Plot of all the energies
	i = 0
	all_plot = plt.figure()
	leg = []
	
	#initial lattice
	lat = main.lattice(N, temp_list[0], field=h, typ=typ)

	#iterate through temps, Temp in J_e/k_b (set = 1)
	for temp in temp_list :
		
		#Preserve spin state, regen random numbers and reset temp 
		lat.new_run(temp)

		print()
		i += 1
		print("Running for temp =", ("%.2f" % temp), ":", i, "/", len(temp_list))
		print()
		
		#Progress the lattice n times 
		lines = lat.progress_n(n, graphs=graphs)

		main.print_lines_to_file("T_" + ("%.2f" % temp) + "_n_" + str(n) + "_N_" + str(N) + ".txt" , lines)
		plt.rcParams.update({'font.size': 18})
		plt.savefig("T_" + ("%.2f" % temp) + ".png")

		plt.figure(all_plot.number)
		add_lines_to_plot(lines)
		leg.append("T = " + str(temp))

		final_state = lat

	plt.legend(leg)
	plt.rcParams.update({'font.size': 22})
	plt.show()
