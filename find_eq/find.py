import sys, random, time, copy
import numpy as np 
from matplotlib import pyplot as plt

sys.path.insert(0, '../..')
sys.path.insert(0, '..')
import main

#Run for 6000 steps on an 100x100 lattice for 2 indep lattices
#Find when spins 

def add_to_plot(rows, energy = True) :
	"""Plots energy (true) or spin from rows output on time steps
	Adds to current figure"""
	
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

	spins = rows[:,2]


	if energy :
		plt.ylabel("Lattice Energy / Joules/(exchange_energy*h_bar^2)")
		plt.plot(steps, energies)
		#plt.figure()
		#plt.plot(steps[1:], log_energies[1:])
	else :
		plt.ylabel("Avg spin /site")
		plt.plot(steps, abs(np.array(spins)))

def find_eq(n_lattice, N, T, h, n_max, typ, graphs) :
	"""Runs two lattices side by side to check for the number of steps required to reach equillibrium"""
	
	lattices = []
	lines =[]
	leg = []
	spin_fig = plt.figure()

	for i in range(n_lattice):

		lattice = main.lattice(N, T, field = h, exch_energy = main.J_e, mag_moment = main.mu_e)
		
		lattices.append(lattice)
	

		lines.append(lattice.progress_n(n_max, graphs=graphs))
		plt.title("Lattice" + str(i))
		plt.rcParams.update({'font.size': 16})

		plt.figure(spin_fig.number)
		add_to_plot(lines[i], energy = False)
		leg.append("Lattice" +str(i))
	

	plt.title("Evolution of Spin for N = " + str(N) + ", T = " + str(T) + ", h = " + str(h))
	plt.rcParams.update({'font.size': 22})
	plt.legend(leg)
	
	return lattices

if __name__ == "__main__" :

	find_eq(2, 50, 0.5, 0, 600, 'r', 6)
	plt.rcParams.update({'font.size': 22})
	plt.show()