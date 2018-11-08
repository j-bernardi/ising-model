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

	#Plot all Ns
	Ns = [200]
	
	#type of lattice (r andom, a+ a- ligned)
	typ = 'a+' 
	#field
	h = 0

	#temp_list = [2, 1]
	temp_list = [x for x in np.linspace(0.5, 2.1, num=5)]
	temp_list += [y for y in np.linspace(2.2, 2.8, num = 10)]
	temp_list += [x for x in np.linspace(2.9, 5, num=5)]

	#steps - whole lattice progressions - use data after j (e.g. j steps to eq)
	n = 100
	j = 60
	
	#graphing
	#number of graphs to display
	graphs = 0
	
	#init plots
	en_plot = plt.figure()
	plt.xlabel("T / J/(exch_en*kb)")
	plt.ylabel("Energy / (J/exchange_en)")
	plt.title("Energy on temp")
	
	mag_plot = plt.figure()
	plt.xlabel("T / J/(exch_en*kb)")
	plt.ylabel("mean spin / site")
	plt.title("mean mag on temp")

	e_leg, m_leg = [], []

	#iterate through temps, Temp in J_e/k_b (set = 1)
	for N in Ns :
		
		#initial lattice
		lat = main.lattice(N, temp_list[0], field=h, typ=typ)

		print("*** Run for N =", N, "***")
		
		#Adjust for lattice size
		i = 0
		if N < 100 :
			n = 100
			j = 60 
		else :
			n = 30
			j =15
		
		#For this N
		ens, mags = [], []
		stdev_ens, stdev_mags = [], []

		#Run for each temp and average the last points
		for temp in temp_list :

			these_ens, these_mags = [], []
			#Preserve spin state, regen random numbers and reset temp 
			lat.new_run(temp)
	
			print()
			i += 1
			print("Running for temp =", ("%.2f" % temp), ":", i, "/", len(temp_list))
			print()
			
			#Progress the lattice n times 
			lines = lat.progress_n(n, graphs=graphs)
	
			for l in lines :
				splt = l.split(",")
				l_en, l_mag = float(splt[1].strip()), float(splt[2].strip())
				
				these_ens.append(l_en/N ** 2)
				these_mags.append(abs(l_mag))
	
			print("Using", len(these_ens) - j, "datapoints")
			
			ens.append(get_mean(these_ens[j:]))
			stdev_ens.append(get_stdev(these_ens[j:]))
			mags.append(get_mean(these_mags[j:]))
			stdev_mags.append(get_stdev(these_mags[j:]))
	
		#Write this temp to file
		with open("N_" + str(N) + "_n_" +str(n) + ".txt", 'w') as f :
			
			for i in range(len(temp_list)) :
	
				f.write(str(temp_list[i]) + ", " + str(ens[i]) + ", " + str(mags[i]))

		#plot for this N
		plt.figure(en_plot.number)
		plt.errorbar(temp_list, ens, yerr=stdev_ens)
		e_leg.append("N_" +str(N) +"_n_" +str(n))
		
	
		plt.figure(mag_plot.number)
		plt.errorbar(temp_list, mags, yerr=stdev_mags)
		m_leg.append("N_" +str(N) +"_n_" +str(n))
	
	
	plt.figure(en_plot.number)
	plt.legend(e_leg)

	plt.figure(mag_plot.number)
	
	ys = (1 -(np.sinh(2/np.array(temp_list)))**-4)**(1/8)
	plt.plot(temp_list, ys)
	
	#plt.axvline(2/np.log(1+np.sqrt(2)))

	m_leg.append("Theoretical_Mag")
	
	plt.legend(m_leg)
	
	plt.show()
