import sys, random, time, copy
import numpy as np 
import scipy.constants as const
import matplotlib as mpl
from matplotlib import pyplot as plt

"""Main script for the ising model simulation
	Run this script with parameters (at bottom of file) to progress a single lattice
	Selected number of times, temperature, field etc
	Prints the results to a heatmap figure to observe, amongst other functionality"""

mu_e = 1 #const.physical_constants["electron mag. mom."][0] #magnetic moment of an electron
h_bar = 1 #const.hbar #Used throughout - check whether used correctly though.
J_e = 1 #spin exchange energy
k_b = 1 #const.k # boltzman constant

class lattice_site :
	"""A class for storing a spin site with its interactions with its neighbours"""

	def __init__(self, spin = 0, nearest_neighbours=[]) :
		"""Initialises a lattice site with random spin (unless specified as +/-1)"""
		
		self.nearest_neighbours = nearest_neighbours

		if spin == 1 or spin == -1  :
			
			self.spin = spin

		else :
			if spin != 0 :
				print("Warning: spin != 0, read " +str(spin), end="\r")

			random_number = np.random.random_sample()
			
			if random_number < 0.5 :
				self.spin = -1
			else :
				self.spin = +1
	
	def flip_spin(self, alt_energy, temp, field=0, exch_energy = J_e, mag_moment = mu_e) : 
		"""If favourable, spin is flipped and energy/spin changes returned
		Calculates alternative energy and progresses if < 0 or boltz prob
		Implementation of metropolis algorithm"""

		flipped, new_alt, dE, dMag = False, alt_energy, 0, 0

		#If not favourable, flip with boltz probability. Else flip (if favourable)
		if alt_energy <= 0 \
			or (alt_energy > 0 and np.random.random_sample() <  np.exp((-1. * alt_energy) /(k_b*temp))) :

			#Update spin, return flag, change in spin+energy
			self.spin *= -1
			dMag = -2 * self.spin
			dE = copy.copy(alt_energy)
			flipped = True
			new_alt = -2 * self.calculate_energy(field=0, exch_energy = J_e, mag_moment = mu_e) 

			#TODO - not updating neighbours!

		return flipped, new_alt, dE, dMag

	def calculate_energy(self, field=0., exch_energy = J_e, mag_moment = mu_e) : 
		"""Calculate and return the energy contribution to the lattice by this spin site
		Dependent on field and nearest neighbour spin interactions."""

		#find the spin interactions
		spin_interaction = 0 

		for neighbour in self.nearest_neighbours :

			spin_interaction += -1. * exch_energy * neighbour.spin * self.spin * h_bar ** 2 

		#Find the field contribution
		field_contribution = -1. * self.spin*h_bar * mag_moment * field 

		return spin_interaction + field_contribution

	def get_spin(self) :
		"""Return the spin of the site"""

		return self.spin

	def set_neighbours(self, neighbour_list) :
		"""Set the nearest neighbours"""

		self.nearest_neighbours = neighbour_list

	def set_location(self, row_pos, col_pos) :
		"""Sets location of the lattice site"""

		self.row_pos = row_pos
		self.col_pos = col_pos

	def get_location(self) :

		return self.row_pos, self.col_pos

class lattice :
	"""Stores an NxN lattice_array of 'lattice_site's and has a number of useful functons pertaining to a whole lattice"""
	
	def __init__(self, N, typ = "r") :
		"""Make an NxN lattice of lattice_site s
		type corresponds to initial spin assignment (r for random)
		returns an NxN dictionary of lattice sites keyed by [x][y] integer location"""

		self.size = N

		#Determine type (random, aligned)
		if typ == "r" :
			spin = 0
		elif typ == "a+" :
			spin = 1
		elif typ == "a-" :
			spin = -1
		else :
			print("Lattice type", typ, "not understood in initialisation.")
			print("Try r, a+, a-.")
			sys.exit()

		#Make the initial arrays
		init_array = np.zeros((N,N)) + spin
		lattice = np.empty((N,N), dtype=lattice_site)

		#Vectorise initialisation - spin argument is held in init_array 
		v_lattice_site = np.vectorize(lattice_site)
		lattice[:,:] = v_lattice_site(init_array)

		#update the object
		self.lattice_array = lattice

		#Set the neighbours, locations and set the lattice dictionary
		self.set_all_neighbours()

		#Set the net lattice energy and spin
		self.net_energy = self.get_net_energy()
		self.net_spin = self.get_net_spin()

		#Set arrays of the current deltas in mag and energy if spin flipped
		v_e_calc = np.vectorize(lattice_site.calculate_energy)
		self.current_alt_Es = -2 * v_e_calc(self.lattice_array)

	def set_all_neighbours(self) :
		"""Sets all the lattice_array's neighbour lists to the appropriate neighbours 
		for every lattice site in the lattice. Also sets their location"""

		N = self.size

		for row in range(N) :
			for col in range(N) :

				#periodic x conditions
				next_row = (row + 1) % self.size
				next_col = (col + 1) % self.size
				prev_row = (row - 1) % self.size
				prev_col = (col - 1) % self.size
				
				neighbours = [self.lattice_array[prev_row, col], self.lattice_array[next_row, col], self.lattice_array[row, prev_col], self.lattice_array[row, next_col]]
				
				self.lattice_array[row, col].set_neighbours(neighbours)
				self.lattice_array[row, col].set_location(row, col)

		return self.lattice_array

	def __repr__(self, loc = False) : 
		"""Prints the grid to a command line - useful for testing
		Prints with locations (true) or without (false)"""
		
		lines = ""
		
		for row in range(self.size) :
			
			line = ""
			
			for col in range(self.size) :
				
				location = str(row) + ", " + str(col)
				spin = self.lattice_array[row, col].spin
				
				if not loc :
					location = ""

				if spin < 0 : 
					line += str(location) + "+"
				else : 
					line += str(location) + "-"
		
			lines += line + "\n"
		
		print(lines)
	
	def produce_spin_list(self) :
		"""Returns the lattice spins as a list of lines"""

		lines = ""
		for row in range(self.size) :
			line = ""
			for col in range(self.size) :
				if self.lattice_array[row, col].spin == -1 :
					line += "+ "
				else : 
					line +=  "- "
			lines += line + "\n"

		return(lines)

	def add_to_plots(self, ax) :
		""""Plots the matrix onto an axis of a subplot, passed as an argument"""
		
		#get a matrix of +/- 1 values 
		v_spin = np.vectorize(lattice_site.get_spin)
		row_col_matrix = v_spin(self.lattice_array)

		#ax.set_xlabel("Col positions", labelpad=-3.5)
		#ax.set_ylabel("Row positions")

		#Heatmap customisations
		cmap= mpl.colors.ListedColormap(['k', 'w'])
		bounds = [-1, 0, 1]
		norm= mpl.colors.BoundaryNorm(bounds, cmap.N)
		heatmap = ax.imshow(row_col_matrix, cmap=cmap, interpolation='none', norm=norm)

		return heatmap

	def graph_locations(self) :
		"""Plots the grid as a heatmap of +/- lattice locations"""

		#get a matrix of +/- 1 values 
		v_spin = np.vectorize(lattice_site.get_spin)
		row_col_matrix = v_spin(self.lattice_array)

		#initialise plot
		fig, ax = plt.subplots()
		
		ax.set_xlabel("Col positions")
		ax.set_ylabel("Row positions")

		#Heatmap customisations
		cmap= mpl.colors.ListedColormap(['k', 'w'])
		bounds = [-1, 0, 1]
		norm= mpl.colors.BoundaryNorm(bounds, cmap.N)
		heatmap = ax.imshow(row_col_matrix, cmap=cmap, interpolation='none', norm=norm)
		
		plt.colorbar(heatmap, ticks=[-1,1])
		
		return ax

	def progress(self, T, field=0, exch_energy = J_e, mag_moment = mu_e) :
		"""Progresses the lattice by stepping through a random row, col
		Then use lattice_site check spin flip method. Updates the net energy and spin also"""

		#print("Beforest", self.current_alt_Es)

		v_prog = np.vectorize(lattice_site.flip_spin)
		
		#print("Before", self.current_alt_Es)

		#Get matrix of flip_bools, energy changes, Magnetisation changes
		flips, new_alts, dEs, dMs = v_prog(self.lattice_array, self.current_alt_Es, T, field=field, exch_energy=exch_energy, mag_moment=mag_moment)

		#Total lattice changes
		DE = dEs.sum()
		DM = dMs.sum()

		#Update the lattice state
		self.current_alt_Es = new_alts

		self.net_energy += DE
		self.net_spin += DM 

		#print("Before again", self.current_alt_Es)

		#calc percent of flips
		tot, flip_count = 0, 0
		for row in range(len(flips)) :
			for col in range(len(flips)) :
				tot += 1 
				if flips[row, col] :
					flip_count += 1
		flip_perc = int(round(100 * flip_count/tot))
		
		return flip_perc

	def progress_n(self, n, T, field=0, exch_energy = J_e, mag_moment = mu_e, graphs = 0, silent = False) :
		"""Progresses the process by n steps- graphs initial, every 'display'th and final
		returns list of lines of steps for later graphing
		step, energy, spin"""

		lns = []

		#Set subplot parameters
		if graphs == 0 : 
			pass
		
		else :
			#Make even
			if graphs % 2 != 0 :
				graphs -= 1

			#Make the figure
			if graphs == 2 :
				fig, axarr = plt.subplots(1, 2)
			else : 
				fig, axarr = plt.subplots(2, int(graphs/2))
			
			#Start after initial, set other plot positions
			plot_row, plot_col = 0, 1
			mid_graphs = (np.round(np.linspace(0, n, num = graphs))[1:]).astype(int)

			#Plot initial
			self.add_to_plots(axarr[0,0])
			axarr[0,0].set_title("Initial (0th) grid, T = " + str(T))
		
		#Print the initial state
		if not silent :
			print("initial net lattice energy, spin:", self.net_energy, self.net_spin)

		#Iterate + progress for n steps
		for i in range(n) :

			#Graph if required (first - e.g. i = 0 -> no progressions)
			if i in mid_graphs :
				
				self.add_to_plots(axarr[plot_row, plot_col])
				axarr[plot_row, plot_col].set_title(str(i) + "th step, T = " + str(T))

				plot_col += 1
				
				if plot_col >= graphs/2 :
					plot_row += 1
					plot_col = 0

			#Add to output lines
			lns.append( str(i) + ", " +  str(self.net_energy) + ", " + str(self.net_spin))

			#progress and append info
			flip_perc = self.progress(T)

			#Track progress
			if not silent :
				print(i + 1, "/", n, "progress, flipped", flip_perc, "% of sites",end="\r" ) 
	
		#Append the final state
		lns.append( str(i+1) + ", " +  str(self.net_energy) + ", " + str(self.net_spin))

		#Plot the final state
		if graphs != 0 :
			
			heatmap = self.add_to_plots(axarr[plot_row, plot_col])
			axarr[plot_row, plot_col].set_title("Final (" + str(i+1) + "th) state, T = " + str(T))
			
			#plt.colorbar(heatmap, ticks=[-1,1])

		#Print the final state
		if not silent :
			print()
			print("Final (" + str(i+1) + "th) net lattice energy, spin:", self.net_energy, self.net_spin)

		return lns

	def get_net_energy(self) :
		"""Return the net magnetic energy in the lattice"""
		
		v_calc = np.vectorize(lattice_site.calculate_energy)
		
		#Grid of energy contributions for each site
		energies = v_calc(self.lattice_array)

		return energies.sum()

	def get_net_spin(self) :
		"""Return the net integer spin in the lattice"""
		
		v_spin = np.vectorize(lattice_site.get_spin)
		spins = v_spin(self.lattice_array)

		return spins.sum()

def print_lines_to_file(filename, lines) :
	"""Prints the lines passed to a file"""

	print("Writing lines to file", filename, "in cwd.")
	with open(filename, 'w') as f :

		for row in lines :
			f.write(row + "\n")

	print("Success")

def read_in_from_file(filename) :

	print("Reading from", filename, "in cwd.") 

	lines = [] 

	with open(filename, 'r') as f :

		for line in f :
			lines.append(line.strip())

	print("Success.")

	return lines

def plot_one(rows, energy = True) :
	"""Plots energy (true) or spin from rows output"""

	plt.figure()
	
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
		plt.plot(steps, energies, 'rx')
		#plt.figure()
		#plt.plot(steps[1:], log_energies[1:])
	else :
		plt.ylabel("Net integer spin")
		plt.plot(steps, spins, 'b*')

def plot_both(rows) :
	"""Takes the step, energy, spin strings and plots"""

	ax1 = (plt.figure()).add_subplot(111)
	ax2 = ax1.twinx()
	
	ax1.set_title("Evolution of net lattice energy and spin with time.")
	ax1.set_xlabel("Time steps taken")

	for i in range(len(rows)):
		
		splt = rows[i].split(", ")

		for el in range(len(splt)) :
			splt[el] = float(splt[el])

		rows[i] = splt

	rows = np.array(rows)

	steps = rows[:,0]
	energies = rows[:,1]
	spins = rows[:,2]

	ax1.set_ylabel("Lattice Energy / Joules/(exchange_energy*h_bar^2)")
	ax1.tick_params('y', colors='b')
	ax1.plot(steps, energies, 'rx')


	ax2.set_ylabel("Net integer spin")
	ax2.plot(steps, spins, 'b*')

if __name__ == "__main__" :

	################################################	
	############# runtime parameters ###############
	
	#size of lattice NxN
	N = 1000
	
	#type of lattice (r andom, a+ a- ligned)
	typ = 'r' 
	
	#field
	h = 0 

	#Temp in J_e/k_b (set = 1)
	temp = 1 

	#steps - whole lattice progressions
	n = 3
	
	#graphing

	#number of graphs to display
	graphs = 8
	################################################

	#Init lattice
	print("Initialising...")
	init_start = time.clock()

	lattice = lattice(N, typ=typ)

	init_end = time.clock() 

	print()
	print(init_end-init_start, "seconds to initialise")
	print()

	#Progress the lattice n times 
	lines = lattice.progress_n(n, temp, graphs=graphs, field = h)

	end = time.clock()

	print()
	print(end - init_end, "seconds to progress", n, "times.", N,"x",N, "grid")
	print((end-init_end)/n,"seconds per step.")
	print()

	plt.show()