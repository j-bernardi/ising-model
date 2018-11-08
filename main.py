import sys, random, time, copy
import numpy as np 
import scipy.constants as const
import matplotlib as mpl
from matplotlib import pyplot as plt

"""Main script for the ising model simulation
	Run this script with parameters (at bottom of file) to progress a single lattice
	Selected number of times, temperature, field etc
	Prints the results to a heatmap figure to observe, amongst other functionality"""

mu = 1 #Mag permeability of free space
mu_e = 1 #const.physical_constants["electron mag. mom."][0] #magnetic moment of an electron
h_bar = 1 #const.hbar
J_e = 1 #spin exchange energy
k_b = 1 #const.k # boltzman constant

class lattice_site :
	"""A class for storing a spin site with its interactions with its neighbours"""

	def __init__(self, row, col, spin = 0, nearest_neighbours=[]) :
		"""Initialises a lattice site with random spin (unless specified as +/-1)"""
		
		self.row, self.col = row, col

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

		self.alt_energy = -1
	
	def flip_spin(self, temp, log_rand_num_kbt, field, exch_energy, mag_moment) : 
		"""If favourable, spin is flipped and energy/spin changes returned
		Calculates alternative energy and progresses if < 0 or boltz prob
		Implementation of metropolis algorithm"""

		flipped, dE, dMag = False, 0, 0

		#If not favourable, flip with boltz probability. Else flip (if favourable)
		if self.alt_energy <= 0 \
			or (self.alt_energy > 0 and log_rand_num_kbt <  -1. * self.alt_energy ):

			dMag = -2 * self.spin
			self.spin *= -1
			dE = self.alt_energy
			flipped = True

		return flipped, dE, dMag

	def calculate_energy(self, field, exch_energy, mag_moment) : 
		"""Calculate and return the energy contribution to the lattice by this spin site
		Dependent on field and nearest neighbour spin interactions."""

		#find the spin interactions
		spin_interaction = 0 

		for neighbour in self.nearest_neighbours :

			spin_interaction += -1. * exch_energy * neighbour.spin * self.spin * h_bar ** 2 

		#Find the field contribution
		field_contribution = -1. * self.spin * h_bar * mag_moment * field 

		#print("field cont", field_contribution)
		#print("field", field)

		return spin_interaction + field_contribution

	def calculate_alt_energy(self, field, exch_energy, mag_moment):

		self.alt_energy = \
		 -2 * self.calculate_energy(field= field, exch_energy = exch_energy, mag_moment = mag_moment)

	def get_spin(self) :
		"""Return the spin of the site"""

		return self.spin

	def get_neighbour_locs(self, N) :

		prev_row = (self.row - 1) % N
		next_row = (self.row + 1) % N

		prev_col = (self.col - 1) % N
		next_col = (self.col + 1) % N

		return prev_row, next_row, prev_col, next_col

	def set_neighbours(self, neighbour_list) :
		"""Set the nearest neighbours"""

		self.nearest_neighbours = neighbour_list

class lattice :
	"""Stores an NxN lattice_array of 'lattice_site's and has a number of useful functons pertaining to a whole lattice"""
	
	def __init__(self, N, T, field = 0, exch_energy = J_e, mag_moment = mu_e, typ = "r") :
		"""Make an NxN lattice of lattice_site s
		type corresponds to initial spin assignment (r for random)
		returns an NxN dictionary of lattice sites keyed by [x][y] integer location"""

		#Lattice parameters
		self.size = N
		self.T = T
		self.field = field
		self.exch_energy = exch_energy
		self.mag_moment = mag_moment

		#Set type of array/spin
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

		#Make and set the lattice array
		self.lattice_array = np.empty((N,N), dtype=lattice_site)

		for row in range(N) :
			for col in range(N) :
				self.lattice_array[row, col] = lattice_site(row, col, spin=spin)

		#Set the neighbours of the lattice sites
		self.set_all_neighbours()

		#Save the current alternative energies (e.g. if spin flipped)
		for row in range(N) :
			for col in range(N) :
				self.lattice_array[row, col].calculate_alt_energy(field= field, exch_energy = exch_energy, mag_moment = mag_moment)

		#Calculate set of random numbers once only for each lattice
		self.rand_array = np.log(np.random.rand(N**2))* k_b * self.T

		#Find initial properties
		self.net_energy = self.get_net_energy()
		self.net_spin = self.get_net_spin()	

	def new_run(self, T, field = 0, exch_energy = J_e, mag_moment = mu_e) :
		"""Uses the same lattice but updates the temp, field and random numbers"""

		self.T = T
		self.field = field
		self.exch_energy = exch_energy
		self.mag_moment = mag_moment

		#Calculate new set of random numbers for this run
		self.rand_array = np.log(np.random.rand(self.size**2))* k_b * self.T

	def permute_rand_array(self) :
		#Iterate through the random list, permute each by 1

		np.random.shuffle(self.rand_array)
	
	def set_all_neighbours(self) :
		"""Sets all the lattice_array's neighbour lists to the appropriate neighbours 
		for every lattice site in the lattice"""

		for row in range(self.size) :
			for col in range(self.size) :
				
				#Get locations
				prev_row, next_row, prev_col, next_col = self.lattice_array[row, col].get_neighbour_locs(self.size)
				
				#Find neighbours
				neighbours = [self.lattice_array[prev_row, col], self.lattice_array[next_row, col], self.lattice_array[row, prev_col], self.lattice_array[row, next_col]]
				
				#set neighbours
				self.lattice_array[row, col].set_neighbours(neighbours)

	def progress_n(self, n, graphs = 0) :
		"""Progresses the process by n steps- graphs initial, every 'display'th and final
		returns list of lines of steps for later graphing
		step, energy, spin"""

		lns = []
		
		#List of locations to progress
		locs = self.get_all_locations()

		plot_row, plot_col, mid_graphs, axarr = self.initialise_figure(n, graphs)
		
		#Print the initial state
		print("initial net lattice energy, spin:", self.net_energy, self.net_spin)

		#Iterate + progress for n steps
		for i in range(n) :

			#Add plot to graphs if in list
			if i in mid_graphs :
				self.add_to_plots(axarr[plot_row, plot_col])
				axarr[plot_row, plot_col].set_title(str(i) + "th step, T = " + ("%.2f" % self.T))
				plot_col += 1
				if plot_col >= graphs/2 :
					plot_row += 1
					plot_col = 0

			#Add to output lines
			lns.append( str(i) + ", " +  str(self.net_energy) + ", " + str(self.net_spin/(self.size **2)))

			#Shuffle locations and progress
			random.shuffle(locs)
			flip_perc = self.progress(locs)

			#Track progress
			print(i + 1, "/", n, "progress. Flipped", flip_perc, "% of sites", end="\r" ) 
	
		
		#Append the final state
		lns.append( str(i+1) + ", " +  str(self.net_energy) + ", " + str(self.net_spin/(self.size ** 2)))

		#Plot the final state
		if graphs != 0 :
			heatmap = self.add_to_plots(axarr[plot_row, plot_col])
			axarr[plot_row, plot_col].set_title("Final (" + str(i+1) + "th) state, T = " + ("%.2f" % self.T))
			#plt.colorbar(heatmap, ticks=[-1,1])

		#Print the final state
		print()
		print("Final (" + str(i+1) + "th) net lattice energy, spin:", self.net_energy, self.net_spin)

		return lns

	def progress(self, locs) :
		"""Progresses the lattice by stepping through a random row, col
		Then use lattice_site check spin flip method. Updates the net energy and spin also"""

		#print("locs", locs)
		flip_count, tot = 0, 0

		#Update each grid cell in a random order
		for loc in locs :

			row, col = loc[0], loc[1]
			
			flipped, dE, dM = \
				self.lattice_array[row, col].flip_spin(\
				 self.T, self.rand_array[row*self.size + col], field=self.field, exch_energy = self.exch_energy, mag_moment = self.mag_moment)
			
			tot += 1	
			
			#Detect a spin flip
			if flipped :
				flip_count += 1

				#Update the alt_energy array for cell and neighbours
				prev_row, next_row, prev_col, next_col = \
				 self.lattice_array[row,col].get_neighbour_locs(self.size)
				
				self.lattice_array[row, col].calculate_alt_energy(field= self.field, exch_energy = self.exch_energy, mag_moment = self.mag_moment)
				self.lattice_array[prev_row, col].calculate_alt_energy(field= self.field, exch_energy = self.exch_energy, mag_moment = self.mag_moment)
				self.lattice_array[next_row, col].calculate_alt_energy(field= self.field, exch_energy = self.exch_energy, mag_moment = self.mag_moment)
				self.lattice_array[row, prev_col].calculate_alt_energy(field= self.field, exch_energy = self.exch_energy, mag_moment = self.mag_moment)
				self.lattice_array[row, next_col].calculate_alt_energy(field= self.field, exch_energy = self.exch_energy, mag_moment = self.mag_moment)

			self.net_energy += dE
			self.net_spin += dM

		#print(self.rand_array)
		self.permute_rand_array()

		flip_perc = int(round(100 * flip_count/tot))

		return flip_perc

	def get_all_locations(self):
		"""Return a list of random locations in the grid"""

		locs = []

		#Shuffle the order of the rows and cols
		rows, cols = list(range(self.size)), list(range(self.size))
		
		for row in rows :
			for col in cols :
				locs.append((row, col))

		return locs

	def initialise_figure(self, n, graphs) :
		"""Initialises the plot space for plotting of the outputs
		Defaults to plotting 1x2 or 2x(g/2) graphs"""

		if graphs == 0 :
			return 0, 0, [], []

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
		axarr[0,0].set_title("Initial (0th) grid, T = " + ("%.2f" % self.T))

		return plot_row, plot_col, mid_graphs, axarr

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

	def get_net_energy(self) :
		"""Return the net magnetic energy in the lattice"""
		
		#Define vectorisation
		v_calc = np.vectorize(lattice_site.calculate_energy)

		#Grid of energy contributions for each site
		energies = v_calc(self.lattice_array, self.field, self.exch_energy, self.mag_moment)

		return energies.sum()

	def get_net_spin(self) :
		"""Return the spin in the lattice,avg per site"""
		
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
	"""Plots energy (true) or spin from rows output on time steps"""

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
	"""Takes energy, spin strings and plots on steps on a two-axis figure"""

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
	N = 50
	
	#type of lattice (r andom, a+ a- ligned), Temp in J_e/k_b (set = 1)
	typ = 'r' 
	temp = 2.5

	#field
	h = 0 

	#steps - whole lattice progressions
	n = 50
	
	#graphing

	#number of graphs to display
	graphs = 0
	################################################

	#Init lattice
	print("Initialising...")
	init_start = time.clock()

	lattice = lattice(N, temp, field = h, exch_energy = J_e, mag_moment = mu_e, typ=typ)

	init_end = time.clock() 

	print()
	print(init_end-init_start, "seconds to initialise")
	print()
	
	#Progress the lattice n times 
	lines = lattice.progress_n(n, graphs=graphs)

	end = time.clock()

	print()
	print(end - init_end, "seconds to progress", n, "times.", N,"x",N, "grid")
	print((end-init_end)/n,"seconds per step.")
	print()

	plt.show()
	