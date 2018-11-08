import sys, random, time
import numpy as np 
import scipy.constants as const
import matplotlib as mpl
from matplotlib import pyplot as plt

mu_e = 1 #const.physical_constants["electron mag. mom."][0] #magnetic moment of an electron
h_bar = 1 #const.hbar #Used throughout - check whether used correctly though.
J_e = 1 #spin exchange energy
k_b = 1 #const.k # boltzman constant

class lattice_site :

	def __init__(self,location, spin, nearest_neighbours) :
		"""Initialised with an integer vector (lattice site), spin 1, -1 (up, down) 
		List of the 4 nearest neighbours (assuming periodic BCs)"""
		
		self.location = location
		self.spin = spin
		self.nearest_neighbours = nearest_neighbours

	def flip_spin(self, temp, field=0, exch_energy = J_e, mag_moment = mu_e) : 
		"""Calculates energy required to flip (ed TD favourability)
		If favourable, spin is flipped, else done so only with boltmann probability"""

		flipped, dE, dMag = False, 0, 0 

		#Change in energy if it were to flip TODO: Check this is correct
		alt_energy = -2 * self.calculate_energy(field= field, exch_energy = exch_energy, mag_moment = mag_moment)

		#If not favourable, flip with boltz probability. Else flip (if favourable)
		if alt_energy <= 0 \
			or (alt_energy > 0 and np.random.random_sample() <  np.exp((-1. * alt_energy) /(k_b*temp))) :

			self.spin *= -1
			dMag = -2 * self.spin
			dE = alt_energy
			flipped = True

		return flipped, dE, dMag

	def calculate_energy(self, field=0., exch_energy = J_e, mag_moment = mu_e) : 

		spin_interaction = 0 

		#print("spin, loc", self.spin, self.location)

		for neighbour in self.nearest_neighbours :

			#print("neighbour loc, spin", neighbour.location, neighbour.spin)
			spin_interaction += -1. * neighbour.spin * self.spin * h_bar ** 2 

		#sys.exit()

		#setting h_bar = 1 here... TODO - decide if I want to do this
		spin_contribution = exch_energy * spin_interaction

		"""TODO: Check this is the correct calc. should be spin square?"""
		field_contribution = -1. * ((self.spin*h_bar) ** 2) * mag_moment * field 

		return spin_contribution + field_contribution

	def set_neighbours(self, neighbour_list) :
		"""Set the nearest neighbours"""

		self.nearest_neighbours = neighbour_list

class lattice :
	"""Stores an NxN lattice_dict of lattice_site s and has a number of useful functons pertaining to a whole lattice"""
	
	def __init__(self, N, type="r") :
		"""Make an NxN lattice of lattice_site s
		type corresponds to initial spin assignment (r for random)
		returns an NxN dictionary of lattice sites keyed by [x][y] integer location"""

		self.size = N

		#dictionary[row][col] e.g. constant time access for a known location.
		sites = {}

		#Create the lattice sites with initial spins
		for row in range(N) :
			sites[row] = {}
			for col in range(N) :
				
				#make a (pseudo) random spin
				random_number = np.random.random_sample()

				if random_number < 0.5 :
					random_spin = -1
				else :
					random_spin = +1
				
				#don't know the neighbours yet - but assign objects and spins
				sites[row][col] = lattice_site(np.array([row,col]), random_spin, [])

		#Set the neighbours and set the lattice dictionary
		sites = self.set_all_neighbours(sites)
		self.lattice_dict = sites

		self.net_energy = self.get_net_energy()
		self.net_spin = self.get_net_spin()

	def set_all_neighbours(self, temp_dict) :
		"""Takes a temporary dictionary, returns updated one"""

		N = self.size

		for row in range(N) :
			for col in range(N) :

				#periodic x conditions
				if row == 0 :
					next_row = row + 1
					prev_row = N - 1
				
				elif row == N-1 :
					next_row = 0
					prev_row = row - 1

				else :
					next_row = row + 1
					prev_row = row - 1 
				
				#Periodic y conditions
				if col == 0 :
					next_col = col + 1
					prev_col = N - 1
				
				elif col == N-1 :
					next_col = 0
					prev_col = col - 1

				else :
					next_col = col + 1
					prev_col = col - 1
				
				neighbours = [temp_dict[prev_row][col], temp_dict[next_row][col], temp_dict[row][prev_col], temp_dict[row][next_col]]
				
				temp_dict[row][col].set_neighbours(neighbours)

		return temp_dict

	def __repr__(self) : 

		lines = ""
		for row in range(self.size) :
			line = ""
			for col in range(self.size) :
				if self.lattice_dict[row][col].spin == -1 :
					line += "+ "
				else : 
					line +=  "- "
			lines += line + "\n"

		return(lines)

	def print_locations(self) :
		
		lines = ""
		
		for row in range(N) :
			
			line = ""
			
			for col in range(N) :
				
				location = lattice.lattice_dict[row][col].location
				spin = lattice.lattice_dict[row][col].spin
				
				if spin < 0 : 
					line += str(location) + "+"
				else : 
					line += str(location) + "-"
		
			lines += line + "\n"
		
		print(lines)

	def graph_locations(self) :

		row_col_matrix = np.zeros((len(self.lattice_dict), len(self.lattice_dict[0])))

		for row in range(len(self.lattice_dict)) :
			for col in range(len(self.lattice_dict[row])) : 
				row_col_matrix[row][col] = self.lattice_dict[row][col].spin

		fig, ax = plt.subplots()
		ax.set_xlabel("Col positions")
		ax.set_ylabel("Row positions")
		cmap= mpl.colors.ListedColormap(['k', 'w'])
		bounds = [-1, 0, 1]
		norm= mpl.colors.BoundaryNorm(bounds, cmap.N)
		heatmap = ax.imshow(row_col_matrix, cmap=cmap, interpolation='none', norm=norm)
		plt.colorbar(heatmap, ticks=[-1,1])
		return ax

	def progress(self, T, field=0, exch_energy = J_e, mag_moment = mu_e) :
		"""Progresses the lattice by stepping through a random row, col
		Then use lattice_site check spin flip method"""
		
		"""TODO - vectorise?"""

		rows = np.array(range(self.size) )
		cols = np.array(range(self.size) )
		
		random.shuffle(rows)
		random.shuffle(cols)

		"""Goes along a whole row still"""
		for row in rows :
			for col in cols :
				flipped, dE, dMag = self.lattice_dict[row][col].flip_spin(T)

				self.net_spin += dMag
				self.net_energy += dE

	def progress_n(self, n, T, field=0, exch_energy = J_e, mag_moment = mu_e, display=0):
		"""Progresses the process by n steps- print every 'display'th and final"""

		lns = []

		if display == 0: 
			display = n + 1

		print("initial net lattice energy, spin:", self.net_energy, self.net_spin)

		for i in range(n) :

			lns.append( str(self.net_energy) + ", " +  str(self.net_spin))

			lattice.progress(temp)
			
			print(i + 1, "/", n, "progress",end="\r" )

			if i % display == 0 :
				ax = lattice.graph_locations()
				ax.set_title(str(i) + "th step")

		print()

		ax = lattice.graph_locations()
		ax.set_title("Final (" + str(i) + "th) step")
		print(str(i) + "th net lattice energy, spin:", self.net_energy, self.net_spin)

		return lns

	def get_net_spin(self) :
		"""Return the net integer spin"""
		tot = 0 
		for row in range(len(self.lattice_dict)) :
			for col in range(len(self.lattice_dict[row])) :
				tot += self.lattice_dict[row][col].spin
		return tot

	def get_net_energy(self) :
		"""Return the net magnetic energy"""
		tot = 0
		for row in range(len(self.lattice_dict)) :
			for col in range(len(self.lattice_dict[row])) :
				tot += self.lattice_dict[row][col].calculate_energy()
		return tot

if __name__ == "__main__" :

	#conditions 
	N = 500 #size of lattice NxN
	typ = 'r' #type of lattice (r andom, a ligned)
	n = 1 #number of progressions (whole-lattice)
	d = 10 #display graph evert dth time
	h = 0
	temp = 1 #kelvin

	start = time.clock()
	
	lattice = lattice(N, type = typ)

	end = time.clock() 

	print(end-start, "seconds to initialise")

	lines = lattice.progress_n(n, temp, display=d, field = h)

	endend = time.clock()

	print(endend - end, "seconds to progress", n, "times.", N,"x",N, "grid")

	plt.show()