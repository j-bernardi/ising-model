import sys, random, time, copy
import numpy as np 
from matplotlib import pyplot as plt

sys.path.insert(0, '..')
import main


#Correlation time stuff
#t counts number of states

#<Q> = lim(T->large) 1/T sum 1-T of Q(t)

#Find timescale at which u(t) and u(t+dt) decreases
#Equillibrium time is a good estimate (a few?)

def get_mean(lst) :

	return sum(lst)/len(lst)

def get_stdev(lst) :

	mean = get_mean(lst)
	var = 0
	for l in lst :
		var += (mean - l) ** 2

if __name__ == "__main__" :

	################################################	
	############# runtime parameters ###############
	
	#size of lattice NxN
	N_list = [50]#, 50, 70, 100]
	
	#type of lattice (r andom, a+ a- ligned), Temp in J_e/k_b (set = 1)
	temp_list = [1.5, 2.0, 2.25, 2.4, 2.8, 2.9]

	#field
	h = 0 

	#Number of steps to take between Ns to get to eq
	n_eq_N = 400
	
	#Number of steps to take between Ts to get to eq
	n_eq_T = 400

	#Number of steps to avg over
	n = 300
	
	#number of graphs to display
	graphs = 0
	################################################

	tau_list = [x for x in np.linspace(1, n)]

	for N in N_list :
		
		print("\n*** FOR N = ", N, "***\n")

		tau_plot, mag_plot = plt.figure(), plt.figure()
		leg, m_leg = [], []
		temp_track = 0
		
		#Get the new lattice to equillibrium - check and manually accept
		if temp_list[0] > 2.27 :
			typ = 'r'
		else :
			typ = 'a+'

		lattice = main.lattice(N, temp_list[0], field = h, exch_energy = main.J_e, mag_moment = main.mu_e, typ=typ)
		print("Progessing to reach Eq for init")
		lattice.progress_n(n_eq_N, 0)

		for temp in temp_list :
			
			temp_track += 1
			print("\nTemp =", temp, "-", temp_track, "/", len(temp_list))
			
			#Init lattice
			lattice.new_run(temp)

			if temp != temp_list[0] :
	
				lattice.progress_n(n_eq_T, 0)
					
			
			#Progress the lattice n times 
			print("Progressing for correlation data")
			lines = lattice.progress_n(n, graphs=graphs)
		
			these_ens, these_mags = [], []
		
			for l in lines :
				splt = l.split(",")
				l_en, l_mag = float(splt[1].strip()), float(splt[2].strip())
				
				these_ens.append(l_en/N ** 2)
				these_mags.append(l_mag)
		
			#Time average the magnetisations with a time gap
			
			autocorrs, autocovs = [], []

			for tau in tau_list :
				
				T, autocov, autocov_0 = 0, 0, 0

				for t in range(n-int(tau)) :
					
					mean_mags = get_mean(these_mags[:n-int(tau)])
					
					T += 1
					
					autocov += (these_mags[t] - mean_mags)*(these_mags[t + int(tau)]-mean_mags)
				
				autocovs.append(autocov)
				autocorrs.append(autocov/autocovs[0])

			plt.figure(tau_plot.number)
			plt.plot(tau_list, autocorrs)
			leg.append(("T = %.2f" % temp))

			plt.figure(mag_plot.number)
			plt.plot(list(range(len(these_mags))), these_mags)
			m_leg.append(("T = %.2f" % temp))
		
		plt.figure(tau_plot.number)
		plt.xlabel("tau")
		plt.ylabel("Autocorrelation")		
		plt.title("Autocorrelation on tau, N = "+ str(N))
		plt.legend(leg)

		plt.figure(mag_plot.number)
		plt.title("Mags")
		plt.legend(m_leg)

	plt.show()