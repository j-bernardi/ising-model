A folder to hold experiment 1:
	H = 0 
	Vary temperature
	Run to equillibrium - start random for high T (5), run for temps down to 0.1 in 0.1 steps, use previous final state as next starting state


exp 1.1
	Outputs: 
		Heatmap of the evolution
		Graph of the evolution of energy and spin (1 axis)
		For each temperature...
		Finds how many runs to reach EQ

exp1.2.py
	Outputs: run of net energy and  mag per temperature. 
	Runs for many n, smooths out points and plots with error bars by taking last n-j points

exp1.2_N_list.py
	Plots multiple on the graph for varying N, varying temperature with each

	Should allude to reasonable n and Ns

	Good plots come out - good for finding a rough estimate of T_C 

	Plotted an N = 200 one 
