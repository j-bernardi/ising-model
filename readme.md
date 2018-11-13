# Ising Model

Current:

	main.py contains the lattice and spin-site objects, and all operations pertaining to updating the lattice.
	
	Enter individual experiments to manipulate the lattice in ways that are helpful for probing the physics.

Old iterations of the main script:

	simple.py
		revert so simple array storage
		iterate through rows and cols
	
		Keeps simple for coding
		Also hoping lack of vectorisation and messing about with interpreted stuff will speed up
	
		Finally should make forward calculation of dE simpler- try this after
	
		CONCERN : 
			Updates spins on the fly- doesn't hold lattice fixed while iterating

	simple_persist.py
		Keeps alternate energy values- some hit in computing these each time
		Better to keep a matrix which is only updated sparsely

	simple_persist_random.py
		Keeps a list of random numbers - only generated once. Shuffle after each iteration

		Large overhead on initialisation
		Significantly fastest! Okay N**2 >> number_of_steps
	forward_attempt_fail.py
		
		Fails as it updates every time and goes checkerboard... not sure why
	
	
	main.py
	
		Recalculates energy at each stage
		Also updates spins on-the-fly, e.g. 
	
	forward_main.py
		
		Stores an alt_energy grid and only recalculates when necessary
	
	
	old_main_t.py
		
		Updates spin as you go through the matrix- this seems like an issue as it updates half way through a time step and affects the neighbour choices
