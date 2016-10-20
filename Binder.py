#!/usr/bin/env python

from __future__ import division #necessary to perform division of ints effectively cast into floats before division

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand

import IsingClass #module written to carry out the required calculations (metropolis, neighbours, energy, magnetization, etc.)


def main():
    steps = 3000 #steps for equilibrium
    temp_steps = 20 #step size for increase in temperature [J/k]

    kT = np.linspace(1.0,4.0,temp_steps) #temperature values

    #arrays for future plotting
    magnetization = np.zeros(temp_steps) 
    mag_square = np.zeros(temp_steps) 
    mag_four = np.zeros(temp_steps)
    binder = np.zeros(temp_steps)

    for l in range(len(kT)):
        temperature = kT[l]
    
        #initialise observables to zero
        M = 0
        M_sq = 0
    
        M_fth = 0
    
        #reset to random initial configuration
        lattice = model.spin_config()  
    
        #set out to reach equilibrium
        for t in range(steps):
            model.metropolis(lattice, temperature)
        
        #collect statistics
        stat_steps = int(steps/2) 
        for k in range(stat_steps):  
        
            model.metropolis(lattice, temperature) 

            mag = model.mag_per_spin_b(lattice)
        
            M += mag
            M_sq += mag**2
            M_fth += mag**4
        
        #observables, after equilibrium is achieved, normalised by sweeps
        magnetization[l] = abs(M)/stat_steps 
        mag_square[l] = M_sq/stat_steps
        mag_four[l] = M_fth/stat_steps


    for l in range(len(kT)):
        binder[l] = 1 - mag_four[l]/(3*(mag_square[l]**2))
        
    plt.plot(kT, binder, label = "N = %d"%N)
    plt.legend(loc = "best")

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.title("Binder cumulant $U_4$ versus temperature")
plt.xlabel("$T$ $[J/k_B]$")
plt.ylabel("$U_4$")


plt.axvline(2.269, color ='r',linestyle='dashed', label="$T_c^{theory} = 2.269$ $[J/k_B]$")
plt.legend()
plt.grid()


N_vals = [8, 10, 25, 32, 50, 64] #side sizes
J = 1.0 #coupling in the hamiltonian
h = 0.0 #external field
   

for N in N_vals:   
    
    model = IsingClass.Ising(N,J,h) #constructor
    main()

plt.show()




