#!/usr/bin/env python

from __future__ import division #necessary to perform division of ints effectively cast into floats before division

import sys
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand

import IsingClass #module written to carry out the required calculations (metropolis, neighbours, energy, magnetization, etc.)


def main():
    steps = 2000 #steps for equilibrium
    temp_steps = 20 #step size for increase in temperature [J/k]

    kT = np.linspace(1.0,4.0,temp_steps) #temperature values

    #arrays for future plotting
    magnetization = np.zeros(temp_steps) 

    for l in range(len(kT)):
        temperature = kT[l]

        #initialise observables to zero
        M = 0

        #reset to random initial configuration
        lattice = model.spin_config()  

        #set out to reach equilibrium
        for t in range(steps):
            model.metropolis(lattice, temperature)

        #collect statistics
        stat_steps = int(steps/2) 
        for k in range(stat_steps):  

            model.metropolis(lattice, temperature) 

            mag = model.mag_per_spin(lattice)

            M += mag

        #observables, after equilibrium is achieved, normalised by sweeps
        magnetization[l] = abs(M)/stat_steps 

    plt.plot(kT, magnetization, '--o', label = "N = %i"%N)
    plt.legend()    


plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.title("Average spin site magentization versus temperature")
plt.xlabel("$T$ $[J/k_B]$")
plt.ylabel("$|<M>|$ $[\mu]$")
plt.grid()


J = 1.0 #coupling constant
h = 0.0 #external field

#vary lattice side size
N_vals = []
args = len(sys.argv)

for i in range(2, args,1):
    N_vals.append(int(sys.argv[i]))


for value in N_vals:
    N = value
    model = IsingClass.Ising(N,J,h) #constructor
    main()

#see if theory matches results (Onsager, 1944)
plt.axvline(2.269, color = "deeppink", label="$T_c^{theory} = 2.269$ $[J/k_B]$")
plt.legend()

#plot all curves in one graph
plt.show()



