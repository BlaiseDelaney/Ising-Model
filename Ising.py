#!/usr/bin/env python

'''
__author__ = Blaise Delaney, JS TP

To be used running the executable Ising.sh. For the given dimensions of
the isotropic lattice the expectation values of physical observables are 
plotted against temperature. The code imports a specific module depending 
on the dimensions of the system. 
'''

from __future__ import division #necessary to perform division of ints effectively cast into floats before division

import sys
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand

import timeit

start = timeit.default_timer()


#select module based on dimensions (isotropic lattice)
dimension = str(sys.argv[1])
tmin = float(sys.argv[2])
tmax = float(sys.argv[3])
N = int(sys.argv[4]) #side size
J = float(sys.argv[5]) #safety, correct calculations (subject to machine rounding error)
h = float(sys.argv[6]) #safety, correct calculations (subject to machine rounding error)


if dimension == "1D":
    import OneClass
    model = OneClass.Ising(N,J,h) #constructor

if dimension == "2D":
    import IsingClass
    model = IsingClass.Ising(N,J,h) #constructor

if dimension == "3D":
    import CubeClass
    model = CubeClass.Ising(N,J,h) #constructor



#main 
steps = 2000 #steps for equilibrium
temp_steps = 20 #step size for increase in temperature [J/k]

kT = np.linspace(tmin,tmax,temp_steps) #temperature values


#arrays for future plotting
magnetization = np.zeros(temp_steps) 
Energy = np.zeros(temp_steps)
susceptibility = np.zeros(temp_steps)
specific_heat = np.zeros(temp_steps)
mag_square = np.zeros(temp_steps) 
e_square = np.zeros(temp_steps) 

for l in range(len(kT)):
    temperature = kT[l]

    #initialise temp observables to zero
    M = 0
    E = 0
    M_sq = 0
    E_sq = 0

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
        mag2 = model.mag2_per_spin(lattice)
        engy = model.energy(lattice)
        engy2 = model.energy2(lattice)

        M += mag
        E += engy
        M_sq += mag2
        E_sq += engy2

    #observables, at desired temperature, normalised by sweeps
    magnetization[l] = abs(M)/stat_steps 
    Energy[l] = E/stat_steps 
    mag_square[l] = M_sq/stat_steps
    e_square[l] = E_sq/stat_steps

#perform manipulation directly on array instead of loop
susceptibility = (mag_square - magnetization**2)/kT
specific_heat = (e_square - Energy**2)/(kT**2)


#plot observables
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
f = plt.figure(figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k'); 


plt.subplot(2, 2, 1 );
plt.grid()
plt.plot(kT, magnetization, 'k-s')
plt.xlabel("$T$ $[J/k_B]$")
plt.ylabel("$|<M>|$ $[\mu]$")

plt.subplot(2, 2, 2 );
plt.grid()
plt.plot(kT, Energy, 'k-s')
plt.xlabel("$T$ $[J/k_B]$")
plt.ylabel("$U$ $[J]$")

plt.subplot(2, 2, 3 );
plt.grid()
plt.plot(kT, susceptibility, 'k-s')
plt.xlabel("$T$ $[J/k_B]$")
plt.ylabel("$\chi$ $[\mu/ k_B]$")

plt.subplot(2, 2, 4 );
plt.grid()
plt.plot(kT, specific_heat, 'k-s')
plt.xlabel("$T$ $[J/k_B$]")
plt.ylabel("$C_V$ $[J/k_B^2]$")


#run time
stop = timeit.default_timer()
print "Run time = ", stop - start, "seconds." #http://stackoverflow.com/questions/5622976/how-do-you-calculate-program-run-time-in-python
plt.show()

