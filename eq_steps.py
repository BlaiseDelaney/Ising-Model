#!/usr/bin/env python 
"""
Python script to analyse the dependence on relaxation time (=steps to reach equilibrium in ferromagnetic regime) 
on grid size, for 2D square lattice. Usage ./eq_steps.py [N values separated by space]

__author__ = Blaise Delaney, JSTP
"""

import numpy as np
import matplotlib as plt
import sys #argv 
import sys
import IsingClass #module to perform metropolis etc


def eq(temp,steps):
  for t in range(steps):
        model.metropolis(lattice, kT)
        if np.sum(lattice) in (-(N**2), N**2): #complete spin allignment: plus/minus NxN
            temp.append(t)
            break
  

#parameters for known behaviour
J = 1.0
h = 0.0
N = int(sys.argv[2])


steps = 3000 #max number of steps for equilibrium
kT = 0.3 #expect full ferrmomagnetic behaviour (complete spin allignment)

model = IsingClass.Ising(N,J,h) #constructor


t_steps = []

weight = 5


for w in range(weight):
  lattice = model.spin_config() #random spin configuration
  for t in range(steps):
        model.metropolis(lattice, kT)
        if np.sum(lattice) in (-(N**2), N**2): #complete spin allignment: plus/minus NxN
            t_steps.append(t)
            break

print "N = ", N, ".", " Steps to full spin alignment = ", np.sum(t_steps)/len(t_steps)
