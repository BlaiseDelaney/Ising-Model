from __future__ import division #necessary to perform division of ints effectively cast into floats before division

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand
from random import randrange
#Ising model basic class
class Ising(object):

    def __init__(self,N,J,h):
        self.Nx = N
        self.J = J
        self.h = h

    def spin_config(self):
        lattice = rand.choice([1,-1],[self.Nx]) 
        return lattice
    
    #plot the spin configuration in the lattice 
    def lattice_plt(self, lattice):
        plt.matshow(lattice, cmap=plt.cm.gray)
        #LaTeX
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        plt.title("Spin configuration in the lattice")
        plt.axis('off')
        plt.show()
    
    def nb_sum(self,lattice,i):
        return lattice[(i+1) % self.Nx] + lattice[(i-1) % self.Nx]
    
    #delta_energy with periodic BC imposed using the modulo operator, per spin site
    def dEnergy(self, lattice,i): 
        nbs = self.nb_sum(lattice, i)
        return 2*(self.J*lattice[i]*nbs + self.h*lattice[i])#include external field
    
    #metropolis with periodic BC
    def metropolis(self, lattice, kT):
        for i in range(self.Nx):
            delta_E = self.dEnergy(lattice,i)
            site = lattice[i]
            if delta_E<0:
                site *= -1 #energy is lowered, accept
            elif np.exp(-delta_E/kT) > rand.rand(): #COMMMENT ON RAND(); uniform distribution in [0,1)
                site *= -1 #accept according to Boltzman dist.
            lattice[i] = site
        return lattice
    
   
    
    #energy per site
    def energy(self, lattice):
        erg = 0
        for i in range(self.Nx):
            Snb = self.nb_sum(lattice,i)
            erg += -self.J*lattice[i]*Snb - self.h*lattice[i]
        norm_energy = erg/(2*(self.Nx)) #divide by twice number of spins to avoid overcounting
        return norm_energy
    
    #energy squared per site                       
    def energy2(self,lattice):
        erg2 = 0
        for i in range(self.Nx):
            #sum of neighbouring spins
            sum_nb2 = self.nb_sum(lattice,i)
            #energy per site (not normalised)
            erg2 += (-self.J*lattice[i]*sum_nb2 - self.h*lattice[i])**2
        norm_energy2 = erg2/(4*self.Nx) 
        return norm_energy2
    
    #magnetization per spin
    def mag_per_spin(self,lattice):
        Mag = 0
        for i in range(self.Nx):
            Mag += lattice[i]
        return Mag/float(self.Nx)
                       
    #magnetization squared per spin
    def mag2_per_spin(self,lattice):
        Mag2 = 0
        for i in range(self.Nx):
            Mag2 += lattice[i]**2
        return Mag2/float(self.Nx)
    
    
    
    
    
        
 

   




