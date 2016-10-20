from __future__ import division #necessary to perform division of ints effectively cast into floats before division

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand
from random import randrange
#Ising model basic class
class Ising(object):
    
    #NOTE: decouple the N's. self.Nx and self.Ny to be specified explicitly to allow for application to 1D by setting Ny=1 or Nx=1

    def __init__(self,N,J,h):
        self.Nx = N
        self.Ny = N
        self.J = J
        self.h = h
    
    def spin_config(self):
        lattice = rand.choice([1,-1],[self.Nx,self.Ny]) 
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
    
    def nb_sum(self,lattice,i,j):
        return lattice[(i+1) % self.Nx][j] + lattice[(i-1) % self.Nx][j] + lattice[i][(j+1) % self.Ny] +  lattice[i][(j-1) % self.Ny]
    
    #delta_energy with periodic BC imposed using the modulo operator, per spin site
    def dEnergy(self, lattice,i,j): 
        nbs = self.nb_sum(lattice, i, j)
        return 2*(self.J*lattice[i][j]*nbs + self.h*lattice[i][j])#include external field
    
    #metropolis with periodic BC
    def metropolis(self, lattice, kT):
        for i in range(self.Nx):
            for j in range(self.Ny):
                delta_E = self.dEnergy(lattice,i ,j)
                site = lattice[i][j]
                if delta_E<0:
                    site *= -1 #energy is lowered, accept
                elif np.exp(-delta_E/kT) > rand.rand(): #COMMMENT ON RAND(); uniform distribution in [0,1)
                    site *= -1 #accept according to Boltzman dist.
                lattice[i][j] = site
        return lattice
    
    #this routine samples the lattice randomly instead of spanning it in an orderly fashion.
    #When comparing very large lattice sizes, for a sufficient large amount of Monte Carlo chain steps, 
    #the explicit dependence of self.metropolis on self.NxN can be minimised,
    #in order to optimise the overall Ising Model runtime.
    def metropolis_randsample(self, lattice, kT):
        i = randrange(self.Nx)
        j = randrange(self.Ny)
        delta_E = self.dEnergy(lattice,i ,j)
        site = lattice[i][j]
        if delta_E<0:
            site *= -1 #energy is lowered, accept
        elif np.exp(-delta_E/kT) > rand.rand(): #COMMMENT ON RAND(); uniform distribution in [0,1)
            site *= -1 #accept according to Boltzman dist.
        lattice[i][j] = site
        return lattice
    
    #energy per site
    def energy(self, lattice):
        erg = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                Snb = self.nb_sum(lattice,i,j)
                erg += -self.J*lattice[i][j]*Snb - self.h*lattice[i][j]
        norm_energy = erg/(2*(self.Nx*self.Ny)) #divide by twice number of spins to avoid overcounting
        return norm_energy
    
    #energy squared per site                       
    def energy2(self, lattice):
        erg2 = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                #sum of neighbouring spins
                sum_nb2 = self.nb_sum(lattice,i,j)
                #energy per site (not normalised)
                erg2 += (-self.J*lattice[i][j]*sum_nb2 - self.h*lattice[i][j])**2
        norm_energy2 = erg2/(4*self.Nx*self.Ny) 
        return norm_energy2
    
    #magnetization per spin
    def mag_per_spin(self,lattice):
        Mag = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                Mag += lattice[i][j]
        return Mag/float(self.Nx*self.Ny)
                       
    #magnetization squared per spin
    def mag2_per_spin(self,lattice):
        Mag2 = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                Mag2 += lattice[i][j]**2
        return Mag2/float(self.Nx*self.Ny)
    
    #magnetization per spin - BINDER patch
    def mag_per_spin_b(self,lattice):
        Mag = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                Mag += lattice[i][j]
        return Mag
                       
    #magnetization squared per spin - BINDER patch
    def mag2_per_spin_b(self,lattice):
        Mag2 = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                Mag2 += lattice[i][j]**2
        return Mag2
    
    #quantity needed for the Binder cumulant 
    def mag4_per_spin(self,lattice):
        Mag4 = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                Mag4 += lattice[i][j]**4
        return Mag4
    
    #look at ground state spin energy with squared (1D line) lattice
    def square_energy_ground_state(self):
        #REFERENCE: http://physics.stackexchange.com/questions/133005/how-to-calculate-the-ground-state-energy-for-the-ising-model
        
        #square lattice
        lattice_groud_state = np.ones((self.Nx,self.Ny), dtype=np.int) #int, save memory - ground state is given by all spins collinear.
        #given the symmetry, simply choose all spins up

        lattice_nb_sum = []

        for i in range(N):
            for j in range(N):
                nbs_sum = lattice_groud_state[i][j]*self.nb_sum(lattice_ground_state,i,j,self.Nx, self.Ny) #projection of spin site onto sum of all edges (topological equivalence of graph and its dual in contribution to energy, Statistical Physics II)
                lattice_nb_sum.append(nbs_sum) #system ground state energy x 2

        ground_energy = (np.sum(lattice_nb_sum)*(-self.J))/(2.0*N*N) #spins can be up or down; energy per site

        E_ground = np.int(ground_energy)
        print "Ground state energy = %dNself.J, for N the total number of spin sites"%E_ground 
        
 

   



