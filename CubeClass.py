from __future__ import division #necessary to perform division of ints effectively cast into floats before division

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand

#Ising model basic class
class Ising(object):

    def __init__(self,N,J,h):
        self.Nx = N
        self.Ny = N
        self.Nz = N
        self.J = J
        self.h = h

    def spin_config(self):
        lattice = rand.choice([1,-1],[self.Nx,self.Ny,self.Nz]) 
        return lattice
    
    def nb_sum(self,lattice,i,j,k):
        return lattice[(i+1) % self.Nx][j][k] + lattice[(i-1) % self.Nx][j][k] + lattice[i][(j+1) % self.Ny][k] +  lattice[i][(j-1) % self.Ny][k] + lattice[i][j][(k+1) % self.Nz] + lattice[i][j][(k-1) % self.Nz]
    
    #delta_energy with periodic BC imposed using the modulo operator, per spin site
    def dEnergy(self, lattice,i,j,k): 
        nbs = self.nb_sum(lattice,i,j,k)
        return 2*(self.J*lattice[i][j][k]*nbs + self.h*lattice[i][j][k])#include external field
    
    #metropolis with periodic BC
    def metropolis(self, lattice, kT):
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    delta_E = self.dEnergy(lattice,i,j,k)
                    site = lattice[i][j][k]
                    if delta_E<0:
                        site *= -1 #energy is lowered, accept
                    elif np.exp(-delta_E/kT) > rand.rand(): #COMMMENT ON RAND(); uniform distribution in [0,1)
                        site *= -1 #accept according to Boltzman dist.
                    lattice[i][j][k] = site
        return lattice
    
    #energy per site
    def energy(self, lattice):
        erg = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    Snb = self.nb_sum(lattice,i,j,k)
                    erg += -self.J*lattice[i][j][k]*Snb - self.h*lattice[i][j][k] 
            norm_energy = erg/(2*(self.Nx*self.Ny*self.Nz)) #divide by twice number of spins to avoid overcounting
        return norm_energy
    
    #energy squared per site                       
    def energy2(self, lattice):
        erg2 = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    #sum of neighbouring spins
                    sum_nb2 = self.nb_sum(lattice,i,j,k)
                    #energy per site (not normalised)
                    erg2 += (-self.J*lattice[i][j][k]*sum_nb2 - self.h*lattice[i][j][k] )**2
        norm_energy2 = erg2/(4*self.Nx*self.Ny*self.Nz) 
        return norm_energy2
    
    #magnetization per spin
    def mag_per_spin(self,lattice):
        Mag = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    Mag += lattice[i][j][k]
        return Mag/float(self.Nx*self.Ny*self.Nz)
                       
    #magnetization squared per spin
    def mag2_per_spin(self,lattice):
        Mag2 = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    Mag2 += lattice[i][j][k]**2
        return Mag2/float(self.Nx*self.Ny*self.Nz)
    
   
