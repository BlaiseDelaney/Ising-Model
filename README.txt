These scripts were written to investigate the Ising Model using the Metropolis Algorithm by Blaise Delaney, JS TP.

Use:

    - [MAIN SCRIPT] run executable Ising.sh to run the analysis on the n-dimensional (n = 1, 2, 3) isotropic lattice. The user will be prompted to input max and min temperatures, N, J and h. The plots of the absolute magnetisation per site, energy, susceptibility and specific hear capacity will be produced.

    - Binder.py: run this script to obtain the plot of the Binder cumulant U4 for set lattice side sizes.

    - Use magnetisation.sh to obtain the graph of the absolute magnetisation vs temperature for given values of J. These will be inputted manually by the user.

    - magnetisation_N.sh: same as magnetization.sh, but for chosen values of lattice side size N instead of J.

    - eq_steps.sh: the user will be prompted to input the dimensions of the lattice side. The script will return the number of iterations necessary for full spin allignement to be achieved at T = 0.3 J/k using the metropolis algorithm.

    - OneClass.py, IsingClass.py and CubeClass.py are modules written to carry out the analysis in 1D, 2D, and 3D respectively. Used by the main script Ising.sh.

    - Ising2D.ipynb is the IPython notebook where the original code for the Ising model was developed. Replicates the results of Ising.sh for the 2D case and allows for a visualisation of the lattice as a function of temperature gradient. This allows to view the transition from a ferromagnetic to the paramagnetic phase, J > 0.