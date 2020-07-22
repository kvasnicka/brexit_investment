#This is the main file. See Readme for how to run the program.

#Activate project environment and install all dependent packages.
#(same versions as the ones used in development)
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Revise
#use Revise for development only. Then we do not need to restart Julia to load updated packages and scripts, which saves a lot of time.
#Also need to use includet instead of include.
#For release versions either just comment this out and switch to standard includes.

#Load necessary packages
using Optim, Parameters, QuantEcon, BenchmarkTools
using .Threads #so we don't have to write Threads.@threads every time

#Load local modules from files (using includet for Revise)
includet("./src/brexDefs.jl")
using .brexDefs

#Run start-up script (see the file for details of what it does)
#Most importantly it loads parameters from file and collects these in array par
#Some outputs are N_S (number of parametrisations/equilibria), N_th (num of threads), par (array of parameters).
includet("./src/startup.jl")


#####################################################
#=
The main body of the program follows. The algorithm proceeds in the following stages:

(0) Prep work. Generate objects used throughout such as objects for storing output.

(1) Compute stationary equilibria for the initial state (pre-referendum) and each of the long-run post-Brexit states.

(2) Compute transition paths.

(3) Compute various statistics, plot results.

=#

#(0) Prep work

#Set up storage for stationary equilibria and transition paths

#(this is just an allocation of memory, it assumes that the number of grid points etc. are the same for all parametrisations, values from par[1] are used).

SE = fill(stat_equil(N_kh = par[1].N_kh,N_z = par[1].N_z),N_S);

#Transition paths (TP). This is a matrix of stationary equilibria struct. Calling it a stationary equilibrium is not accurate but it is essentially a path of distribution, value functions, etc. over time, where in each period we store the same objects that would be stored in a stationary equilibrium.

#(It is assumed that the length of all transition paths is the same and equal to the value in par[1].T_max)
TP = fill(stat_equil(N_kh = par[1].N_kh,N_z = par[1].N_z),par[1].T_max,N_S);


#=
To do:


- Run some benchmarks for the utility and other functions... Maybe it's easier in Jupyter

- Fill in more parameters in brexDefs (from the paper)

- then continue with the skeleton construction.
    - choose types for value function, policy function, etc.

- Bit for initialisation (get an initial guess of stationary equilibirum) [Later on - loading from file should be implemented].

- Bit for solving for stationary equilibrium.

- Then implement bit for shooting algorithm.

- Analysis of results, saving data, etc.

=#
