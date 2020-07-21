#=
This is the main program file for the Brexit investment paper.

Use: 'julia --parfile parameterfile' runs the program with parameters contained in file './parameters/parameterfile'.

File parameters/example.jl is an example parameter file, used by the command 'julia main.jl --parfile example.jl', and is also used as a default if the program is called without specifying the parameter file name.

The role of each parameter is described in module BrexDefs (in BrexDefs.jl), in the definition of struct pars, together with their default values (these are used when the program is called without arguments).

=#

#Activate project environment and install all dependent packages.
#This should need to be done only once
import Pkg

#Activate and instantiate to make sure the same version of packages are used which worked during development (running without instantiate uses the currently installed version of packages which may be better but is also risky).
Pkg.activate(".")
Pkg.instantiate()

#Load necessary packages
using Optim, Parameters, QuantEcon
using .Threads #so we don't have to write Threads.@threads every time
const N_th = Threads.nthreads()

#Include files containing modules and load them
include("./src/brexDefs.jl")
using .brexDefs

#Run start-up script (see the file for details of what it saves)
#most importantly it loads parameters from file and collects these in array par
include("./src/startup.jl")


#####################################################
#=
The main body of the program follows. The algorithm proceeds in the following stages:

(0) Prep work. Generate objects used throughout such as objects for storing output.

(1) Compute stationary equilibria for the initial state (pre-referendum) and each of the long-run post-Brexit states.

(2) Compute transition paths.

(3) Compute various statistics, plot results.

=#

#(0) Prep work


#Create a copy of parameters for each of the realisations of Brexit (an arbitrary number of parameters can differ - in the baseline it is just tariffs τ). Store these in an array.

#(There is some unnecessary copying but it makes the code easier to generalise if we want to consider different values of more parameters than just tarrifs in different steady states. For this reason struct pars should be kept small - no large arrays etc.)

#In the baseline element with index 1 is the pre-Brexit state, 2 is soft Brexit, 3 is Hard Brexit.



#Construct stationary equilibria structs, which also allocated memory

#change - just loop over the elements of par[i]
#Change these: put this all in an array SE[i]

#SE_eu = stat_equil(N_kh = par.N_kh,N_z = par.N_z) #EU
#SE_soft = stat_equil(N_kh = par.N_kh,N_z = par.N_z) #soft Brexit
#SE_hard = stat_equil(N_kh = par.N_kh,N_z = par.N_z) #hard Brexit


#=
To do:
1) Array of parameters structs - generalise construction for a given number of states (using a for loop) and fix the issue with immutable structs.

2)
Transition paths
choose the storage - probably a composite type: array of stat_equil (indices for time period and aggregate shock realisation)
- then it should be called just equil or something else (it's a lot of objects together, not exactly aggregate state).

3) Fill in the remaining values in stat_equil and par struct (production function, etc.)

 =#
