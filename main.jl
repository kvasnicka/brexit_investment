#This is the main file. See Readme for how to run the program.

#Activate project environment and install all dependent packages.
#(same versions as the ones used in development)
import Pkg
Pkg.activate(".")
Pkg.instantiate()

#Load necessary packages
using Optim, Parameters, QuantEcon, BenchmarkTools
using .Threads #so we don't have to write Threads.@threads every time

#Load local modules from files (using includet if Revise is loaded)
#This is convenient during development - any changes in the modules will be reflected in the program without having to restart Julia.
#Revise should be configured in the julia startup file (see Pkg Revise documentation for details.) If it's not loaded we use standard include.
if isdefined(Revise,:includet)
    includet("./src/brexDefs.jl")
    includet("./src/brexMod.jl")
else
    include("./src/brexDefs.jl")
    include("./src/brexMod.jl")
end
using .brexDefs
using .brexMod

#Run start-up script (see the file for details of what it does)
#Most importantly it loads parameters from file and collects these in array par
#Some outputs are N_S (number of parametrisations/equilibria), N_th (num of threads), par (array of parameters).
include("./src/startup.jl")

#=
The main body of the program follows. The algorithm proceeds in the following stages:

(0) Prep work. Generate objects used throughout such as objects for storing output.

(1) Compute stationary equilibria for the initial state (pre-referendum) and each of the long-run post-Brexit states.

(2) Compute transition paths.

(3) Compute various statistics, plot results, save results
=#



#*************(0) Prep work**********************

#Set up storage for stationary equilibria and transition paths

#(It assumes that the number of grid points etc. are the same for all parametrisations, values from par[1] are used.).

SE = fill(stat_equil(N_kh = par[1].N_kh,N_z = par[1].N_z),N_S);

#Transition paths (TP). This is a matrix of stationary equilibria struct. Calling it a stationary equilibrium is not accurate but it is essentially a path of distribution, value functions, etc. over time, where in each period we store the same objects that would be stored in a stationary equilibrium.

#(It is assumed that the length of all transition paths is the same and equal to the value in par[1].T_max)
TP = fill(stat_equil(N_kh = par[1].N_kh,N_z = par[1].N_z),par[1].T_max,N_S);

#************(1) Stationary equilibrium ***************


#************(2) Transition paths************************



########Saving results#################
#Implement the following:
#Unless save_results = false
#Upon running the program, a folder is generated in subfolder results.
#Name of the folder is yymmdd_hhmmss_parfile and it contains (to begin with)
#A copy of the parameter file.
#A copy of the results (SE to begin with)
#
#Change .gitignore to make sure these are not uploaded to github.
