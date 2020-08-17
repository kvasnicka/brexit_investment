#This is the main file. See Readme.md for how to run the program.

#It is assumed that the working directory is the directory where main.jl is located and contains file Project.toml.

#Activate project environment and install all dependent packages.
#(same versions as the ones used in development)
import Pkg
Pkg.activate(".")
Pkg.instantiate()
#Load necessary packages
using Optim, Parameters, QuantEcon, BenchmarkTools,JLD
using .Threads #so we don't have to write Threads.@threads every time

#Import the local module brexit_investment which contains all new functions, data types, etc. used in the project
import brexit_investment
using brexit_investment


#Run start-up script (see the file for details of what it does)
#Most importantly it loads parameters from file and collects these in array par.
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

#=
#Set up storage for stationary equilibria and transition paths

#SE - array of stationary equilibria (one for each element of par_diff).

#TP -transition paths
 This is a matrix of stationary equilibria struct. Calling it a stationary equilibrium is not accurate but it is essentially a path of distribution, value functions, etc. over time, where in each period we store the same objects that would be stored in a stationary equilibrium.

#It is assumed that:
1) the number of grid points etc. are the same for all parametrisations, values from par[1] are used.
2) the length of all transition paths is the same and equal to the value in par[1].T_max)
=#

if ((@isdefined loadData) && (@isdefined loadFolder) && loadData)
    #loading SE, TP from files
    println("Loading initial guess from folder results/$loadFolder.")
    SE = load("results/$loadFolder/SE.jld","SE")
    TP = load("results/$loadFolder/TP.jld","TP")
else
    #Initialisation with default values
    SE = fill(stat_equil(N_kh = par[1].N_kh,N_z = par[1].N_z),N_S);
    TP = fill(stat_equil(N_kh = par[1].N_kh,N_z = par[1].N_z),par[1].T_max,N_S);
end


#************(1) Stationary equilibrium ***************
#Compute stationary equilibria for each set of parameters.

println("Computing Stationary Equilibria...\n")
for i=1:N_S
println("Equilibrium $i out of $N_S")
    #SE[i] is initial guess
    SE[i] = SE_compute!(par,SE[i])
end

#************(2) Transition paths************************


#************(3) Statistics, Plots, Saving results*****************


#Saving results
#(function is defined in brexTools.jl, see there for details)
saveAll(foldername,SE,TP)
