#This is the main file. See Readme.md for how to run the program.

#It is assumed that the working directory is the directory where main.jl is located and contains file Project.toml.

#Activate project environment and install all dependent packages.
#(same versions as the ones used in development)

#Comment this out during development (debugger is a lot faster this way)
import Pkg

#Can comment this out during debugging to save time
Pkg.activate(".")
Pkg.instantiate()

#Load necessary packages
using Optim, Parameters, QuantEcon, BenchmarkTools,JLD, Interpolations
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
    dataloaded = true
else
    #Initialisation with default values

    #Initial guess for value function 0: this leads to a corner solution in the first iteration (choose capital equal to the first point on the grid). This is usually fine but may not be, depending on the optimisation subroutine that we choose.
    #Ultimately, the initial guess does not matter (checked numerically)
    SE = fill(stat_equil(N_kh = par[1].N_kh,N_k = par[1].N_k,N_z = par[1].N_z),N_S);

    #Alternative - initial guess for value function is log of capital. This avoids the corner solution in the first iteration.
    #SE = fill(stat_equil(N_kh = par[1].N_kh,N_k = par[1].N_k,N_z = par[1].N_z,V=log.(fillV(par[1].N_k,par[1].N_z,par[1].k_gr))),N_S);

    TP = fill(stat_equil(N_kh = par[1].N_kh,N_k = par[1].N_k,N_z = par[1].N_z),par[1].T_max,N_S);
    dataloaded = false
end


#JIT compilation acceleration (first call the functions with tiny arguments so the compiler precompiles functions). Subsequent callSEs are fast.
partiny = pars(N_z = 2,N_k = 2,VFI_maxiter=1,SE_maxiter=1)
SEtiny = stat_equil(N_z=2,N_k=2,N_kh=2)
SE_compute!(partiny,SEtiny)

#************(1) Stationary equilibrium ***************
#Compute stationary equilibria for each set of parameters.
println("Computing Stationary Equilibria...\n")
for i=1:N_S
println("________________________________
Equilibrium $i out of $N_S")
    #SE[i] is initial guess in the first iteration. In the following iteration, the initial guess is the stationary equilibrium computed in the first call (SE[1]), unless we loaded the guess from file (which should be better)
    if(!dataloaded)
        SE[i] = SE_compute!(par[i],SE[i])
    else
        SE[i] = SE_compute!(par[i],SE[1])
    end

    print("Only computing the first stationary equilibrium during development./n")
    break
end

#************(2) Transition paths************************

#As initial guess for value function, we should use the value functions from stationary equilibria.

#************(3) Statistics, Plots*****************

#Generate a few illustrative plots

#(The historgram plot as a sanity check)
#range for plots (grid points indices)
a1 = 1
a2 = 50
zind = 5 #index of shock for which we plot the distribution
zval = par[1].shock_mc.state_values[zind]
p = plot(par[1].k_gr[a1:a2],SE[1].μ[a1:a2,1],title = "μ", label = "z = $zval",lw = 2)
display(p)

#Plot of the policy function and value function
a1 = 1
a2 = par[1].N_k
zind = 5 #index of shock for which we plot the distribution
zval = par[1].shock_mc.state_values[zind]
pV = plot(par[1].k_gr[a1:a2],SE[1].V[a1:a2,1],title = "V", label = "z = $zval",lw = 2)
display(pV)

#Also plot ξc:

#Saving results
#(function saveAll is defined in brexTools.jl, see there for details)
saveAll(foldername,SE,TP)

println("
To do:


- fix issue with Uc (premultiplication of everything leads to divergence of the value function if Uc != 1.0. Only some parts should be premultiplied, excluding the continuation value!)

- check stopping rule for policy function (so far we just iterate many times so convergence is guaranteed but it is unnecessarily expensive)


- h decreasing in A. This is counterintuitive. This also explains why h is lower for higher z (individual productivity) since this enters the firm's problem the same way as A. Need to investigate this futher.


- possible issue: check that the Tauchen approximation is correct (par[1].shock_mc.state_values[5]=2.0, not 1.0).

- baseline.jl: contains experimental value of A.

 - check boundary conditions at lower end of the grid (extrapolation penalty?, no depreciation, generous bound on labour supply)

")
