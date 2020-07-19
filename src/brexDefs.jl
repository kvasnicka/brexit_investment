#=
File BrexDefs.jl contains module brexDefs, which is the module defining data
types used in the project
=#

module brexDefs
#This module file contains definitions of data structures used in the project
using Parameters #so we can use with_kw macro.
using QuantEcon #for Markov chains

export pars

#=
Struct pars contains all model parameters and their default values.
It also includes other objects which are pre-computed at the outset depending on the value of various parameters (for example grids, Markov chains, ...)

Example of use:
par = pars(β = 0.96) constructs a parameter struct with all values default
except for the discount factor.

Then we can access the values by par.β or par.σ directly. We can also unpack
the parameters into variables for direct access using for example
α,β = @unpack par


!!! Declare all fields as a concrete type (for JIT performance)
For example Int64, not Integer.
 =#

@with_kw struct pars
    ################ Households ##########################
    #Utility function (CRRA, if σ = 1 it is log)
    σ::Float64 = 1
    u = σ == 1 ? x -> log(x) : x -> (x^(1 - σ) - 1) / (1 - σ)

    #discount factor
    β::Float64 = 0.95

    ############### Firms ############################
    #Production function parameters


    ############## Shocks ############################

    ############## What to do #######################
    simonly::Bool = false
    #If simonly == true then the program loads a previously saved solution
    #and simulates the economy (and generates plots).

    ############## Solution algorithm parameters ###############
    #Grid for capital in individual firm's problem
    N_k::Int32 = 100
    k_min::Float64 = 0.01 #should be > 0 (for stability)
    k_max::Float64 = 20 #needs to be checked ex-post

    #Grid for capital of individual firm
    k_gr::StepRangeLen = range(k_min,k_max,length = N_k)


    #Number of capital grid points in the histogram of firm's distribution.
    #The same grid boundaries as in the individual firms' problem are used, (finer grid can be a good idea)
    N_kh::Int32 = N_k
    k_gr_hist::StepRangeLen = range(k_min,k_max,length = N_kh)

    #Productivity shock process is either AR(1) approximated by Tauchen
    N_z::Int32 = 9 #number of shock realisations in approximation
    AR1_μ::Float64 = 1.0 #AR(1) mean
    AR1_ρ::Float64 = 0.5 #AR(1) autocorrelation
    AR1_σ::Float64 = 0.1 #std deviation of innovation
    AR1_stdev::Float64 = 2 #number of standard deviations to be covered by the grid

    #shock process is a Markov chain type (defined in QuantEcon)
    shock_mc::MarkovChain = tauchen(N_z,AR1_ρ,AR1_σ,AR1_μ,AR1_stdev)
    #shock_mc.p is the transition matrix, shock_mc.state_values are the actual realisations.

    #Tmax is the number of periods after which we assume that the model reaches
    #the new stationary distribution. An acceptable value needs to be found
    #experimentally (and depends on other parameters of the model)
    T_max::Int32 = 100

    #VFI_maxiter is the maximum number of iterations in VFI algorithm (solving the
    #individual firm's problem).
    VFI_maxiter::Int64 = 500
    VFI_howard::Bool = true #if true Howard's acceleration algoritm will be used
    VFI_howard_c::Int64 = 20 #Maximisation is performed in iterations
    #(1,1+VFI_howard_c,1+2*VFI_howard_c,...)

end


end
