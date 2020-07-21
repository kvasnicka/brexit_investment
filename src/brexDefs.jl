#=
File BrexDefs.jl contains module brexDefs, which is the module defining data types used in the project
=#

module brexDefs
#This module file contains definitions of data structures used in the project
using Parameters #so we can use with_kw macro.
using QuantEcon #for Markov chains

export pars,stat_equil

#=
Struct pars contains all model parameters and their default values.
It also includes other objects which are pre-computed at the outset depending on the value of various parameters (for example grids, Markov chains, ...)

Example of use:
par = pars(β = 0.96) constructs a parameter struct with all values default
except for the discount factor.

Then we can access the values by par.β or par.σ directly. We can also unpack
the parameters into variables for direct access using for example
α,β = @unpack par

#Need to be careful with abstract types. For example TFL is abstract float. But Integer is NOT a subtype of AbstractFloat, so if we for example use k_min::TFL = 1 we get an error because we are trying to assign Int64 constant to a type which does not contain integers.
Later on - can replace AbstractFloat with something more general like Real (more robust but potentially performance could suffer due to JIT issues), or manually perform conversions.
 =#

@with_kw struct pars{TFL<:AbstractFloat,TINT<:Integer}

    ################ Households ##########################
    #Utility function (CRRA, if σ = 1 it is log)
    σ::TFL = 1.0
    u = σ == 1.0 ? x -> log(x) : x -> (x^(1 - σ) - 1) / (1 - σ)

    β::TFL = 0.95 #discount factor

    τ::TFL = 0.0 #tarrif

    ############### Firms ############################
    #Production function parameters
    A::TFL = 1.0 #Total Factor Productivity

    ############## Shocks ############################

    ############## What to do #######################
    simonly::Bool = false
    #If simonly == true then the program loads a previously saved solution
    #and simulates the economy (and generates plots).

    ############## Solution algorithm parameters ###############
    #Grid for capital in individual firm's problem
    N_k::TINT = 100
    k_min::TFL = 0.01 #should be > 0 (for stability)
    k_max::TFL = 20.0 #needs to be checked ex-post

    #Grid for capital of individual firm
    k_gr::StepRangeLen = range(k_min,k_max,length = N_k)

    #Number of capital grid points in the histogram of firm's distribution.
    #The same grid boundaries as in the individual firms' problem are used, (finer grid can be a good idea)
    N_kh::TINT = N_k
    k_gr_hist::StepRangeLen = range(k_min,k_max,length = N_kh)

    #Productivity shock process is either AR(1) approximated by Tauchen
    N_z::TINT = 9 #number of shock realisations in approximation
    AR1_μ::TFL = 1.0 #AR(1) mean
    AR1_ρ::TFL = 0.5 #AR(1) autocorrelation
    AR1_σ::TFL = 0.1 #std deviation of innovation
    AR1_stdev::TINT = 2 #number of standard deviations to be covered by the grid (has to be an integer)

    #shock process is a Markov chain type (defined in QuantEcon)
    shock_mc::MarkovChain = tauchen(N_z,AR1_ρ,AR1_σ,AR1_μ,AR1_stdev)
    #shock_mc.p is the transition matrix, shock_mc.state_values are the actual realisations.

    #Tmax is the number of periods after which we assume that the model reaches he new stationary distribution. [assumed to be the same for all shock realisations]
    T_max::TINT = 100

    #VFI_maxiter is the maximum number of iterations in VFI algorithm (solving the
    #individual firm's problem).
    VFI_maxiter::TINT = 500
    VFI_howard::Bool = true #if true Howard's acceleration algoritm will be used
    VFI_howard_c::TINT = 20 #Maximisation is performed in iterations
    #(1,1+VFI_howard_c,1+2*VFI_howard_c,...)

end

@with_kw mutable struct stat_equil{TFL<:AbstractFloat,TINT<:Integer}
    N_kh::TINT #Number of grid points for capital and shock realisation, no default, must be supplied
    N_z::TINT

    #Distribution of firms (first index corresponds to each value of capital, second index corresponds to shock realisation)
    μ::Array{TFL,2} = zeros(N_kh,N_z)

    #Value function

    #Policy function

    #Prices


end

#To be defined:
#mutable struct for transition paths. Will contain essentially the same objects as the mutable struct stat_equil, but with two extra indices (one for time period, and one for each realisation of the Brexit uncertainty shock [indexed 1,2,3])

end
