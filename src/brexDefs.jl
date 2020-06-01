#=
File BrexDefs.jl contains module brexDefs, which is the module defining data
types used in the project
=#

module BrexDefs
#This module file contains definitions of data structures used in the project
using Parameters

export pars

#=
Struct pars contains all model parameters and their default values

Example of use:
par = pars(β = 0.96) constructs a parameter struct with all values default
except for the discount factor.

Then we can access the values by par.β or par.σ directly. We can also unpack
the parameters into variables for direct access using for example
α,β = @unpack par

 =#

@with_kw struct pars
    ################ Households ##########################
    #Utility function (CRRA, if σ = 1 it is log)
    σ = 1
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
    N_k::Int = 100
    k_min::Int = 0
    k_max::Int = 20

    #Number of capital grid points in the histogram of firm's distribution.
    #The same grid boundaries as in the individual firms' problem are used.
    N_kh::Int = 100

    #Productivity shock process


    #Tmax is the number of periods after which we assume that the model reaches
    #the new stationary distribution. An acceptable value needs to be found
    #experimentally (and depends on other parameters of the model)
    T_max::Int = 100


    #VFI_maxiter is the maximum number of iterations in VFI algorithm (solving the
    #individual firm's problem).
    VFI_maxiter::Int = 500
    VFI_howard::Bool = true #if true Howard's acceleration algoritm will be used
    VFI_howard_c::Int = 20 #Maximisation is performed in iterations
    #(1,1+VFI_howard_c,1+2*VFI_howard_c,...)

end



end
