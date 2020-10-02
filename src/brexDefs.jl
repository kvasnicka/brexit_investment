#=
BrexDefs.jl contains data types used in the project, and some functions working on the types such as checking parameters
=#


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

 Note: Every anonymous function needs to have its own parametric type (TF,TF2,...). For example see the declaration of u::TF1 and y::TF2
 Otherwise the program crashes due to type mismatch when the parameters struct is instantiated.

 Also, every anonymous function needs to have a parametric type, otherwise
 performance suffers due to issues in JIT complilation (see Performance Tips
 in Julia documentation for details).
 =#

 @with_kw struct pars{TF<:AbstractFloat,TI<:Integer,TF1<:Function,TF2<:Function}
     ################ Households ##########################
     #Utility function (CRRA, if σ = 1 it is log)
     σ::TF = 1.0
     χ::TF = 1.0 #disutility of work (linear)
     β::TF = 0.95 #discount factor

     #Very important to use parametric type TF for the anonymous function,
     #Otherwise JIT compilation often does not work properly (huge performance drop)
#     u::TF1 = σ == 1.0 ? (c,n) -> log(c) - χ*n : (c,n) -> (c^(1 - σ) - 1) / (1 - σ) - χ*n
    #We assume a log linear utility function as Bloom et al, BT and others (simplifies computation a lot as real wage is then directly related to MUc and we need to iterate over one fewer price)
     u::TF1 =  (c,n) -> log(c) - χ*n

     τ::TF = 0.0 #tarrif

     ############### Firms ############################
     #Production function parameters
     #Intermediate goods firms
     A::TF = 1.0 #Total Factor Productivity
     α::TF = 0.3
     ν::TF = 0.6
     #Production function - Cobb-Douglas with decreasing returns.
     #A is TFP which is a constant so is not an input argument.
     y::TF2 = (z,k,n) -> A*z*k^α*n^ν

     #Nlim is the maximum labour supply (for some calibrations, we can have very high labour supply if there is almost no capital and productivity is low)
     Nmax::TF = 3.0

     #Production of final goods
     ϵ::TF = 0.5 #elasticity of substitution
     ω::TF = 0.5


     ############## Shocks ############################
     #Idiosyncratic productivity shock
     #Productivity shock process is either AR(1) approximated by Tauchen
     N_z::TI = 9 #number of shock realisations in approximation
     AR1_μ::TF = 1.0 #AR(1) mean
     AR1_ρ::TF = 0.5 #AR(1) autocorrelation
     AR1_σ::TF = 0.1 #std deviation of innovation
     AR1_stdev::TI = 2 #number of standard deviations to be covered by the grid (has to be an integer)

     #shock process is a Markov chain type (defined in QuantEcon)
     shock_mc::MarkovChain = tauchen(N_z,AR1_ρ,AR1_σ,AR1_μ,AR1_stdev)
     #shock_mc.p is the transition matrix, shock_mc.state_values are the actual realisations.

     #Aggregate shock
     tb::TI = 3 #period when Brexit happens (default = 3)
     PS::Array{TF,1} = [0.0,0.5,0.5] #Probability of each Brexit outcome (aggregate shock S realisation)

     ########## What the program should do################
     simonly::Bool = false
     #If simonly == true then the program loads a previously saved solution
     #and simulates the economy (and generates plots).

     ############## Solution algorithm parameters ###############
     #Grid for capital in individual firm's problem
     N_k::TI = 100
     k_min::TF = 0.01 #should be > 0 (for stability)
     k_max::TF = 20.0 #needs to be checked ex-post

     #Grid for capital of individual firm
     k_gr::StepRangeLen = range(k_min,k_max,length = N_k)
     #Warning: With cubic spline interpolation (see Vint_mode below) we need to use equispaced grids, with linear interpolation any grid is fine.

     #Number of capital grid points in the histogram of firm's distribution.
     #The same grid boundaries as in the individual firms' problem are used, (finer grid can be a good idea)
     N_kh::TI = N_k
     k_gr_hist::StepRangeLen = range(k_min,k_max,length = N_kh)

     #Vint_mode determines what interpolation type is used in value function interpolation. So far implemented values are 1 (linear interpolation) and 2 (cubic spline). 1 should be more robust, 2 could result in performance gains if stable.
     Vint_mode::TI = 1

     #Tmax is the number of periods after which we assume that the model reaches the new stationary distribution. It is the total number of periods, not the number of periods after Brexit happens.
     T_max::TI = 100

     #SE_maxiter is the maximum number of iterations in finding stationary equilibrium
     SE_maxiter::TI = 200

     #VFI_maxiter is the maximum number of iterations in the VFI algorithm - updates of the value function using the updated policy function.
     #If Howard accelaration is used (VFI_howard = k>1), then maximisation is performed only every k-th iteration. This can speed things up quite a bit if the maximisation is a relatively expensive step in the computation.
     VFI_maxiter::TI = 500
     VFI_howard::TI = 10 #Default value is 1, a value of around 20 should be reasonable.

 end

#mutable struct stat_equil contains everything that describes a stationary equilibrium: distribution of firms, prices, value and policy functions, etc.
#!!!Warning. If this is changed we also need to change the copy function for the struct defined below, otherwise copying it will result in errors.
@with_kw mutable struct stat_equil{TF<:AbstractFloat,TI<:Integer}
    #Number of grid points for capital and shock realisation
    N_k::TI=100 #policy and value function grid points
    N_kh::TI=100 #histogram grid points
    N_z::TI=9

    #Distribution of firms (first index corresponds to each value of capital, second index corresponds to shock realisation)
    μ::Array{TF,2} = zeros(N_kh,N_z)

    #Value function V(k,z) (beginning of period before ξ observed)
    V::Array{TF,2} = zeros(N_k,N_z)

    #Policy functions and value functions:

    #h is desired level of investment if there were no adjustment costs (does not depend on adjustment costs OR current capital stock).
    #This is because the adjustment costs paid on investment adjustment do not depend on the current level of capital (a simplificaiton common in the literature). (for possible future generalisations, just turn this into a matrix like ξc)
    #E is the associated value function used for computing cutoff cost
    h::Array{TF,1} = zeros(N_z)
    E::Array{TF,1} = zeros(N_z)

    #cutoff adjustment costs (it is optimal to invest if the realised adjustment cost is less than this).
    ξc::Array{TF,2} = zeros(N_k,N_z)

    #Labour demand policy function
    N::Array{TF,2} = zeros(N_k,N_z)

    #Prices
    Uc::TF = 1.0 #marginal utility of consumption
    Q::TF = 1.0 #Real exchange rate
    pd::TF = 1.0 #Relative price of domestic tradeable goods
    w::TF = 1.0 #real wage (only saved for convenience, it is actually a function of Uc).
end

#Adding method for copying structs
Base.copy(s::stat_equil) = stat_equil(N_kh = s.N_kh, N_z = s.N_z, μ = s.μ, V = s.V, h = s.h, ξc = s.ξc, Uc = s.Uc, Q = s.Q, pd = s.pd,w = s.w,N=s.N)

#This function performs checks of parameters
function check_par(par,N_S)
    if par.σ <= 0
        error("σ must be positive")
    end
    if par.σ != 1.0
        error("σ must be equal to 1.0 (only log utility implemented so far)")
    end
    if par.T_max <= par.tb
        error("T_max must be greater than t_brex")
    end
    if par.k_min > par.k_max
        error("k_min must be greater than k_max")
    end
    if length(par.PS) != N_S
        error("Parameter PS and parr_diff set in parameter file must have the same length")
    end
end
