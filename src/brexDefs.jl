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

     #Note: It is VERY important to use parametric type TF for the anonymous function,
     #Otherwise JIT compilation does not work properly (huge performance drop)
     #u::TF1 = σ == 1.0 ? (c,n) -> log(c) - χ*n : (c,n) -> (c^(1 - σ) - 1) / (1 - σ) - χ*n

     #We assume a log linear utility function as Bloom et al, BT and others (simplifies computation a lot as real wage is then directly related to MUc and we need to iterate over one fewer price)
     u::TF1 =  (c,n) -> log(c) - χ*n

     ############### Firms ############################
     #Production function parameters
     #Intermediate goods firms
     A::TF = 1.0 #Total Factor Productivity (WARNING: THIS IS currently OVERWRITTEN BY THE EXAMPLE baseline.jl parameter file)

     α::TF = 0.3
     ν::TF = 0.6
     #Production function - Cobb-Douglas with decreasing returns.
     #A is TFP which is a constant in the baseline but it could vary later (effect of Brexit on productivity)
     y::TF2 = (A,z,k,n) -> A*z*k^α*n^ν

     #depreciation rate
     δ::TF = 0.05

     ξbar = 1.0 #maximum possible realisation of the adjustment cost (ξ is U[0,ξbar])

     #Nmax is the maximum labour demanded by firms and supplied by households. So far the bound is put at a very high level (effectively unrestricted) to avoid non-differentiabilities, and we check ex post that it is not unreasonably high.
     Nmax::TF = 1000000.0

     #Production of final goods
     ϵ::TF = 0.5 #elasticity of substitution
     ω::TF = 0.5


     ############## Shocks ############################
     #Idiosyncratic productivity shock
     #Productivity shock process is either AR(1) approximated by Tauchen
     N_z::TI = 9 #number of shock realisations in approximation
     #(high N_z results in smooth distributions but increases computation time linearly)
     AR1_μ::TF = 1.0 #AR(1) mean
     AR1_ρ::TF = 0.5 #AR(1) autocorrelation
     AR1_σ::TF = 0.1 #std deviation of innovation
     AR1_stdev::TI = 2 #number of standard deviations to be covered by the grid (has to be an integer)

     #shock process is a Markov chain type (defined in QuantEcon)
     shock_mc::MarkovChain = tauchen(N_z,AR1_ρ,AR1_σ,AR1_μ,AR1_stdev)
     #shock_mc.p is the transition matrix, shock_mc.state_values are the actual realisations.

     #Aggregate shock
     tb::TI = 3 #period when Brexit happens (default = 3)
     PS::Array{TF,1} = [0.2,0.5,0.3] #Probability of each Brexit outcome (aggregate shock S realisation)
     #(no Brexit, deal/soft Brexit, no deal/hard Brexit)

     τ::TF = 0.0 #tarrif in baseline (no Brexit)
     τD::TF = 0.05 #tariff in "deal" Brexit
     τND::TF = 0.2 #tariff in "no deal" Brexit

     ########## What the program should do################
     simonly::Bool = false
     #If simonly == true then the program loads a previously saved solution
     #and simulates the economy (and generates plots).

     ############## Solution algorithm parameters ###############
     #Grid for capital in individual firm's problem
     N_k::TI = 100 #temporarily low in development
     k_min::TF = 0.01 #should be > 0 (for stability)
     k_max::TF = 20.0 #needs to be checked ex-post

     #Grid for capital of individual firm
     k_gr::StepRangeLen = range(k_min,k_max,length = N_k)
     #Warning: With cubic spline interpolation (see Vint_mode below) we need to use equispaced grids, with linear interpolation any grid is fine.

     #Number of capital grid points in the histogram of firm's distribution.
     #The same grid boundaries as in the individual firms' problem are used, (finer grid can be a good idea)
     N_kh::TI = N_k
     k_gr_hist::StepRangeLen = range(k_min,k_max,length = N_kh)

     #Vint_mode determines what interpolation type is used in value function interpolation. So far implemented values are 1 (linear interpolation) and 2 (cubic spline).
     #WARNING: With linear interpolation we maximise a nondifferentiable function and the optimal choice of capital always lies on the grid (because we maximise piecewise linear function minus k'). We need a VERY fine grid, otherwise we get a poor approximation
     #Cubic spline performs vastly better!)
     Vint_mode::TI = 2

     #Tmax is the number of periods after which we assume that the model reaches the new stationary distribution. It is the total number of periods, not the number of periods after Brexit happens.
     T_max::TI = 100

     #SE_maxiter is the maximum number of iterations in finding stationary equilibrium
     #This should be a fairly large number (maybe 1-5k), temporarily low during development
     SE_maxiter::TI = 10

     #VFI_maxiter is the maximum number of iterations in the VFI algorithm - updates of the value function using the updated policy function.
     #If Howard accelaration is used (VFI_howard = k>1), then maximisation is performed only every k-th iteration. This can speed things up quite a bit if the policy improvement step is relatively expensive.
     VFI_maxiter::TI = 1000 #just for development - a greater value should be set (and a stopping criterion used)
     VFI_howard::TI = 10 #Default value is 1, a value of around 20 should be reasonable.
 end

#mutable struct stat_equil contains everything that describes a stationary equilibrium: distribution of firms, prices, value and policy functions, etc.
#!!!Warning. If this is changed we also need to change the copy function for the struct defined below, otherwise copying it will result in errors.
@with_kw mutable struct stat_equil{TF<:AbstractFloat,TI<:Integer}
    #Number of grid points for capital and shock realisation
    #(these are only default values - will be overwritten if these values are passed to the constructor when SE objects are generated)
    N_k::TI=100 #policy and value function grid points
    N_kh::TI=100 #histogram grid points
    N_z::TI=9


    #Distribution of firms (first index corresponds to each value of capital, second index corresponds to shock realisation)
    μ::Array{TF,2} = ones(N_kh,N_z)/(N_kh*N_z)
    #default value is a naive guess of uniform distribution

    #Value function V(k,z) (beginning of period before ξ observed)
    V::Array{TF,2} = zeros(N_k,N_z)
    #Warning: zeros are a really bad initial guess because they will result in a corner solution in the first iteration - and with some unconstrained optimisation algorithms this could even result in a negative capital choice.
    #So some better value should be used when initialising the SE object. See function fillV below for details.

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
#!!!For matrices we need to use copy, otherwise the function just creates a pointer, and the data is not copied.
Base.copy(s::stat_equil) = stat_equil(N_kh = s.N_kh, N_z = s.N_z, μ = copy(s.μ), V = copy(s.V), h = copy(s.h), ξc = copy(s.ξc), Uc = s.Uc, Q = s.Q, pd = s.pd,w = s.w,N=copy(s.N))


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

#=function fill_V returns a value function, where the value is current capital, independent on z. This can then be used to construct some reasonable initial guess of the value function.

It should be the case that the value function is concave and increasing, with slope initially greater than 1 (log(k) is a reasonable)
(this guarantees that we will get an interior optimum in the first iteration, if the upper bound on capital is large enough)
=#
function fillV(N_k,N_z,k_gr)
    V = zeros(N_k,N_z)
    for i=1:N_k
        V[i,:] .= k_gr[i]
    end
    return V
end
