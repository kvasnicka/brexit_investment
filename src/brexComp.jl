#brexComp.jl contains functions for computing stationary equilibrium and transition paths.

#Function SE_compute! computes stationary equilibrium for given parameters par.
#SE is the initial guess and is overwritten.
#
function SE_compute!(par,SE)
#=Computation of stationary equilibrium proceeds in the following steps

(1) Given prices, solve the firm's problem
(2) Given solution of the firm's problem and prices, find a stationary distribution (simulation)
(3) Find the excess demands (check market-clearing conditions)
(4) Update prices given excess demands.

Steps (1) - (4) are repeated until convergence of prices or the maximum number of iterations is reached.

In principle, we just need an excess demand function (ED(prices) and could use some standard root finding algorithm, possibly using numerical derivatives. The above is standard for univariate (one price) problems in heterogenous agents models and exploits the structure of the problem. With three prices it is unclear whether it still makes sense to use this algorithm and update multiple prices in step (4), or whether it is better to use a standard root finding procedure. To begin with we follow the standard procedure but later on can switch to standard root finding to find prices such that ED(prices) = [0,0,0] approximately.

=#

#Iterate over prices.
for SEind = 1:par.SE_maxiter

#Given prices, find the excess demands.
#(The goal is to find prices such that all excess demands are close to zero)
#Later on we will use some algorithm to find a vector of prices such that the excess demands are approximately zero (market-clearing).
#In one dimension a bisection algorithm would be a possible approach, in multiple dimensions we will need to do something more sophisticated.

#Get excess demands and new stationary equilibria
excess,SEn = ED(par,SE)

#Stopping rule - to be implemented. At this stage of development break after the first loop.
println("Stopping rule in SE_compute not implemented yet, only one iteration performed.")
break

end

end #of SE_compute!


#Function ED computes the excess demands for given prices. SEg is the initial guess for a stationary equilibrium (contains value function, policy function, prices).
#Function ED returns a pair: excess, SEn
#where SEn is the new candidate for stationary equilibrium (prices, policy functions, etc.). This is better than simply overwriting the initial guess, especially when using a generic root finding algorithm (for some vector of prices off-equilibrium the policy functions and stationary distribution might be very far from the equilibrium one, and may be completely non-sensical so we might run into convergence problems).

#The prices for which we compute excess demands are part of initial guess SEg.
#Later on we can either write a wrapper function of prices only, or include prices as explicit argument here (if using a generic root-finding algorithm).
function ED(par,SEg)
    excess = [0.0,0.0,0.0] #Initialise excess demands
    #Compute the real wage (using households FOC and log linear utility) before copying new SE candidate (in any case, prices in SEn should not be used anywhere as they were not updated)
    SEg.w = par.χ/SEg.Uc


    SEn = copy(SEg) #Initialise the new stationary equilibrium
    #(this might be a bit wasteful since most of the values are overwritten later - maybe using an empty constructor is faster).

    #Solve the firm's problem given prices (these are contained in struct SEg along with the initial guess of value and policy functions), saving the new policy functions etc. in SEn (new stationary equilibrium candidate)
    firm_solve!(par,SEn,SEg)


    #Given the distribution of firms and policy functions, get excess demands.


#Return excess demands and the new candidate for stationary equilibrium
return excess,SEn

end


#Function firm_solve! solves the firm's problem given prices.
#SEn is the new stationary equilibrium (candidate) where firm's value function, policy function, etc. will be saved.
#SEg is the initial guess which shall not be changed by the function call in any way.
#The function uses Caretesian indexing implemented in Julia - see https://julialang.org/blog/2016/02/iteration/
function firm_solve!(par,SEn,SEg)

    #Generate iterator (this takes no space and is fast)
    V_ind = CartesianIndices((1:par.N_kh,1:par.N_k))

    #Value function iteration
    #All of the below (maximisation step and update will be put inside a loop, iterated until convergence).


    #parallel loop over grid points (commented out during development for debugging purposes)
    #@threads for i in eachindex(V_ind)
    for i in eachindex(V_ind)
        #V_ind[i] contains the Cartesian index for the i-th element of the matrix, V_ind[i][j] the index for the j-th state which corresponds to the grid point.

        #index of capital is V_ind[i][1]
        #Index of shock realisation is V_ind[i][2]

        #current capital and shock realisation
        k_ind = V_ind[i][1]
        z_ind = V_ind[i][2]

        k = par.k_gr[k_ind]

    end


end
