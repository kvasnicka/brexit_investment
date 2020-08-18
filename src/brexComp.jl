#brexComp.jl contains functions for computing stationary equilibrium and transition paths.

#Function SE_compute! computes stationary equilibrium for given parameters par.
#SE is the initial guess and is overwritten.
#
function SE_compute!(par,SE)
#=Computation of stationary equilibrium proceeds in the following steps

(1) Given prices, solve the firm's problem
(2) Given solution of the firm's problem and prices, find a stationary distirbution (simulation)
(3) Find the excess demands (check market-clearing conditions)
(4) Update prices given excess demands.

Steps (1) - (4) are repeated until convergence of prices or the maximum number of iterations is reached.

In principle, we we just need an excess demand function (ED(prices) and could use some standard root finding algorithm, possibly using numerical derivatives. The above is standard for univariate (one price) problems in heterogenous agents models and exploits the structure of the problem. With three prices it is unclear whether it still makes sense to use this algorithm and update multiple prices in step (4), or whether it is better to use a standard root finding procedure. To begin with we follow the standard procedure but later on can switch to standard root finding to find prices such that ED(prices) = [0,0,0] approximately.

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
break

end

end #of SE_compute!


#Function ED computes the excess demands for given prices. SEg is the initial guess for a stationary equilibrium (contains value function, policy function, prices).
#Function ED returns a pair: excess, SEn
#where SEn is the new candidate for stationary equilibrium (prices, policy functions, etc.). This is better than simply overwriting the initial guess, especially when using a generic root finding algorithm (for some vector of prices off-equilibrium the policy functions and stationary distribution might be very far from the equilibrium one, and may be completely non-sensical so we might run into convergence problems).

#The prices for which we compute excess demands are part of initial guess SEg.
#Later on we can either write a wrapper function of prices only, or include prices as explicit argument here (if using a generic root-finding algorithm).

#!This is actually quite expensive.
function ED(par,SEg)
    excess = [0.0,0.0,0.0] #excess demands placeholder
    SEn = copy(SEg) #Initialise the new stationary equilibrium

#Solve the firms' problem given prices


#Return excess demands and the new candidate for stationary equilibrium
return excess,SEn

end


#Function firm_solve! solves the firm's problem given prices.
