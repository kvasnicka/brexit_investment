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
    #Compute the optimal labour supply (this is a static problem, depends only on prices and not on the continuation value).
    compute_N!(par,SEn)

    #Value function iteration
    for VFIind = 1:par.VFI_maxiter

        #If this is the first iteration, or if the iteration index is a multiple of VFI_howard, update the policy function.
        if (mod(VFIind,par.VFI_howard) == 0 || VFIind == 1)
            #Sen contains the current policy function (and value function, prices, etc.) and the policy function part will be overwritten.
            #poll_diff is the stopping criterion for policy function.
            pol_diff = update_pol!(par,SEn)
        end

        #Update the value function (again, SEn contains the updated policy function, etc.).
        #V_diff is the stopping criterion value for value function
        V_diff = update_V!(par,SEn)

        #Check stopping criteria (for policy function or value function)
        #(to do - so far no stopping criterion, maximum number of iterations always performed)

    end

end #firm_solve!


#Function compute_N! computes the labour supply for each point in the state space (and saves it in SE.N)
function compute_N!(par,SE)
    V_ind = CartesianIndices((1:par.N_k,1:par.N_z))
    #debug -no threads
    @threads for i in eachindex(V_ind)
        k = par.k_gr[V_ind[i][1]]
        z = par.shock_mc.state_values[V_ind[i][2]]
        SE.N[i] = min((SE.w/(par.A*z*k^par.α*par.ν))^(1/(1-par.ν)),par.Nmax)
    end
end

#Function update_pol! updates the policy function (and overwrites the previous one). It returns value corresponding to the stopping criterion.
function update_pol!(par,SE)
    #Generate iterator (this takes almost no memory and is fast)
    V_ind = CartesianIndices((1:par.N_k,1:par.N_z))

    #First compute the unrestricted optimal capital level. This does not depend on the current level of capital, or adjustment costs, so we only loop over the current productivity.

    #Construct interpolator for the ex ante value function
    Vint = get_Vint(par,SE.V)

    #Objective function - expectation of the ex ante value function.
    #Need a function which takes, as an input, current state (z_ind), transition matrix, some value of capital kpr (next-period capital)


    #The number of productivity shock realisations tends to be quite small, so there might not even be a performance gain from multi-threading.
    for z_ind=1:par.N_z
        #Compute optimal level of capital in the absence of adjustment costs
        #This will be saved in policy function h, value saved in E

        #First write the objective function (so I can get this value of any level of k'. Then run a maximisation routine.)

        #Implement maximisation.

        #As a precaution only overwrite the previous policy if the return is strictly greater than the previous return.

    end

    #parallel loop over grid points
    @threads for i in eachindex(V_ind)
        #V_ind[i] contains the Cartesian index for the i-th element of the matrix, V_ind[i][j] the index for the j-th state which corresponds to the grid point.

        #index of capital is V_ind[i][1]
        #Index of shock realisation is V_ind[i][2]

        #current capital and shock realisation
        k_ind = V_ind[i][1]
        z_ind = V_ind[i][2]
        k = par.k_gr[k_ind]

        #Compute the optimal adjustment cost threshold.

    end


    #placeholder for stopping criterion
    return 0.0
end

#Function get_Vint constructs an interpolator for the ex ante value function.
#It returns object Vint, where Vint[z_ind](k) gives the value for idiosyncratic productivity with index z_ind and value of capital k.
#(The object is an array of interpolators, one for each possible shock realisation)
#Inputs are the parameter struct (which contains capital grid and number of grid points) and the value function on the grid.
function get_Vint(par,V)
    if par.Vint_mode == 1 #linear
        Vint = fill(LinearInterpolation(par.k_gr,V[1:par.N_k,1]),par.N_z)
        for zind = 2:par.N_z
            Vint[zind] = LinearInterpolation(par.k_gr,V[1:par.N_k,zind])
        end
    elseif par.Vint_mode == 2 #cubic spline
        Vint = fill(CubicSplineInterpolation(par.k_gr,V[1:par.N_k,1]),par.N_z)
        for zind = 2:par.N_z
            Vint[zind] = CubicSplineInterpolation(par.k_gr,V[1:par.N_k,zind])
        end
    else
        error("Unsupported value of Vint_mode given.")
    end

    return Vint
end

#Expected ex ante value function. kpr is next-period capital level, z_ind the index of current shock realisation,P is the transition matrix, Vint is the value function interpolant obtained using function get_Vint
function EV(kpr,z_ind,P,N_z,Vint)
    EV = 0.0
    for zpr = 1:N_z #cycle over next-period shock realisations, P[z_ind,zpr] is the transition probability
        EV += P[z_ind,zpr]*Vint[zpr](kpr)
        print(P[z_ind,zpr])
    end
    return EV
end

#Function update_V! updates the value function (and overwrites the previous one). It returns value corresponding to the stopping criterion.
function update_V!(par,SE)

    #placeholder for stopping criterion
    return 0.0
end
