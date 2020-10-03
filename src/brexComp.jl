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
        #V_diff is the stopping criterion value for the value function
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

    #Not using @threads in development. Try it later to see what performance difference it makes.
    for z_ind=1:par.N_z
        #Compute optimal level of capital in the absence of adjustment costs
        #This will be saved in policy function h, value saved in E

        SE.h[z_ind],SE.E[z_ind] = find_KU(z_ind,par.shock_mc.p,par.N_z,Vint,par.k_min,par.k_max)

        #Warning! - the continuation value E does not include the left-over capital (1-δ)*k. So it is not exactly Vadj as in our model notation, but more like E in KT2008 paper. In the numerical implementation it is better to add this later (otherwise we would have to keep the value function E for all k)
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


        #value of waiting (next period capital is depreciated current capital)
        #(If the depreciated capital falls below the lowest grid point then truncate it to avoid extrapolation issues)
        Vwait = EV(max(par.k_min,(1-par.δ)*k),z_ind,par.shock_mc.p,par.N_z,Vint)
        #Vadj = E + (1-δ)k (E does not contains the policy)
        Vadj = SE.E[z_ind] + (1-par.δ)*k

        SE.ξc[i] = (Vadj - Vwait)/(SE.w*SE.pd)

        #debug
        if(SE.ξc[i]<0.0)
            error("ξc should be non-negative. Check the code.")
        end

    end


    #placeholder for stopping criterion
    return 0.0
end

#Function get_Vint constructs an interpolator for the ex ante value function.
#It returns object Vint, where Vint[z_ind](k) gives the value for idiosyncratic productivity with index z_ind and value of capital k.
#(The object is an array of interpolators, one for each possible shock realisation)
#Inputs are the parameter struct (which contains capital grid and number of grid points) and the value function on the grid.
#Note: Currently, interpolant returns an out of bounds error if extrapolation is attempted.
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
    end
    return EV
end

#Function find_Ku find the vlaue of capital which maximises the ex ante value function (Ku for Capital unconstrained by adjustment costs)
#At this stage it uses a simple search algorithm withing bounds of the capital grid. The issue is that this does not use any initial guess so is not very efficient and may not be very stable either in case there are local optima (due to limitations of the Optim.jl package where univariate bounded optimisation does not use an initial guess). If this is a performance bottleneck, or if local instability is a problem, a way forward is to (1) write a wrapper function which allows evaluation of EV outside of capital grid (nearest neighbour plus a steep convex penalty function of distance from grid boundary), (2) use an unconstrained optimisation algorithm which uses an initial guess - policy function from previous iteration in VFI - and should converge faster (3) check if the optimal value falls withing boundaries and if not, use the bounded optimisation search algorithm.
function find_KU(z_ind,P,N_z,Vint,k_min,k_max)
    #This uses default values for algorithm settings, does not use an initial guess!
    #note -EV (because it's a minimisation function, we want to maximise)
    #Also note that this does not include the (1-δ)k term which does not depend on the choice of next-period capital and is added to the value manually outside of this function.
    res = optimize(kpr -> -(-kpr + EV(kpr,z_ind,P,N_z,Vint)),k_min,k_max);

    #Return a pair - the optimal capital choice, and the associated maximum value
    return Optim.minimizer(res),-Optim.minimum(res) #(-1 again due to maxmin)

    #Comment: If there are issues with convergence etc. implement some robustness,exception handling, etc.
    #robustness  - compare the solution with the evaluated initial guess. If the initial guess is better, then return that.
end

#Function update_V! updates the value function (and overwrites the previous one). It returns value corresponding to the stopping criterion.
function update_V!(par,SE)

    #placeholder for stopping criterion
    return 0.0
end
