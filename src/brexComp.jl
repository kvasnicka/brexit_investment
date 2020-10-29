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

#Note: When looking for stationary equilibrium we can just normalise Uc = 1 (that only matters when we have aggregate shocks so the marginal utility of consumption is different in different states)

=#

#Iterate over prices.
for SEind = 1:par.SE_maxiter

#Given prices, find the excess demands.
#(The goal is to find prices such that all excess demands are close to zero)
#Later on we will use some algorithm to find a vector of prices such that the excess demands are approximately zero (market-clearing).
#In one dimension a bisection algorithm would be a possible approach, in multiple dimensions we will need to do something more sophisticated.

#So far we do not loop over prices to find a statoinary equilibrium with 0 excess demands, just compute the excess demands (and the new stationary equilibrium candidate - which will contain updated value functions and policy functions, but not updated prices).

#Get excess demands and new stationary equilibrium
excess,SEn = ED(par,SE)


#Stopping rule - to be implemented. At this stage of development break after the first loop.
println("Only one iteration in SE_compute! preformed, updating prices not implemented yet.\n")

#Return the new stationary equilibrium.
return SEn

break

end

end #of SE_compute!


#Function ED computes the excess demands for given prices. SEg is the initial guess for a stationary equilibrium (contains value function, policy function, prices).
#Function ED returns a pair: excess, SEn
#where SEn is the new candidate for stationary equilibrium (prices, policy functions, etc.). This is safer than simply overwriting the initial guess, especially when using a generic root finding algorithm (for some vector of prices off-equilibrium the policy functions and stationary distribution might be very far from the equilibrium one, and may be completely non-sensical so we might run into convergence problems).

#The prices for which we compute excess demands are part of initial guess SEg.
#Later on we can either write a wrapper function of prices only, or include prices as explicit argument here (if using a generic root-finding algorithm).
function ED(par,SEg)
    SEn = copy(SEg) #Initialise the new stationary equilibrium (do not overwrite the guess)

    excess = [0.0,0.0,0.0] #Initialise excess demands
    #Compute the real wage (using households FOC and log linear utility) before copying new SE candidate (in any case, prices in SEn should not be used anywhere as they were not updated)
    SEg.w = par.χ/SEg.Uc

    #Solve the firm's problem given prices (these are contained in struct SEg along with the initial guess of value and policy functions), saving the new policy functions etc. in SEn (new stationary equilibrium candidate)
    firm_solve!(par,SEn,SEg)

    #Find the stationary distribution of firms
    stationary_μ!(par,SEn)


    #To do: compute excess demands using market clearing conditions

#Return excess demands and the new candidate for stationary equilibrium
return excess,SEn

end


function stationary_μ!(par,SE)
    #Function stationary_μ! finds the stationary distribution μ using a simulation method.
    #Note: We use a simulation method as opposed to constructing a transition matrix and looking for a stationary distribution directly because the major part of the algorithm focuses on transition paths, where the distribution has to be updated by applying the policy function only once. In this case, it is not efficient to construct the matrix (costly) and then apply matrix multiplication, we might as well directly perform the simulation step.


    #Construct also interpolator for the adjustment cost threshold ξc (it is needed only if we are using a finer histogram than the policy function grid)
    #ξcint = get_ξcint(par,SE.ξc)
    #Pass this to the function

    #Pre-compute the closest points in the grid for both the optimal choice and for the depreciated capital, and the associated weights. This saves a lot of time compared to finding these points repeatedly (costly grid search).
    hclose,hclosew,kdclose,kdclosew = prepClosePoints(SE,par);

    #Move the below to a function: sim_step()

    for i = 1:10 #iterations to update the distribution: So far just do a few. later, do a lot and implement a stopping rule.

    #Update the distribution
    μn = update_μ(par,SE,hclose,hclosew,kdclose,kdclosew)

    #Update the distribution in SE (in place). We can also check convergence criteria.
    SE.μ = copy(μn)

    end

end #stationary_μ!

#=
Function update_μ takes a given distribution of firms and returns the next-period distribution of firms μn by applying the firm's policy functions for capital and the law of motion for exogenous shocks. The current distribution and policy functions are saved in struct SE.

Note: we need to precompute and pass the indices and weights of closest grid points (hclose,hclosew,kdclose,kdclosew).
=#
function update_μ(par,SE,hclose,hclosew,kdclose,kdclosew)
    μn = zeros(par.N_kh,par.N_z) #initialise the distribution
    μ_ind = CartesianIndices((1:par.N_kh,1:par.N_z)) #get indices for looping

    #frequently used - shock process transition matrix and state values.
    P = par.shock_mc.p #frequently used
    zvals = par.shock_mc.state_values

    #!!! Issue: multithreading does not work here because different points are trying to write into the same memory address (because multiple points assign mass to the same index). For now just do it serially. But in the future might need to write something more sophisticated.

    #cycle over all points in the historgram
    #For now use at least @simd instead of threads. Once we benchmark it on a super-computer, we can see what difference @simd makes, and if it makes sense to write a proper parallel implementation (or try to fix @threads by using some lock tricks - see Julia documentation on multi-threading).
    #We can use simd because the order of the loop does not matter, and because this is the inner-most loop.
    @simd for i in eachindex(μ_ind)
        #index of capital is μ_ind[i][1]
        #Index of shock realisation is μ_ind[i][2]

        #current capital and shock realisation
        k_ind = μ_ind[i][1]
        z_ind = μ_ind[i][2]
        k = par.k_gr_hist[k_ind] #histogram grid, can be finer than the grid used to solve the firms' problem.
        z = par.shock_mc.state_values[μ_ind[i][2]]

        #=Depending on the adjustment cost realisation ξ, the firm will (A) adjust and set its capital equal to the optimal unadjusted capital (h(Z)), or (B) not adjust and set capital to (1-δ)k.

        The share of firms that adjust is G(ξc) = ξc/ξbar, where ξc is the cutoff point for the current state, given by policy function ξc(k,z). All these firms choose that same capital stock. The share of firms that do not adjust is 1-G(ξc), and all these firms let capital depreciate.

        Because the choices of capital in the policy function h and the depreciated capital do not fall exactly on the histogram grid, we need to find the two closest points in the grid, and split the measure of firms between the two points (using a linear weighting which puts more weight on the closer grid point). Because the choices are fixed, we have precomputed the indices of grid points and the weights.
        =#

        #(A) firms that adjust
        if par.N_k == par.N_kh
            #Same grid for policy function and histogram.
            G = SE.ξc[k_ind,z_ind]/par.ξbar #Share of firms that adjust
        else
            #different grids (need to use interpolation to get the cutoff point)
            error("Different grids not supported yet. Generalise function prepClosePoints.")
        end

        #h = SE.h[z_ind] #capital choice for firm that adjust - just for debugging (we already precomputed the closest points and weights)
        #μn[hclose[z_ind,1],:] is the row of the distribution (associated with the capital point on the grid) which is updated.
        #SE.μ[k_ind,z_ind] is the mass to be distributed between all points
        #G is the share of firms that adjust

        #The closest grid points and weights:
        h1 = hclose[z_ind,1]
        h2 = hclose[z_ind,2]
        w1 = hclosew[z_ind,1]
        w2 = hclosew[z_ind,2]

        #Add mass to the first close point (index hclose[z_ind,1], weight hclosew[z_ind,1])
        μn[hclose[z_ind,1],:] +=  SE.μ[k_ind,z_ind]*G*hclosew[z_ind,1]*P[z_ind,:]
        #2nd point:
        μn[hclose[z_ind,2],:] +=  SE.μ[k_ind,z_ind]*G*hclosew[z_ind,2]*P[z_ind,:]


        #(B) firms that do not adjust
        #The same as above but the points and weights are in kdepclose, and the indexing is in terms of capital, because the "policy" depends on capital now (1-δ)k
        μn[kdclose[k_ind,1],:] +=  SE.μ[k_ind,z_ind]*(1.0-G)*kdclosew[k_ind,1]*P[z_ind,:]
        #2nd point:
        μn[kdclose[k_ind,2],:] +=  SE.μ[k_ind,z_ind]*(1.0-G)*kdclosew[k_ind,2]*P[z_ind,:]
    end

    #debug: print sum of the distribution
    #a = sum(μn)
    #print("debug in update_μ. sum of μ mass = $a \n")

    return μn
end

#Function firm_solve! solves the firm's problem given prices.
#SEn is the new stationary equilibrium (candidate) where firm's value function, policy function, etc. will be saved.
#SEg is the initial guess which shall not be changed by the function call in any way.
#The function uses Caretesian indexing implemented in Julia - see https://julialang.org/blog/2016/02/iteration/
function firm_solve!(par,SEn,SEg)
    #Compute the optimal labour supply (this is a static problem, depends only on prices and not on the continuation value).
    compute_N!(par,SEn)

    debug = true #write messages about convergence if true

    #Value function iteration
    for VFIind = 1:par.VFI_maxiter
        Vint = get_Vint(par,SEn.V) #construct interpolator for the current value function (which is SE.V)

        #If this is the first iteration, or if the iteration index is a multiple of VFI_howard, update the policy function.
        if (mod(VFIind,par.VFI_howard) == 0 || VFIind == 1)
            #SEn contains the current policy function (and value function, prices, etc.) and the policy function part will be overwritten.
            #poll_diff is the stopping criterion for policy function.
            pol_diff = update_pol!(par,SEn,Vint)

            #debug: report the difference in poll_diff

            if debug
            a = pol_diff[1]
            b = pol_diff[2]
            print("Debug in firm_solve!. Iteration $VFIind. AAD policy difference: h = $a, ξc = $b \n")
            end
        else
            #We do not update the policy but we still have to update the value of choosing the policy h!
            for z_ind = 1:par.N_z
                #Value of optimal adjustment is just the expected discounted ex ante value of choosing the optimal capital level h, minus the chosen capital level.
                SEn.E[z_ind] = par.β* EV(SEn.h[z_ind],1,par.shock_mc.p,par.N_z,Vint) - SEn.h[z_ind]
            end
        end

        #Get the new updated value function
        V_new = get_V0(par,SEn,Vint)

        #Check stopping criteria for the value function (difference b/w Vnew and SEn.V)

        if debug
        #Absolute relative difference at each grid point:
        ARD = (abs.(V_new-SEn.V)./(abs.(SEn.V) .+ 0.001))
        MARD = maximum(ARD) #maximum relative absolute deviation
        AARD = mean(ARD) #average relative absolute deviation

        print("Value function MARD = $MARD, AARD = $AARD \n")
        end


        #The VFI converges monotonically. Right now, the stopping rule is not checked and the maximum number of iterations is performed.

        #Update the value function
        SEn.V = copy(V_new)

    end #VFI for loop
end #firm_solve!


#Function compute_N! computes the labour supply for each point in the state space (and saves it in SE.N)
function compute_N!(par,SE)
    V_ind = CartesianIndices((1:par.N_k,1:par.N_z))
    @threads for i in eachindex(V_ind)
        k = par.k_gr[V_ind[i][1]]
        z = par.shock_mc.state_values[V_ind[i][2]]
        SE.N[i] = min((SE.w/(par.A*z*k^par.α*par.ν))^(1/(1-par.ν)),par.Nmax)
    end
end

#Function update_pol! updates the policy function (and overwrites the previous one). It returns value corresponding to the stopping criterion.
function update_pol!(par,SE,Vint)
    #Generate iterator (this takes almost no memory and is fast)
    V_ind = CartesianIndices((1:par.N_k,1:par.N_z))

    #First compute the unrestricted optimal capital level. This does not depend on the current level of capital, or adjustment costs, so we only loop over the current productivity.

    poldiff = [0.0,0.0] #Vector for reporting difference in the policy function (average absolute deviation in the capital and policy)

    #Get a copy of the policy functions so we can compute difference
    hcopy = copy(SE.h)
    ξccopy = copy(SE.ξc)



    @threads for z_ind=1:par.N_z
        #Compute optimal level of capital in the absence of adjustment costs
        #This will be saved in policy function h, value saved in E

        SE.h[z_ind],SE.E[z_ind] = find_KU(z_ind,par.β,par.shock_mc.p,par.N_z,Vint,par.k_min,par.k_max)

        #Warning! - the continuation value E does not include the left-over capital (1-δ)*k. So it is not exactly Vadj as in our model notation, but more like E in KT2008 paper. In the numerical implementation it is better to add this later (otherwise we would have to keep the value function E for all k)
    end


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

        SE.ξc[i] = min((Vadj - Vwait)/(SE.w*SE.pd),par.ξbar)
        #truncate so the adjustment threshold is always in [0,ξbar], where ξbar is the maximum adjustment cost (this is important so that we get values of the cdf of ξ at ξc in [0,1] later)

        #For development only
        #

        #If the adjustment cost is negative, it means that the optimal continuation value (without paying adjustment costs) is less than the value of waiting, which should never be the case. It happens only rarely (small numerical errors).
        if(SE.ξc[i]<0.0)
            SE.ξc[i]=0.0
        end
    end


    #Return the maximum absolute deviation (relative for capital, absolute for ξc - because this is frequently close to 0)
    poldiff = [maximum(abs.(hcopy-SE.h)./(abs.(SE.h).+0.0001 )),maximum(abs.(ξccopy-SE.ξc))]

    #placeholder for stopping criterion
    return poldiff
end

#Function get_V0 uses the policy functions and the value function contained in SE to return a new updated ex ante value function. Unlike update_pol it does not overwrite the original value function in place (due to synchronisation issues and checking convergence). A few things such as generating the iterator, interpolator, and looping, are the same as in updat_pol! and are explained there in more detail.
function get_V0(par,SE,Vint)
    Vnew = similar(SE.V) #initialise new value function

    V_ind = CartesianIndices((1:par.N_k,1:par.N_z))

    #parallel loop over grid points
    @threads for i in eachindex(V_ind)
        #V_ind[i] contains the Cartesian index for the i-th element of the matrix, V_ind[i][j] the index for the j-th state which corresponds to the grid point.

        #index of capital is V_ind[i][1]
        #Index of shock realisation is V_ind[i][2]

        #current capital and shock realisation
        k_ind = V_ind[i][1]
        z_ind = V_ind[i][2]
        k = par.k_gr[k_ind]
        z = par.shock_mc.state_values[z_ind]

        #Get the updated ex ante value:
        #current return (this does not depend on ξ so no need to take expectation)
        Vnew[i] = (par.y(par.A,z,k,SE.N[i]))*SE.pd

        #=
        Now use the facts that ξ is U[0,ξbar]. For all ξ>ξc(k,z), the firm does not adjust, in which case the ex post value V1(k,z,ξ) does not depend on ξ, and we get that part of the expectation trivially (just a constant times G(ξc), where G is the distribution function of ξ (G(ξ) = ξ/ξbar).

        For ξ < ξc, the ex ante value does depend on ξ, but only through the adjustment cost term AC(ξ) which is a linear function so we can take an analytical expectation of this term.


        G = G(ξbar) is the share of firms that adjust their investment.
        =#
        G = SE.ξc[i]/par.ξbar
        #value if not adjusting * share of firms not adjusting (prob of not adjusting)
        Vnew[i] +=  (1-G)*EV(max(par.k_min,(1-par.δ)*k),z_ind,par.shock_mc.p,par.N_z,Vint)

        #value if adjusting excluding the terms that depends on ξ * share of firms adjusting + expeted AC(ξ) (conditional expectation)

        #For adjusting firms, do not forget that we need to addd the (1-δ)k term to E to get Vadj as in the paper.
        Vadj = SE.E[z_ind] + (1-par.δ)*k
        Vnew[i] += G*Vadj

        #Subtract the expectation of adjustment costs paid, which is E(ξ|ξ<ξc) = ξc^2/(2*ξbar).
        Vnew[i] -= (SE.ξc[i]^2)/(2.0*par.ξbar) * SE.w*SE.pd

        #Multiply the value function by the marginal utility (since the value function is in terms of utils. Everything is premultiplied by this - current return, continuation value, adjustment costs).
        Vnew[i] *= SE.Uc

    end


    return Vnew
end

#Function Kpol(k,δ,ξ,ξc,h) (pol for policy) takes as input the current capital stock, depreciation, adjustment cost ξ, threshold ξc, and the optimal capital in the absence of adjustment costs h, returns the actual value of next-period capital
function Kpol(k,δ,ξ,ξc,h,k_min)
    if ξ < ξc
        return h
    else
        return max((1-δ)*k,k_min) #truncate to avoid extrapolation issues at the low end of the grid (if (1-δ)k would fall off the grid, such as at the first grid point, we are effectively letting the firm keep their current capital without paying the adjustment cost)
    end
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

#Function ξcint generates interpolator for the capital adjustment threshold.
#Same setting for the interpolation is used as for the value function.
function get_ξcint(par,ξc)
    if par.Vint_mode == 1 #linear
        ξcint = fill(LinearInterpolation(par.k_gr,ξc[1:par.N_k,1]),par.N_z)
        for zind = 2:par.N_z
            ξcint[zind] = LinearInterpolation(par.k_gr,ξc[1:par.N_k,zind])
        end
    elseif par.Vint_mode == 2 #cubic spline
        ξcint = fill(CubicSplineInterpolation(par.k_gr,ξc[1:par.N_k,1]),par.N_z)
        for zind = 2:par.N_z
            ξcint[zind] = CubicSplineInterpolation(par.k_gr,ξc[1:par.N_k,zind])
        end
    else
        error("Unsupported value of Vint_mode.")
    end
    return ξcint
end

#Expected ex ante value function. kpr is next-period capital level, z_ind the index of current shock realisation,P is the transition matrix, Vint is the value function interpolant obtained using function get_Vint
function EV(kpr,z_ind,P,N_z,Vint)
    EV = 0.0
    for zpr = 1:N_z #cycle over next-period shock realisations, P[z_ind,zpr] is the transition probability
        EV += P[z_ind,zpr]*Vint[zpr](kpr)
    end
    return EV
end

#Function find_Ku find the value of capital which maximises the ex ante value function (Ku for Capital unconstrained by adjustment costs)
#At this stage it uses a simple search algorithm withing bounds of the capital grid. The issue is that this does not use any initial guess so is not very efficient and may not be very stable either in case there are local optima (due to limitations of the Optim.jl package where univariate bounded optimisation does not use an initial guess). If this is a performance bottleneck, or if local instability is a problem, a way forward is to (1) write a wrapper function which allows evaluation of EV outside of capital grid (nearest neighbour plus a steep convex penalty function of distance from grid boundary), (2) use an unconstrained optimisation algorithm which uses an initial guess - policy function from previous iteration in VFI - and should converge faster (3) check if the optimal value falls withing boundaries and if not, use the bounded optimisation search algorithm.
function find_KU(z_ind,β,P,N_z,Vint,k_min,k_max)
    #This uses default values for algorithm settings, does not use an initial guess!
    #note -EV (because it's a minimisation function, we want to maximise)
    #Also note that this does not include the (1-δ)k term which does not depend on the choice of next-period capital and is added to the value manually outside of this function.
    res = optimize(kpr -> -(-kpr + β*EV(kpr,z_ind,P,N_z,Vint)),k_min,k_max);

    #Return a pair - the optimal capital choice, and the associated maximum value
    return Optim.minimizer(res),-Optim.minimum(res) #(-1 again due to maxmin)

    #Comment: If there are issues with convergence etc. implement some robustness,exception handling, etc.
    #robustness  - compare the solution with the evaluated initial guess. If the initial guess is better, then return that.
end

#Function findclosest2 finds the two points in a grid X which are closest to point x.
#If point x is exactly equal to some grid point, then both returned points will be equal to x.
function findclosest2(x,X)
    absdist = abs.(X .- x)

    #value and index of the closest element
    (minval, i) = findmin(absdist)

    if x == X[i]
       return [i,i] #exact match
    end

    if i == 1
       return [1,2] #first point the closset, so the 2 closest are 1 and 2
    elseif i == length(X)
       return [length(X)-1,length(X)]
    end

    #The most common case, we are not at the edge of the grid, so we can take a look at the left
    #and the right neighbouring point to the minimum distance one and pick the closer one.
    if(abs(x-X[i-1]) < abs(x-X[i+1]))
       return [i-1,i]
    else
        return [i,i+1]
    end
end



function prepClosePoints(SE,par)
    #=
    Function prepClosePoints finds the closest two points in the histogram for:
    (1) the capital choices contained in the policy function (h)
    (2) depreciated capital choice (1-δ)k
    It also finds the weights on these two points.
    =#

    hclose = fill(0,par.N_z,2) #fill so it's integer type
    hclosew = zeros(par.N_z,2) #weight
    for z_ind = 1:par.N_z
        hclose[z_ind,:] = findclosest2(SE.h[z_ind],par.k_gr_hist)
        #(SE.h[z_ind] is the optimal capital level at this grid point)

        #generate weights on the two points (convex combination)
        k1 = par.k_gr_hist[hclose[z_ind,1]]
        k2 = par.k_gr_hist[hclose[z_ind,2]]

        #If both points are the same, it means that the policy implies a corner solution.
        #We need to treat this case separately to avoid division by 0. Assign equal weights.
        if(k1 == k2)
            hclosew[z_ind,1:2] .= 0.5
        else
            hclosew[z_ind,1] = (k2 - SE.h[z_ind])/(k2-k1)
            hclosew[z_ind,2] = 1-hclosew[z_ind,1]
        end
    end

    kdclose = fill(0,par.N_k,2)
    kdclosew = zeros(par.N_k,2)
    for k_ind = 1:par.N_k
        k = par.k_gr_hist[k_ind] #current stock
        kd = max((1-par.δ)*k,par.k_min) #depreciated capital truncated from below (so we don't fall off the grid)
        #find the closest gridpoint to next-period depreciated capital (1-δ)k
        kdclose[k_ind,:] = findclosest2(kd,par.k_gr_hist)

        k1 = par.k_gr_hist[kdclose[k_ind,1]]
        k2 = par.k_gr_hist[kdclose[k_ind,2]]

        if(k1 == k2)
            kdclosew[k_ind,1:2] .= 0.5
        else
            kdclosew[k_ind,1] = (k2 - kd)/(k2-k1)
            kdclosew[k_ind,2] = 1-kdclosew[k_ind,1]
        end
    end
    return hclose,hclosew,kdclose,kdclosew
end
