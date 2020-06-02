#=
This file contains module brexPar.jl, which exports functions used for parallelisation.
=#
module brexPar

#In the past this module contained manual functions for transforming linear indices into Cartesian indices.
#However, this was redundant because these functions are directly implemented in Julia.
#There are a lot of other advantages - see https://julialang.org/blog/2016/02/iteration/

#=
We can for example create an iterator:
V_ind = CartesianIndices((1:N,1:N)).

#This takes minimal storage space so it's no issue creating an iterator for a huge matrix.

Then we can just loop
for i in eachindex(ValueFunction)
    V_ind[i] #this contains Cartesian index corresponding to an element of the value function

    V_ind[i][j] #index for j-th state

end

Then V_ind[i] will return Cartesian index,and
=#

#In a future release this module may be removed, or it may create iterators (globally available everywhere)

end
