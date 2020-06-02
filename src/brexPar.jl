#=
This file contains module brexPar.jl, which exports functions used for parallelisation.
=#
module brexPar

export V_ind_unfold

#Function V_ind_unfold takes the one-dimensional (folded) index (ind) and unfolds it into index for capital and shock in the firms' value function problem.
#This is returned as a 2-dimensional array
#Note: This is for the start-of-period value function V(k,z;agg. state).
function V_ind_unfold(ind,N_k,N_z)
    #Index of capital grid point
    i_k = convert(Int, ceil(ind/N_z));
    #Index for capital
    i_z = convert(Int, floor(mod(ind-0.05, N_z))+1);

    return [i_k,i_z]
end


end
