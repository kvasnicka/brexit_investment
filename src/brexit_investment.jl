module brexit_investment

using Revise

using Optim, Parameters, QuantEcon, BenchmarkTools
using .Threads #so we don't have to write Threads.@threads every time


#Type definitions and their functions (brexDefs.jl)
export pars,stat_equil,check_par

#Functions from brexComp
export SE_compute

include("brexDefs.jl")
include("brexComp.jl")


end # module
