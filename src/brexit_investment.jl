module brexit_investment

using Revise

using Optim, Parameters, QuantEcon, BenchmarkTools,Dates,JLD
using .Threads #so we don't have to write Threads.@threads every time


#Type definitions and their functions (brexDefs.jl)
export pars,stat_equil,check_par

#Functions from brexComp.jl
export SE_compute!

#Export functions from brexTools.jl
export genFolderName,prep0,saveAll,saveInFile,loadFromFile


### Include source files ###
include("brexDefs.jl")
include("brexComp.jl")
include("brexTools.jl")


end # module
