#=
This is the main program file for the Brexit investment paper.

Use: 'julia --parfile parameterfile' runs the program with parameters contained
in file './parameters/parameterfile'. File parameters/example.jl is an example
parameter file, used by the command 'julia main.jl --parfile example.jl'

The role of each parameter is described in module BrexDefs (in BrexDefs.jl),
in the definition of struct pars, together with their default values
(these are used when the program is called without arguments).

The parameters describe not only the economic environment, but also what the
program should do (for example solve the model or simulate only), and various
parameters related to the solution algorithm.
=#

#Activate project environment and install all dependent packages.
#This should need to be done only once
import Pkg
Pkg.activate(".")
Pkg.instantiate()

#Load necessary packages
using Optim, Parameters
#using Parameters
N_th = Threads.nthreads()

#Include file containing module BrexDefs (definition of data types) and load it
include("./src/brexDefs.jl")
using .BrexDefs

#Include file containing module BrexPar (parallelisation tools) and load it
include("./src/brexPar.jl")
using .BrexPar

#Get command line arguments, save them in parsed_args named tuple (global variable)
#parfile is the name of the parameter file
include("./src/commandLineArgs.jl")


println("****** Brexit investment model ******")
println("_____________________________________")
println("Number of available threads: $N_th")
if N_th == 1
    println("Consider increasing the number of threads by 'export JULIA_NUM_THREADS=#'")
end


#Load parameters from the given file unless no file was given
#(in this case use the constructor pars() with no arguments which results in default values)
if length(parfile) == 0
    par = pars()
    println("Default values of parameters used.")
else
    include("./parameters/$parfile")
    println("Parameters loaded from file ./parameters/$parfile")
end

#####################################################
#=
The main body of the program follows. The algorithm proceeds in the following
stages:

(0) Initialise variables containing output

(1) Compute stationary equilibria for the initial state (pre-referendum) and
each of the long-run post-Brexit states.

(2) Compute transition paths by a shooting algorithm.

(3) Compute various statistics, plot results.

=#
