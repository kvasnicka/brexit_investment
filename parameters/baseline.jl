#=
This is an example parameter file.

Parameters not set here will be set to their default values. For default values and the meaning of parameters see the definition of struct
type pars in module BrexDefs.

The script needs to define the following variables: par_comm, par_diff, PS, t_brex (see below for their example values)

!Warning: The types of constants uses in this script must agree with the abstract types given in the struct definition. This mainly means that parameters which are supposed to be floats need to be given as Float constants (for example k_min = 1.0 is fine, k_min = 1 is not, because the program is trying to save an Integer variable into an Float type - no automatic conversion is performed in the current version of the program)

There are 2 groups of parameters:
    (1) Parameters specific to a stationary equilibrium (aggregate state)
        - set using par_diff
    (2) Parameters non-specific to stationary equilibrium
        - set using par_comm
=#

#(1) Parameters specific to a stationary equilibrium (aggregate state):

#Example: Set the discount factor and preferences in all states [everything else will be equal to default values]
par_comm = "
β = 0.9,
σ = 2.0
"
#=
Example2: Leaving everything at default but changing the probability of different Brexit outcomes (first number is probability of no Brexit, probabilities have to sum to 1)
par_comm = "
PS = [0.0,0.5,0.5]
"
=#

#(Warning: putting a comma after the last entry results in an error!)

#(2) Parameters non-specific to stationary equilibrium:
#The number of elements of par_diff determines the number of states everywhere in the program.
#The first state corresponds to pre-brexit state (initial stationary equilibrium)

#Example: Tarrifs before Brexit are 0, in soft Brexit they are 0.05, in hard Brexit they are 0.1. TFP is also possibly affected, but no change in the baseline..
par_diff = [
    "τ = 0.0,A=1.0",
    "τ = 0.05,A=1.0",
    "τ = 0.1,A=1.0"
]

#Warning: Setting the values of some variables to be different between Brexit states will break the program - usually these are parameters related to the solution technique (see brexDefs.jl for details). In general, changing values of parameters which have economic interpretation should be fine.

#=If loadData = true, the program reads initial guesses for the Stationary equilibrium and transition paths from folder results/loadFolder

This can be very useful when the computation takes a long time, or when convergence is problematic (then we can use as initial guess stationary equilibria and transition paths obtained for not too different set of parameters)

These parameters are optional, if not defined here the program works fine.

Warning: If the folder does not contain the data or contains wrong data (for example data obtained for a different number of grid points or length of transition path) the program will crash or may result in unpredictable behaviour!
=#

loadData = false
loadFolder = "20200810_213249_example"
