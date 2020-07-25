#=
This is an example parameter file.

Parameters not set here will be set to their default values. For default values and the meaning of parameters see the definition of struct
type pars in module BrexDefs.

The script needs to define the following variables: par_comm, par_diff, PS, t_brex (see below for their example values)

!Warning: The types of constants uses in this script must agree with the abstract types given in the struct definition. This mainly means that parameters which are supposed to be floats need to be given as Float constants (for example k_min = 1.0 is fine, k_min = 1 is not, because the program is trying to save an Integer variable into an Float type - no automatic conversion is performed in the current version of the program)

=#

#Parameters which differ from the default values but are common to all aggregate states are in par_comm

#Example: Set the discount factor and preferences in all states [everything else will be equal to default values]
par_comm = "β = 0.9,σ=2.0"

#Parameters which are different between aggregate states are set in par_diff.
#The number of elements of par_diff determines the number of states everywhere in the program.
#The first state corresponds to pre-brexit state (initial stationary equilibrium)

#Example: Tarrifs before Brexit are 0, in soft Brexit they are 0.05, in hard Brexit they are 0.1,
#TFP is also possibly affected, but no change in the baseline..

par_diff = [
    "τ = 0.0,A=1.0",
    "τ = 0.05,A=1.0",
    "τ = 0.1,A=1.0"
]


#Probabilities for each aggregate state when Brexit is resolved.
#The dimension must be the same as the number of different parametrisations set in par_diff
const PS = fill(1/length(par_diff),length(par_diff))

#Another example - probability 0 of No Brexit, and 50/50 soft/hard:
#const PS = [0.0,0.5,0.5]

#Time of Brexit (period numbering starts at 1, so t_brex=3 for example means 2 years of Pre-Brexit state before the aggregate uncertainty is realised)
const t_brex = 3
