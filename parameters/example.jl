#=
This is an example parameter file.

A parameter file must construct a struct named par of type pars().
Various arguments in the constructor call change values of different parameters.

The example here changes the value of discount factor to 0.99 and leaves all other
parameters at their default value.

For default values and the meaning of parameters see the definition of struct
type pars in module BrexDefs.

=#
par = pars(β = 0.99)
