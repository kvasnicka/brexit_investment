#=
This is an example parameter file.

A parameter file must construct a struct named par of type pars().

Various arguments in the constructor call change values of different parameters.

The example here changes the value of discount factor to 0.99 and leaves all other parameters at their default value.

For default values and the meaning of parameters see the definition of struct
type pars in module BrexDefs.

!Warning: The types of constants uses in this script must agree with the abstract types given int he struct definition. This mainly means that parameters which are supposed to be floats need to be given as Float constants (for example k_min = 1.0 is fine, k_min = 1 is not, because the program is trying to save an Integer variable into an Float type - no automatic conversion is performed)

=#
par = pars(β = 0.99)
