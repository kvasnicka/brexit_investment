#Activate project environment and install all dependent packages. This should need to be run only once.
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Optim

#Parameters (so far hard-coded here)
