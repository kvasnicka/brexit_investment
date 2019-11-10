#Activate project environment and install all dependent packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Optim
println("Packages loaded...")
