#=
startup.jl is a startup script which does several things:

1) Parses command line arguments
2) Loads parameters from given parameter file and constructs an array of parameters
3) Displays some messages

After running the script, the following variables are available:

const N_th::Int64 - the number of threads
const N_S::Int64 - number of aggregate states (parametrisations)
parfile::String - name of the parameter file
par::Array{pars{Float64,Int64},1} - parameters for each steady state
=#

#Get command line arguments, save them in parsed_args named tuple (global variable). Also save parfile - the name of the parameter file
using ArgParse #for parsing command line arguments.
#Load command line arguments. Move this to an included file later:
#Set up table of arguments passed by command line - see doc of ArgParse for details
s = ArgParseSettings()
@add_arg_table! s begin
    "--parfile"
        help = "parameter file name"
        arg_type = String
        default = ""
    "--numworkers"
        help = "number of workers"
        arg_type = Int
        default = 1
    "--slurm"
        help = "Use ClusterManagers for SLURM. 1 yes, 0 no (default)."
        arg_type = Int
        default = 0
end
#now get the parsed arguments
parsed_args = parse_args(ARGS, s)

#Copy the often used variables from the dictionary
parfile = parsed_args["parfile"]

#Intro message
const N_th = Threads.nthreads()
println("**************************************")
println("****** Brexit investment model *******")
println("**************************************\n")

println("Number of available threads: $N_th\n")
if N_th == 1
    println("Consider increasing the number of threads by 'export JULIA_NUM_THREADS=#'")
end

#Load parameters from file
if length(parfile) == 0
    println("No parameter file given. Using default parameter file baseline.jl.")
    parfile="baseline.jl"
end
include("../parameters/$parfile") #directory above
println("Parameters loaded from file parameters/$parfile\n")

#Number of aggregate states of the economy/parametrisations
const N_S = length(par_diff)

#Generate array par of parameters for different long-run equilibria
#using values of parameters from the parameter file.
expr = "par=["
for i=1:N_S
    global expr
   expr*="pars("
   if length(par_comm)>0 expr *= par_comm*"," end
   if length(par_diff)>0 expr *= par_diff[i] end
    expr*="),"
end
expr*="]"

#parse the string as a command and evaluate
expr2 = Meta.parse(expr)
eval(expr2)
expr = nothing;expr2 = nothing

#Check parameters (element-wise)
check_par.(par,N_S)

#Initialise results folder:
foldername = "results/"*genFolderName(parfile)
#Create results folder if it does not already exist
if !isdir("results")
    mkdir("results")
end
#Create the results/foldername folder. If it already exists, an error will occur
#(This could happen if the program is executed multiple times within a single second so the folder names are the same).
mkdir(foldername)
#Copy patemeter file in the results folder (for reproducibility)
cp("parameters/$parfile","$foldername/$parfile")

#Print a message about parameters:
println("*********Summary of parameters:*************")
println("All baseline parameters (pre-Brexit):")
@show par[1] #par[1] is the baseline (pre-Brexit)
println("_______________________")
println("Parameters that differ between stationary equilibria:")
@show par_diff
println(" ")
println("********************************************")
