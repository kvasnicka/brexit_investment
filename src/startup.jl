#=
startup.jl is a startup script which does several things:

1) Parses command line arguments
2) Loads parameters from given parameter file and constructs an array of parameters
3) Displays some initial messages

After running the script, the following variables are available:

N_th::Int64 - the number of threads
N_S::Int64 - number of aggregate states (parametrisations)

parfile::String - name of the parameter file


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
println("****** Brexit investment model ******")
println("_____________________________________")
println("Number of available threads: $N_th")
if N_th == 1
    println("Consider increasing the number of threads by 'export JULIA_NUM_THREADS=#'")
end

#Load parameters from file
if length(parfile) == 0
    println("No parameter file given. Using default parameter file example.jl.")
    parfile="example.jl"
end
include("../parameters/$parfile") #directory above
println("Parameters loaded from file parameters/$parfile")

#Number of aggregate states of the economy/parametrisations
const N_S = length(par_diff)

#Generate array par of parameters for different long-run equilibria.
par = fill(pars(),N_S)
#Generate different parameters using the strings above
for i=1:N_S
   expr = "par[$i] = pars("
   if length(par_comm)>0
        expr = expr*par_comm*","
    end
   if length(par_diff)>0
        expr = expr*par_diff[i]
    end
    expr = expr*");"
    #parse the string as a command and evaluate
    ex1 = Meta.parse(expr)
    eval(ex1)
end
