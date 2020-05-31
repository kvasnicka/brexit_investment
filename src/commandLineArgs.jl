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
