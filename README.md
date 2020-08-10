# brexit_investment
Code for solving the model of investment response to Brexit in SOE

File main.jl is the main script which solves the model, plots results, compute statistics, etc.

Use: 'julia main.jl --parfile parameterfile' runs the program with parameters contained in file 'parameters/parameterfile'.

File parameters/example.jl is an example parameter file, used by the command 'julia main.jl --parfile example.jl', and is also used as default if the program is called without specifying the parameter file.

The role of each parameter is described in module BrexDefs (in BrexDefs.jl), in the definition of struct pars, together with their default values (these are used when the program is called without arguments).

Other files:

Folder parameters contains various parameter configurations (see file parameters/baseline.jl for the baseline parameter file which can be used as template).

Folder src contains source files for various modules. The convention is that one file contains one module, and each module can contain multiple related functions, definitions of data types, etc.

Description of various files follows:

src/brexDefs.jl contains definitions of data types.
src/startup.jl is a startup script
