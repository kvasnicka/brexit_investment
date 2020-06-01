# brexit_investment
Code for solving the model of investment response to Brexit in SOE

This is very much work in progress.

File main.jl is the main script which will solve the model, plot results, compute statistics, etc. See the comments inside the file for details.

Folder parameters contains various parameter configurations (see file parameters/example.jl for an example parameter file).

Folder src contains source files for various modules. The convention is that one file contains one module, and each module can contain multiple related functions, definitions of data types, etc.

Description of various data types follows:

src/brexDefs.jl contains definitions of data types.
src/brexPar.jl contains functions used in parallelisation.
