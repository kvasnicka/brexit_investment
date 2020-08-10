#brexTools.jl contains functions for input/output (saving a reading results), and various other tools.

#Function gen_filename generates the name of folder where results are saved. This is based on time of script execution, and the name of the parameter file.
function genFolderName(parfile)
    #If parfile ends with ".jl", do not use that part in generating folder name
    #(this does not handle other file types, but parameter files should be .jl and
    #it is not a problem in Unix if the folder name contains a dot)
    if parfile[end-2:end] == ".jl"
        parfileloc = parfile[1:end-3]
    else
        parfileloc = parfile
    end

    #Get year, month, day, hour, minute, second as strings, prepend zeros as necessary
    #(for correct sorting when some numbers are less than 10)
    dat = Dates.now()
    y = prep0(year(dat),4)
    m = prep0(month(dat),2)
    d = prep0(day(dat),2)
    h = prep0(hour(dat),2)
    mi = prep0(minute(dat),2)
    s = prep0(second(dat),2)

    #Construct foldername and return it.
    #format is yymmdd_hhmmss_parfile
    foldername = y*m*d*"_"*h*mi*s*"_"*parfileloc
    return foldername
end

#Function prep0 converts number num to a string, and adds leading zeroes such that the length of the string is at least minlength
function prep0(num,minlength)
    num_str = string(num)
    #prepend zeroes
    for i=1:(minlength-length(num_str))
       num_str = "0"*num_str
    end
    return num_str
end

#Function saveAll saves stationary equilibria and transition paths in JLD (binary) format
function saveAll(foldername,SE,TP)
    #For saving complex data types use JLD format (need using JLD)
    #The file needs to have suffix .jld so the type is correctly inferred.
    save("$foldername/SE.jld","SE",SE)
    save("$foldername/TP.jld","TP",TP)

    #Below: saving data in text format (to be used for simulation if desired)
    #saveInFile(V,"test.out",foldername,text = false,replace = true)
    #loadFromFile!(V,"test.out",foldername,text = false)
end

#Functions saveInFile can be used to manually save data in binary or text format, not using JLD format. This can be useful if we want to visually inspect the output of the program (such as simulation results), or import it in some other program for further analysis.

#Function saveInFile saves variable x in folder foldername
#By default data is saved in binary form (this is usually faster).
#If text = true, the data is saved in text form so we can read it easily
function saveInFile(x,filename,foldername;text = false,replace = true)
    fullpath = "$foldername/$filename"

    #If the file exists and replace = true (default), delete it and continue
    #If the file exists and replace = false, write a warning message and exit the function
    if isfile(fullpath)
       if replace
            rm(fullpath)
        else
           println("Warning: File $fullpath already exists but replace = false, so is not replaced.")
            return nothing
        end
    end
    io = open(fullpath,"w")
    if text #text output
            show(io,x)
        else #binary output
            write(io,x)
    end
    close(io)
end

#Function loadFromFile loads a given variable from file
#text has to be the same as when the variable was saved
function loadFromFile!(x,filename,foldername;text = false)
    fullpath = "$foldername/$filename"

    #If text = true then the original data was saved in text form using function show,
    #Implementation of reading text data not done yet
    #A way to do this is to use package DelimitedFiles to save data and then load them
    #(low priority change)
    if text
       error("text = true not yet implemented in function loadFromFile.")
    end
    io = open(fullpath,"r")
            read!(io,x)
    close(io)
end
