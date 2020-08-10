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
