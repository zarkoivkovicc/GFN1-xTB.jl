module PeriodicTable
export anumber,amass, asymbol
atomicnumber = Dict{String,Int64}("H" => 1,
                    "C" => 6,
                    "N" => 7,
                    "O" => 8)
atomicmass = Dict{String,Float64}("H" => 1.008,
"C" => 12.011,
"N" => 14.007,
"O" => 15.999)
atomicsymbol = Dict{Float64,String}(1 => "H",
6 => "C",
7 => "N",
8 => "O")
function anumber(element::AbstractString)::Int64
    return atomicnumber[element]
end
function amass(element::AbstractString)::Float64
    return atomicmass[element]
function asymbol(elnumber::Int64)::String
    return atomicsymbol[elnumber]
end
end
end