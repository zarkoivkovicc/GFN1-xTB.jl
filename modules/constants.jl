module PeriodicTable
export anumber
atomicnumber = Dict{String,Int64}("H" => 1,
    "C" => 6,
    "N" => 7,
    "O" => 8)
function anumber(element::AbstractString)::Int64
    return atomicnumber[element]
end
end