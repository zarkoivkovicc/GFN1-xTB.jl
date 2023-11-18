module Parsers
include("constants.jl")
using .PeriodicTable
export parsexyz, loadparams, binpath
const global binpath::String = dirname(abspath(PROGRAM_FILE))

"""
    parsexyz(filename)

Extracts number of elements, elements list and atomic coordinates (in Angstroms) from an xyz file.

Input: filename::String \n
Output: natoms::Int64, elements::Array{String}(natoms), coordinates::Array{Float64}(natoms,3)

"""
function parsexyz(filename::String)
    xyzfile = readlines(filename)
    natoms = parse(Int64, xyzfile[1])
    elements = Array{String}(undef, natoms)
    coordinates = Matrix{Float64}(undef, 3, natoms)
    for atom in 1:natoms
        line = split(strip(replace(xyzfile[2+atom], r"\t{1,}|\s{2,}" => ' ')))
        elements[atom] = line[1]
        coordinates[:, atom] = parse.(Float64, line[2:end])
    end
    return natoms, elements, coordinates
end
"""
    This internal function that returns list of cleaned lines. If only serves to prepare parameters for parsing.
"""
function readparams(filename::String)
    rawdata = String[]
    lines = eachline(filename)
    for line in lines
        endline = occursin('#', line) ? findfirst('#', line) - 1 : length(line)
        cleanline = strip(replace(line[1:endline], r"\t{1,}|\s{2,}" => ' ')) #RegEx for multiple spaces or multiple tabs
        if !isempty(cleanline)
            push!(rawdata, cleanline)
        end
    end
    return rawdata
end

function parseline(type::Type, line::AbstractString)
    return tryparse.(type, split(line))
end

function parsematrix(dim1::Int64, dim2::Int64, rawmatrix::Array{String})
    matrix = Array{Float64}(undef, dim1, dim2)
    for i in 1:dim1
        matrix[i, :] = parseline(Float64, rawmatrix[i])
    end
    return matrix
end

function parseparams(rawdata::Array{String})
    BORH_TO_Å = parse(Float64, rawdata[1])
    AU_TO_EV = parse(Float64, rawdata[2])
    ntypes = parse(Int64, rawdata[3]) # number of atom types in paramteres file
    indx = Dict{String,Int64}(split(rawdata[4]) .=> 1:ntypes) #index of elements for parameters
    #Zero order repulsion energy parameters
    k_f = parse(Float64, rawdata[5])
    α = parseline(Float64, rawdata[6])
    z_eff = parseline(Float64, rawdata[7])

    #Dispersion parameters
    a1, a2, s6, s8, k_cn, k_l = parseline(Float64, rawdata[8])
    q_a = parseline(Float64, rawdata[9])
    unsc_radii = parseline(Float64, rawdata[10])
    sc_radii = 4 / 3 / BORH_TO_Å .* unsc_radii
    max_numrefcn = parse(Int64, rawdata[11])
    numrefcn = parseline(Int64, rawdata[12])
    refcn = Matrix{Float64}(undef,max_numrefcn,max_numrefcn)
    for i in 1:ntypes
        refcn[1:numrefcn[i],i] = parseline(Float64, rawdata[12+i])
    end
    # 12 + ntypes of lines have been used
    global linenum = 12 + ntypes
    c_ref = Matrix{Matrix{Float64}}(undef,ntypes,ntypes)

    for _ in 1:ntypes*(ntypes+1)/2
        temp = parseline(Int64, rawdata[linenum+1])[1:2]
        index = temp[1],temp[2]
        c_ref[index[1],index[2]] = parsematrix(numrefcn[index[1]], numrefcn[index[2]],
            rawdata[linenum+2:linenum+2+numrefcn[index[1]]])
        global linenum += numrefcn[index[1]] + 1
    end

    # Basis set parameters
    k_cn_l = parseline(Float64, rawdata[end]) # must be parsed here for efficient allocation of matrices
    nshtyp::Int64 = length(k_cn_l) # number of shell types
    shtyp = [Vector{Int64}() for _ = 1:ntypes] 
    nprim = [Vector{Int64}() for _ = 1:ntypes]
    ζ = [Vector{Vector{Float64}}() for _ = 1:ntypes]
    d = [Vector{Vector{Float64}}() for _ = 1:ntypes]
    linenum += 1 # Go to the next line
    global totnsh = 0 # total number of shells
    while true
        index, sh, np = parseline(Int64, rawdata[linenum])
        if np == length(parseline(Float64, rawdata[linenum+1])) == length(parseline(Float64, rawdata[linenum+2]))
            push!(nprim[index],np)
            push!(shtyp[index],sh)
            push!(ζ[index],parseline(Float64, rawdata[linenum+1]))
            push!(d[index],parseline(Float64, rawdata[linenum+2]))
        else
            break
        end
        global linenum += 3
        global totnsh += 1
    end
    k_ab = Dict{Tuple{Int64,Int64,Int64,Int64},Float64}()
    while true
        temp = split(rawdata[linenum])
        if length(temp) == 5
            index1 = Tuple(parse.(Int64, temp[1:4]))
            index2 = (index1[3],index1[4],index1[1],index1[2])
            k_ab[index1] = parse(Float64, temp[5])
            k_ab[index2] = parse(Float64, temp[5])
        else
            break
        end
        global linenum += 1
    end

    nk_ll = parse(Int64, rawdata[linenum]) # number of k_ll parameters
    k_ll = Matrix{Float64}(undef,nshtyp,nshtyp)
    for i in 1:nk_ll
        temp = split(rawdata[linenum+i])
        index = parse.(Int64, temp[1:2])
        k_ll[index[1],index[2]] = parse(Float64, temp[3])
    end
    linenum += nk_ll + 1
    electroneg = parseline(Float64, rawdata[linenum])
    k_en = parse(Float64, rawdata[linenum+1])
    r_cov = parseline(Float64, rawdata[linenum+2]) ./ BORH_TO_Å
    γ = parseline(Float64, rawdata[linenum+3])
    linenum += 3
    std_sh_pop = Matrix{Int64}(undef,ntypes,nshtyp)
    η_al = Matrix{Float64}(undef,ntypes,nshtyp)
    h_al = Matrix{Float64}(undef,ntypes,nshtyp)
    k_poli = Matrix{Float64}(undef,ntypes,nshtyp)
    for i in 1:totnsh
        a, l, temp... = split(rawdata[linenum+i])
        index = parse(Int64, a), parse(Int64, l)
        std_sh_pop[index[1],index[2]] = parse(Int64, temp[1])
        η_al[index[1],index[2]], h_al[index[1],index[2]], k_poli[index[1],index[2]] = parse.(Float64, temp[2:end])
    end
    h_al = h_al ./ AU_TO_EV
    return BORH_TO_Å, AU_TO_EV, indx,
    k_f, α, z_eff, a1, a2, s6, s8, k_cn, k_l, q_a, sc_radii, numrefcn, refcn, c_ref,
    nprim, shtyp, ζ, d, k_ab, nk_ll, k_ll,
    electroneg, k_en, r_cov, γ, std_sh_pop, η_al, h_al, k_poli, k_cn_l
end

"""
Reads data from the parameters file and outputs parsed parameters.

Input: filename::String (default: "parameters/parameters.dat")
Output: parameters::Tuple
"""
function loadparams(filename="$binpath/parameters/parameters_best.dat"::String)
    rawdata = readparams(filename)
    return parseparams(rawdata)
end
end