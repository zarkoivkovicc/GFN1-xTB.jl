module Methods
using LinearAlgebra: norm
export pairwise_distance, E_rep, E_disp, gen_basisfun
function pairwise_distance(natoms :: Int64, coords :: Matrix{Float64})
    r = Matrix{Float64}(undef,natoms,natoms) #pairwise_distance
    for i in 1:natoms, j in i:natoms
        r[j,i]  = norm(coords[:,j]-coords[:,i])
    end
    return r
end
function E_rep(natoms :: Int64, id :: Vector{Int64}, r :: Matrix{Float64},
    α :: Vector{Float64}, z_eff :: Vector{Float64}, k_f :: Float64)
    E_rep ::Float64 = 0
    for i in 1:natoms, j in i+1:natoms
        E_rep +=
        z_eff[id[i]]z_eff[id[j]]/r[j,i] * ℯ^(-√(α[id[i]]α[id[j]])r[j,i]^k_f)
    end
    return E_rep
end
function E_disp(natoms :: Int64, id :: Vector{Int64}, r :: Matrix{Float64},
    a1 :: Float64, a2 :: Float64, s6 :: Float64, s8 :: Float64, k_cn :: Float64, k_l :: Float64,
    q_a :: Vector{Float64}, sc_radii :: Vector{Float64}, numrefcn :: Vector{Int64},
    refcn :: Matrix{Float64}, c_ref :: Matrix{Matrix{Float64}})
    function f6_damp(r::Float64,a1::Float64,a2::Float64,q_a::Float64,q_b::Float64)
        return r^6/
        (r^6+(a1*(9q_a*q_b)^0.25+a2)^6)
    end
    function f8_damp(r::Float64,a1::Float64,a2::Float64,q_a::Float64,q_b::Float64)
        return r^8/
        (r^8+(a1*(9q_a*q_b)^0.25+a2)^8)
    end
    function cn(a::Int64,r::Matrix{Float64})
        res::Float64 = 0
        for i in 1:a-1
            res += (1+ ℯ^(-k_cn*((sc_radii[id[a]]+sc_radii[id[i]])/r[a,i]-1)))^-1
        end
        for i in a+1:natoms
            res += (1+ ℯ^(-k_cn*((sc_radii[id[a]]+sc_radii[id[i]])/r[i,a]-1)))^-1
        end
        return res
    end
    function l(k_l::Float64,cn_a::Float64,cn_a_ref::Float64)
        return ℯ^(-k_l*(cn_a - cn_a_ref)^2)
    end
    function c6(a::Int64,b::Int64,cn_a::Float64,cn_b, refcn :: Matrix{Float64})
        nom::Float64, denom::Float64 = 0, 0
        if a > b
            a, b = b, a
        end
        for i in 1:numrefcn[a], j in 1:numrefcn[b] # HERE TAKE CARE THE ORDER!!!
            l_ij = l(k_l,cn_a,refcn[i,a])*l(k_l,cn_b,refcn[j,b])
            nom += c_ref[a,b][i,j]*l_ij
            denom += l_ij
        end
        return nom/denom
    end
    E_disp :: Float64 = 0
    for i in 1:natoms, j in i+1:natoms
        E_disp += -s6*c6(id[j],id[i],cn(j,r),cn(i,r),refcn)/r[j,i]^6*f6_damp(r[j,i],a1,a2,q_a[id[j]],q_a[id[i]])-
        s8*3*√(q_a[id[j]]q_a[id[i]])*c6(id[j],id[i],cn(j,r),cn(i,r),refcn)/r[j,i]^8*f8_damp(r[j,i],a1,a2,q_a[id[j]],q_a[id[i]])
    end
    return E_disp
end
function gen_basisfun(natoms:: Int64, id :: Vector{Int64}, shtyp:: Vector{Vector{Int64}})
    basis_functions = Vector{Vector{Int64}}()
    for i in 1:natoms, sh in shtyp[id[i]]
        if sh == 1 || sh == 2
            push!(basis_functions,[i,sh,0])
        elseif sh ==3
            push!(basis_functions,[i,sh,1])
            push!(basis_functions,[i,sh,2])
            push!(basis_functions,[i,sh,3])
        end
    end
    return basis_functions
end
function nbasisf(id::Vector{Int64},shtyp :: Vector{Vector{Int64}})
    n = 0
    nbasis_per_type :: Vector{Int64} = [1,1,3] # 1 s, 1 s', 3 p orbitals
    for i in id, j in shtyp[i]
        n += nbasis_per_type[j]
    end
    return n
end
end