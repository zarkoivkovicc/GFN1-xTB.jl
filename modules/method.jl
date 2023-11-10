module Methods
export zero_ord_rep_energy, E_disp
function E_rep(natoms :: Int64, elnum :: Vector{Int64}, r :: Matrix{Float64},
    α :: Dict{Int64,Float64}, z_eff :: Dict{Int64,Float64}, k_f :: Float64)
    E_rep ::Float64 = 0
    for i in 1:natoms, j in i+1:natoms
        E_rep +=
        z_eff[elnum[i]]z_eff[elnum[j]]/r[j,i] * ℯ^(-√(α[elnum[i]]α[elnum[j]])r[j,i]^k_f)
    end
    return E_rep
end
function E_disp(natoms :: Int64, elnum :: Vector{Int64}, r :: Matrix{Float64},
    a1 :: Float64, a2 :: Float64, s6 :: Float64, s8 :: Float64, k_cn :: Float64, k_l :: Float64,
    q_a :: Dict{Int64,Float64}, sc_radii :: Dict{Int64,Float64}, numrefcn :: Dict{Int64,Int64},
    refcn :: Dict{Int64,Array{Float64}}, c_ref :: Dict{Tuple{Int64,Int64},Matrix{Float64}})

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
            res += (1+ ℯ^(-k_cn*((sc_radii[elnum[a]]+sc_radii[elnum[i]])/r[a,i]-1)))^-1
        end
        for i in a+1:natoms
            res += (1+ ℯ^(-k_cn*((sc_radii[elnum[a]]+sc_radii[elnum[i]])/r[i,a]-1)))^-1
        end
        return res
    end
    function l(k_l::Float64,cn_a::Float64,cn_a_ref::Float64)
        return ℯ^(-k_l*(cn_a - cn_a_ref)^2)
    end
    function c6(a::Int64,b::Int64,cn_a::Float64,cn_b)
        nom::Float64, denom::Float64 = 0, 0
        for i in 1:numrefcn[a], j in 1:numrefcn[b]
            l_ij = l(k_l,cn_a,refcn[a][i])*l(k_l,cn_b,refcn[b][j])
            nom += c_ref[(a,b)][i,j]*l_ij
            denom += l_ij
        end
        return nom/denom
    end
    E_disp :: Float64 = 0
    for i in 1:natoms, j in i+1:natoms
        E_disp += -s6*c6(elnum[j],elnum[i],cn(j,r),cn(i,r))/r[j,i]^6*f6_damp(r[j,i],a1,a2,q_a[elnum[j]],q_a[elnum[i]])-
        s8*3*√(q_a[elnum[j]]q_a[elnum[i]])*c6(elnum[j],elnum[i],cn(j,r),cn(i,r))/r[j,i]^8*f8_damp(r[j,i],a1,a2,q_a[elnum[j]],q_a[elnum[i]])
    end
    return E_disp
end
end