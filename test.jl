include("io.jl")
include("constants.jl")
using LinearAlgebra: norm
using .Parsers, .PeriodicTable
using Profile
natoms, elementsym, coordinates = parsexyz("$binpath/tests/data/methanol.xyz")
elnum = anumber.(elementsym)

# Load parameters
BORH_TO_Å, EV_TO_AU, 
k_f, α, z_eff, a1, a2, s6, s8, k_cn, k_l, q_a, sc_radii, numrefcn, refcn, c_ref,
nprim, ζ, d, k_ab, nk_ll, k_ll,
electroneg, k_en, r_cov, γ, std_sh_pop, η_a, h_al, k_poli, k_cn_l = loadparams()

coordinates  ./= BORH_TO_Å # Convert to angstroms
r = Matrix{Float64}(undef,natoms,natoms) #pairwise_distance
for i in 1:natoms, j in i:natoms
    r[j,i]  = norm(coordinates[:,j]-coordinates[:,i])
end
# zero_ord_energy
E_0 ::Float64 = 0
for i in 1:natoms, j in i+1:natoms
        global E_0 +=
        z_eff[elnum[i]]z_eff[elnum[j]]/r[j,i] * ℯ^(-√(α[elnum[i]]α[elnum[j]])r[j,i]^k_f)
end
# zero_ord_disp_energy
E_disp::Float64 = 0
function f6_damp(r::Float64,a1::Float64,a2::Float64,q_a::Float64,q_b::Float64)
    return r^6/
    (r^6+(a1*(9q_a*q_b)^0.25+a2)^6)
end
function f8_damp(r::Float64,a1::Float64,a2::Float64,q_a::Float64,q_b::Float64)
    return r^8/
    (r^8+(a1*(9q_a*q_b)^0.25+a2)^8)
end
function cn(a::Int64,r::Matrix{Float64})
    cn::Float64 = 0
    for i in 1:a-1
        cn += (1+ ℯ^(-k_cn*((sc_radii[elnum[a]]+sc_radii[elnum[i]])/r[a,i]-1)))^-1
    end
    for i in a+1:natoms
        cn += (1+ ℯ^(-k_cn*((sc_radii[elnum[a]]+sc_radii[elnum[i]])/r[i,a]-1)))^-1
    end
    return cn
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
for i in 1:natoms, j in i+1:natoms
    global E_disp += -s6*c6(elnum[j],elnum[i],cn(j,r),cn(i,r))/r[j,i]^6*f6_damp(r[j,i],a1,a2,q_a[elnum[j]],q_a[elnum[i]])-
    s8*3*√(q_a[elnum[j]]q_a[elnum[i]])*c6(elnum[j],elnum[i],cn(j,r),cn(i,r))/r[j,i]^8*f8_damp(r[j,i],a1,a2,q_a[elnum[j]],q_a[elnum[i]])
end
print(E_disp)


