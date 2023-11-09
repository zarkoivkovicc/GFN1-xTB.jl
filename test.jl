include("io.jl")
include("constants.jl")
using LinearAlgebra: norm
using .Parsers, .PeriodicTable
using Profile
natoms, elements, coordinates = parsexyz("$binpath/tests/data/methanol.xyz")
elnum = anumber.(elements)

# Load parameters
BORH_TO_Å, EV_TO_AU, 
k_f, α, z_eff, a1, a2, s6, s8, k_cn, k_l, sc_radii, refcn, c_ref,
nprim, ζ, d, k_ab, nk_ll, k_ll,
electroneg, k_en, r_cov, γ, std_sh_pop, η_a, h_al, k_poli, k_cn_l = loadparams()

coordinates  ./= BORH_TO_Å # Convert to angstroms
r = Matrix{Float64}(undef,natoms,natoms) #pairwise_distance
for i in 1:natoms, j in i:natoms
    r[j,i]  = norm(coordinates[j,:]-coordinates[i,:])
end

zero_ord_rep_en ::Float64 = 0
for i in 1:natoms, j in i+1:natoms
        global zero_ord_rep_en +=
        z_eff[elnum[i]]z_eff[elnum[j]]/r[j,i] * ℯ^(-√(α[elnum[i]]α[elnum[j]])r[j,i]^k_f)
end
print(zero_ord_rep_en)

