include("modules/io.jl")
include("modules/constants.jl")
include("modules/method.jl")
using LinearAlgebra: norm
using .Parsers, .PeriodicTable, .Methods

natoms, elementsym, coordinates = parsexyz("$binpath/tests/data/methanol.xyz")

# Load parameters
BORH_TO_Å, EV_TO_AU, indx,
k_f, α, z_eff, a1, a2, s6, s8, k_cn, k_l, q_a, sc_radii, numrefcn, refcn, c_ref,
nprim, shtyp, ζ, d, k_ab, nk_ll, k_ll,
electroneg, k_en, r_cov, γ, std_sh_pop, η_a, h_al, k_poli, k_cn_l = loadparams()
id = map( x -> get(indx,x,-1),elementsym) # Each atom gets id of its type
coordinates  ./= BORH_TO_Å # Convert to angstroms
r = pairwise_distance(natoms,coordinates)
basis_fun = gen_basisfun(natoms,id,shtyp)
print(basis_fun)
print("Repulsion energy: \n")
print(E_rep(natoms, id, r, α, z_eff, k_f),'\n')
print("Dispersion energy: \n")
print(E_disp(natoms, id , r, a1, a2, s6, s8, k_cn, k_l, q_a, sc_radii, numrefcn, refcn, c_ref))




