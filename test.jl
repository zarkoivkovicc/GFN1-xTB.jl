include("modules/io.jl")
include("modules/constants.jl")
include("modules/method.jl")
#import Pkg
#Pkg.add(["ArgParse", "Memoization"],io=devnull)

using ArgParse, Memoization
using LinearAlgebra: norm, eigen, Hermitian, Diagonal
using Printf
using .Parsers, .PeriodicTable, .Methods

natoms, elementsym, coordinates = parsexyz("$binpath/tests/data/secbutylamine.xyz")

# Load parameters
BORH_TO_Å, AU_TO_EV, indx,
k_f, α, z_eff, a1, a2, s6, s8, k_cn, k_l, q_a, sc_radii, numrefcn, refcn, c_ref,
nprim, shtyp, ζ, d, k_ab, nk_ll, k_ll,
electroneg, k_en, r_cov, γ, std_sh_pop, η_a, h_al, k_poli, k_cn_l = loadparams()
id = map( x -> get(indx,x,-1),elementsym) # Each atom gets id of its type
coordinates  ./= BORH_TO_Å # Convert to angstroms
r = pairwise_distance(natoms,coordinates)
basis_fun = gen_basisfun(natoms,id,shtyp) # basis function is represented as (atom,sh_type,orientation)
S = overlapmatrix(basis_fun,id,shtyp,coordinates,ζ,d)
S_sqrt_inv = S_inv_sqrt(S)
H_0 = gen_H_0(basis_fun,natoms,r,id,k_ab,k_ll,electroneg,k_en,r_cov,h_al,k_poli,k_cn_l,sc_radii,k_cn,S)

print(H_0[9,4])
print("Repulsion energy: \n")
print(E_rep(natoms, id, r, α, z_eff, k_f),'\n')
print("Dispersion energy: \n")
print(E_disp(natoms, id , r, a1, a2, s6, s8, k_cn, k_l, q_a, sc_radii, numrefcn, refcn, c_ref))




