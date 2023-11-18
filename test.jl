include("modules/io.jl")
include("modules/constants.jl")
include("modules/method.jl")
#import Pkg
#Pkg.add(["ArgParse", "Memoization"],io=devnull)

using ArgParse, Memoization
using LinearAlgebra: norm, eigen, Hermitian, Diagonal, ⋅, Transpose
using Printf
using .Parsers, .PeriodicTable, .Methods

natoms, elementsym, coordinates = parsexyz("$binpath/tests/data/methanol.xyz")

# Load parameters
BORH_TO_Å, AU_TO_EV, indx,
k_f, α, z_eff, a1, a2, s6, s8, k_cn, k_l, q_a, sc_radii, numrefcn, refcn, c_ref,
nprim, shtyp, ζ, d, k_ab, nk_ll, k_ll,
electroneg, k_en, r_cov, Γ, std_sh_pop, η_al, h_al, k_poli, k_cn_l = loadparams()
id = map( x -> get(indx,x,-1),elementsym) # Each atom gets id of its type
coordinates  ./= BORH_TO_Å # Convert to angstroms
r = get_distances(natoms,coordinates)
basis_fun, shells = get_basisfun_shells(natoms,id,shtyp)
# basis function is represented as (atom,sh_type,orientation,sh_index)
# shell of an atom is represented as (atom,sh_type)
S = get_S(basis_fun,id,shtyp,coordinates,ζ,d)
S_sqrt_inv = Inv_sqrt(S)
H_0 = get_H_0(basis_fun,natoms,r,id,k_ab,k_ll,electroneg,k_en,r_cov,h_al,k_poli,k_cn_l,sc_radii,k_cn,S)
γ_shpairs = get_γ_abllp(shells,r,id,η_al)
C = Matrix{Float64}(undef,size(S)...)
F = Matrix{Float64}(undef,size(S)...)
Fp = Matrix{Float64}(undef,size(S)...)
atomic_charges:: Vector{Float64} = zeros(natoms)
shell_charges:: Vector{Float64} = zeros(length(shells))
get_F!(F,shell_charges,atomic_charges,γ_shpairs,Γ,H_0,S,basis_fun,id)
Fp = S_sqrt_inv'*H_0*S_sqrt_inv


#print("Repulsion energy: \n")
#print(E_rep(natoms, id, r, α, z_eff, k_f),'\n')
#print("Dispersion energy: \n")
#print(E_disp(natoms, id , r, a1, a2, s6, s8, k_cn, k_l, q_a, sc_radii, numrefcn, refcn, c_ref))




