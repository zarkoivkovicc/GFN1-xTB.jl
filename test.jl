include("modules/io.jl")
include("modules/constants.jl")
include("modules/method.jl")
#import Pkg
#Pkg.add(["ArgParse", "Memoization"],io=devnull)

using ArgParse, Memoization
using LinearAlgebra: norm, eigen, Hermitian, Diagonal, ⋅, Transpose
using Printf
using .Parsers, .PeriodicTable, .Methods
charge::Int64 = 0
λ_damp::Float64 = 0.4
natoms, elementsym, coordinates = parsexyz("$binpath/tests/data/methanol.xyz")

# Load parameters
const BORH_TO_Å, AU_TO_EV, indx,
k_f, α, z_eff, a1, a2, s6, s8, k_cn, k_l, q_a, sc_radii, numrefcn, refcn, c_ref,
nprim, shtyp, ζ, d, k_ab, nk_ll, k_ll,
electroneg, k_en, r_cov, Γ, std_sh_pop, η_al, h_al, k_poli, k_cn_l = loadparams()
id = map( x -> get(indx,x,-1),elementsym) # Each atom gets id of its type
coordinates  ./= BORH_TO_Å # Convert to angstroms
r = get_distances(natoms,coordinates)
basis_fun, shells, shells_bf = get_basisfun_shells(natoms,id,shtyp)
# basis function is represented as (atom,sh_type,orientation,sh_index)
# shell of an atom is represented as (atom,sh_type)
nelec = get_nelec(std_sh_pop,shells,id)
S = get_S(basis_fun,id,shtyp,coordinates,ζ,d)
S_sqrt_inv = Inv_sqrt(S)
H_0 = get_H_0(basis_fun,natoms,r,id,k_ab,k_ll,electroneg,k_en,r_cov,h_al,k_poli,k_cn_l,sc_radii,k_cn,S)
γ_shpairs = get_γ_abllp(shells,r,id,η_al)
F::Matrix{Float64} = zeros(size(S)...)
P::Matrix{Float64} = zeros(size(S)...)
atomic_charges:: Vector{Float64} = zeros(natoms)
shell_charges:: Vector{Float64} = zeros(length(shells))
Ee::Float64 = 0
# take care which charges you use for which energy!
for cycle in 1:50
    print("Cycle number: $cycle \n")
    get_F!(F,shell_charges,atomic_charges,γ_shpairs,Γ,H_0,S,basis_fun,id)
    global E_orb, C, perm = get_eigen_from_F(F,S_sqrt_inv)
    ΔP:: Float64 = get_P!(P,nelec,C[:,perm])
    E1::Float64 = sum(P .* H_0)
    E2::Float64 = 0
    for i in 1:length(shells), j in 1:length(shells)
        E2 += shell_charges[i]*shell_charges[j]*γ_shpairs[i,j]
    end
    E2 /= 2
    E3::Float64 = 0
    for i in 1:natoms
        E3 += atomic_charges[i]^3*Γ[id[i]]
    end
    E3 /= 3
    E_prev::Float64 = Ee
    global Ee = E1 + E2 + E3
    ΔE::Float64 = abs(Ee-E_prev)
    get_shell_charges!(shell_charges,id,std_sh_pop,S,P,shells,basis_fun)
    get_atomic_charges!(atomic_charges,shell_charges,shells)
    print("E1 = $E1, E2 = $E2, E3 = $E3, Ee = $Ee \n")
    print("ΔP = $ΔP, ΔE = $ΔE \n")
    if ΔP <= 1e-4 && ΔE <= 1e-7
        break
    end
end
#display(S)
#display(H_0)
#display(S_sqrt_inv)
#display(F)
#display(P)
#display(atomic_charges)
#display(E_orb[perm])
#display(C[:,perm])
print("Repulsion energy: \n")
print(get_E_rep(natoms, id, r, α, z_eff, k_f),'\n')
print("Dispersion energy: \n")
print(get_E_disp(natoms, id , r, a1, a2, s6, s8, k_cn, k_l, q_a, sc_radii, numrefcn, refcn, c_ref))