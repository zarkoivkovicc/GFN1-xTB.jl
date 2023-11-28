module GFN1_xTB
include("io.jl")
include("method.jl")
using .Parsers, .QuantumChem
using Printf
export main
function main(output:: String,xyz_file:: String, verbose::Bool,parameters::String,
    maxiter::Int64, λ_damp::Float64, charge::Int64)
    # Parse information from xyz file
    natoms, elementsym, coordinates = parsexyz(xyz_file)
    # Parse information from parameters file
    BORH_TO_Å, AU_TO_EV, indx,
    k_f, α, z_eff, a1, a2, s6, s8, k_cn, k_l, q_a, sc_radii, numrefcn, refcn, c_ref,
    nprim, shtyp, ζ, d, k_ab, nk_ll, k_ll,
    electroneg, k_en, r_cov, Γ, std_sh_pop, η_al, h_al, k_poli, k_cn_l = loadparams("$binpath/parameters/$parameters")
    # Each atom will be assigned an id, the index of the atom type in all paramters
    id = map( x -> get(indx,x,-1),elementsym)
    coordinates  ./= BORH_TO_Å # Convert to angstroms
    r = get_distances(natoms,coordinates)
    E_0_rep::Float64 = get_E_rep(natoms,id,r,α,z_eff,k_f)
    E_0_disp::Float64 = get_E_disp(natoms,id,r,a1,a2,s6,s8,k_cn,k_l,q_a,sc_radii,numrefcn,refcn,c_ref)
    # Generate all basis functions and shells. Check docs of get_basisfun_shells for more info
    basis_fun, shells = get_basisfun_shells(natoms,id,shtyp)
    nelec = get_nelec(std_sh_pop,shells,id) - charge
    S = get_S(basis_fun,id,shtyp,coordinates,ζ,d)
    S_sqrt_inv = Inv_sqrt(S)
    H_0 = get_H_0(basis_fun,natoms,r,id,k_ab,k_ll,electroneg,k_en,r_cov,h_al,k_poli,k_cn_l,sc_radii,k_cn,S)
    γ_shpairs = get_γ_abllp(shells,r,id,η_al)
    F::Matrix{Float64} = zeros(size(S)...)
    P::Matrix{Float64} = zeros(size(S)...)
    atomic_charges:: Vector{Float64} = zeros(natoms)
    shell_charges:: Vector{Float64} = zeros(length(shells))
    Ee::Float64 = E1::Float64 = E2::Float64 = E3::Float64 = E_prev::Float64 = ΔE::Float64 = 0
    for _ in 1:maxiter
        get_F!(F,shell_charges,atomic_charges,γ_shpairs,Γ,H_0,S,basis_fun,id)
        E_orb, C, perm = get_eigen_from_F(F,S_sqrt_inv)
        ΔP = get_P!(P,nelec,C[:,perm])
        new_shell_charges = get_shell_charges(id,std_sh_pop,S,P,shells,basis_fun)
        new_atomic_charges = get_atomic_charges(natoms,new_shell_charges,shells)
        E1 = sum(P .* H_0)
        E2 = 0
        for i in 1:length(shells), j in 1:length(shells)
            E2 += new_shell_charges[i]*new_shell_charges[j]*γ_shpairs[i,j]
        end
        E2 /= 2
        E3 = 0
        for i in 1:natoms
            E3 += new_atomic_charges[i]^3*Γ[id[i]]
        end
        E3 /= 3
        damp_charges!(atomic_charges,shell_charges,new_atomic_charges,new_shell_charges,λ_damp,1e-3)
        E_prev = Ee
        Ee = E1 + E2 + E3
        ΔE = abs(Ee-E_prev)
        @printf("Ee = %.8f, E1 = %.6f, E2 = %.6f, E3 = %.6f \n",Ee,E1,E2,E3)
        @printf("ΔP = %.2e, ΔE = %.2e \n",ΔP,ΔE)
        if ΔP <= 1e-4 && ΔE <= 1e-7
            break
        end
    end
end
end
