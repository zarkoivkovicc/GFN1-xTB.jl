module GFN1_xTB
include("io.jl")
include("method.jl")
using .Parsers, .QuantumChem, .Output
using Printf
export main
function main(out:: String,xyz_file:: String, verbose::Bool,parameters::String,
maxiter::Int64, λ_damp::Float64, charge::Int64)
appendf(out,"THE GFN1-XTB TIGHT BINDING PROGRAM STARTED ")
appendf(out,"Reference paper: A Robust and Accurate Tight-Binding Quantum
Chemical Method for Structures, Vibrational Frequencies, and Noncovalent Interactions of Large
Molecular Systems Parametrized for All spd-Block Elements (Z = 1 - 86)”, S. Grimme, C.
Bannwarth and P. Shushkov, J. Chem. Theory Comput. 2017, 13, 1989 - 2009, https://doi.
org/10.1021/acs.jctc.7b00118.")
# Parse information from parameters file
BORH_TO_Å, AU_TO_EV, indx,
k_f, α, z_eff, a1, a2, s6, s8, k_cn, k_l, q_a, sc_radii, numrefcn, refcn, c_ref,
nprim, shtyp, ζ, d, k_ab, nk_ll, k_ll,
electroneg, k_en, r_cov, Γ, std_sh_pop, η_al, h_al, k_poli, k_cn_l = loadparams("$binpath/parameters/$parameters")
appendf(out,"SUCCESFULLY LOADED PARAMETERS FROM A FILE")
# Parse information from xyz file
natoms, elementsym, coordinates = parsexyz(xyz_file)
appendf(out,"SUCCESFULLY LOADED MOLECULE FROM A FILE")
# Each atom will be assigned an id, the index of the atom type in all paramters
id = map( x -> get(indx,x,-1),elementsym)
coordinates  ./= BORH_TO_Å # Convert to angstroms
r = get_distances(natoms,coordinates)
E_0_rep::Float64 = get_E_rep(natoms,id,r,α,z_eff,k_f)
appendf(out,"SUCCESFULLY CALCULATED ZERO-ORDER REPULSION ENERGY")
if verbose == true
    open(out,"a") do file
        @printf(file,"E_rep: %.8f \n",E_0_rep)
    end
end
E_0_disp::Float64 = get_E_disp(natoms,id,r,a1,a2,s6,s8,k_cn,k_l,q_a,sc_radii,numrefcn,refcn,c_ref)
appendf(out,"SUCCESFULLY CALCULATED ZERO-ORDER DISPERSION ENERGY")
if verbose == true
    open(out,"a") do file
        @printf(file,"E_disp: %.8f \n",E_0_disp)
    end
end
# Generate all basis functions and shells. Check docs of get_basisfun_shells for more info
basis_fun, shells = get_basisfun_shells(natoms,id,shtyp)
nbasisfun::Int64 = length(basis_fun); nshells::Int64 = length(shells)
nelec = get_nelec(std_sh_pop,shells,id) - charge
appendf(out,"SUCCESFULLY GENERATED $nbasisfun BASIS FUNCTIONS AND $nshells SHELLS ")
appendf(out,"TOTAL NUMBER OF ELECTRONS IS $nelec ")
S = get_S(basis_fun,id,shtyp,coordinates,ζ,d)
S_sqrt_inv = Inv_sqrt(S)
if verbose == true
    appendf(out,"SUCCESFULLY CALCULATED OVERLAP MATRIX ")
    appendf(out,S,"Overlap matrix S: ","%9.6f")
    appendf(out,"SUCCESFULLY CALCULATED SQUARE ROOT OF INVERSE OF OVERLAP MATRIX S^-1/2: ")
    appendf(out,S_sqrt_inv,"Overlap matrix S: ","%9.6f")
end
H_0 = get_H_0(basis_fun,natoms,r,id,k_ab,k_ll,electroneg,k_en,r_cov,h_al,k_poli,k_cn_l,sc_radii,k_cn,S)
if verbose == true
    appendf(out,"SUCCESFULLY CALCULATED ZERO-ORDER HAMILTONIAN MATRIX ")
    appendf(out,H_0,"Zero-order Hamiltonian H_0: ","%9.6f")
end
γ_shpairs = get_γ_abllp(shells,r,id,η_al)
if verbose == true
    appendf(out,"SUCCESFULLY CALCULATED MATRIX OF COULOMB-LIKE FACTORS γ ")
    appendf(out,γ_shpairs,"Matrix of Coulomb factors γ: ","%9.6f")
end
F::Matrix{Float64} = zeros(size(S)...)
P::Matrix{Float64} = zeros(size(S)...)
E_orb::Vector{Float64} = zeros(length(basis_fun))
perm = Vector{Int64}(undef,length(basis_fun))
C::Matrix{Float64} = zeros(size(S)...)
atomic_charges:: Vector{Float64} = zeros(natoms)
shell_charges:: Vector{Float64} = zeros(length(shells))
Ee::Float64 = E1::Float64 = E2::Float64 = E3::Float64 = E_prev::Float64 = E::Float64 = 0 
ΔE::Float64 = ΔP::Float64 = Δt::Float64 = ΔT::Float64 = 0
nsteps::Int64 = 0
appendf(out,"SUCCESFULLY CREATED EVERYTHING, ENTERING SCF...!")
open(out,"a") do file
    println(file," N         Ee            E1          E2          E3          ΔE          ΔP         t      t_tot")
    println(file,"---------------------------------------------------------------------------------------------------")
end
for step in 1:maxiter
    Δt = @elapsed begin
        get_F!(F,shell_charges,atomic_charges,γ_shpairs,Γ,H_0,S,basis_fun,id)
        E_orb, C, perm = get_eigen_from_F(F,S_sqrt_inv)
        ΔP = get_P!(P,nelec,C[:,perm])
        new_shell_charges = get_shell_charges(id,std_sh_pop,S,P,shells,basis_fun)
        global new_atomic_charges = get_atomic_charges(natoms,new_shell_charges,shells)
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
        #display(E_orb[perm])
        #display(C[:,perm])
        #display(P)
        #display(atomic_charges)
        E_prev = Ee
        Ee = E1 + E2 + E3
        ΔE = abs(Ee-E_prev)
    end
    ΔT += Δt
    appendf(out,step,Ee,E1,E2,E3,ΔE,ΔP,Δt,ΔT)
    if ΔP <= 1e-4 && ΔE <= 1e-7
        nsteps += step
        break
    end
end
if nsteps != 0
    appendf(out,"SCF SUCCESFULLY CONVERGED AFTER $nsteps ITERATIONS! ")
else
    appendf(out,"WARNING: SCF DID NOT CONVERGE AFTER $nsteps ITERATIONS ")
end
appendf(out,E_orb[perm],"Orbital energies","%9.6f")
appendf(out,C[:,perm],"%9.6f")
appendf(out,new_atomic_charges,"Atomic charges","%9.6f")
E = Ee + E_0_rep + E_0_disp
open(out,"a") do file
    @printf(file,"FINAL SINGLE POINT ENERGY: %.8f \n",E)
    @printf(file,"Energy contributions   Ee: %.8f    E_rep: %.8f    E_disp: %.8f \n",Ee,E_0_rep,E_0_disp)
end
appendf(out,"PROGRAM TERMINATED NORMALLY")
end

end
