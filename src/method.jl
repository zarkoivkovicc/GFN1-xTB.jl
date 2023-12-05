module QuantumChem
using LinearAlgebra: norm, Hermitian, Diagonal, eigen, ⋅
using  Memoization: @memoize
using ..Parsers: Parameters
export get_distances, get_E_rep, get_E_disp, get_basisfun_shells, get_S, get_eigen_from_F,
get_nelec, get_P!,Inv_sqrt, get_H_0, get_γ_abllp, get_F!, get_shell_charges,
get_atomic_charges, damp_charges!
"""
    get_distances(natoms::Int64, coords::Matrix{Float64})
    
    Returns the pairwise distance matrix. Only top upper-right part of the matrix is actually computed.

"""
function get_distances(natoms::Int64, coords::Matrix{Float64})
    r = Matrix{Float64}(undef,natoms,natoms) #pairwise_distance
    for i in 1:natoms, j in i:natoms
        r[j,i]  = norm(coords[:,j]-coords[:,i])
    end
    return r
end
"""
    get_nelec(shells::Vector{Vector{Int64}},id::Vector{Int64},params::Parameters)

    Returns the number of electrons considered in this tight-binding model assuming neutrality.
"""
function get_nelec(shells::Vector{Vector{Int64}},id::Vector{Int64},params::Parameters)
    (; std_sh_pop) = params
    nelec::Int64 = 0
    for sh in shells
        nelec += std_sh_pop[id[sh[1]],sh[2]]
    end
    return nelec
end
"""
    get_E_rep(natoms::Int64, id::Vector{Int64}, r::Matrix{Float64}, params::Parameters)

    Returns zero-order repulsion energy.
"""
function get_E_rep(natoms::Int64, id::Vector{Int64}, r::Matrix{Float64}, params::Parameters)
    (; α, z_eff, k_f) = params
    E_rep::Float64 = 0
    for i in 1:natoms, j in i+1:natoms
        E_rep +=
        z_eff[id[i]]z_eff[id[j]]/r[j,i] * ℯ^(-√(α[id[i]]α[id[j]])r[j,i]^k_f)
    end
    return E_rep
end
"""
    cn(a::Int64,natoms::Int64, id::Vector{Int64}, r::Matrix{Float64}, params::Parameters)

    Returns the coordination number of a-th atom in the list of atoms.
"""
function cn(a::Int64,natoms::Int64, id::Vector{Int64}, r::Matrix{Float64}, params::Parameters)
    (; sc_radii, k_cn) = params
    res::Float64 = 0
    for i in 1:a-1
        res += (1+ ℯ^(-k_cn*((sc_radii[id[a]]+sc_radii[id[i]])/r[a,i]-1)))^-1
    end
    for i in a+1:natoms
        res += (1+ ℯ^(-k_cn*((sc_radii[id[a]]+sc_radii[id[i]])/r[i,a]-1)))^-1
    end
    return res
end
"""
    get_E_disp(natoms::Int64, id::Vector{Int64}, r::Matrix{Float64},params::Parameters)

    Returns the dispersion energy.
"""
function get_E_disp(natoms::Int64, id::Vector{Int64}, r::Matrix{Float64},params::Parameters)
    (; a1, a2, s6, s8, k_l, q_a, numrefcn, refcn, c_ref) = params
    function f6_damp(r::Float64,a1::Float64,a2::Float64,q_a::Float64,q_b::Float64)
        return r^6/
        (r^6+(a1*(9q_a*q_b)^0.25+a2)^6)
    end
    function f8_damp(r::Float64,a1::Float64,a2::Float64,q_a::Float64,q_b::Float64)
        return r^8/
        (r^8+(a1*(9q_a*q_b)^0.25+a2)^8)
    end
    function l(k_l::Float64,cn_a::Float64,cn_a_ref::Float64)
        return ℯ^(-k_l*(cn_a - cn_a_ref)^2)
    end
    function c6(a::Int64,b::Int64,cn_a::Float64,cn_b, refcn::Matrix{Float64})
        nom::Float64, denom::Float64 = 0, 0
        if a > b
            a, b = b, a
        end
        for i in 1:numrefcn[a], j in 1:numrefcn[b]
            l_ij = l(k_l,cn_a,refcn[i,a])*l(k_l,cn_b,refcn[j,b])
            nom += c_ref[a,b][i,j]*l_ij
            denom += l_ij
        end
        return nom/denom
    end
    E_disp::Float64 = 0
    for i in 1:natoms, j in i+1:natoms
        E_disp += -s6*c6(id[j],id[i],cn(j,natoms,id,r,params),cn(i,natoms,id,r,params),refcn)/
        r[j,i]^6*f6_damp(r[j,i],a1,a2,q_a[id[j]],q_a[id[i]])
        E_disp += -s8*3*√(q_a[id[j]]q_a[id[i]])*c6(id[j],id[i],cn(j,natoms,id,r,params),cn(i,natoms,id,r,params),refcn)/
        r[j,i]^8*f8_damp(r[j,i],a1,a2,q_a[id[j]],q_a[id[i]])
    end
    return E_disp
end
"""
    get_basisfun_shells(natoms::Int64, id::Vector{Int64}, params::Parameters)
    
    Returns the list of basis functions and list of shells.
    Each basis function is three component vector:
        1. number of atom where the basis function is centered
        2. type of shell the basis function belongs to
        3. type and orientation of basis function (0:s 1:px 2:py 3:pz)
        4. number of shell the basis function belongs to
    Each shell is represented as two component vector:
        1. number of atom the shell is centered
        2. type of shell
"""
function get_basisfun_shells(natoms::Int64, id::Vector{Int64}, params::Parameters)
    (; shtyp) = params
    basis_functions = Vector{Vector{Int64}}()
    shells = Vector{Vector{Int64}}()
    indxsh::Int64 = 1
    indxbf::Int64 = 1
    for i in 1:natoms
        for sh in shtyp[id[i]]
            if sh == 1 || sh == 2
                push!(basis_functions,[i,sh,0,indxsh])
                indxbf += 1
            elseif sh ==3
                push!(basis_functions,[i,sh,1,indxsh])
                push!(basis_functions,[i,sh,2,indxsh])
                push!(basis_functions,[i,sh,3,indxsh])
                indxbf += 3
            end
            push!(shells,[i,sh])
            indxsh +=1
        end
    end
    return basis_functions, shells
end
"""
    get_S(basis_fun::Vector{Vector{Int64}},id::Vector{Int64},coord::Matrix{Float64},params::Parameters)

    Returns the overlap matrix.

"""
function get_S(basis_fun::Vector{Vector{Int64}},id::Vector{Int64},coord::Matrix{Float64},params::Parameters)

    """
    ∫_bf_bf(bf1::Vector{Int64},bf2::Vector{Int64},id::Vector{Int64},coord::Matrix{Float64},params::Parameters)

        Returns the integral between two basis functions.
    """
function ∫_bf_bf(bf1::Vector{Int64},bf2::Vector{Int64},id::Vector{Int64},coord::Matrix{Float64},params::Parameters)
    (; shtyp, ζ, d) = params
    """
        ∫s_s(ζ::Float64, ξ::Float64,r1::Vector{Float64},r2::Vector{Float64})

        Returns the integral between two primites of s type. This function is cached.
    """
    @memoize function ∫s_s(ζ::Float64, ξ::Float64,r1::Vector{Float64},r2::Vector{Float64})
                return ℯ^(-ξ*norm(r1-r2)^2)*(π/ζ)^1.5
    end
    """
        ∫prim_prim(ζ1::Float64, ζ2::Float64,r1::Vector{Float64},r2::Vector{Float64},
            typ1::Int64, typ2::Int64)

        Returns the integral between two primitves.
    """
    function ∫prim_prim(ζ1::Float64, ζ2::Float64,r1::Vector{Float64},r2::Vector{Float64},
                        typ1::Int64, typ2::Int64)
                    ζ = ζ1 + ζ2; ξ = ζ1*ζ2/(ζ); r_p = (ζ1*r1 + ζ2*r2)/ζ
                    s_s = ∫s_s(ζ, ξ, r1, r2)
                    if typ1 == typ2 == 0
                        return s_s
                    elseif typ1 == 0
                        return (r_p-r2)[typ2]*s_s
                    elseif typ2 == 0
                        return (r_p-r1)[typ1]*s_s
                    elseif typ1 == typ2
                        return (r_p-r1)[typ1]*(r_p-r2)[typ2]*s_s + s_s/(2*ζ)
                    else
                        return (r_p-r1)[typ1]*(r_p-r2)[typ2]*s_s
                    end
                end
        integral::Float64 = 0
        indx1 = findfirst(x -> x == bf1[2],shtyp[id[bf1[1]]])
        indx2 = findfirst(x -> x == bf2[2],shtyp[id[bf2[1]]])
        ζ1_ = ζ[id[bf1[1]]][indx1]; ζ2_ = ζ[id[bf2[1]]][indx2]
        d1_ = d[id[bf1[1]]][indx1]; d2_ = d[id[bf2[1]]][indx2]
        for (ζ1,d1) in zip(ζ1_,d1_), (ζ2,d2) in zip(ζ2_,d2_)
            integral += d1*d2*∫prim_prim(ζ1,ζ2,coord[:,bf1[1]],coord[:,bf2[1]],bf1[3],bf2[3])
        end
        return integral
    end
    dim = length(basis_fun)
    S = Matrix{Float64}(undef,dim,dim)
    for i in 1:dim, j in i:dim
    S[j,i] = S[i,j] = ∫_bf_bf(basis_fun[i],basis_fun[j],id,coord,params)
    end
    return Hermitian(S)
end
"""
    Inv_sqrt(S::Hermitian{Float64})
    Returns S^-1/2 of the matrix S
"""
function Inv_sqrt(S::Hermitian{Float64})
    λ, L = eigen(S)
    return L * Diagonal(1 ./ sqrt.(λ)) * L'
end

"""
    get_H_0(basis_fun::Vector{Vector{Int64}},natoms::Int64, r::Matrix{Float64}, id::Vector{Int64},
    S::Hermitian{Float64}, params::Parameters)

    Returns zero-order hamlitonian.
"""
function get_H_0(basis_fun::Vector{Vector{Int64}},natoms::Int64, r::Matrix{Float64}, id::Vector{Int64},
    S::Hermitian{Float64}, params::Parameters)
    """
    zeroh(bf1::Vector{Int64},bf2::Vector{Int64},natoms::Int64, r::Matrix{Float64}, id::Vector{Int64}, S::Float64,
    params::Parameters)

        Returns the zero-order hamlitonian element between two basis functions.
    """
    function zeroh(bf1::Vector{Int64},bf2::Vector{Int64},natoms::Int64, r::Matrix{Float64}, id::Vector{Int64}, S::Float64,
        params::Parameters)
        (; k_ab, k_ll, electroneg, k_en, r_cov, h_al, k_poli, k_cn_l, sc_radii, k_cn) = params
        function Π(bf1::Vector{Int64},bf2::Vector{Int64},r::Matrix{Float64}, id::Vector{Int64},
            k_poli::Matrix{Float64},r_cov::Vector{Float64})
            return (1+k_poli[id[bf1[1]],bf1[2]]*(r[bf2[1],bf1[1]]/(r_cov[id[bf1[1]]]+r_cov[id[bf2[1]]]))^0.5)*
            (1+k_poli[id[bf2[1]],bf2[2]]*(r[bf2[1],bf1[1]]/(r_cov[id[bf1[1]]]+r_cov[id[bf2[1]]]))^0.5)
        end

        K = get(k_ab,(id[bf1[1]],bf1[2],id[bf2[1]],bf2[2]),1.0)
        indx1, indx2 = min(bf1[2],bf2[2]), max(bf1[2],bf2[2])
        k_llp = k_ll[indx1,indx2]
        h_al1 = h_al[id[bf1[1]],bf1[2]]*(1+k_cn_l[bf1[2]]*cn(bf1[1],natoms,id,r,params))
        h_al2 = h_al[id[bf2[1]],bf2[2]]*(1+k_cn_l[bf2[2]]*cn(bf2[1],natoms,id,r,params))
        if bf1[2] == 2 || bf2[2] == 2
            off_diag_scal::Float64 = 1
        else
            off_diag_scal = 1.0 + k_en*(electroneg[id[bf1[1]]] - electroneg[id[bf2[1]]])^2
        end
        if bf1[1] == bf2[1]
            if bf1 == bf2
                return h_al1
            else
                return 0
            end
        else
            return K*k_llp*(h_al1+h_al2)/2*S*off_diag_scal*Π(bf1,bf2,r,id,k_poli,r_cov)
        end
    end
    dim = length(basis_fun)
    H_0 = Matrix{Float64}(undef,dim,dim)
    for i in 1:dim, j in i:dim
        H_0[j,i] = H_0[i,j] = zeroh(basis_fun[i],basis_fun[j],natoms,r,id,S[i,j],params)
    end
    return Hermitian(H_0)
end
"""
    get_γ_abllp(shells::Vector{Vector{Int64}}, r::Matrix{Float64},id::Vector{Int64},params::Parameters)

    Returns the matrix of coulomb-like law prefactors γ
"""
function get_γ_abllp(shells::Vector{Vector{Int64}}, r::Matrix{Float64},id::Vector{Int64},params::Parameters)
    (; η_al) = params
    function γ_abllp(sh1::Vector{Int64}, sh2::Vector{Int64}, r::Matrix{Float64},id::Vector{Int64},
        η_al::Matrix{Float64},k_g::Float64 = 2.0)
        η1 = η_al[id[sh1[1]],sh1[2]]; η2 = η_al[id[sh2[1]],sh2[2]]
        return (r[sh2[1],sh1[1]]^k_g+(0.5*(1/η1 + 1/η2))^k_g)^(-1/k_g)
    end
    dim = length(shells)
    res = Matrix{Float64}(undef,dim,dim)
    for i in 1:dim, j in i:dim
        res[j,i] = res[i,j] = γ_abllp(shells[i],shells[j],r,id,η_al)
    end
    return res
end

"""
    get_F!(F::Matrix{Float64}, shell_charges::Vector{Float64}, atomic_charges::Vector{Float64}, γ::Matrix{Float64},
    H_0::Hermitian{Float64}, S::Hermitian{Float64}, basis_fun::Vector{Vector{Int64}}, id::Vector{Int64}, params::Parameters)

    Calculates (in-place) fock matrix F.
"""
function get_F!(F::Matrix{Float64}, shell_charges::Vector{Float64}, atomic_charges::Vector{Float64}, γ::Matrix{Float64},
    H_0::Hermitian{Float64}, S::Hermitian{Float64}, basis_fun::Vector{Vector{Int64}}, id::Vector{Int64}, params::Parameters)
    (; Γ) = params
    """
    F_(μ::Int64, ν::Int64, basis_fun::Vector{Vector{Int64}},shell_charges::Vector{Float64},
        atomic_charges::Vector{Float64},γ::Matrix{Float64},Γ::Vector{Float64},
        H_0::Hermitian{Float64},S::Hermitian{Float64},id::Vector{Int64})
        
        Return fock matrix element F[μ,ν], where μ and ν are basis function numbers
    """
function F_(μ::Int64, ν::Int64, basis_fun::Vector{Vector{Int64}},shell_charges::Vector{Float64},
        atomic_charges::Vector{Float64},γ::Matrix{Float64},Γ::Vector{Float64},
        H_0::Hermitian{Float64},S::Hermitian{Float64},id::Vector{Int64})
        function δϵ_al(sh_indx::Int64,γ::Matrix{Float64}, shell_charges::Vector{Float64})
            return γ[sh_indx,:]⋅shell_charges
        end
        function δΕ_al(a::Int64,id::Vector{Int64},Γ::Vector{Float64}, atomic_charges::Vector{Float64})
            return Γ[id[a]]*atomic_charges[a]^2
        end
        return H_0[μ,ν] - 0.5*S[μ,ν]*(δϵ_al(basis_fun[μ][4],γ,shell_charges)+
                                      δϵ_al(basis_fun[ν][4],γ,shell_charges)+
                                      δΕ_al(basis_fun[μ][1],id,Γ,atomic_charges)+
                                      δΕ_al(basis_fun[ν][1],id,Γ,atomic_charges))
    end
    dim = size(S)
    for i in 1:dim[1], j in 1:dim[1]
        F[i,j] = F_(i,j,basis_fun,shell_charges,atomic_charges,γ,Γ,H_0,S,id)
    end
end

"""
    get_P!(P::Matrix{Float64}, nelec::Int64, C::Matrix{Float64})::Float64

    Calculates desity matrix from coefficient matrix and returns the change in norm of density matrix.
"""
function get_P!(P::Matrix{Float64}, nelec::Int64, C::Matrix{Float64})::Float64
    dim = size(C,1)
    n_occ::Int64 = nelec/2
    norm_old::Float64 = norm(P)
    for i in 1:dim, j in i:dim
        P[i,j] = 0
        P[i,j] += C[i,1:n_occ]⋅C[j,1:n_occ]
        P[j,i] = P[i,j]
    end
    P .*= 2
    norm_new::Float64 = norm(P)
    return abs(norm_new-norm_old)/length(P)
end

"""
    get_eigen_from_F(F::Matrix{Float64},S_sqrt_inv::Matrix{Float64})

    Returns oribal energies, coefficient matrix and perumatation that sorts
    the eigenvectors and eigenvalues in the order of raising energies.
"""
function get_eigen_from_F(F::Matrix{Float64},S_sqrt_inv::Matrix{Float64})

    Fp::Matrix{Float64} = S_sqrt_inv' * F * S_sqrt_inv
    E_orb::Vector{Float64}, Cp::Matrix{Float64} = eigen(Fp,sortby=nothing)
    C::Matrix{Float64} = S_sqrt_inv*Cp
    return E_orb, C, sortperm(E_orb)
end

"""
    get_shell_charges(id::Vector{Int64},S::Hermitian{Float64},P::Matrix{Float64},shells::Vector{Vector{Int64}},
    basis_fun::Vector{Vector{Int64}}, params::Parameters)

    Returns the shell charges vector calculated from the density matrix.
"""
function get_shell_charges(id::Vector{Int64},S::Hermitian{Float64},P::Matrix{Float64},shells::Vector{Vector{Int64}},
    basis_fun::Vector{Vector{Int64}}, params::Parameters) ::Vector{Float64}
    (; std_sh_pop) = params
    nbf = length(basis_fun)
    nsh = length(shells)
    shell_charges::Vector{Float64} = zeros(nsh)
    for i in 1:nbf, j in 1:nbf
        shell_charges[basis_fun[i][4]] -= S[i,j]*P[i,j]
    end
    for i in 1:nsh
        shell_charges[i] += std_sh_pop[id[shells[i][1]],shells[i][2]]
    end
    return shell_charges
end

"""
    get_atomic_charges(natoms::Int64,shell_charges::Vector{Float64},shells::Vector{Vector{Int64}})

    Returns the atomic charges calculated from the shell charges.
"""
function get_atomic_charges(natoms::Int64,shell_charges::Vector{Float64},shells::Vector{Vector{Int64}})::Vector{Float64}
    charges::Vector{Float64} = zeros(natoms)
    for i in eachindex(shell_charges)
        charges[shells[i][1]] += shell_charges[i]
    end
    return charges
end
"""
    damp_charges!(atomic_charges::Vector{Float64},shell_charges::Vector{Float64},
    new_atomic_charges::Vector{Float64},new_shell_charges::Vector{Float64};
    λ_damp::Float64,Δq_max::Float64)

    Updates the charges. Performs damping if necessary.
"""
function damp_charges!(atomic_charges::Vector{Float64},shell_charges::Vector{Float64},
    new_atomic_charges::Vector{Float64},new_shell_charges::Vector{Float64},
    λ_damp::Float64;Δq_max::Float64=1e-3)
    if maximum(abs.(new_shell_charges .- shell_charges)) >= Δq_max
        atomic_charges .= atomic_charges .+ λ_damp.*(new_atomic_charges.-atomic_charges)
        shell_charges .= shell_charges .+ λ_damp.*(new_shell_charges.-shell_charges)
    else
        atomic_charges .= new_atomic_charges
        shell_charges .= new_shell_charges
    end
end
end