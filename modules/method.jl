module Methods
using LinearAlgebra: norm, Hermitian, Diagonal, eigen
using  Memoization: @memoize
export pairwise_distance, E_rep, E_disp, gen_basisfun, overlapmatrix, S_inv_sqrt, gen_H_0
function pairwise_distance(natoms :: Int64, coords :: Matrix{Float64})
    r = Matrix{Float64}(undef,natoms,natoms) #pairwise_distance
    for i in 1:natoms, j in i:natoms
        r[j,i]  = norm(coords[:,j]-coords[:,i])
    end
    return r
end
function E_rep(natoms :: Int64, id :: Vector{Int64}, r :: Matrix{Float64},
    α :: Vector{Float64}, z_eff :: Vector{Float64}, k_f :: Float64)
    E_rep ::Float64 = 0
    for i in 1:natoms, j in i+1:natoms
        E_rep +=
        z_eff[id[i]]z_eff[id[j]]/r[j,i] * ℯ^(-√(α[id[i]]α[id[j]])r[j,i]^k_f)
    end
    return E_rep
end
function cn(a::Int64,natoms::Int64, id :: Vector{Int64}, r::Matrix{Float64},sc_radii :: Vector{Float64},k_cn :: Float64)
    res::Float64 = 0
    for i in 1:a-1
        res += (1+ ℯ^(-k_cn*((sc_radii[id[a]]+sc_radii[id[i]])/r[a,i]-1)))^-1
    end
    for i in a+1:natoms
        res += (1+ ℯ^(-k_cn*((sc_radii[id[a]]+sc_radii[id[i]])/r[i,a]-1)))^-1
    end
    return res
end
function E_disp(natoms :: Int64, id :: Vector{Int64}, r :: Matrix{Float64},
    a1 :: Float64, a2 :: Float64, s6 :: Float64, s8 :: Float64, k_cn :: Float64, k_l :: Float64,
    q_a :: Vector{Float64}, sc_radii :: Vector{Float64}, numrefcn :: Vector{Int64},
    refcn :: Matrix{Float64}, c_ref :: Matrix{Matrix{Float64}})
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
    function c6(a::Int64,b::Int64,cn_a::Float64,cn_b, refcn :: Matrix{Float64})
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
    E_disp :: Float64 = 0
    for i in 1:natoms, j in i+1:natoms
        E_disp += -s6*c6(id[j],id[i],cn(j,natoms,id,r,sc_radii,k_cn),cn(i,natoms,id,r,sc_radii,k_cn),refcn)/
        r[j,i]^6*f6_damp(r[j,i],a1,a2,q_a[id[j]],q_a[id[i]])
        E_disp += -s8*3*√(q_a[id[j]]q_a[id[i]])*c6(id[j],id[i],cn(j,natoms,id,r,sc_radii,k_cn),cn(i,natoms,id,r,sc_radii,k_cn),refcn)/
        r[j,i]^8*f8_damp(r[j,i],a1,a2,q_a[id[j]],q_a[id[i]])
    end
    return E_disp
end
function gen_basisfun(natoms:: Int64, id :: Vector{Int64}, shtyp:: Vector{Vector{Int64}})
    basis_functions = Vector{Vector{Int64}}()
    for i in 1:natoms, sh in shtyp[id[i]]
        if sh == 1 || sh == 2
            push!(basis_functions,[i,sh,0])
        elseif sh ==3
            push!(basis_functions,[i,sh,1])
            push!(basis_functions,[i,sh,2])
            push!(basis_functions,[i,sh,3])
        end
    end
    return basis_functions
end
function overlapmatrix(basis_fun :: Vector{Vector{Int64}},
    id ::Vector{Int64},
    shtyp:: Vector{Vector{Int64}},
    coord :: Matrix{Float64},
    ζ:: Vector{Vector{Vector{Float64}}},
    d:: Vector{Vector{Vector{Float64}}})

    function ∫_bf_bf(bf1 :: Vector{Int64},bf2 :: Vector{Int64},
        id ::Vector{Int64},
        shtyp:: Vector{Vector{Int64}},
        coord :: Matrix{Float64},
        ζ:: Vector{Vector{Vector{Float64}}},
        d:: Vector{Vector{Vector{Float64}}})
            @memoize function ∫s_s(ζ :: Float64, ξ :: Float64,
                r1 :: Vector{Float64},r2 :: Vector{Float64})
                        return ℯ^(-ξ*norm(r1-r2)^2)*(π/ζ)^1.5
            end
            function ∫prim_prim(ζ1 :: Float64, ζ2 :: Float64,
                    r1 :: Vector{Float64},r2 :: Vector{Float64},
                    typ1 :: Int64, typ2 :: Int64)
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
        integral :: Float64 = 0
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
    S[j,i] = S[i,j] = ∫_bf_bf(basis_fun[i],basis_fun[j],id,shtyp,coord,ζ,d)
    end
    return Hermitian(S)
end
function S_inv_sqrt(S :: Hermitian{Float64})
    λ, L = eigen(S)
    return L * Diagonal(1 ./ sqrt.(λ)) * L'
end
function gen_H_0(basis_fun:: Vector{Vector{Int64}},natoms::Int64, r::Matrix{Float64}, id :: Vector{Int64},
    k_ab :: Dict{Tuple{Int64,Int64,Int64,Int64},Float64},k_ll:: Matrix{Float64},
    electroneg :: Vector{Float64}, k_en :: Float64, r_cov :: Vector{Float64},h_al::Matrix{Float64},
    k_poli::Matrix{Float64}, k_cn_l :: Vector{Float64},
    sc_radii::Vector{Float64},k_cn :: Float64, S :: Hermitian{Float64})

    function zeroh(bf1:: Vector{Int64},bf2 :: Vector{Int64},natoms::Int64, r::Matrix{Float64}, id :: Vector{Int64},
        k_ab :: Dict{Tuple{Int64,Int64,Int64,Int64},Float64},
        k_ll:: Matrix{Float64}, electroneg :: Vector{Float64}, k_en :: Float64,
        r_cov :: Vector{Float64},h_al::Matrix{Float64}, k_poli::Matrix{Float64},
        k_cn_l :: Vector{Float64},sc_radii::Vector{Float64},k_cn :: Float64, S :: Float64)
        function Π(bf1:: Vector{Int64},bf2 :: Vector{Int64},r::Matrix{Float64}, id :: Vector{Int64},
            k_poli::Matrix{Float64},r_cov :: Vector{Float64})
            return (1+k_poli[id[bf1[1]],bf1[2]]*(r[bf2[1],bf1[1]]/(r_cov[id[bf1[1]]]+r_cov[id[bf2[1]]]))^0.5)*
            (1+k_poli[id[bf2[1]],bf2[2]]*(r[bf2[1],bf1[1]]/(r_cov[id[bf1[1]]]+r_cov[id[bf2[1]]]))^0.5)
        end
        K = get(k_ab,(id[bf1[1]],bf1[2],bf2[1],bf2[2]),1.0)
        indx1, indx2 = min(bf1[2],bf2[2]), max(bf1[2],bf2[2])
        k_llp = k_ll[indx1,indx2]
        h_al1 = h_al[id[bf1[1]],bf1[2]]*(1+k_cn_l[bf1[2]]*cn(bf1[1],natoms,id,r,sc_radii,k_cn))
        h_al2 = h_al[id[bf2[1]],bf2[2]]*(1+k_cn_l[bf2[2]]*cn(bf2[1],natoms,id,r,sc_radii,k_cn)) 
        if bf1[2] == bf2[2] == 2
            off_diag_scal = 0
        else
            off_diag_scal = 1 + k_en*(electroneg[id[bf1[1]]] - electroneg[id[bf2[1]]])^2
        end
        if bf1[1] == bf2[1]
            if bf1 == bf2
                return h_al1
            else
                return 0
            end
        else
            pii = Π(bf1,bf2,r,id,k_poli,r_cov)
            return K*k_llp*(h_al1+h_al2)/2*S*off_diag_scal*pii
        end
    end

    H_0 = Matrix{Float64}(undef,size(S)...)
    dim = length(basis_fun)
    for i in 1:dim, j in i:dim
        H_0[j,i] = H_0[i,j] = zeroh(basis_fun[i],basis_fun[j],natoms,r,id,k_ab,k_ll,electroneg,
        k_en,r_cov,h_al,k_poli,k_cn_l,sc_radii,k_cn,S[i,j])
    end
    return Hermitian(H_0)
end


end