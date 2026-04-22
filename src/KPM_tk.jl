"""
    KPM_Tn(H_mpo, N, sites; scale=nothing, maxdim=40,
           dmrg_nsweeps, dmrg_maxdim, dmrg_linkdim) -> (Tn_list, scale, center)

Build the list of Chebyshev MPOs `T_n((H−center·I)/scale)` for `n = 0…N`.

## Scale argument
- If `scale` is provided the spectrum is assumed centered at zero and `H/scale`
  is used directly (backward-compatible path).
- If `scale` is `nothing` the spectral bounds are estimated automatically by
  running a short DMRG minimisation of H (ground state, E_min) and −H (maximum,
  E_max).  The scale and center are then:
      center = (E_max + E_min) / 2
      scale  = (E_max − E_min) / 2 × 1.1    ← 10 % buffer on each side

## Return value
Returns `(Tn_list, scale, center)`.  To convert a physical energy ω to the
scaled argument passed to `get_ldos_w_from_Tn`:
    ω_r = (ω − center) / scale  ∈ (−1, 1)
"""
function KPM_Tn(H_mpo::MPO, N::Int, sites;
                scale::Union{Real, Nothing} = nothing,
                maxdim::Int   = 40,
                dmrg_nsweeps::Int  = 5,
                dmrg_maxdim        = [10, 20, 40],
                dmrg_linkdim::Int  = 4,
                cutoff::Real       = 1e-8)

    # ── Spectral bounds ───────────────────────────────────────────────────
    if isnothing(scale)
        println("KPM_Tn: estimating spectral bounds via DMRG…")
        E_min, _ = dmrg_gs(H_mpo, sites;
                           nsweeps      = dmrg_nsweeps,
                           maxdim       = dmrg_maxdim,
                           linkdim_init = dmrg_linkdim,
                           noise        = [1e-6, 1e-7, 0.0],
                           outputlevel  = 0)
        E_max_neg, _ = dmrg_gs((-1.0) * H_mpo, sites;
                                nsweeps      = dmrg_nsweeps,
                                maxdim       = dmrg_maxdim,
                                linkdim_init = dmrg_linkdim,
                                noise        = [1e-6, 1e-7, 0.0],
                                outputlevel  = 0)
        E_max  = -E_max_neg
        center = (E_max + E_min) / 2
        scale  = (E_max - E_min) / 2 * 1.1
        println("  E_min = $(round(E_min; digits=4)),  E_max = $(round(E_max; digits=4))")
        println("  center = $(round(center; digits=4)),  scale = $(round(scale; digits=4))")
    else
        center = 0.0
    end

    # ── Scaled Hamiltonian: (H − center·I) / scale ────────────────────────
    I_mpo   = MPO(sites, "Id")
    Ham_n   = (1 / scale) * +(H_mpo, (-center) * I_mpo; cutoff = cutoff)

    # ── Chebyshev recursion T_0 = I,  T_1 = H_scaled,  T_k = 2H·T_{k-1} − T_{k-2}
    T_k_minus_2 = I_mpo
    T_k_minus_1 = Ham_n
    Tn_list = [T_k_minus_2, T_k_minus_1]

    for k in 3:N+1
        T_k = +(2 * apply(Ham_n, T_k_minus_1; cutoff = cutoff),
                -T_k_minus_2; maxdim = maxdim)
        T_k = ITensorMPS.truncate!(T_k; cutoff = cutoff)
        T_k_minus_2 = T_k_minus_1
        T_k_minus_1 = T_k
        push!(Tn_list, T_k)
        println(ITensorMPS.maxlinkdim(T_k))
    end

    return Tn_list, scale, center
end

function get_density_from_Tn(Tn_list,N;fermi=0,maxdim=40)  

    jackson_kernel = [(N - n) * cos(π * n / N) + sin(π * n / N) / tan(π / N) for n in 0:N-1]

    function G_n(n)
        if n == 1
            return acos(fermi)
        else
            return sin((n-1) * acos(fermi)) / (n-1)
        end
    end

    # Compute electronic density
    A = Tn_list[1] * G_n(1) * jackson_kernel[1] 
    for n in 2:N
        A = +(A,  2 *  Tn_list[n] * G_n(n) * jackson_kernel[n]; maxdim=maxdim)
        A = ITensorMPS.truncate!(A;cutoff=1e-8)
    end
    A /= (π* N)
    
    return  A
end

function get_Green_retarded_from_Tn(Tn_list, N, ω; η=1e-2, maxdim=40)

    # Jackson kernel for damping
    jackson_kernel = [(N - n) * cos(π * n / N) + sin(π * n / N) / tan(π / N)
                      for n in 0:N-1]

    # Chebyshev propagator term
    function G_n(n, ω, η)
        z = ω + 1im*η         # retarded Green’s function: +iη
        θ = acos(z)
        return -2im/(1+ ==(n-1,0)) * exp(-1im * (n-1) * θ) / (sqrt(1 - z^2))
    end

    # Build Green’s function expansion
    G = Tn_list[1] * G_n(1, ω, η) * jackson_kernel[1]
    for n in 2:N
        G = +(G, Tn_list[n] * G_n(n, ω, η) * jackson_kernel[n]; maxdim=maxdim)
        G = ITensorMPS.truncate!(G; cutoff=1e-8)
    end
    G /= N

    return G
end


function get_ldos_w_from_Tn(Tn_list,N,ω;maxdim=40)  

    jackson_kernel = [(N - n) * cos(π * n / N) + sin(π * n / N) / tan(π / N) for n in 0:N-1]

    function G_n(n)
       return cos((n - 1) * acos(ω))/(π* sqrt(1-ω^2))
    end

    # Compute electronic density
    A = Tn_list[1] * G_n(1) * jackson_kernel[1] 
    for n in 2:N
        A = +(A,  2 *  Tn_list[n] * G_n(n) * jackson_kernel[n]; maxdim=maxdim)
        A = ITensorMPS.truncate!(A;cutoff=1e-8)
    end
    A /= (π* N)
    
    return  A
end


function get_PH_from_Tn(Tn_list,N,ω;maxdim=40)  

    jackson_kernel = [(N - n) * cos(π * n / N) + sin(π * n / N) / tan(π / N) for n in 0:N-1]

    function G_n(n)
       return cos((n - 1) * acos(ω))/(π* sqrt(1-ω^2))
    end

    # Compute electronic density
    A = Tn_list[1] * G_n(1) * jackson_kernel[1] 
    for n in 2:N
        A = +(A,  2 *  Tn_list[n] * G_n(n) * jackson_kernel[n]; maxdim=maxdim)
        A = ITensorMPS.truncate!(A;cutoff=1e-8)
    end
    A /= (π* N)
    
    return  A
end


#for getting electron densities
function get_density_quantics(A,L)
    
    xvals = range(0, (2^L - 1); length=2^L)
    f(x) =  1 -  inner(random_mps(sites,to_binary_vector(Int(x),L))',A, random_mps(sites,to_binary_vector(Int(x),L)))
    qtt, ranks, errors = quanticscrossinterpolate(Float64, f,  xvals ; tolerance=1e-8)

    tt = TCI.tensortrain(qtt.tci)
    density_mps = ITensors.MPS(tt;sites)
  
    density_mpo = outer(density_mps',density_mps)
    for i in 1:L
        density_mpo.data[i] =  Quantics._asdiagonal(density_mps.data[i],sites[i])
    end
    
    return qtt,density_mpo,density_mps
end