function MPS_from_number(sites, n::Int)
    N = length(sites)
    number = n - 1
    # Generate N-bit binary string with leading zeros (MSB-first)
    binary_str = string(number, base=2, pad=N)
    # Reverse only if Quantics grid expects LSB-first (uncomment if needed)
    # binary_str = reverse(binary_str)
    init_states = [bit == '1' ? "1" : "0" for bit in binary_str]
    
    new_sites = siteinds("Qubit", N)
    psi = MPS(new_sites, init_states)
    
    for i in 1:N
        psi[i] = replaceind(psi[i], new_sites[i], sites[i])
    end
    
    return psi
end

function get_correlator_quantics(A,L,sites)
    f(x,y) = inner(MPS_from_number(sites, Int(x)),apply(A,MPS_from_number(sites, Int(y))))
    mpo = MPS_to_MPO_Q(quanticsMPSMAT(f, 2^L, 1e-8),L)
    return mpo
end

function checker(mpo, sites, i, j)
    return inner(MPS_from_number(sites, Int(i)),apply(mpo,MPS_from_number(sites, Int(j))))
end

function get_pos_quantics(L,sites)
    f(x) = div(x + 1, 2)
    xvals = range(1, (2^L); length=2^L)
    qtt, ranks, errors = quanticscrossinterpolate(Float64, f,  xvals ; tolerance=1e-8)
    tt = TCI.tensortrain(qtt.tci)
    density_mps = ITensors.MPS(tt;sites)
    #does not have to be this outer product
    density_mpo = outer(density_mps',density_mps) 
    
    for i in 1:L
        density_mpo.data[i] =  Quantics._asdiagonal(density_mps.data[i],sites[i])
    end
    return density_mpo
end


function get_sz_quantics(L,sites)
    f(x) = (-1)^(x + 1)
    xvals = range(1, (2^L); length=2^L)
    println(xvals)
    qtt, ranks, errors = quanticscrossinterpolate(Float64, f,  xvals ; tolerance=1e-8)
    tt = TCI.tensortrain(qtt.tci)
    density_mps = ITensors.MPS(tt;sites)
    #does not have to be this outer product
    density_mpo = outer(density_mps',density_mps) 
    
    for i in 1:L
        density_mpo.data[i] =  Quantics._asdiagonal(density_mps.data[i],sites[i])
    end
    return density_mpo
end

function Cr(i)
    i1 = MPS_from_number(sites, Int(2*i-1))
    i2 = MPS_from_number(sites, Int(2*i))
    f1 = inner(i1,apply(C_op,i1))
    f2 = inner(i2,apply(C_op,i2))
    return (f1 + f2)
end