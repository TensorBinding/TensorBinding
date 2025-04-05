function convert_mps_to_mpo(mps,new_sites)
    #input is old MPS, new_sites are the actual sites of the system
    N = length(mps)
    new_N = N ÷ 2
    new_mpo = MPO(new_N)
    
    for i in 1:new_N
        idx1 = 2i - 1
        idx2 = 2i
        
         
        A = mps[idx1]
        B = mps[idx2]
        combined_T = A * B
        
     
        old_s1 = siteind(mps, idx1)
        old_s2 = siteind(mps, idx2)
  
        new_s_in = new_sites[i]'
        new_s_out = new_sites[i]
         
        new_mpo[i] = replaceinds(combined_T, [old_s1, old_s2] => [new_s_in, new_s_out])
    end
    
    return new_mpo
end

function binary_to_MPS(n, size,sites)
    # Convert to binary string
    binary_str = string(n, base=2)
    
    # Pad the binary string with leading zeros to match the desired size
    padded_binary_str = lpad(binary_str, size, '0')
    
    # Convert the padded string into MPS
    return random_mps(sites,collect(padded_binary_str) |> x -> map(s -> string(s), x))
end
