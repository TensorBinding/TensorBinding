# -*- coding: utf-8 -*-
using ITensors: MPO, MPS, OpSum, expect, inner, siteinds
using ITensors
using ITensorMPS
using LinearAlgebra
using Plots
using QuanticsTCI
import TensorCrossInterpolation as TCI
using Quantics


ITensors.op(::OpName"sigma_plus",::SiteType"Qubit") =
 [0 1
  0 0]

ITensors.op(::OpName"sigma_minus",::SiteType"Qubit") =
 [0 0
  1 0]

#the making of trivial hopping MPO
function kinetic(L, sites) 
 
    kinetic = OpSum()
    for i in 1:L
        os = OpSum()
        os += 1,"sigma_plus",L-(i-1)

        for i in 1:L-i 
            os *=  ("Id",i) 
        end

        for i in L+2-i :L 
            os *=  ("sigma_minus",i) 
        end
        
        kinetic += os
    end
        
    for i in 1:L
        os = OpSum()
        os += 1,"sigma_minus",L-(i-1)

        for i in 1:L-i 
            os *=  ("Id",i) 
        end

        for i in L+2-i :L 
            os *=  ("sigma_plus",i) 
        end
        
        kinetic += os
    end
 
    k_mpo = MPO(kinetic,sites)
    return k_mpo 
end

#the KPM process
function KPM_Tn(H,N)
     
    Ham_n = H/10 
    T_k_minus_2 = Id_op/L
    T_k_minus_1 = Ham_n   
    Tn_list = [T_k_minus_2,T_k_minus_1]

    for k in 1:N
        if k == 1
            T_k = T_k_minus_2
        elseif k == 2
            T_k = T_k_minus_1
        else
            #truncate to increase efficiency
            T_k = +(2 * apply(Ham_n, T_k_minus_1;cutoff=1e-8) , -T_k_minus_2;  maxdim =100) 
            T_k = ITensorMPS.truncate!(T_k;cutoff=1e-8)
            
            T_k_minus_2 = T_k_minus_1 
            T_k_minus_1 =  T_k    
            push!(Tn_list,T_k)
         
        end
    end
    return Tn_list 
end

#function for building MPS vector of atomic site
function to_binary_vector(n, size)
    # Convert to binary string
    binary_str = string(n, base=2)
    
    # Pad the binary string with leading zeros to match the desired size
    padded_binary_str = lpad(binary_str, size, '0')
    
    # Convert the padded string into a vector of strings (each character is a string)
    return collect(padded_binary_str) |> x -> map(s -> string(s), x)
end

#calculation of density matrix
function get_density_from_Tn(Tn_list,N,fermi=0)
       
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
        A = +(A,  2 *  Tn_list[n] * G_n(n) * jackson_kernel[n] ;maxdim=100)
        A = ITensorMPS.truncate!(A;cutoff=1e-8)
    end
    A /= (π* N)
    
    return  A
end

#get the density MPS
function get_density_quantics(A,L )
    
    xvals = range(0, (2^L - 1); length=2^L)
    f(x) =  1 -  inner(random_mps(sites,to_binary_vector(Int(x),L))',A, random_mps(sites,to_binary_vector(Int(x),L)))
    qtt, ranks, errors = quanticscrossinterpolate(Float64, f,  xvals ; tolerance=1e-8)
    tt = TCI.tensortrain(qtt.tci)
    density_mps = ITensors.MPS(tt;sites)
    #does not have to be this outer product
    density_mpo = outer(density_mps',density_mps) 
    
    for i in 1:L
        density_mpo.data[i] =  Quantics._asdiagonal(density_mps.data[i],sites[i])
    end
    
    return qtt,density_mpo,density_mps 
end

#general SCF
function SCF_Hubbard(L,
    sites, 
    u_mpo, 
    max_iter, 
    N,  
    threshold,
    mpo_guess_up,
    mpo_guess_down,
    mps_guess_up,
    mps_guess_down,
    qtt_den_old_up,
    qtt_den_old_down,
    mix)
     
    #initial MF Hamiltonian
    u_mpo_ini_up = apply(u_mpo, mpo_guess_up)
    u_mpo_ini_down = apply(u_mpo, mpo_guess_down)
    k_mpo = kinetic(L, sites) 

    #half-filling
    ham_up = +(k_mpo  , u_mpo_ini_down ,-  0.5 * u_mpo;maxdim=100)
    ham_down = +(k_mpo ,  u_mpo_ini_up   , -  0.5 * u_mpo;maxdim=100)

    conv_error = []

    for i in 1:max_iter
        
        Tn_list_up = KPM_Tn(ham_up,N)
        Tn_list_down = KPM_Tn(ham_down,N)
        
        A_up = get_density_from_Tn(Tn_list_up,N)
        A_down = get_density_from_Tn(Tn_list_down,N)
      
        qtt_den_new_up, den_mpo_up,den_mps_up = get_density_quantics(A_up,L)
        qtt_den_new_down, den_mpo_down,den_mps_down  = get_density_quantics(A_down,L)

        # Compute the difference for both up and down components
        diff_up = norm(den_mps_up - mps_guess_up)/(norm(Id_op)/L)
        diff_down = norm(den_mps_down -  mps_guess_down)/(norm(Id_op)/L)

        max_diff = (diff_up + diff_down)/2
    
        push!(conv_error, max_diff)
  
        # Check convergence
        if max_diff< threshold
            println("SCF converged in $(i) iterations.")
            return qtt_den_new_up,qtt_den_new_down, Tn_list_up,Tn_list_down, conv_error 
        end
        
        #dynamical mixing
        den_mpo_up_new = mix * den_mpo_up + (1-mix) * mpo_guess_up
        den_mpo_down_new = mix * den_mpo_down + (1-mix) * mpo_guess_down
        den_mps_up_new = mix * den_mps_up + (1-mix) * mps_guess_up
        den_mps_down_new = mix * den_mps_down + (1-mix) * mps_guess_down
        
        #update
        u_den_op_up = apply(u_mpo, den_mpo_up_new;cutoff=1e-8)
        u_den_op_down = apply(u_mpo, den_mpo_down_new;cutoff=1e-8)

        ham_up = +(k_mpo ,  u_den_op_down, - 0.5 *u_mpo;maxdim=100)# 
        ham_down = +(k_mpo ,  u_den_op_up ,- 0.5 *u_mpo;maxdim=100)
        
        mpo_guess_up = den_mpo_up_new
        mpo_guess_down = den_mpo_down_new
        mps_guess_up = den_mps_up_new
        mps_guess_down = den_mps_down_new
        
        qtt_den_old_up = qtt_den_new_up
        qtt_den_old_down = qtt_den_new_down
      
    end
    println("SCF can't converged.")
    return qtt_den_new_up,qtt_den_new_down,conv_error
end


