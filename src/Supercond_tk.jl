# Supercond_tk.jl — Spin and Nambu auxiliary degrees of freedom for MPO Hamiltonians
#
# Follows the same prepend-core pattern as twisted_tk.jl: an auxiliary site
# (spin or particle/hole) is prepended to a position-qubit MPO, extending
# it by one site.  Multiple prepends can be chained:
#
#   [nambu_s, spin_s, pos_qubits...]   ← BdG with spin (call prepend_spin first)
#   [spin_s,  pos_qubits...]           ← spin-resolved tight-binding
#   [nambu_s, pos_qubits...]           ← spinless BdG
#
# Operator convention (both spin and Nambu use 2-state 1-indexed basis):
#   spin:  state 1 = ↑,        state 2 = ↓
#   Nambu: state 1 = particle,  state 2 = hole

# ─────────────────────────────────────────────────────────────────
# 0.  General matrix-based prepend  (public; handles complex types)
# ─────────────────────────────────────────────────────────────────

"""
    prepend_core_op(H_mpo, s, mat) -> MPO

Extend `H_mpo` by prepending a single-site operator given by the matrix `mat`
acting on index `s`.  Entry `mat[i,j]` = ⟨i|op|j⟩ (1-indexed basis of `s`).

This is the public, type-aware generalisation of the internal
`_prepend_layer_core` used in twisted_tk.jl.  It correctly handles complex
matrices (e.g., σ_y, τ_y) by allocating the ITensor with the element type of
`mat`.

The returned MPO has site indices `[s; original sites…]`.
"""
function prepend_core_op(H_mpo::MPO, s::Index, mat::AbstractMatrix{T}) where T <: Number
    Lh    = length(H_mpo)
    bond0 = Index(1, "Link,l=0")
    Op    = ITensor(T, s', s, bond0)
    for j in axes(mat, 2), i in axes(mat, 1)
        iszero(mat[i, j]) || (Op[s' => i, s => j, bond0 => 1] = mat[i, j])
    end
    delta0 = ITensor(bond0);  delta0[bond0 => 1] = 1.0
    H1_ext = H_mpo[1] * delta0
    ext    = MPO(Lh + 1)
    ext[1] = Op
    ext[2] = H1_ext
    for k in 3:Lh+1
        ext[k] = H_mpo[k-1]
    end
    return ext
end

# Convenience: accept an untyped matrix by promoting to ComplexF64
prepend_core_op(H::MPO, s::Index, mat::AbstractMatrix) =
    prepend_core_op(H, s, ComplexF64.(mat))


# ─────────────────────────────────────────────────────────────────
# 1.  Spin-½ site index and operators
# ─────────────────────────────────────────────────────────────────

"""
    spin_index() -> Index

Create a dim-2 Index tagged "Spin" (state 1 = ↑, state 2 = ↓).
Pass the result as `spin_s` to all `prepend_spin` calls.
"""
spin_index() = Index(2, "Spin")


# 2×2 spin-½ operator matrices, ComplexF64 throughout for uniformity.
# Basis: |↑⟩ = 1, |↓⟩ = 2.
const _SPIN_OPS = Dict{Symbol, Matrix{ComplexF64}}(
    :Id   => [1   0;  0   1],
    :Pup  => [1   0;  0   0],          # |↑⟩⟨↑|  — spin-up projector
    :Pdn  => [0   0;  0   1],          # |↓⟩⟨↓|  — spin-down projector
    :Sz   => [1/2 0;  0  -1/2],        # S_z = ½σ_z
    :Sp   => [0   1;  0   0],          # S_+ = |↑⟩⟨↓|  (spin-flip ↓→↑)
    :Sm   => [0   0;  1   0],          # S_- = |↓⟩⟨↑|  (spin-flip ↑→↓)
    :Sx   => [0   1/2; 1/2  0],        # S_x = ½σ_x
    :Sy   => [0  -1im/2; 1im/2  0],    # S_y = ½σ_y
    :iSy  => [0   1;  -1   0],         # i·σ_y  — singlet pairing spin factor
    :miSy => [0  -1;   1   0],         # (i·σ_y)† = −i·σ_y
)


"""
    prepend_spin(H_mpo, spin_s, op) -> MPO

Extend `H_mpo` by prepending the spin-½ operator `op` on `spin_s`
(created with `spin_index()`).

`op` may be a `Symbol` naming a built-in operator, or an explicit
2×2 `AbstractMatrix` for custom operators.

| Symbol  | Matrix             | Typical use                           |
|---------|--------------------|---------------------------------------|
| `:Id`   | I₂                 | Spin-degenerate term                  |
| `:Pup`  | diag(1,0)          | Spin-up projector                     |
| `:Pdn`  | diag(0,1)          | Spin-down projector                   |
| `:Sz`   | diag(½,−½)         | Zeeman / exchange field               |
| `:Sp`   | \\|↑⟩⟨↓\\|         | Spin-flip ↓→↑ (SOC, spin-orbit)      |
| `:Sm`   | \\|↓⟩⟨↑\\|         | Spin-flip ↑→↓ (SOC, spin-orbit)      |
| `:Sx`   | ½σ_x              | In-plane exchange                     |
| `:Sy`   | ½σ_y              | In-plane exchange                     |
| `:iSy`  | i·σ_y = [[0,1],[-1,0]] | Singlet pairing spin structure   |
| `:miSy` | −i·σ_y            | h.c. of singlet pairing               |

Basis: state 1 = ↑, state 2 = ↓.
"""
function prepend_spin(H_mpo::MPO, spin_s::Index,
                      op::Union{Symbol, AbstractMatrix})
    mat = op isa Symbol ?
          (haskey(_SPIN_OPS, op) ? _SPIN_OPS[op] :
           error("Unknown spin op :$op.  Known: $(sort(collect(keys(_SPIN_OPS))))")) :
          ComplexF64.(op)
    return prepend_core_op(H_mpo, spin_s, mat)
end


# ─────────────────────────────────────────────────────────────────
# 2.  Nambu (particle–hole) site index and operators
# ─────────────────────────────────────────────────────────────────

"""
    nambu_index() -> Index

Create a dim-2 Index tagged "Nambu" (state 1 = particle, state 2 = hole).
Pass the result as `nambu_s` to all `prepend_nambu` calls.
"""
nambu_index() = Index(2, "Nambu")


# 2×2 Nambu operator matrices, ComplexF64 throughout.
# Basis: |particle⟩ = 1, |hole⟩ = 2.
const _NAMBU_OPS = Dict{Symbol, Matrix{ComplexF64}}(
    :Id => [1   0;  0   1],
    :Pp => [1   0;  0   0],            # |p⟩⟨p|  — particle-sector projector
    :Ph => [0   0;  0   1],            # |h⟩⟨h|  — hole-sector projector
    :tz => [1   0;  0  -1],            # τ_z  — kinetic sign in BdG
    :tx => [0   1;  1   0],            # τ_x  — real pairing (spinless p-wave)
    :ty => [0  -1im; 1im  0],          # τ_y  — imaginary / chiral pairing
    :tp => [0   1;  0   0],            # τ_+ = |p⟩⟨h|  — pairing Δ
    :tm => [0   0;  1   0],            # τ_- = |h⟩⟨p|  — pairing Δ† (h.c.)
)


"""
    prepend_nambu(H_mpo, nambu_s, op) -> MPO

Extend `H_mpo` by prepending the Nambu (particle–hole) operator `op` on
`nambu_s` (created with `nambu_index()`).

`op` may be a `Symbol` naming a built-in operator, or an explicit
2×2 `AbstractMatrix`.

| Symbol | Matrix              | Typical use                          |
|--------|---------------------|--------------------------------------|
| `:Id`  | I₂                  | Particle + hole (no asymmetry)       |
| `:Pp`  | diag(1,0)           | Particle-sector projector            |
| `:Ph`  | diag(0,1)           | Hole-sector projector                |
| `:tz`  | diag(1,−1)          | Kinetic τ_z in BdG                   |
| `:tp`  | \\|p⟩⟨h\\|          | Pairing amplitude Δ                  |
| `:tm`  | \\|h⟩⟨p\\|          | Pairing Δ† (h.c.)                    |
| `:tx`  | σ_x                 | Real pairing (spinless p-wave)       |
| `:ty`  | σ_y                 | Imaginary / chiral pairing           |

Basis: state 1 = particle, state 2 = hole.

Typical spinless BdG assembly:
```julia
H_BdG = prepend_nambu(H_kin,         nambu_s, :tz)   # τ_z ⊗ H_kin
      + prepend_nambu(H_pair,        nambu_s, :tp)   # τ_+ ⊗ Δ
      + prepend_nambu(dag(H_pair),   nambu_s, :tm)   # τ_- ⊗ Δ†
```
"""
function prepend_nambu(H_mpo::MPO, nambu_s::Index,
                       op::Union{Symbol, AbstractMatrix})
    mat = op isa Symbol ?
          (haskey(_NAMBU_OPS, op) ? _NAMBU_OPS[op] :
           error("Unknown Nambu op :$op.  Known: $(sort(collect(keys(_NAMBU_OPS))))")) :
          ComplexF64.(op)
    return prepend_core_op(H_mpo, nambu_s, mat)
end


# ─────────────────────────────────────────────────────────────────
# 3.  Higher-level assemblers
# ─────────────────────────────────────────────────────────────────

"""
    spin_hamiltonian(H_up, H_down, spin_s;
                     H_Zeeman=nothing, cutoff=1e-8) -> MPO

Build a spin-resolved Hamiltonian on `[spin_s; pos_sites…]`:

    H = P_↑ ⊗ H_up  +  P_↓ ⊗ H_down  [+  S_z ⊗ H_Zeeman]

`H_up`, `H_down` are MPOs on the same position sites (they may differ for
spin-orbit coupling or magnetic exchange).  `H_Zeeman` is an optional
position-MPO encoding a local magnetic field `h(x)`; it enters as
`S_z ⊗ H_Zeeman` so spin-↑ gains `+½ h(x)` and spin-↓ gains `−½ h(x)`.
"""
function spin_hamiltonian(H_up::MPO, H_down::MPO, spin_s::Index;
                          H_Zeeman::Union{MPO, Nothing} = nothing,
                          cutoff::Real = 1e-8)
    H = +(prepend_spin(H_up,   spin_s, :Pup),
          prepend_spin(H_down, spin_s, :Pdn); cutoff=cutoff)
    isnothing(H_Zeeman) || (H = +(H, prepend_spin(H_Zeeman, spin_s, :Sz); cutoff=cutoff))
    return H
end


"""
    bdg_hamiltonian(H_kin, H_pair, nambu_s; cutoff=1e-8) -> MPO

Build a **spinless** Bogoliubov–de Gennes Hamiltonian on `[nambu_s; pos_sites…]`:

    H_BdG = τ_z ⊗ H_kin  +  τ_+ ⊗ H_pair  +  τ_- ⊗ dag(H_pair)

`H_kin` is the single-particle kinetic/hopping MPO measured from the chemical
potential (`H_kin = H_tb − μ·I`).  `H_pair` encodes the pairing amplitude
`Δ(i,j)`; for spatially uniform s-wave pairing build it as a diagonal MPO:
```julia
H_pair = qtt_mpo(L, 0:N-1, sites, _ -> Δ)
```
or use `hopping2MPO` for a spatially varying or p-wave pairing.

The result is Hermitian for any `H_kin = H_kin†` and any complex `H_pair`.
"""
function bdg_hamiltonian(H_kin::MPO, H_pair::MPO, nambu_s::Index;
                         cutoff::Real = 1e-8)
    H_pair_adj = swapprime(dag(H_pair), 0, 1)
    return +(+(prepend_nambu(H_kin,         nambu_s, :tz),
               prepend_nambu(H_pair,        nambu_s, :tp); cutoff=cutoff),
               prepend_nambu(H_pair_adj,    nambu_s, :tm); cutoff=cutoff)
end


"""
    bdg_spin_hamiltonian(H_kin_up, H_kin_down, H_pair, spin_s, nambu_s;
                         H_soc=nothing, cutoff=1e-8) -> MPO

Build a **spin-½ singlet** BdG Hamiltonian on `[nambu_s, spin_s; pos_sites…]`.

**Nambu–spin convention**: `Ψ = (c_↑, c_↓, c†_↓, −c†_↑)ᵀ` (standard BCS).

    H_BdG = τ_z⊗P_↑ ⊗ H_kin_up  +  τ_z⊗P_↓ ⊗ H_kin_down
          + τ_+⊗(i·σ_y) ⊗ H_pair  +  τ_-⊗(−i·σ_y) ⊗ H_pair†
          [+ τ_z⊗S_z ⊗ H_soc]

The first two lines are the kinetic energy (allowing spin-dependent fields,
e.g. Zeeman: pass `H_kin_up = H_tb − (μ+h)·I`, `H_kin_down = H_tb − (μ−h)·I`).

The pairing lines implement singlet Cooper-pair creation via the antisymmetric
spin factor `i·σ_y = [[0,1],[−1,0]]`.  For on-site s-wave pairing, build
`H_pair` as a diagonal MPO with `Δ` on the diagonal.

`H_soc` (optional) adds an Ising-type spin-orbit coupling `τ_z⊗S_z⊗H_soc`.

The result is Hermitian for real or complex `H_pair` and any `H_kin_up/dn`.
"""
function bdg_spin_hamiltonian(
    H_kin_up::MPO, H_kin_down::MPO, H_pair::MPO,
    spin_s::Index, nambu_s::Index;
    H_soc::Union{MPO, Nothing} = nothing,
    cutoff::Real = 1e-8,
)
    # τ_z ⊗ P_↑ ⊗ H_kin_up  and  τ_z ⊗ P_↓ ⊗ H_kin_down
    H = +(prepend_nambu(prepend_spin(H_kin_up,   spin_s, :Pup), nambu_s, :tz),
          prepend_nambu(prepend_spin(H_kin_down, spin_s, :Pdn), nambu_s, :tz); cutoff=cutoff)

    # τ_+ ⊗ (i·σ_y) ⊗ Δ  +  τ_- ⊗ (−i·σ_y) ⊗ Δ†
    # i·σ_y = [[0,1],[-1,0]] (real matrix): P_↑ pairs with h-↓, P_↓ pairs with h-↑ (singlet)
    H_pair_adj = swapprime(dag(H_pair), 0, 1)
    H_tp = prepend_nambu(prepend_spin(H_pair,       spin_s, :iSy),  nambu_s, :tp)
    H_tm = prepend_nambu(prepend_spin(H_pair_adj,   spin_s, :miSy), nambu_s, :tm)
    H    = +(+(H, H_tp; cutoff=cutoff), H_tm; cutoff=cutoff)

    # Optional Ising SOC: τ_z ⊗ S_z ⊗ H_soc
    if !isnothing(H_soc)
        H = +(H, prepend_nambu(prepend_spin(H_soc, spin_s, :Sz), nambu_s, :tz); cutoff=cutoff)
    end
    return H
end
