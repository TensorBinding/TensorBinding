#!/usr/bin/env julia
# Generates TensorBinding module dependency diagram as SVG

W, H = 1550, 960

parts = String[]

push!(parts, """<?xml version="1.0" encoding="UTF-8"?>
<svg width="$W" height="$H" xmlns="http://www.w3.org/2000/svg">
<defs>
  <marker id="arr" markerWidth="9" markerHeight="6" refX="9" refY="3" orient="auto">
    <polygon points="0 0,9 3,0 6" fill="#555"/>
  </marker>
  <marker id="arr_ext" markerWidth="9" markerHeight="6" refX="9" refY="3" orient="auto">
    <polygon points="0 0,9 3,0 6" fill="#999"/>
  </marker>
</defs>
<rect width="$W" height="$H" fill="#fafafa"/>
""")

# ── helpers ────────────────────────────────────────────────────────────────────

function box!(parts, cx, cy, w, h, label, fill, stroke; fs=13, italic=false)
    x, y = cx-w/2, cy-h/2
    style = italic ? "font-style='italic'" : ""
    push!(parts, """<rect x="$x" y="$y" width="$w" height="$h" rx="7" fill="$fill" stroke="$stroke" stroke-width="2"/>""")
    push!(parts, """<text x="$cx" y="$(cy+5)" text-anchor="middle" dominant-baseline="middle" font-family="'Courier New',monospace" font-size="$fs" font-weight="bold" $style fill="#111">$label</text>""")
end

function label!(parts, cx, cy, txt; fs=11, color="#444")
    push!(parts, """<text x="$cx" y="$cy" text-anchor="middle" font-family="Arial,sans-serif" font-size="$fs" fill="$color" font-style="italic">$txt</text>""")
end

# straight arrow (bottom of src → top of dst)
function arr!(parts, x1, y1, x2, y2; color="#555", dash="", marker="arr")
    da = dash=="" ? "" : """stroke-dasharray="$dash" """
    push!(parts, """<line x1="$x1" y1="$y1" x2="$x2" y2="$y2" stroke="$color" stroke-width="1.6" marker-end="url(#$marker)" $da/>""")
end

# quadratic bezier arrow
function carr!(parts, x1, y1, x2, y2, bx, by; color="#555", dash="", marker="arr")
    da = dash=="" ? "" : """stroke-dasharray="$dash" """
    push!(parts, """<path d="M $x1 $y1 Q $bx $by $x2 $y2" fill="none" stroke="$color" stroke-width="1.6" marker-end="url(#$marker)" $da/>""")
end

# cubic bezier arrow
function cubic!(parts, x1,y1, cx1,cy1, cx2,cy2, x2,y2; color="#555", dash="", marker="arr")
    da = dash=="" ? "" : """stroke-dasharray="$dash" """
    push!(parts, """<path d="M $x1 $y1 C $cx1 $cy1 $cx2 $cy2 $x2 $y2" fill="none" stroke="$color" stroke-width="1.6" marker-end="url(#$marker)" $da/>""")
end

# ── node registry ──────────────────────────────────────────────────────────────
# name => (cx, cy, w, h)
nodes = Dict{String,NTuple{4,Float64}}()

function defnode!(nodes, name, cx, cy, w, h)
    nodes[name] = (Float64(cx), Float64(cy), Float64(w), Float64(h))
end

# bottom-center of node
bot(n) = (nodes[n][1], nodes[n][2] + nodes[n][4]/2)
# top-center of node
top(n) = (nodes[n][1], nodes[n][2] - nodes[n][4]/2)
# left-center
lft(n) = (nodes[n][1] - nodes[n][3]/2, nodes[n][2])
# right-center
rgt(n) = (nodes[n][1] + nodes[n][3]/2, nodes[n][2])

# ── layout constants ───────────────────────────────────────────────────────────
H_NODE = 38

y_ext   = 85.0
y_found = 200.0
y_struct= 325.0
y_comp  = 455.0
y_anal  = 585.0
y_ext2  = 715.0

# ── LAYER LABELS ──────────────────────────────────────────────────────────────
push!(parts, """<text x="18" y="$(y_ext+5)" text-anchor="start" font-family="Arial,sans-serif" font-size="11" fill="#999" font-style="italic">External</text>""")
push!(parts, """<text x="18" y="$(y_found+5)" text-anchor="start" font-family="Arial,sans-serif" font-size="11" fill="#2060b0" font-style="italic">Foundation</text>""")
push!(parts, """<text x="18" y="$(y_struct+5)" text-anchor="start" font-family="Arial,sans-serif" font-size="11" fill="#007030" font-style="italic">Geometry / System</text>""")
push!(parts, """<text x="18" y="$(y_comp+5)" text-anchor="start" font-family="Arial,sans-serif" font-size="11" fill="#b06000" font-style="italic">Computation</text>""")
push!(parts, """<text x="18" y="$(y_anal+5)" text-anchor="start" font-family="Arial,sans-serif" font-size="11" fill="#5030a0" font-style="italic">Observables</text>""")
push!(parts, """<text x="18" y="$(y_ext2+5)" text-anchor="start" font-family="Arial,sans-serif" font-size="11" fill="#b02050" font-style="italic">Extensions</text>""")

# horizontal separators
for y in [y_ext+25, y_found+25, y_struct+25, y_comp+25, y_anal+25, y_ext2+25]
    push!(parts, """<line x1="145" y1="$y" x2="1520" y2="$y" stroke="#ddd" stroke-width="1"/>""")
end

# ── EXTERNAL PACKAGES ─────────────────────────────────────────────────────────
defnode!(nodes, "ext_itensors",  280, y_ext, 240, 36)
defnode!(nodes, "ext_quantics",  730, y_ext, 330, 36)
defnode!(nodes, "ext_other",    1210, y_ext, 220, 36)

box!(parts,  280, y_ext, 240, 36, "ITensors · ITensorMPS · NDTensors", "#f0f0f0","#aaa"; fs=11)
box!(parts,  730, y_ext, 330, 36, "QuanticsTCI · TCI · Quantics · QuanticsGrids", "#f0f0f0","#aaa"; fs=11)
box!(parts, 1210, y_ext, 220, 36, "LinearAlgebra · Plots · PyCall", "#f0f0f0","#aaa"; fs=11)

# ── FOUNDATION ────────────────────────────────────────────────────────────────
defnode!(nodes, "utils",       310, y_found, 130, H_NODE)
defnode!(nodes, "hamiltonian", 840, y_found, 155, H_NODE)

box!(parts, 310, y_found, 130, H_NODE, "utils.jl",       "#cce4ff","#2060b0")
box!(parts, 840, y_found, 155, H_NODE, "Hamiltonian.jl", "#cce4ff","#2060b0")

# ── GEOMETRY / SYSTEM ─────────────────────────────────────────────────────────
defnode!(nodes, "lattice2d", 220, y_struct, 150, H_NODE)
defnode!(nodes, "tbsystem",  760, y_struct, 140, H_NODE)

box!(parts, 220, y_struct, 150, H_NODE, "2D_lattice.jl", "#c8ffd8","#007030")
box!(parts, 760, y_struct, 140, H_NODE, "TBSystem.jl",   "#c8ffd8","#007030")

# ── COMPUTATION ───────────────────────────────────────────────────────────────
defnode!(nodes, "kpm",    360, y_comp, 120, H_NODE)
defnode!(nodes, "qft",    640, y_comp, 120, H_NODE)
defnode!(nodes, "timeev", 910, y_comp, 130, H_NODE)
defnode!(nodes, "dmrg",  1170, y_comp, 120, H_NODE)

box!(parts,  360, y_comp, 120, H_NODE, "KPM_tk.jl",    "#ffd8a0","#b06000")
box!(parts,  640, y_comp, 120, H_NODE, "QFT_tk.jl",    "#ffd8a0","#b06000")
box!(parts,  910, y_comp, 130, H_NODE, "Timeev_tk.jl", "#ffd8a0","#b06000")
box!(parts, 1170, y_comp, 120, H_NODE, "DMRG_tk.jl",   "#ffd8a0","#b06000")

# ── OBSERVABLES / ANALYSIS ────────────────────────────────────────────────────
defnode!(nodes, "topology",     130, y_anal, 155, H_NODE)
defnode!(nodes, "purification", 390, y_anal, 170, H_NODE)
defnode!(nodes, "meanfi",       680, y_anal, 140, H_NODE)
defnode!(nodes, "rpa",          940, y_anal, 110, H_NODE)

box!(parts,  130, y_anal, 155, H_NODE, "Topology_tk.jl",     "#e0ccff","#5030a0")
box!(parts,  390, y_anal, 170, H_NODE, "Purification_tk.jl", "#e0ccff","#5030a0")
box!(parts,  680, y_anal, 140, H_NODE, "Meanfi_tk.jl",       "#e0ccff","#5030a0")
box!(parts,  940, y_anal, 110, H_NODE, "RPA_tk.jl",          "#e0ccff","#5030a0")

# ── EXTENSIONS ────────────────────────────────────────────────────────────────
defnode!(nodes, "supercond", 180, y_ext2, 160, H_NODE)
defnode!(nodes, "twisted",   580, y_ext2, 140, H_NODE)
defnode!(nodes, "bilayer",  1040, y_ext2, 140, H_NODE)

box!(parts,  180, y_ext2, 160, H_NODE, "Supercond_tk.jl", "#ffc8d8","#b02050")
box!(parts,  580, y_ext2, 140, H_NODE, "twisted_tk.jl",   "#ffc8d8","#b02050")
box!(parts, 1040, y_ext2, 140, H_NODE, "bilayer_tk.jl",   "#ffc8d8","#b02050")

# ── EDGES ─────────────────────────────────────────────────────────────────────

# External → internal (dashed, gray)
# ext_quantics → hamiltonian
let (x1,y1)=bot("ext_quantics"), (x2,y2)=top("hamiltonian")
    carr!(parts, x1,y1, x2,y2, (x1+x2)/2+40, (y1+y2)/2; color="#aaa", dash="5,3", marker="arr_ext")
end
# ext_quantics → QFT (long dashed)
let (x1,y1)=bot("ext_quantics"), (x2,y2)=top("qft")
    cubic!(parts, x1,y1, x1,y1+60, x2,y2-60, x2,y2; color="#bbb", dash="5,3", marker="arr_ext")
end
# ext_itensors → timeev
let (x1,y1)=bot("ext_itensors"), (x2,y2)=top("timeev")
    cubic!(parts, x1,y1, x1,y1+80, x2,y2-80, x2,y2; color="#bbb", dash="5,3", marker="arr_ext")
end
# ext_itensors → dmrg
let (x1,y1)=bot("ext_itensors"), (x2,y2)=top("dmrg")
    cubic!(parts, x1,y1, x1+200,y1+150, x2-200,y2-150, x2,y2; color="#bbb", dash="5,3", marker="arr_ext")
end

# utils → hamiltonian
let (x1,y1)=rgt("utils"), (x2,y2)=lft("hamiltonian")
    arr!(parts, x1,y1, x2,y2; color="#2060b0")
end
# utils → 2D_lattice
let (x1,y1)=bot("utils"), (x2,y2)=top("lattice2d")
    carr!(parts, x1,y1, x2,y2, x1-40,(y1+y2)/2; color="#2060b0")
end
# utils → tbsystem
let (x1,y1)=bot("utils"), (x2,y2)=top("tbsystem")
    carr!(parts, x1,y1, x2,y2, (x1+x2)/2,(y1+y2)/2; color="#2060b0")
end
# utils → supercond (long, left side)
let (x1,y1)=lft("utils"), (x2,y2)=top("supercond")
    cubic!(parts, x1,y1, x1-80,y1+100, x2-60,y2-100, x2,y2; color="#2060b0")
end
# utils → topology
let (x1,y1)=lft("utils"), (x2,y2)=top("topology")
    cubic!(parts, x1,y1, x1-60,y1+150, x2-30,y2-80, x2,y2; color="#2060b0")
end

# Hamiltonian → 2D_lattice
let (x1,y1)=bot("hamiltonian"), (x2,y2)=top("lattice2d")
    cubic!(parts, x1,y1, x1-100,y1+40, x2+60,y2-40, x2,y2; color="#2060b0")
end
# Hamiltonian → TBSystem
let (x1,y1)=bot("hamiltonian"), (x2,y2)=top("tbsystem")
    carr!(parts, x1,y1, x2,y2, (x1+x2)/2+20,(y1+y2)/2; color="#2060b0")
end
# Hamiltonian → Timeev
let (x1,y1)=bot("hamiltonian"), (x2,y2)=top("timeev")
    carr!(parts, x1,y1, x2,y2, (x1+x2)/2+20,(y1+y2)/2; color="#2060b0")
end
# Hamiltonian → twisted (long right curve)
let (x1,y1)=bot("hamiltonian"), (x2,y2)=top("twisted")
    cubic!(parts, x1,y1, x1,y1+120, x2,y2-120, x2,y2; color="#2060b0")
end

# 2D_lattice → twisted
let (x1,y1)=bot("lattice2d"), (x2,y2)=top("twisted")
    cubic!(parts, x1,y1, x1+60,y1+100, x2-80,y2-100, x2,y2; color="#007030")
end
# 2D_lattice → bilayer (via long curve)
let (x1,y1)=rgt("lattice2d"), (x2,y2)=top("bilayer")
    cubic!(parts, x1,y1, x1+200,y1+50, x2+30,y2-100, x2,y2; color="#007030")
end

# TBSystem → KPM
let (x1,y1)=bot("tbsystem"), (x2,y2)=top("kpm")
    cubic!(parts, x1,y1, x1-80,y1+40, x2+60,y2-40, x2,y2; color="#007030")
end
# TBSystem → QFT
let (x1,y1)=bot("tbsystem"), (x2,y2)=top("qft")
    carr!(parts, x1,y1, x2,y2, (x1+x2)/2-10,(y1+y2)/2; color="#007030")
end
# TBSystem → twisted
let (x1,y1)=bot("tbsystem"), (x2,y2)=top("twisted")
    carr!(parts, x1,y1, x2,y2, (x1+x2)/2+40,(y1+y2)/2+20; color="#007030")
end
# TBSystem → Purification
let (x1,y1)=bot("tbsystem"), (x2,y2)=top("purification")
    cubic!(parts, x1,y1, x1-60,y1+60, x2+40,y2-60, x2,y2; color="#007030")
end
# TBSystem → Meanfi
let (x1,y1)=bot("tbsystem"), (x2,y2)=top("meanfi")
    carr!(parts, x1,y1, x2,y2, (x1+x2)/2,(y1+y2)/2; color="#007030")
end

# KPM → Topology
let (x1,y1)=bot("kpm"), (x2,y2)=top("topology")
    cubic!(parts, x1,y1, x1-40,y1+40, x2+30,y2-40, x2,y2; color="#b06000")
end
# KPM → Purification
let (x1,y1)=bot("kpm"), (x2,y2)=top("purification")
    carr!(parts, x1,y1, x2,y2, (x1+x2)/2-10,(y1+y2)/2; color="#b06000")
end
# KPM → Meanfi
let (x1,y1)=bot("kpm"), (x2,y2)=top("meanfi")
    cubic!(parts, x1,y1, x1+80,y1+40, x2-80,y2-40, x2,y2; color="#b06000")
end
# KPM → RPA
let (x1,y1)=rgt("kpm"), (x2,y2)=top("rpa")
    cubic!(parts, x1,y1, x1+160,y1+10, x2+10,y2-80, x2,y2; color="#b06000")
end

# twisted → bilayer
let (x1,y1)=rgt("twisted"), (x2,y2)=lft("bilayer")
    arr!(parts, x1,y1, x2,y2; color="#b02050")
end

# Hamiltonian → bilayer (right long curve)
let (x1,y1)=rgt("hamiltonian"), (x2,y2)=top("bilayer")
    cubic!(parts, x1,y1, x1+120,y1+80, x2+60,y2-120, x2,y2; color="#2060b0")
end

# ── TITLE ─────────────────────────────────────────────────────────────────────
push!(parts, """<text x="775" y="35" text-anchor="middle" font-family="Arial,sans-serif" font-size="21" font-weight="bold" fill="#222">TensorBinding.jl — Module Dependency Map</text>""")

# ── LEGEND ────────────────────────────────────────────────────────────────────
lx, ly = 1350, 830
push!(parts, """<rect x="$lx" y="$ly" width="170" height="110" rx="6" fill="white" stroke="#ccc" stroke-width="1"/>""")
push!(parts, """<text x="$(lx+85)" y="$(ly+18)" text-anchor="middle" font-family="Arial,sans-serif" font-size="11" font-weight="bold" fill="#333">Legend</text>""")

for (i,(fill,stroke,lbl)) in enumerate([
    ("#cce4ff","#2060b0","Foundation"),
    ("#c8ffd8","#007030","Geometry / System"),
    ("#ffd8a0","#b06000","Computation"),
    ("#e0ccff","#5030a0","Observables"),
    ("#ffc8d8","#b02050","Extensions"),
])
    ry = ly + 24 + (i-1)*16
    push!(parts, """<rect x="$(lx+8)" y="$ry" width="14" height="12" rx="2" fill="$fill" stroke="$stroke" stroke-width="1.5"/>""")
    push!(parts, """<text x="$(lx+27)" y="$(ry+10)" font-family="Arial,sans-serif" font-size="10" fill="#333">$lbl</text>""")
end

# dashed = external dep note
push!(parts, """<line x1="$(lx+8)" y1="$(ly+106)" x2="$(lx+22)" y2="$(ly+106)" stroke="#999" stroke-width="1.5" stroke-dasharray="4,2"/>""")
push!(parts, """<text x="$(lx+27)" y="$(ly+110)" font-family="Arial,sans-serif" font-size="10" fill="#666">External dependency</text>""")

# ── CLOSE ─────────────────────────────────────────────────────────────────────
push!(parts, "</svg>")

open("TensorBinding_DepMap.svg", "w") do f
    write(f, join(parts, "\n"))
end
println("Saved TensorBinding_DepMap.svg")
