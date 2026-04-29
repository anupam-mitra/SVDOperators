import ITensors
import ITensorMPS

# Physical indices for two qubits
sa = ITensorMPS.siteind("S=1/2")
sb = ITensorMPS.siteind("S=1/2")

# For every output of a layer, we need primed indices
sa1 = ITensors.prime(sa)
sa2 = ITensors.prime(sa1)
sa3 = ITensors.prime(sa2)
sa4 = ITensors.prime(sa3)

sb1 = ITensors.prime(sb)
sb2 = ITensors.prime(sb1)
sb3 = ITensors.prime(sb2)
sb4 = ITensors.prime(sb3)

function show_singular_values(label, gate, left_inds)
    println(label)
    u, s, v = ITensors.svd(gate, left_inds)
    s_matrix = Array(s, ITensors.inds(s)...)
    singular_values = [s_matrix[i, i] for i in axes(s_matrix, 1)]
    @show singular_values
    println()
end

function rxx_gate(theta, sa, sb)
    c = cos(theta / 2)
    s = sin(theta / 2)
    matrix = ComplexF64[
        c 0 0 -im * s
        0 c -im * s 0
        0 -im * s c 0
        -im * s 0 0 c
    ]
    return ITensors.op(matrix, sa, sb)
end

gate::ITensors.ITensor = 
    ITensorMPS.op("T", sa) * ITensorMPS.op("I", sb) * 
    ITensorMPS.op("CX", sa1, sb1) * 
    ITensorMPS.op("I", sa2) * ITensorMPS.op("T", sb2) * 
    ITensorMPS.op("CX", sb3, sa3)

show_singular_values("T I CXab I T CXba", gate, [sa, sa4])

cz = ITensorMPS.op("CZ", sa, sb)
show_singular_values("CZ", cz, [sa, sa1])

iswap = ITensorMPS.op("iSWAP", sa, sb)
show_singular_values("iSWAP", iswap, [sa, sa1])

productGate = 
    ITensorMPS.op("T", sa) * ITensorMPS.op("S", sb)
show_singular_values("Product gate", productGate, [sa, sa1])

swap = ITensorMPS.op("SWAP", sa, sb)
show_singular_values("SWAP", swap, [sa, sa1])

for (label, theta) in [
    ("pi", pi),
    ("pi/2", pi / 2),
    ("pi/3", pi / 3),
    ("pi/4", pi / 4),
    ("pi/8", pi / 8),
]
    show_singular_values("RXX(" * label * ")", rxx_gate(theta, sa, sb), [sa, sa1])
end
