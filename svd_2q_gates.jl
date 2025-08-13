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

gate::ITensors.ITensor = 
    ITensorMPS.op("T", sa) * ITensorMPS.op("I", sb) * 
    ITensorMPS.op("CX", sa1, sb1) * 
    ITensorMPS.op("I", sa2) * ITensorMPS.op("T", sb2) * 
    ITensorMPS.op("CX", sb3, sa3)

println("T I CXab I T CXba")
u, s, v = ITensors.svd(gate, [sa, sa4])
@show s

println("CZ")
cz = ITensorMPS.op("CZ", sa, sb) 
u, s, v = ITensors.svd(cz, [sa, sa1])
@show s

println("iSWAP")
iswap = ITensorMPS.op("iSWAP", sa, sb) 
u, s, v = ITensors.svd(iswap, [sa, sa1])
@show s

println("Product gate")
productGate = 
    ITensorMPS.op("T", sa) * ITensorMPS.op("S", sb)
u, s, v = ITensors.svd(productGate, [sa, sa1])
@show s
