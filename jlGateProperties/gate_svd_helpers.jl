function parse_dimension(args)
    if isempty(args)
        return 2
    elseif length(args) == 1
        d = tryparse(Int, args[1])
        d === nothing && error("expected integer local dimension d, got $(args[1])")
        d < 2 && error("expected local dimension d >= 2, got $d")
        return d
    end

    error("usage: julia --project=. jlGateProperties/svd_2q_gates.jl [d]")
end

function make_two_site_indices(d)
    if d == 2
        sa = ITensorMPS.siteind("S=1/2")
        sb = ITensorMPS.siteind("S=1/2")
    else
        sa = ITensorMPS.siteind("Qudit"; dim=d)
        sb = ITensorMPS.siteind("Qudit"; dim=d)
    end

    # Each gate layer acts on a new prime level of the same two sites.
    sa1 = ITensors.prime(sa)
    sa2 = ITensors.prime(sa1)
    sa3 = ITensors.prime(sa2)
    sa4 = ITensors.prime(sa3)

    sb1 = ITensors.prime(sb)
    sb2 = ITensors.prime(sb1)
    sb3 = ITensors.prime(sb2)
    sb4 = ITensors.prime(sb3)

    return (; sa, sb, sa1, sa2, sa3, sa4, sb1, sb2, sb3, sb4)
end

computational_basis_index(a, b, d) = a * d + b + 1

function gate_singular_values(gate, left_inds)
    _, s, _ = ITensors.svd(gate, left_inds)
    s_matrix = Array(s, ITensors.inds(s)...)
    return [s_matrix[i, i] for i in axes(s_matrix, 1)]
end

function print_singular_values(label, singular_values)
    println(label)
    @show singular_values
    println()
end

function composite_gate(indices)
    return ITensorMPS.op("T", indices.sa) * ITensorMPS.op("I", indices.sb) *
           ITensorMPS.op("CX", indices.sa1, indices.sb1) *
           ITensorMPS.op("I", indices.sa2) * ITensorMPS.op("T", indices.sb2) *
           ITensorMPS.op("CX", indices.sb3, indices.sa3)
end

function product_gate(indices)
    return ITensorMPS.op("T", indices.sa) * ITensorMPS.op("S", indices.sb)
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

function sum_gate_matrix(d)
    matrix = zeros(ComplexF64, d^2, d^2)
    for a in 0:(d - 1), b in 0:(d - 1)
        input_index = computational_basis_index(a, b, d)
        output_index = computational_basis_index(a, mod(a + b, d), d)
        matrix[output_index, input_index] = 1.0
    end
    return matrix
end

function cz_gate_matrix(d)
    matrix = zeros(ComplexF64, d^2, d^2)
    ω = cis(2 * π / d)
    for a in 0:(d - 1), b in 0:(d - 1)
        index = computational_basis_index(a, b, d)
        matrix[index, index] = ω^(a * b)
    end
    return matrix
end

function swap_gate_matrix(d)
    matrix = zeros(ComplexF64, d^2, d^2)
    for a in 0:(d - 1), b in 0:(d - 1)
        input_index = computational_basis_index(a, b, d)
        output_index = computational_basis_index(b, a, d)
        matrix[output_index, input_index] = 1.0
    end
    return matrix
end

function sum_gate(d, sa, sb)
    return ITensors.op(sum_gate_matrix(d), sa, sb)
end

function cz_gate(d, sa, sb)
    return ITensors.op(cz_gate_matrix(d), sa, sb)
end

function swap_gate(d, sa, sb)
    return ITensors.op(swap_gate_matrix(d), sa, sb)
end

function qubit_fixed_gate_cases(indices)
    return [
        ("T I CXab I T CXba", composite_gate(indices), [indices.sa, indices.sa4]),
        ("CZ", ITensorMPS.op("CZ", indices.sa, indices.sb), [indices.sa, indices.sa1]),
        ("iSWAP", ITensorMPS.op("iSWAP", indices.sa, indices.sb), [indices.sa, indices.sa1]),
        ("Product gate", product_gate(indices), [indices.sa, indices.sa1]),
        ("SWAP", ITensorMPS.op("SWAP", indices.sa, indices.sb), [indices.sa, indices.sa1]),
    ]
end

function qudit_fixed_gate_cases(d, indices)
    return [
        ("SUM (CX_d)", sum_gate(d, indices.sa, indices.sb), [indices.sa, indices.sa1]),
        ("CZ_d", cz_gate(d, indices.sa, indices.sb), [indices.sa, indices.sa1]),
        ("SWAP_d", swap_gate(d, indices.sa, indices.sb), [indices.sa, indices.sa1]),
    ]
end

function rxx_sweep_cases(indices)
    return [
        ("RXX(pi)", rxx_gate(pi, indices.sa, indices.sb), [indices.sa, indices.sa1]),
        ("RXX(pi/2)", rxx_gate(pi / 2, indices.sa, indices.sb), [indices.sa, indices.sa1]),
        ("RXX(pi/3)", rxx_gate(pi / 3, indices.sa, indices.sb), [indices.sa, indices.sa1]),
        ("RXX(pi/4)", rxx_gate(pi / 4, indices.sa, indices.sb), [indices.sa, indices.sa1]),
        ("RXX(pi/8)", rxx_gate(pi / 8, indices.sa, indices.sb), [indices.sa, indices.sa1]),
    ]
end

function run_cases(cases)
    for (label, gate, left_inds) in cases
        print_singular_values(label, gate_singular_values(gate, left_inds))
    end
end

function main(d)
    indices = make_two_site_indices(d)

    if d == 2
        run_cases(qubit_fixed_gate_cases(indices))
        run_cases(rxx_sweep_cases(indices))
        return nothing
    end

    println("d = $d")
    println("Using standard qudit gate generalizations only.")
    println()
    run_cases(qudit_fixed_gate_cases(d, indices))
    return nothing
end
