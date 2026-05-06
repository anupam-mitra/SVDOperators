# SVDOperators

This repository compares operator Schmidt singular values for two-site gates in two implementations:

- Julia with `ITensors` / `ITensorMPS`
- Python with `NumPy` and Qiskit

Both implementations support:

- `d = 2`: qubit-specific gate comparisons
- `d > 2`: standard qudit generalizations only

There is no automated test suite or CI. The main verification workflow is to rerun the scripts and inspect the printed singular values.

## Repository Layout

- `jlGateProperties/svd_2q_gates.jl`: Julia entrypoint
- `jlGateProperties/gate_svd_helpers.jl`: Julia helper functions and gate definitions
- `pyGateProperties/Gate_Schmidt_Decomposition.py`: Python script and Jupytext source of truth
- `pyGateProperties/Gate_Schmidt_Decomposition.ipynb`: notebook mirror of the Python script
- `Project.toml`: Julia dependencies
- `requirements.txt`: Python dependencies
- `pyproject.toml`: Python project metadata and Jupytext config

## Prerequisites

- Julia installed and available on `PATH`
- Python 3.10+ installed and available on `PATH`

## Julia Workflow

### First-time setup

Install Julia dependencies from the repo root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The project environment is defined at the repo root by `Project.toml`.

### Run for qubits (`d = 2`)

Run the default qubit workflow:

```bash
julia --project=. jlGateProperties/svd_2q_gates.jl
```

This prints singular values for:

- `T I CXab I T CXba`
- `CZ`
- `iSWAP`
- `Product gate`
- `SWAP`
- `RXX(pi)`
- `RXX(pi/2)`
- `RXX(pi/3)`
- `RXX(pi/4)`
- `RXX(pi/8)`

### Run for qudits (`d > 2`)

Pass the local dimension as the first positional argument:

```bash
julia --project=. jlGateProperties/svd_2q_gates.jl 3
julia --project=. jlGateProperties/svd_2q_gates.jl 4
```

For `d > 2`, the Julia workflow intentionally switches to standard qudit gate generalizations only. It prints singular values for:

- `SUM (CX_d)`
- `CZ_d`
- `SWAP_d`

It does not define qudit analogues of the qubit-only cases such as `iSWAP`, `RXX`, `S`, or `T`.

### Julia implementation notes

- `jlGateProperties/svd_2q_gates.jl` is a thin entrypoint that calls `main(parse_dimension(ARGS))`.
- `jlGateProperties/gate_svd_helpers.jl` contains the gate builders, singular-value helper, and the `d = 2` / `d > 2` dispatch.
- For `d = 2`, the composite circuit uses a different tensor partition than the built-in two-qubit gates:
  - composite circuit: `[sa, sa4]`
  - built-in two-site gates: `[sa, sa1]`
- For `d > 2`, Julia uses `Qudit` indices plus explicit matrix-defined gates for `SUM`, `CZ_d`, and `SWAP_d`.

## Python Workflow

### Environment setup

Create and activate a virtual environment from the repo root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Optional: install the local Python project metadata too:

```bash
source .venv/bin/activate
pip install -e .
```

### Run for qubits (`d = 2`)

Run the default qubit workflow:

```bash
source .venv/bin/activate
python pyGateProperties/Gate_Schmidt_Decomposition.py
```

This prints singular values for the same qubit cases as Julia:

- `T I CXab I T CXba`
- `CZ`
- `iSWAP`
- `Product gate`
- `SWAP`
- `RXX(pi)`
- `RXX(pi/2)`
- `RXX(pi/3)`
- `RXX(pi/4)`
- `RXX(pi/8)`

### Run for qudits (`d > 2`)

Pass the local dimension as the first positional argument:

```bash
source .venv/bin/activate
python pyGateProperties/Gate_Schmidt_Decomposition.py 3

source .venv/bin/activate
python pyGateProperties/Gate_Schmidt_Decomposition.py 4
```

For `d > 2`, the Python workflow intentionally switches to standard qudit gate generalizations only. It prints singular values for:

- `SUM (CX_d)`
- `CZ_d`
- `SWAP_d`

As in Julia, qudit analogues of `iSWAP`, `RXX`, `S`, and `T` are intentionally omitted.

### Notebook / Jupytext workflow

`pyGateProperties/Gate_Schmidt_Decomposition.py` is a Jupytext percent-format script mirrored to `pyGateProperties/Gate_Schmidt_Decomposition.ipynb`.

If you edit the Python implementation, resync the notebook with:

```bash
source .venv/bin/activate
jupytext --sync pyGateProperties/Gate_Schmidt_Decomposition.py
```

If you want a clean notebook regenerated from the script, use:

```bash
source .venv/bin/activate
jupytext --to ipynb pyGateProperties/Gate_Schmidt_Decomposition.py -o pyGateProperties/Gate_Schmidt_Decomposition.ipynb
```

### Python implementation notes

- For `d = 2`, Python uses Qiskit circuits and `RXXGate(theta)`.
- For `d > 2`, Python uses explicit `NumPy` matrices rather than Qiskit qudit APIs.
- The Schmidt bipartition is computed by reshaping the operator as `(d, d, d, d)`, transposing with `(0, 2, 1, 3)`, then reshaping back to `(d * d, d * d)` before running SVD.

## Comparing Julia and Python

The two implementations are intended to match numerically.

- Julia computes singular values through `ITensors.svd(...)` with explicit tensor index partitions.
- Python computes singular values by reshaping the operator matrix into a rank-4 tensor and applying `numpy.linalg.svdvals(...)`.

For qubits, both implementations should agree on the qubit gate set and the `RXX(theta)` sweep.

For qudits, both implementations should agree on the standard-only cases:

- `SUM (CX_d)`
- `CZ_d`
- `SWAP_d`

## Verification Workflow

### Qubit cross-check

```bash
julia --project=. jlGateProperties/svd_2q_gates.jl

source .venv/bin/activate
python pyGateProperties/Gate_Schmidt_Decomposition.py
```

### Qudit cross-check

```bash
julia --project=. jlGateProperties/svd_2q_gates.jl 3

source .venv/bin/activate
python pyGateProperties/Gate_Schmidt_Decomposition.py 3
```

You can repeat the same pattern for `d = 4` or higher.

### Expected qualitative checks

- `SWAP` / `SWAP_d`: all singular values should be `1`
- `SUM (CX_d)`: the first `d` singular values should be near `sqrt(d)`, with the rest near `0`
- `CZ_d`: the same Schmidt spectrum pattern as `SUM (CX_d)`

## Limitations

- No automated test suite
- No CI workflow
- Julia dependency resolution is intentionally unpinned in git because `Manifest*.toml` is ignored
- Qudit support is standard-only and matrix-based in both languages
- Non-standard qudit analogues of `RXX`, `iSWAP`, `S`, `T`, and the qubit composite circuit are not implemented

## Developer Notes

- Re-run the relevant script after changes:
  - Julia: `julia --project=. jlGateProperties/svd_2q_gates.jl [d]`
  - Python: `source .venv/bin/activate && python pyGateProperties/Gate_Schmidt_Decomposition.py [d]`
- Keep the Python `.py` and `.ipynb` files in sync through Jupytext.
- In Julia, do not collapse or reuse the prime-chain indices (`sa1..sa4`, `sb1..sb4`) unless you also rework the tensor wiring.
- In Python, do not change the transpose order `(0, 2, 1, 3)` unless you also intend to change the Schmidt bipartition convention.
