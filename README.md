# Invariants code

Symbolic computation, in Python with SymPy, of the coefficient matrix that links
the invariants of an *n*-partite system to its correlation functions on equality
patterns. The matrix is written in two reduced bases:

- the **orbit basis**, one basis vector per orbit of S₃ acting on equality patterns;
- the **irrep basis**, adapted to the decomposition R[Π] = U ⊕ W into the trivial
  and the non-trivial part of the representation. In this basis the coefficient
  matrix is block diagonal.

This notebook accompanies the article: Third-Order Quantum Information from Distributed Randomized Measurements.

## Contents

| File | Description |
|---|---|
| `Invariants_code.ipynb` | Jupyter notebook with the full construction |

## Requirements

- Python 3
- [SymPy](https://www.sympy.org)
- Jupyter (Notebook or Lab)

```bash
pip install sympy jupyter
```

Tested with Python 3.13 and SymPy 1.14.

## Usage

Open the notebook, run all cells **from the top, in order**, then call the
functions with the number of parties `n`:

```python
n = 2

M(n)          # full coefficient matrix  (5**n x 5**n)
M_orb(n)      # coefficient matrix in the orbit basis
M_irrep(n)    # coefficient matrix in the irrep basis (block diagonal)
```

## Mathematical setup

**Single party.**

- G = S₃ = {E, T₁, T₂, T₃, C₁, C₂}: the identity, three transpositions and two 3-cycles.
- Π = {D, T₁, T₂, T₃, Eq}: the five equality patterns (all distinct, the three
  patterns with exactly one pair equal, all equal). `gam[g]` is the 5×5 permutation
  matrix by which g acts on R[Π].
- `PU = (1/6) Σ_g gam[g]` projects onto the invariant subspace U, and `PW = I − PU`
  projects onto its orthogonal complement W.
- `S` (5×6) maps R[S₃] to R[Π]. `W` is the Weingarten matrix with parameters e, t, c,
  and s = e + 3t + 2c.
- `Cmap` (5×6) merges the two 3-cycles C₁ and C₂, because the orientation of a cycle
  cannot be distinguished. Since the C₁ and C₂ columns of S·W coincide,
  S·W = SW · Cmap, where `SW` is S·W with its last column removed.

**n parties.**

$$y = \Big(\bigotimes_{i=1}^n S W_i\Big)\,x = M(n)\,z,\qquad z = C_n\,x,\qquad M(n) = \bigotimes_{i=1}^n SW_i$$

where each party *i* has its own symbols e_i, t_i, c_i, s_i.

**Orbit basis.** The columns of `R_orb(n)` are the orbit vectors
r_O = Σ_{π∈O} e_π. They are orthogonal, R_orbᵀR_orb = diag(|O|), and
L_orb = diag(1/|O|) R_orbᵀ is a left inverse. In orbit coordinates
(y_orb = L_orb y, z_orb = L_orb z):

$$y_{orb} = M_{orb}\,z_{orb},\qquad M_{orb} = L_{orb}\,M(n)\,R_{orb}$$

**Irrep basis.** `R_irrep(n)` collects the invariant vectors of every block of
(U ⊕ W)^{⊗n}, where each party's factor is either U or W, in a fixed order of the
parties. With Q = L_orb R_irrep (so that z_orb = Q z_irrep):

$$M_{irrep} = Q^{-1} M_{orb}\, Q$$

The entries of `M_irrep` are factored, and e_i + 3t_i + 2c_i is replaced by s_i.

## Main functions

| Function | Output | Size |
|---|---|---|
| `x_full(n)` | symbolic vector of invariants x | 6ⁿ |
| `Cn(n)` | Cmap^{⊗n} | 5ⁿ × 6ⁿ |
| `z(n)` | reduced invariants z = Cn(n)·x | 5ⁿ |
| `y_patt(n)` | symbolic vector of correlation functions y | 5ⁿ |
| `M(n)` | SW₁ ⊗ … ⊗ SWₙ | 5ⁿ × 5ⁿ |
| `P(n)` | projector onto the S₃-invariant subspace | 5ⁿ × 5ⁿ |
| `R_orb(n)`, `L_orb(n)` | orbit basis and its left inverse | 5ⁿ × N, N × 5ⁿ |
| `R_irrep(n)` | irrep-adapted basis of the invariant subspace | 5ⁿ × N |
| `Q(n)` | change of basis, irrep → orbit coordinates | N × N |
| `M_orb(n)` | coefficient matrix in the orbit basis | N × N |
| `M_irrep(n)` | coefficient matrix in the irrep basis | N × N |

N is the number of orbits, N = (5ⁿ + 3·3ⁿ + 2·2ⁿ)/6: N = 3, 10, 37, 150 for
n = 1, 2, 3, 4.

## Conventions

- Tensor products follow the Kronecker convention: the first party is the slowest
  index, the same order as `itertools.product`.
- The columns of `R_orb` are ordered by the lexicographically first index tuple of
  each orbit, with respect to the basis order (D, T₁, T₂, T₃, Eq).
- The columns of `R_irrep` are grouped by the number k of parties in W
  (k = 0, …, n) and, within each k, by the positions of those parties.
- In the irrep basis, blocks with exactly one party in W contain no invariants.
  For example, `M_irrep(2)` splits into blocks of size 9 (U⊗U) and 1 (W⊗W), and
  `M_irrep(3)` into blocks of size 27, 3, 3, 3 and 1.

## Notes

- **Performance.** All matrices are built as dense symbolic matrices, so the cost
  grows quickly with n. `M_irrep(3)` takes a few seconds; n ≥ 4 is considerably slower.
- **Name clashes.** `from sympy import *` is used, and some names in the notebook
  (`E`, `S`, `Eq`, `Q`) replace SymPy's objects of the same name. Avoid naming
  variables `gamma`, `C` or other SymPy names. If you get unexpected errors,
  restart the kernel and run all cells from the top.

## Author
Corrado Cosimo Mattiaci, Institute of Physics, Faculty of Physics, Astronomy and Informatics, 
Nicolaus Copernicus University, Grudziądzka 5/7, 87-100 Toruń, Poland
corrmatt@doktorant.umk.pl

## License

[e.g. MIT]
