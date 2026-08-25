# Certified phase-sweep amplification for sign changes of `psi(x)-x`

This archive verifies the new finite computation used in the accompanying
manuscript.  Starting from the public `N=21` certificate of
Grzeskowiak--Kaczorowski--Pankowski--Radziejewski (GKPR), it certifies

    2*kappa > 2789/100000 > 1/36,

and hence, by Kaczorowski's reduction,

    liminf V(T)/log T > gamma_0/pi + 2789/100000
                      > gamma_0/pi + 1/36.

## Logical scope

The proof is modular.

* The published GKPR paper and its Mathematica/Arb notebook certify the
  `N=21` contour, tail domination, Diophantine erosion, and the inclusion of
  the tabulated product cells in the eroded admissible phase set.
* The verifier in this archive checks every additional deduction introduced by
  the phase sweep: the restricted slice, the exact directed convolution, the
  seam loss, the common outer width, and the final rational inequalities.
* It does **not** independently reconstruct the upstream contour, tail, zeta
  zero, or lattice proof.  Exact copies of the imported data are fixed by
  SHA-256 hashes.

No Monte Carlo, Sobol sampling, floating-point volume estimate, root search, or
heuristic lattice reduction occurs in the new trusted calculation.

## Certified parameters and inequalities

The phase sweep uses

    M = 24000,  h = 1/(2M) = 1/48000,
    delta = 132737/1000000,
    epsilon_0 = 332/10^15.

For coordinates `2,...,21`, let

    P_n(X) = sum_{k=1}^{16000} (x_{n,k}-x_{n,k-1}) X^k.

The required coefficient sum is the sum of degrees at most `15999` in
`P_2 ... P_21`.  Each input and each intermediate coefficient is rounded
downward at scale `2^B`.  Because all coefficients are nonnegative, the result
is a rigorous lower bound.  Kronecker blocks have width 256 bits; the verifier
checks that all fixed-point coefficients are less than one, so

    16001 * 2^(2B) < 2^256

prevents carries between blocks for `B=80` and `B=112`.

The exact output integers are

    I_80  = 16808857862779264,
    I_112 = 72193494803811619816853673,

where the swept phase mass is at least `I_B / 2^(B-20)`.  The seam loss is at
most

    13944/10^15,

and the outer width is bounded by

    W_0 = 132737/10^6 + (23999/24000)*(44000/48209).

The `B=80` output already proves

    2*(I_80/2^60 - 13944/10^15)/W_0 > 2789/100000.

The exact positive margin is

    54622583536897606889634853
    --------------------------------
    20779650200857873481728000000000.

The coarser mass bound `1457/100000` independently gives the clean target
`1/36`, with exact positive margin

    8244591184361251 / 85045417930687500000.

## Files

* `verify_phase_sweep.c` - C11/GMP verifier.
* `independent_check.py` - independently written Python big-integer
  implementation of the parser and convolution.
* `upstream_data/x1.txt`, ..., `x21.txt` - exact rational phase endpoints.
* `upstream_data/zetaZeros22.txt` - upstream interval data used for the coarse
  ordinate bounds.
* `upstream_data/Verification.nb` - upstream rigorous Mathematica notebook.
* `certificate.json` - machine-readable parameters and exact outputs.
* `results/` - reference verifier outputs.
* `SHA256SUMS` - hashes of all proof inputs, sources, and reference outputs.

## Build and run

Requirements: a C11 compiler, GNU Make, GMP, and Python 3.

Run the C verifier and compiler/precision cross-checks with

    ./run_verifier.sh

Run the independent exact-integer implementation with

    ./run_crosscheck.sh

or separately with

    python3 independent_check.py 80
    python3 independent_check.py 112

Every successful run ends with `RESULT=PASS`.  No network access is required.

## Upstream provenance

The upstream data are from

    https://github.com/maciej-radziejewski/sign-changes

and are used under the scope described in `UPSTREAM.md`.  The local copies in
this archive, rather than mutable remote files, are the inputs to the verifier.
