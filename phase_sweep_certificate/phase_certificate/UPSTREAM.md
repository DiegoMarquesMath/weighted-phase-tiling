# Upstream certificate used by the phase sweep

The new proof imports the following statements from the published GKPR paper
and its public `Verification.nb` notebook.

1. The phase set associated with the published `N=21` contour is admissible:
   every certified orbit hit produces a zero of the complete normalized
   `k`-function.
2. The product cells encoded by `x1.txt`, ..., `x21.txt` lie in the published
   Diophantine erosion of that admissible phase set.
3. The zeta-zero intervals, tail bound, contour inequalities, and lattice
   inequalities used to prove (1) and (2) were checked in the upstream
   Mathematica/Arb certificate.

The present C and Python verifiers do not duplicate those analytic checks.
They verify all new consequences of the imported cells: the first-coordinate
slice, the 24,000 phase translates, the conservative fundamental-cube seam
loss, the fixed outer-width majorant, the exact volume convolution, and the
final rational bounds.

Upstream URL:

    https://github.com/maciej-radziejewski/sign-changes

The exact local files used by this proof are fixed by `SHA256SUMS`.
