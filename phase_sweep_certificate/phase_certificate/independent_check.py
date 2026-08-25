#!/usr/bin/env python3
"""Independent exact-integer cross-check for the directed convolution.

This implementation shares no C parser or convolution code with
verify_phase_sweep.c.  It reads the exact rational endpoint files, performs
Kronecker-substitution convolution using Python's arbitrary-precision integers,
and checks the exact final targets.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys

L = 16000
BLOCK_BITS = 256
BLOCK_MASK = (1 << BLOCK_BITS) - 1
EXPECTED = {
    80: 16808857862779264,
    112: 72193494803811619816853673,
}


def read_poly(root: Path, n: int, scale_bits: int) -> int:
    scale = 1 << scale_bits
    packed = 0
    prev = Fraction(0)
    lines = (root / f"x{n}.txt").read_text(encoding="ascii").splitlines()
    if len(lines) != L:
        raise ValueError(f"x{n}.txt has {len(lines)} entries, expected {L}")
    for k, line in enumerate(lines, 1):
        cur = Fraction(line)
        if cur < prev:
            raise ValueError(f"x{n}.txt decreases at entry {k}")
        diff = cur - prev
        coeff = (diff.numerator * scale) // diff.denominator
        if not 0 <= coeff < scale:
            raise ValueError(f"coefficient out of range in x{n}.txt at {k}")
        packed |= coeff << (BLOCK_BITS * k)
        prev = cur
    return packed


def multiply_and_round(a: int, b: int, scale_bits: int) -> int:
    scale = 1 << scale_bits
    raw_product = a * b
    packed = 0
    for k in range(L + 1):
        raw_coeff = (raw_product >> (BLOCK_BITS * k)) & BLOCK_MASK
        coeff = raw_coeff >> scale_bits
        if coeff >= scale:
            raise ValueError(f"intermediate coefficient at degree {k} is not < 1")
        packed |= coeff << (BLOCK_BITS * k)
    return packed


def main() -> None:
    bits = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    if bits not in EXPECTED:
        raise SystemExit("usage: independent_check.py {80|112}")
    root = Path(__file__).resolve().parent / "upstream_data"

    product = read_poly(root, 2, bits)
    for n in range(3, 22):
        product = multiply_and_round(product, read_poly(root, n, bits), bits)

    numerator = sum((product >> (BLOCK_BITS * k)) & BLOCK_MASK for k in range(L))
    if numerator != EXPECTED[bits]:
        raise AssertionError(
            f"unexpected I_{bits}: {numerator}; expected {EXPECTED[bits]}"
        )

    phase_mass = Fraction(numerator, 1 << (bits - 20))
    seam = Fraction(13944, 10**15)
    width = Fraction(132737, 10**6) + Fraction(23999, 24000) * Fraction(44000, 48209)
    gain = 2 * (phase_mass - seam) / width
    if gain <= Fraction(2789, 100000):
        raise AssertionError("strong target 2789/100000 was not certified")
    if gain <= Fraction(1, 36):
        raise AssertionError("clean target 1/36 was not certified")

    print(f"B={bits}")
    print(f"I_B={numerator}")
    print(f"phase_mass={phase_mass}")
    print(f"2*kappa_lower={gain}")
    print(f"margin_over_2789/100000={gain - Fraction(2789, 100000)}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
