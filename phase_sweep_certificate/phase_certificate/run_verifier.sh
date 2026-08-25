#!/usr/bin/env sh
set -eu
trap 'make clean >/dev/null 2>&1 || true' EXIT
sha256sum -c SHA256SUMS
make clean
make check
make crosscheck
if command -v clang >/dev/null 2>&1; then
    make compiler-crosscheck
fi
