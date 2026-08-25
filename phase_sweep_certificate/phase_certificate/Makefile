CC ?= cc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra
LDLIBS ?= -lgmp

all: verify_phase_sweep

verify_phase_sweep: verify_phase_sweep.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

check: verify_phase_sweep
	./verify_phase_sweep upstream_data

crosscheck:
	$(CC) $(CFLAGS) -DSCALE_BITS=80 -o verify_phase_sweep_80 verify_phase_sweep.c $(LDLIBS)
	./verify_phase_sweep_80 upstream_data
	$(CC) $(CFLAGS) -DSCALE_BITS=112 -o verify_phase_sweep_112 verify_phase_sweep.c $(LDLIBS)
	./verify_phase_sweep_112 upstream_data

python-crosscheck:
	python3 independent_check.py 80
	python3 independent_check.py 112

compiler-crosscheck:
	clang $(CFLAGS) -DSCALE_BITS=112 -o verify_phase_sweep_clang verify_phase_sweep.c $(LDLIBS)
	./verify_phase_sweep_clang upstream_data

clean:
	rm -f verify_phase_sweep verify_phase_sweep_80 verify_phase_sweep_112 \
	      verify_phase_sweep_clang verify_phase_sweep_san
