#define _GNU_SOURCE
#include <stdio.h>
#include <gmp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define L 16000
#define N 21
#ifndef SCALE_BITS
#define SCALE_BITS 112
#endif
#define BLOCK_WORDS 4

#if 2 * SCALE_BITS + 15 >= 256
#error "SCALE_BITS too large for the 256-bit Kronecker blocks"
#endif

typedef struct { uint64_t lo, hi; } u128s;

static void die(const char *msg)
{
    fprintf(stderr, "VERIFIER ERROR: %s\n", msg);
    exit(EXIT_FAILURE);
}

static int equal128(u128s a, u128s b)
{
    return a.lo == b.lo && a.hi == b.hi;
}

static int cmp128(u128s a, u128s b)
{
    if (a.hi < b.hi) return -1;
    if (a.hi > b.hi) return 1;
    if (a.lo < b.lo) return -1;
    if (a.lo > b.lo) return 1;
    return 0;
}

static u128s pow2_128(unsigned e)
{
    u128s r = {0, 0};
    if (e < 64) r.lo = UINT64_C(1) << e;
    else if (e < 128) r.hi = UINT64_C(1) << (e - 64);
    else die("power does not fit 128 bits");
    return r;
}

static u128s mpz_to_u128(const mpz_t z)
{
    if (mpz_sgn(z) < 0 || mpz_sizeinbase(z, 2) > 128)
        die("integer does not fit 128 bits");
    uint64_t w[2] = {0, 0};
    size_t count = 0;
    mpz_export(w, &count, -1, sizeof(uint64_t), 0, 0, z);
    u128s r = {w[0], w[1]};
    return r;
}

static void u128_to_mpz(mpz_t z, u128s a)
{
    const uint64_t w[2] = {a.lo, a.hi};
    mpz_import(z, 2, -1, sizeof(uint64_t), 0, 0, w);
}

/* Divide a 256-bit base digit by 2^SCALE_BITS, rounding downward. */
static u128s shr_block(const uint64_t *w)
{
    const unsigned limb = SCALE_BITS / 64;
    const unsigned off = SCALE_BITS % 64;
    const uint64_t a0 = limb < 4 ? w[limb] : 0;
    const uint64_t a1 = limb + 1 < 4 ? w[limb + 1] : 0;
    const uint64_t a2 = limb + 2 < 4 ? w[limb + 2] : 0;
    u128s r;
    if (off == 0) {
        r.lo = a0;
        r.hi = a1;
    } else {
        r.lo = (a0 >> off) | (a1 << (64 - off));
        r.hi = (a1 >> off) | (a2 << (64 - off));
    }
    return r;
}

/* Encode sum a[k] X^k at X=2^256. */
static void pack_poly(mpz_t z, const u128s *a)
{
    const u128s one = pow2_128(SCALE_BITS);
    const size_t nw = (size_t)(L + 1) * BLOCK_WORDS;
    uint64_t *words = calloc(nw, sizeof(uint64_t));
    if (!words) die("allocation failure while packing a polynomial");
    for (int k = 0; k <= L; ++k) {
        if (cmp128(a[k], one) >= 0)
            die("a fixed-point polynomial coefficient is not below one");
        words[(size_t)k * BLOCK_WORDS] = a[k].lo;
        words[(size_t)k * BLOCK_WORDS + 1] = a[k].hi;
    }
    mpz_import(z, nw, -1, sizeof(uint64_t), 0, 0, words);
    free(words);
}

/*
   Exact nonnegative convolution followed coefficientwise by floor(/2^B).
   Each input coefficient is < 2^B.  A raw convolution coefficient is below
   (L+1)2^(2B) < 2^256, so Kronecker blocks cannot carry into one another.
   Therefore out/2^B is a coefficientwise lower bound for the product of the
   represented fixed-point polynomials.
*/
static void mul_fixed(u128s *out, const u128s *a, const u128s *b)
{
    mpz_t A, B, P;
    mpz_inits(A, B, P, NULL);
    pack_poly(A, a);
    pack_poly(B, b);
    mpz_mul(P, A, B);

    const size_t nw = (size_t)(2 * L + 1) * BLOCK_WORDS;
    uint64_t *words = calloc(nw, sizeof(uint64_t));
    if (!words) die("allocation failure while unpacking a product");
    size_t got = 0;
    mpz_export(words, &got, -1, sizeof(uint64_t), 0, 0, P);
    if (got > nw) die("unexpected product length");

    for (int k = 0; k <= L; ++k)
        out[k] = shr_block(words + (size_t)k * BLOCK_WORDS);

    const u128s one = pow2_128(SCALE_BITS);
    for (int k = 0; k <= L; ++k)
        if (cmp128(out[k], one) >= 0)
            die("a convolved fixed-point coefficient is not below one");

    free(words);
    mpz_clears(A, B, P, NULL);
}

static void convolution_self_test(void)
{
    u128s *a = calloc(L + 1, sizeof(*a));
    u128s *b = calloc(L + 1, sizeof(*b));
    u128s *c = calloc(L + 1, sizeof(*c));
    if (!a || !b || !c) die("allocation failure in self-test");

    a[1] = pow2_128(SCALE_BITS - 2); /* 1/4 */
    a[2] = pow2_128(SCALE_BITS - 3); /* 1/8 */
    b[1] = pow2_128(SCALE_BITS - 1); /* 1/2 */
    b[3] = pow2_128(SCALE_BITS - 4); /* 1/16 */
    mul_fixed(c, a, b);

    if (!equal128(c[2], pow2_128(SCALE_BITS - 3)) ||
        !equal128(c[3], pow2_128(SCALE_BITS - 4)) ||
        !equal128(c[4], pow2_128(SCALE_BITS - 6)) ||
        !equal128(c[5], pow2_128(SCALE_BITS - 7)))
        die("Kronecker-convolution self-test failed");

    free(a); free(b); free(c);
}

static void floor_scaled_q(u128s *r, const mpq_t q)
{
    mpz_t t, z;
    mpz_inits(t, z, NULL);
    mpz_mul_2exp(t, mpq_numref(q), SCALE_BITS);
    mpz_fdiv_q(z, t, mpq_denref(q));
    *r = mpz_to_u128(z);
    mpz_clears(t, z, NULL);
}

static void read_increment_poly(const char *path, u128s *a)
{
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(EXIT_FAILURE); }

    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    mpq_t prev, cur, diff;
    mpq_inits(prev, cur, diff, NULL);
    mpq_set_ui(prev, 0, 1);
    memset(a, 0, (L + 1) * sizeof(*a));

    int k = 1;
    while (k <= L && (len = getline(&line, &cap, f)) != -1) {
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
            line[--len] = '\0';
        if (mpq_set_str(cur, line, 10) != 0)
            die("malformed rational in an x-file");
        mpq_canonicalize(cur);
        if (mpq_cmp(cur, prev) < 0)
            die("an x-file is not nondecreasing");
        mpq_sub(diff, cur, prev);
        floor_scaled_q(&a[k], diff);
        mpq_set(prev, cur);
        ++k;
    }
    if (k != L + 1)
        die("an x-file does not contain 16000 entries");

    free(line);
    fclose(f);
    mpq_clears(prev, cur, diff, NULL);
}

static void parse_centered_interval(const char *line, mpq_t lo, mpq_t hi)
{
    const char *p = strchr(line, '[');
    if (!p) die("malformed CenteredInterval");
    ++p;
    const char *star = strstr(p, "* 2^");
    if (!star) die("missing midpoint exponent");
    const char *mend = star;
    while (mend > p && mend[-1] == ' ') --mend;
    char *ms = strndup(p, (size_t)(mend - p));
    if (!ms) die("allocation failure while parsing a midpoint");

    char *endp;
    const long em = strtol(star + 4, &endp, 10);
    const char *comma = strchr(endp, ',');
    if (!comma) die("missing interval comma");
    ++comma;
    while (*comma == ' ') ++comma;
    const char *starr = strstr(comma, "* 2^");
    if (!starr) die("missing radius exponent");
    const char *rend = starr;
    while (rend > comma && rend[-1] == ' ') --rend;
    char *rs = strndup(comma, (size_t)(rend - comma));
    if (!rs) die("allocation failure while parsing a radius");
    const long er = strtol(starr + 4, &endp, 10);

    mpz_t m, r;
    mpz_inits(m, r, NULL);
    if (mpz_set_str(m, ms, 10) != 0 || mpz_set_str(r, rs, 10) != 0)
        die("malformed ball integer");
    free(ms); free(rs);

    mpq_t qm, qr;
    mpq_inits(qm, qr, NULL);
    mpq_set_z(qm, m);
    mpq_set_z(qr, r);
    if (em < 0) mpz_mul_2exp(mpq_denref(qm), mpq_denref(qm), (unsigned)(-em));
    else mpz_mul_2exp(mpq_numref(qm), mpq_numref(qm), (unsigned)em);
    if (er < 0) mpz_mul_2exp(mpq_denref(qr), mpq_denref(qr), (unsigned)(-er));
    else mpz_mul_2exp(mpq_numref(qr), mpq_numref(qr), (unsigned)er);
    mpq_canonicalize(qm);
    mpq_canonicalize(qr);
    mpq_sub(lo, qm, qr);
    mpq_add(hi, qm, qr);

    mpq_clears(qm, qr, NULL);
    mpz_clears(m, r, NULL);
}

int main(int argc, char **argv)
{
    const char *dir = argc > 1 ? argv[1] : "upstream_data";
    const unsigned long M = 24000UL;

    convolution_self_test();

    u128s *cur = calloc(L + 1, sizeof(*cur));
    u128s *p = calloc(L + 1, sizeof(*p));
    u128s *tmp = calloc(L + 1, sizeof(*tmp));
    if (!cur || !p || !tmp) die("allocation failure for coefficient arrays");

    char path[1024];
    for (int n = 2; n <= 21; ++n) {
        snprintf(path, sizeof(path), "%s/x%d.txt", dir, n);
        read_increment_poly(path, p);
        if (n == 2) {
            memcpy(cur, p, (L + 1) * sizeof(*cur));
        } else {
            mul_fixed(tmp, cur, p);
            u128s *swap = cur; cur = tmp; tmp = swap;
            memset(tmp, 0, (L + 1) * sizeof(*tmp));
        }
    }

    /* C is a fixed-point lower bound for the rest-volume through degree L-1. */
    mpz_t Cint, z;
    mpz_inits(Cint, z, NULL);
    mpz_set_ui(Cint, 0);
    for (int k = 0; k <= L - 1; ++k) {
        u128_to_mpz(z, cur[k]);
        mpz_add(Cint, Cint, z);
    }

    /* Check h=1/(2M) <= x_{1,1}, so only the first x_1-cell is used. */
    snprintf(path, sizeof(path), "%s/x1.txt", dir);
    FILE *fx = fopen(path, "r");
    if (!fx) { perror(path); exit(EXIT_FAILURE); }
    char *line = NULL;
    size_t cap = 0;
    if (getline(&line, &cap, fx) == -1) die("empty x1 file");
    fclose(fx);
    mpq_t x11, h;
    mpq_inits(x11, h, NULL);
    if (mpq_set_str(x11, line, 10) != 0) die("malformed x1 first entry");
    mpq_canonicalize(x11);
    free(line);
    mpq_set_ui(h, 1, 2 * M);
    if (mpq_cmp(h, x11) >= 0) die("failed to prove h < x_{1,1}");

    /* Verify the coarse ordinate bounds directly from the certified balls. */
    snprintf(path, sizeof(path), "%s/zetaZeros22.txt", dir);
    FILE *fz = fopen(path, "r");
    if (!fz) { perror(path); exit(EXIT_FAILURE); }
    char *l0 = NULL, *l1 = NULL;
    size_t c0 = 0, c1 = 0;
    if (getline(&l0, &c0, fz) == -1 || getline(&l1, &c1, fz) == -1)
        die("zetaZeros22 has fewer than two entries");
    fclose(fz);

    mpq_t g0lo, g0hi, g1lo, g1hi, test;
    mpq_inits(g0lo, g0hi, g1lo, g1hi, test, NULL);
    parse_centered_interval(l0, g0lo, g0hi);
    parse_centered_interval(l1, g1lo, g1hi);
    free(l0); free(l1);
    mpq_set_ui(test, 14135, 1000);
    if (mpq_cmp(g0hi, test) >= 0) die("failed to prove gamma_0 < 14.135");
    mpq_set_ui(test, 21022, 1000);
    if (mpq_cmp(g1lo, test) <= 0) die("failed to prove gamma_1 > 21.022");

    /*
       mass = 2^20 C.  The factor is M * 2^21 * (1/(2M)).
       W < delta + ((M-1)/M) * (44/7)/(6.887), using pi < 22/7.
    */
    mpq_t mass, mass_target, seam_loss, safe_mass, safe_mass_target,
          W, delta, fac, reciprocal_theta_upper,
          gain, coarse_gain, target, margin, strong_target, strong_margin;
    mpq_inits(mass, mass_target, seam_loss, safe_mass, safe_mass_target,
              W, delta, fac, reciprocal_theta_upper,
              gain, coarse_gain, target, margin, strong_target, strong_margin,
              NULL);

    mpq_set_z(mass, Cint);
    mpz_mul_2exp(mpq_numref(mass), mpq_numref(mass), 20);
    mpz_mul_2exp(mpq_denref(mass), mpq_denref(mass), SCALE_BITS);
    mpq_canonicalize(mass);
    mpq_set_ui(mass_target, 1457, 100000); /* 0.01457 */
    if (mpq_cmp(mass, mass_target) <= 0)
        die("certified phase mass is not greater than 0.01457");

    /* Conservative loss at the seam of the fundamental cube:
       2*N*epsilon_0 with N=21 and epsilon_0=332/10^15. */
    if (mpq_set_str(seam_loss, "13944/1000000000000000", 10) != 0)
        die("failed to initialize the seam-loss rational");
    mpq_canonicalize(seam_loss);
    mpq_sub(safe_mass, mass, seam_loss);
    mpq_sub(safe_mass_target, mass_target, seam_loss);
    if (mpq_sgn(safe_mass_target) <= 0)
        die("seam correction exhausts the coarse phase mass");

    mpq_set_ui(delta, 132737, 1000000);
    mpq_set_ui(fac, M - 1, M);
    mpq_set_ui(reciprocal_theta_upper, 44000, 48209);
    mpq_mul(W, fac, reciprocal_theta_upper);
    mpq_add(W, W, delta);

    mpq_mul_2exp(gain, safe_mass, 1);
    mpq_div(gain, gain, W);

    /* A compact exact claim used in the paper: replace the computed mass by
       1457/100000 and subtract the same exact seam correction. */
    mpq_mul_2exp(coarse_gain, safe_mass_target, 1);
    mpq_div(coarse_gain, coarse_gain, W);
    mpq_set_ui(target, 1, 36);
    mpq_set_ui(strong_target, 2789, 100000);
    if (mpq_cmp(gain, strong_target) <= 0)
        die("computed lower bound is not greater than 2789/100000");
    if (mpq_cmp(gain, target) <= 0)
        die("computed lower bound is not greater than 1/36");
    if (mpq_cmp(coarse_gain, target) <= 0)
        die("coarse paper inequality is not greater than 1/36");
    mpq_sub(strong_margin, gain, strong_target);
    mpq_sub(margin, coarse_gain, target);

    printf("weighted phase-sweep certificate\n");
    printf("fixed-point bits       : %d\n", SCALE_BITS);
    printf("number of phase charts: %lu\n", M);
    printf("h                      : 1/%lu\n", 2 * M);
    printf("C_lower                : %.17g\n", mpq_get_d(mass) / (double)(1UL << 20));
    printf("phase_mass_lower       : %.17g\n", mpq_get_d(mass));
    printf("seam_loss_upper        : %.17g\n", mpq_get_d(seam_loss));
    printf("safe_phase_mass_lower  : %.17g\n", mpq_get_d(safe_mass));
    printf("outer_width_upper      : %.17g\n", mpq_get_d(W));
    printf("2*kappa lower          : %.17g\n", mpq_get_d(gain));
    printf("coarse 2*kappa lower   : %.17g\n", mpq_get_d(coarse_gain));
    printf("strong target          : 2789/100000 = %.17g\n", 2789.0 / 100000.0);
    printf("strong exact margin    : ");
    mpq_out_str(stdout, 10, strong_margin);
    printf("\n");
    printf("clean target           : 1/36 = %.17g\n", 1.0 / 36.0);
    printf("coarse exact margin    : ");
    mpq_out_str(stdout, 10, margin);
    printf("\n");
    printf("RESULT                 : PASS\n");

    mpq_clears(mass, mass_target, seam_loss, safe_mass, safe_mass_target,
               W, delta, fac, reciprocal_theta_upper,
               gain, coarse_gain, target, margin, strong_target, strong_margin,
               g0lo, g0hi, g1lo, g1hi, test, x11, h, NULL);
    mpz_clears(Cint, z, NULL);
    free(cur); free(p); free(tmp);
    return EXIT_SUCCESS;
}
