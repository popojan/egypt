#include <flint/arb.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <stdlib.h>

void sqrtr(fmpq_t out, fmpq_t err_inv, fmpz_t nn, long k, long m) {
    fmpq_t xm1y, xm1;
    fmpz_t n, nnm;
    fmpq_init(xm1y);
    fmpq_init(xm1);
    fmpz_init(n);
    fmpz_init(nnm);
    fmpz_set(nnm, nn);
    if (m > 0)
    {
        fmpz_set_ui(n, m);
        fmpz_mul_ui(n, n, 2);
        fmpz_add_ui(n, n, 1);
        fmpz_pow_fmpz(nnm, nn, n);
    } else {
        fmpz_set(nnm, nn);
    }
    fmpz_sqrt(n, nnm);
    fmpz_add_ui(n, n, 1);
    fmpq_set_fmpz_frac(xm1y, nnm, n);
    fmpq_set_fmpz(xm1, nnm);
    fmpq_add_fmpz(xm1, xm1,nnm);
    fmpz_mul(n, n, n);
    fmpz_sub(n, n, nnm);
    fmpq_div_fmpz(xm1, xm1, n);
    fmpq_set_fmpz(xm1y, n);
    fmpz_sqrt(n, nnm);
    fmpz_add_ui(n, n, 1);
    fmpz_mul_ui(n, n, 2);
    fmpq_div_fmpz(xm1y, xm1y, n);

    fmpq_t a, b, c, d;
    fmpq_init(a);
    fmpq_init(c);
    fmpq_init(d);
    fmpq_one(c);
    fmpz_t one;
    fmpz_init(one);
    fmpq_init(b);
    fmpq_zero(a);
    fmpq_set(b, xm1);
    fmpq_set(c, xm1);
    fmpq_add_ui(c, c, 1);
    fmpq_div(c, b, c);
    fmpq_add(a, a, b);
    fmpq_add(a, a, c);
    for(ulong i = 2; i <= k; ++i) {
        // a(k+1) = a(k)^2 / (a(k-1)*(1 + a(k)))
        fmpq_add_ui(d, c, 1);
        fmpq_mul(d, d, b);
        fmpq_swap(b, c);
        fmpq_mul(c, b, b);
        fmpq_div(c, c, d);
        fmpq_add(a, a, c);
    }
    if (m  > 0) {
        fmpz_pow_ui(n, nn, m);
        fmpq_div_fmpz(a, a, n);
    }
    fmpq_inv(err_inv, c);
    fmpq_mul(err_inv, err_inv, xm1y);
    fmpq_mul(out, a, xm1y);

    fmpq_clear(a);
    fmpq_clear(c);
    fmpz_clear(one);
    fmpq_clear(b);
    fmpq_clear(d);
    fmpq_clear(xm1y);
    fmpq_clear(xm1);
    fmpz_clear(n);
    fmpz_clear(nnm);
}

int main(int argc, char *argv[])
{
    if(argc < 2) {
        flint_printf("Usage: ./sqrt <n> [<k> [<m>[<r>]]]"
                     "\n\tn ... number to take square root of"
                     "\n\tk ... egyptian summand count => precision (default 10)"
                     "\n\tm ... extra exponent to speedup convergence (default 0)"
                     "\n\tr ... number of repetitions to time (default 1)"
                     "\n");
        return 0;
    }
    fmpz_t x0i, sqrtx0i;
    fmpz_init(x0i);
    fmpz_init(sqrtx0i);
    arb_t x0, a, b;
    arb_init(x0);
    arb_init(a);
    arb_init(b);
    fmpz_set_str(x0i, argv[1], 10);

    if(fmpz_is_square(x0i)){
        fmpz_sqrt(sqrtx0i, x0i);
        fmpz_print(sqrtx0i);
        flint_printf("\t1\n");
        fmpz_fprint(stderr,x0i);
        flint_fprintf(stderr, " is perfect square of ");
        fmpz_fprint(stderr,sqrtx0i);
        flint_fprintf(stderr, "\n");
        return 1;
    }

    fmpq_t err;
    fmpq_init(err);
    fmpq_t xx;
    fmpq_init(xx);

    long n = 10, m = 0;
    if (argc > 2) {
        n = atol(argv[2]);
        if(n < 2) n = 2;
    }

    if (argc > 3) {
        m = atol(argv[3]);
        if(m < 0) m = 0;
    }
    //pell_rational(x, y, x0i);
    //fmpz_sub_ui(x, x, 1);
    //sqrtr(xx, err, x0i, n, m);

    int rep = 1;
    if(argc > 4) {
        rep = atol(argv[4]);
    }
    if(rep < 1) {
        return 0;
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    for(ulong i = 0; i< rep; ++i) {
        sqrtr(xx, err, x0i, n, m);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    long delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;

    const long LOW_PREC = 16;

    arb_t tmp;
    arb_init(tmp);
    arb_set_fmpq(a, err, LOW_PREC);

    arb_log(a, a, LOW_PREC);

    arb_set_fmpz(x0, x0i);

    fmpq_fprint(stderr, xx);
    flint_fprintf(stderr, "\n");

    flint_fprintf(stderr, "pell sqrt\t% 10d microseconds per %d repetition(s)\t", delta_us, rep);

    arb_t x1;
    arb_init(x1);
    arb_set_fmpz(x0, x0i);
    arb_set(x1, x0);

    fmpz_sqrt(sqrtx0i, x0i);
    arb_set_fmpz(tmp, sqrtx0i);
    arb_log(tmp, tmp, LOW_PREC);
    arb_add(a, a, tmp, LOW_PREC);
    arb_log_ui(tmp, 10, LOW_PREC);
    arb_div(a, a, tmp, LOW_PREC);

    ulong prec = (long)(arf_get_d(arb_midref(a),0))+1;
    arb_mul(a, a, tmp, LOW_PREC);
    arb_log_ui(tmp, 2, LOW_PREC);
    arb_div(a, a, tmp, LOW_PREC);
    long bin_prec = (long)(arf_get_d(arb_midref(a),0))+1;
    if (bin_prec < 2) bin_prec = 2;
    if (prec < 2) prec = 2;

    arb_set_fmpq(a, xx, prec);

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    for(long i = 0; i < rep; ++i) {
        arb_sqrt(x0, x1, bin_prec);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;

    arb_fprintn(stderr,a, prec-1, 0);
    flint_fprintf(stderr, "\n");

    flint_fprintf(stderr," arb sqrt\t% 10d microseconds per %d repetition(s)\t", delta_us, rep);

    arb_fprintn(stderr,x0, prec - 1, 0);
    flint_fprintf(stderr, "\n");

    flint_cleanup();
    arb_clear(x0);
    arb_clear(x1);
    fmpq_clear(err);
    arb_clear(tmp);
    arb_clear(a);
    arb_clear(b);
    fmpz_clear(sqrtx0i);
    return 0;
}
