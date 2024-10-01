#include <flint/arb.h>
#include <flint/fmpq.h>
#include <time.h>
#include <stdlib.h>

void sqrt3(fmpq_t out, fmpq_t err_inv, fmpz_t x, fmpz_t y, int n, int debug) {
    fmpq_t a;
    fmpq_t c;
    fmpq_init(a);
    fmpq_init(c);
    fmpq_one(c);
    fmpz_t one;
    fmpz_t b;
    fmpz_t d;
    fmpz_init(one);
    fmpz_init(b);
    fmpz_init(d);
    fmpq_one(a);
    if(debug) {
        fmpz_print(x);
        flint_printf("\t");
        fmpz_print(y);
        flint_printf("\n");
        flint_printf("*\n1\t1\n");
    }
    for(ulong j = 1; j <= n; ++j) {
        fmpz_one(one);
        for(ulong i = 1; i <= j; ++i) {
            fmpz_fac_ui(b, j+i);
            fmpz_fac_ui(d, j-i);
            fmpz_divexact(b, b, d);
            fmpz_pow_ui(d, x, i);
            fmpz_mul(b, b, d);
            fmpz_mul_2exp(b, b, i-1);
            fmpz_fac_ui(d, i+i);
            fmpz_divexact(b, b, d);
            fmpz_add(one, one, b);
        }
        if(j==n) {
            fmpq_set_fmpz_frac(err_inv, y, x);
            //fmpq_add_ui(err_inv, err_inv, 1);
            fmpq_mul_fmpz(err_inv, err_inv, one);
        }
        if(debug) {
            flint_printf("1\t");
            fmpz_print(one);
            flint_printf("\n");
        }
        fmpq_set_fmpz(c, one);
        fmpq_inv(c, c);
        fmpq_add(a, a, c);
    }
    fmpq_mul_fmpz(a, a, x);
    fmpq_div_fmpz(a, a, y);
    fmpq_set(out, a);

    fmpq_clear(a);
    fmpq_clear(c);
    fmpz_clear(one);
    fmpz_clear(b);
    fmpz_clear(d);
}
void sqrt3_rec(fmpq_t out, fmpq_t err_inv, fmpz_t x, fmpz_t y, int n, int debug) {
    fmpq_t a, b, c, d;
    fmpq_init(a);
    fmpq_init(c);
    fmpq_init(d);
    fmpq_one(c);
    fmpz_t one;
    fmpz_init(one);
    fmpq_init(b);
    fmpq_zero(a);
    fmpq_set_fmpz(b, x);
    fmpq_set_fmpz(c, x);
    fmpq_add_ui(c, c, 1);
    fmpq_div(c, b, c);
    fmpq_add(a, a, b);
    fmpq_add(a, a, c);
    if(debug) {
        fmpq_set_fmpz_frac(d, x, y);
        fmpq_canonicalise(d);
        fmpz_print(fmpq_numerator_ptr(d));
        flint_printf("\t");
        fmpz_print(fmpq_denominator_ptr(d));
        flint_printf("\n");
        flint_printf("*\n1\t1\n");
    }
    for(ulong i = 2; i <= n; ++i) {
        // a(k+1) = a(k)^2 / (a(k-1)*(1 + a(k)))
        fmpq_add_ui(d, c, 1);
        fmpq_mul(d, d, b);
        fmpq_swap(b, c);
        fmpq_mul(c, b, b);
        fmpq_div(c, c, d);
        fmpq_add(a, a, c);
        if(debug) {
            flint_printf("1\t");
            fmpz_print(fmpq_denominator_ptr(c));
            flint_printf("\n");
        }
    }
    fmpq_inv(err_inv, c);
    fmpq_mul_fmpz(err_inv, err_inv, y);
    fmpq_div_fmpz(out, a, y);

    fmpq_clear(a);
    fmpq_clear(c);
    fmpz_clear(one);
    fmpq_clear(b);
    fmpq_clear(d);
}
//courtesy of Pim Spelier: https://github.com/pimsp/PWS_chakravala
void pell(fmpz_t a_out, fmpz_t b_out, fmpz_t d) {
    fmpz_t a, b, k, m_0, m_1, p_1, m_2, p_2, m, tmp;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(k);
    fmpz_init(m_0);
    fmpz_init(m_1);
    fmpz_init(p_1);
    fmpz_init(m_2);
    fmpz_init(p_2);
    fmpz_init(m);
    fmpz_init(tmp);

    //k = a**2-d*b**2
    fmpz_sqrt(a, d);
    fmpz_one(b);
    fmpz_add(a, a, b);
    fmpz_mul(k, a, a);
    fmpz_mul(m, b, b);
    fmpz_mul(m, m, d);
    fmpz_sub(k, k, m);


    while (!fmpz_is_one(k)) {
        fmpz_invmod(m_0, b, k);
        fmpz_mul(m_0, m_0, a);
        fmpz_neg(m_0, m_0);
        fmpz_mod(m_0, m_0, k);

        //m_1 = (int(d**0.5+0.5)//k)*k + m_0
        fmpz_sqrt(m_1, d);
        fmpz_add(m_1, m_1, d);
        fmpz_sqrt(m_1, m_1);

        fmpz_add(m_0, m_0, m_1);
        fmpz_fdiv_r(m_1, m_1, k);
        fmpz_sub(m_1, m_0, m_1);

        fmpz_mul(p_1, m_1, m_1);
        //p_1 = m_1**2-d
        fmpz_sub(p_1, p_1, d);

        fmpz_zero(m);
        if (fmpz_cmp(p_1, m) > 0) {
            fmpz_abs(m, k);
            fmpz_sub(m_2, m_1, m);
        } else {
            fmpz_abs(m, k);
            fmpz_add(m_2, m_1, m);
        }
        fmpz_mul(p_2, m_2, m_2);
        fmpz_sub(p_2, p_2, d);
        if (fmpz_cmpabs(p_1, p_2) < 0) {
            fmpz_set(m, m_1);
        }
        else {
            fmpz_set(m, m_2);
        }
        //a = (a*m+d*b)//abs(k)
        fmpz_abs(p_2, k);
        fmpz_fmma(tmp, a,m,d,b);
        fmpz_divexact(tmp, tmp, p_2);

        //b = (a+b*m)//abs(k)
        fmpz_mul(b, b, m);
        fmpz_add(b, b, a);
        fmpz_divexact(b, b, p_2);

        //k = (m**2-d)//k
        fmpz_mul(p_2, m, m);
        fmpz_sub(p_2, p_2, d);
        fmpz_divexact(k, p_2, k);

        fmpz_set(a, tmp);
    }
    fmpz_set(a_out, a);
    fmpz_set(b_out, b);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(k);
    fmpz_clear(m_0);
    fmpz_clear(m_1);
    fmpz_clear(p_1);
    fmpz_clear(m_2);
    fmpz_clear(p_2);
    fmpz_clear(m);
    fmpz_clear(tmp);
}

int main(int argc, char *argv[])
{
    if(argc < 2) {
        flint_printf("Usage: ./sqrt <n> [<k> [<r>]]"
                     "\n\tn ... number to take square root of"
                     "\n\tk ... egyptian summand count => precision (default 10)"
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

    fmpz_t x, y;
    fmpq_t err;
    fmpz_init(x);
    fmpz_init(y);
    fmpq_init(err);
    fmpz_one(y);
    fmpq_t xx;
    fmpq_init(xx);

    int n = 10;
    if (argc > 2) {
        n = atol(argv[2]);
        if(n < 2) n = 2;
    }
    pell(x, y, x0i);
    fmpz_sub_ui(x, x, 1);
    sqrt3_rec(xx, err, x, y, n, 1);

    int rep = 1;
    if(argc > 3) {
        rep = atol(argv[3]);
    }
    if(rep < 1) {
        return 0;
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    pell(x, y, x0i);
    fmpz_sub_ui(x, x, 1);

    for(ulong i = 0; i< rep; ++i) {
        sqrt3_rec(xx, err, x, y, n, 0);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    fmpz_add_ui(x, x, 1);
    flint_fprintf(stderr,"pell sqrt\t");
    fmpz_fprint(stderr,x);
    flint_fprintf(stderr,"^2 - ");
    fmpz_fprint(stderr,x0i);
    flint_fprintf(stderr," * ");
    fmpz_fprint(stderr,y);
    flint_fprintf(stderr,"^2 == 1\t");

    long delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;

    const long LOW_PREC = 16;

    arb_t tmp;
    arb_init(tmp);
    arb_set_fmpq(a, err, LOW_PREC);

    arb_log(a, a, LOW_PREC);

    arb_set_fmpz(x0, x0i);
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
    arb_set_fmpq(a, xx, bin_prec);

    fmpq_fprint(stderr, xx);
    flint_fprintf(stderr, "\n");

    flint_fprintf(stderr, "pell sqrt\t% 10d microseconds per %d repetition(s)\t", delta_us, rep);

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    for(long i = 0; i < rep; ++i) {
        arb_set_fmpz(x0, x0i);
        arb_sqrt(x0, x0, bin_prec);
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
    fmpz_clear(x);
    fmpz_clear(y);
    fmpq_clear(err);
    arb_clear(tmp);
    arb_clear(a);
    arb_clear(b);
    fmpz_clear(sqrtx0i);
    return 0;
}
