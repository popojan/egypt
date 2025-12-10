# Egyptian Fractions

Fast algorithm for representing rational numbers as egyptian fractions.

## Properties
* suitable for very large inputs
* returns rather small denominators

## Usage

```
Egyptian Fractions

Usage: egypt [OPTIONS] [NUMERATOR] [DENOMINATOR]

Arguments:
  [NUMERATOR]    [default: 1]
  [DENOMINATOR]  [default: 1]

Options:
  -r, --reverse        Reverse merge strategy
  -m, --merge          Extra O(n^2) merge step possibly reducing number of terms
      --raw            Output minimal number of raw quadruplets (aka symbolic sums)
      --bisect         Output raw quadruplets bisected according to --limit
  -s, --silent         No output
      --batch          Batch mode (expects numerator and denominator on each line of stdin)
  -l, --limit <LIMIT>  Maximum number of terms for breaking large symbolic sums [default: 8]
  -h, --help           Print help
  -V, --version        Print version
```

## Performance
```
$ time ./egypt -s '2 9689 ^ 1 -' '2 9941 ^ 1 -'

real    0m0.345s
user    0m0.296s
sys     0m0.047s
```

```
$ time ./egypt -s 162259276829213363391578010288127 170141183460469231731687303715884105727

real    0m0.002s
user    0m0.001s
sys     0m0.001s
```

```
$ time ./egypt --limit 2 999999 1000000
1       2
1       4
1       8
1       16
1       33
1       50
1       83
1       12450
1       32912
1       49368
1       90387
1       285716
1       571432
1       1999996
1       685610442
1       2057423870
1       11904714285
1       83332333335
1       499999000000

real    0m0.002s
user    0m0.000s
sys     0m0.002s
```

---

## Examples

### 7 / 19

* Wolfram|Alpha
  * 1 / 3 + 1 / 29 + 1 / 1653   
* `egypt --merge --limit 19 7 19`
  * 1 / 3 + 1 / 33 + 1 / 209    

### 2023 / 2024
* Wolfram|Alpha
  * 1 / 2 + 1 / 3 + 1 / 7 + 1 / 43 + 1 / 16768 + 1 / 766160103 + 1 / 978335504948790912
* `egypt --merge --limit 2023 2023 2024`
  * 1 / 2 + 1 / 3 + 1 / 8 + 1 / 33 + 1 / 92
* `egypt --reverse --merge --limit 2023 2023 2024`
  * 1 / 2 + 1 / 3 + 1 / 7 + 1 / 43 + 1 / 18447 + 1 / 184184
* `egypt --limit 2 2023 2024`
    *   1 / 2 + 1 / 4 + 1 / 8 + 1 / 11 + 1 / 33 + 1 / 674 + 1 / 899 + 1 / 2442 + 1 / 4044 + 1 / 24938 + 1 / 2046264 + 1 / 2423704

## Irrational / Transcendental Numbers

Supports RPN expressions with constants: `pi`, `e`, `phi` (golden ratio), `sqrt2`, `gamma` (Euler-Mascheroni).

```bash
# Pi/4 as Egypt fractions (raw symbolic tuples)
$ egypt pi 4 --raw -p 64
1	1	1	3
4	5	1	1
9	14	1	15
219	452	1	72
32763	33215	1	9
331698	364913	1	17
6535219	6900132	1	2
20335483	27235615	1	5
```

Each tuple `(u, v, i, j)` represents: `sum_{k=i}^{j} 1/((u-v+vk)(u+vk))`

The `-p` flag controls precision in bits (default 256). Higher precision = more CF terms = more tuples.

```bash
# Golden ratio
$ egypt phi 1 --raw -p 64 | head -5
1	0	0	0
1	1	1	1
1	2	1	1
2	3	1	1
3	5	1	1
```

Note: For irrationals, output represents a finite-precision rational approximation.
Early tuples are stable convergents of the true constant; later ones may diverge.

### Pell Equation Solver

Post-process Egypt output to find Pell solutions (p² - D·q² = ±1):

```bash
# sqrt(13): fundamental solution 649² - 13·180² = 1
$ egypt "13 sqrt" 1 --raw -p 64 | python3 scripts/pell.py 13
q	p	norm
1	3	-4
1	4	3
2	7	-3
3	11	4
5	18	-1
...
180	649	1
# Fundamental solution (norm=1): p=649, q=180

# Cattle problem (D=4729494): 41-digit solution in 22ms
$ time egypt "4729494 sqrt" 1 --raw -p 2048 | python3 scripts/pell.py 4729494 | tail -1
50549485234315033074477819735540408986340	109931986732829734979866232821433543901088049	1
real	0m0.022s
```

## Note

> * returns rather small denominators
>
When using legacy configuration `egypt --merge --limit <LIMIT> <NUMERATOR> <DENOMINATOR>`, where `LIMIT >= DENOMINATOR - 1`,
largest denominator factor should not be greater than original denominator. Fast default limit is however `2`,
which means that *bisecting* large symbolic sums can introduce bigger denominators. Moreover, `--limit` argument is itself
limited by `usize`, as opposed to other `BigInt` inputs.

## Relation to Continued Fractions

```
<< "wl/Egypt.wl"

compare[Rational[p_, q_]] := {
  Total /@ Partition[ Differences @ Convergents[ p/q ], 2],
  ReleaseHold @ EgyptianFractions[ p/q , Method -> "Expression"]
}
```

**Theorem**: Egypt values equal paired differences of continued fraction convergents.
This explains the monotonicity property: paired CF differences cancel the alternating sign pattern.

---

## Theory & Related Work

For theoretical background on the symbolic telescoping representation:

- **Paper**: [Egyptian Fractions via Modular Inverse: Symbolic Telescoping Representation](https://github.com/popojan/orbit/blob/main/docs/papers/egyptian-fractions-telescoping.tex)
- **Wolfram implementation**: [`Orbit/Kernel/EgyptianFractions.wl`](https://github.com/popojan/orbit/blob/main/Orbit/Kernel/EgyptianFractions.wl)
- **CF-Egypt Bijection (Proven Dec 2025)**: [Session documentation](https://github.com/popojan/orbit/blob/main/docs/sessions/2025-12-10-cf-egypt-equivalence/README.md)
- **γ-Egypt Simplification**: [Characterization theorems](https://github.com/popojan/orbit/blob/main/docs/sessions/2025-12-10-cf-egypt-equivalence/gamma-egypt-simplification.md)

The symbolic representation `{u, v, i, j}` compresses consecutive unit fractions into
telescoping sums with closed form: `(j-i+1) / ((u-v+vi)(u+vj))`.

This reduces complexity from O(numerator) expanded fractions to O(log denominator) symbolic tuples.

### CF-Egypt Bijection (✅ Proven)

For q = a/b with CF [0; a₁, a₂, ..., aₙ] and convergent denominators {q₀=1, q₁, ..., qₙ=b}:

| Case | u_k | v_k | j_k |
|------|-----|-----|-----|
| Regular (k < ⌈n/2⌉ or n even) | q_{2k-2} | q_{2k-1} | a_{2k} |
| Last tuple, odd CF | q_{n-1} | q_n - q_{n-1} | 1 |

**Key insight:** XGCD quotients ARE CF coefficients, enabling single-pass computation.

---

## Illustrations

* Lissajous Curves: https://mathworld.wolfram.com/LissajousCurve.html

![lissajous1](doc/7_11.png)

![lissajous2](doc/8_11.png)

![lissajous2](doc/53_57_83.png)

## Approximating rational numbers

by different rational numbers with denominator coprime to the original denominator

### 7 / 11
![approx_7_11](doc/approx_7_11.png)

### 7 / 11 errors
![approx_7_11_err](doc/approx_7_11_err.png)
