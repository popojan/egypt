# Divisor-residue descent (`--descent`)

Experimental strategy after the descent step of
[*Short Egyptian fractions*](https://github.com/openai/math/blob/main/preprints/Short-Egyptian-fractions-September-25-2026/Short-Egyptian-fractions-September-25-2026.pdf)
(OpenAI, 2026-09-25; [Lean formalization](https://github.com/openai/math/blob/main/lean/docs/025.md),
not peer reviewed). The paper proves O(log log b) terms for every a/b (Erdős #304). The step used
here is the one in its proof outline and in Section 5 (`descent.tex`, eq. `rescaled-descent`).

For a remainder u/(MQ) and a divisor t of a fixed auxiliary M:

    u/(MQ) = 1/((M/t)·z) + h/(z·M·Q),   z = ⌈Q·t/u⌉,   h = u·z − Q·t = (−Q·t) mod u

The input enters as u = a·M/gcd(b, M), Q = b/gcd(b, M). Each step emits one unit fraction, and the
remainder is h/(MQ') with Q' = zQ and h < u.

The work went in two rounds. The first used M = p_K# with K fixed, as the task suggested, and
the hypothesis below was written before it. The second changed M and K on the evidence of the
first; that is what `--descent` does now.

## Hypothesis (written before implementation)

Model: with M = p_k# there are 2^k divisors. If the residues (−Q·t) mod u behave like random
points of [0, u), the smallest is about u/2^k, so each step removes about k bits from the numerator.
The numerator starts with log2(a) + log2(M) bits, so

    terms ≈ (log2 a + log2 M)/k + 1–2 (tail: h = 0 needs u/gcd(u,Q) to be squarefree and p_k-smooth)

log2(M)/k is 2.9 for k = 8, 3.6 for k = 12, 4.05 for k = 16 and 4.4 for k = 20. This is the fixed
cost of multiplying the numerator by M. So:

1. **Term count, fixed k.** The descent is a greedy step that is about 0.7·k times faster.
   Greedy removes about 1 nat of numerator per step, the descent about k bits. The count stays
   linear in log a, not log log b. With k = 12 I expect about 4–5 terms at b ≈ 10³, 6 at 10⁶, 8 at
   10¹² and 12 at 10³⁰. The CF-based strategies grow at least linearly in the CF length
   (≈ 0.84·ln b tuples), so the descent should win clearly at 10¹² and 10³⁰.
2. **Small inputs.** No gain. For 7/19 and 2023/2024 I expect 4–6 terms. `--merge` already gives
   3 and 5 terms, and the log2(M) overhead dominates there. Smaller k helps here.
3. **Scaling k with b.** With M ≈ b, as in the paper (k ≈ ln b/ln ln b), the formula gives
   terms ≈ 2.9·ln ln b: the paper's O(log log b) in heuristic form. The cost is 2^k candidates per
   step, which is feasible up to b ≈ 10³⁰ (k ≈ 20) and hopeless far beyond it.
4. **Denominators.** Once the chosen t exceeds u/Q, z > 1 and the remainder denominator squares,
   so its bit length roughly doubles every step, as in classical greedy:
   B(n+1) ≈ 2·B(n) + k − ½·log2 M. For k = 12 that is about 4·2^n bits after n terms: roughly
   100–250 bits at b ≈ 10³–10⁶, about 1000 at 10¹² and 10⁴ or more at 10³⁰. Default, `--limit` and
   `--merge` keep the maximum near b^2 to b^8, i.e. tens to hundreds of bits. `--greedy` is
   between the two.
5. **Larger k helps both measures** (fewer terms, so fewer doublings) until the 2^k scan per step
   becomes too slow. The trade-off is runtime, not denominator size.
6. a = b−1 should behave like a random a of the same size. Nothing in the descent uses the
   CF structure that makes b−1 special for the CF-based methods.

Expected verdict: a large term-count win only for large b, bought with doubly exponential
denominators. This matches the paper, which bounds length only.

## Round 1: M = p_K#, fixed K

The first version took `--descent-primes K` (default 12). It handed the remainder to the default
expansion once a denominator passed 65536 bits. Cells are **terms / max denominator bits**. For
random groups (20 samples, b uniform in [10^e, 2·10^e), seed 304) a cell is median (max) / median.

| K | random 10³ | random 10⁶ | random 10¹² | random 10³⁰ |
|---|---|---|---|---|
| 4 | 4.5 / 21 | 6 / 89.5 | 10 / 1215.5 | 38.5 / 137872 (cap) |
| 8 | 4.5 / 21.5 | 6 / 63.5 | 9 / 460.5 | 17 / 46601.5 (cap 8 of 21) |
| 12 | 5 / 22.5 | 7 / 96.5 | 9 / 371.5 | 14 / 8592 |
| 16 | 5 / 22.5 | 8 / 115.5 | 9 / 345.5 | 14 / 4640 |
| 20 | 5 / 22.5 | 9 / 182.5 | 10 / 409 | 14 / 3907 |

Per-step trace (random a/b ≈ 10³⁰, K = 16; numerator bits u → h, then bits of the emitted
denominator):

    u159→158 →156 →154 →151 →141 →127 →114 →100 →86 →67 →52 →39 →23 →10 →0
    d:  2, 2, 4, 5, 9, 19, 32, 63, 107, 199, 383, 747, 1482, 2956, 5884

An independent Python re-implementation reproduces the binary's term counts. No duplicate-term
collision occurred in any run.

### Hypothesis check

| claim | verdict |
|---|---|
| 1. term count, fixed K: linear in log a, about k bits per step, clear win at 10¹² and 10³⁰ | **Confirmed, but 1–2 terms worse than predicted.** K = 12: 5 / 7 / 9 / 14 terms against the predicted 4–5 / 6 / 8 / 12. At 10³⁰ that is 14 terms where `--merge` needs 56.5; at 10¹² it is 9 against 20.5. |
| 2. no gain on small inputs | **Confirmed.** 7/19 gives the same 3 terms (for K ≥ 12, literally the greedy 3, 29, 1653). 2023/2024 gives 6 terms against 5 for `--merge --limit 2023`. |
| 3. scaling K with b gives about 2.9·ln ln b | **Falsified as stated.** It plateaus: 14 terms at both K = 16 and K = 20 for 10³⁰, because the model missed a start-up phase (below). |
| 4. denominators roughly double in bits per step | **Confirmed** (trace above). The start is slower than predicted: 96 rather than about 250 bits at 10⁶, 370 rather than about 1000 at 10¹². |
| 5. larger K helps both measures | **Partly.** True at 10³⁰ up to K = 16. Flat from 16 to 20, and worse for small b. |
| 6. a = b−1 behaves like random a | **Roughly.** It is a little above the random median (K = 12: 11 against 9 at 10¹², 17 against 14 at 10³⁰), but far below every CF method, for which b−1 is the worst case (`default`: 139 terms at 10³⁰). |

**What the model missed.** While the remainder r = u/(MQ) is large, h/u = ⌈x⌉ − x with
x = t/(rM) ∈ (0, 1/r]. The divisors then cover only about 1/r unit intervals, so the residues are
far from uniform, and the first steps remove 1–5 bits instead of K. Full K-bit steps come once
r < 2^−K. Because r squares every step, that takes about log₂ K near-greedy steps:

    terms ≈ log2 K + (log2 a + log2(M/gcd(b, M)))/K + 1

This gives 15 for K = 16 at 10³⁰ (measured 14) and explains the plateau. Asymptotically it is
still O(log log b) when M ≈ b, but only far beyond b ≈ 10³⁰, where 2^K divisors per step become
infeasible (K ≈ 43 for 10¹⁰⁰).

## Round 2: M = lcm(p_K#, b), K from the input, no fallback

Three findings from round 1 changed the design:

- **The gcd(b, M) term.** If M contains the prime factors of b, the log₂ M overhead and the
  denominator factor Q both shrink. This was tested in Python with the same min-h rule, for three
  choices of M: plain M = p_K#; M = p_K#·rad(b), which needs the factorization of b; and
  M = lcm(p_K#, b) with only the divisors d and d·c, where d | p_K# and c = M/p_K#. The last one
  costs one gcd and gives Q = 1 for every b.

  | b, K | p_K# | p_K#·rad(b) | lcm(p_K#, b) |
  |---|---|---|---|
  | 10⁶, K = 4 | 6 / 121 | 6 / 47.5 | 6 / 63.5 |
  | 10¹², K = 12 | 9 / 423.5 | 9 / 174.5 | 9 / 254.5 |
  | 10³⁰, K = 8 | 15.5 / 36324.5 | 13 / 2498.5 | 15 / 12823.5 |
  | 10³⁰, K = 12 | 14 / 8969.5 | 12 / 1350 | 13 / 4767 |
  | 10³⁰, K = 16 | 13 / 4527.5 | 12 / 1237 | 13 / 2175 |

  rad(b) is best but needs factoring, which a general-purpose CLI cannot afford. lcm keeps the
  numerator of the plain version, u = a·p_K#/gcd(b, p_K#), and has 2^(K+1) instead of 2^(K+ω(b))
  divisors. It takes about half of rad's denominator gain at no cost.
- **K.** The optimum of the corrected estimate is K = 4 / 7 / 13 for 10 / 20 / 40-bit numerators,
  close to what was measured. `--descent` uses K = ⌈bits(a)/3⌉, clamped to [4, 16]. The upper bound
  keeps the scan at 2^16 residues per term; K = 20 costs 2–3 s at 10³⁰ and gains nothing.
- **The fallback hurt.** Handing a remainder with a 65k-bit denominator to the default expansion
  gives terms of about b⁸ that size. At K = 16, one random input each, with the fallback / without
  it: 10⁶⁰ gave 57 / 19 terms (872K / 164K of output), and 10⁸⁰ gave 92 / 21 terms (3.5M / 564K).
  The arithmetic is unbounded (rug/GMP); the practical limits are memory and decimal output. Without
  any cap, 10¹⁰⁰ took 27 terms and 17–27 MB in 1–1.5 s with `-x` (7 s in decimal), 10¹²⁰ took
  120–312 MB and 1.2 GB RSS, and 10¹⁴⁰ passed 2³² bits after 95 s and 4.8 GB. Every further term
  doubles the size. So `--descent` stays a pure descent and exits with status 2 once a remainder
  denominator would exceed 2²⁴ bits.

## Implementation

`--descent` (`as_egyptian_fraction_descent` in `src/main.rs`) is opt-in and changes no other output.
All 24 earlier goldens are byte-identical.

- K = ⌈bits(a)/3⌉ ∈ [4, 16], M = lcm(p_K#, b) = p_K#·c. The divisors are d and d·c for d | p_K#.
  The start is u = a·M/b, Q = 1.
- Each step scans r = (Q·t) mod u over all divisors without allocating. The smallest h is r = 0,
  or else the largest r. Ties go to the smaller t, which means a smaller z and a smaller remainder
  denominator. A term 1/((M/t)z) that was already emitted is skipped; all residues are sorted only
  if every smallest-h term is taken.
- Exit status 2 if no unused term exists (never observed) or a remainder denominator would exceed
  2²⁴ bits. The bit test sums the factors' bit lengths, because rug's `significant_bits()` returns
  `u32`.
- `--limit`, `--greedy`, `--merge` and `--reverse` do not apply. With `--raw`, 1/d is printed as the
  one-term tuple (1, d−1, 1, 1), since (u−v+v·1)(u+v·1) = 1·d, and an integer part n as (n, 0, 0, 0).
- Tests: the unit invariant test, the CLI switch matrix, a quarter of the random-rational test,
  M107/M127, a b with repeated small primes, a primorial denominator, the raw mapping, and the
  golden `descent_999999`.

## Results

`scripts/bench_descent.py target/release/egypt 20` checks every row for an exact sum and distinct
unit terms. The cells use the same format as in round 1. README commands: 7/19 `-m -l 19` gives
3 / 8, 2023/2024 `-m -l 2023` gives **5 / 7**, and 2023/2024 `-r -m -l 2023` gives 6 / 18.

| strategy | 7/19 | 2023/2024 | (10^3-1)/10^3 | random, b~10^3 | (10^6-1)/10^6 | random, b~10^6 | (10^12-1)/10^12 | random, b~10^12 | (10^30-1)/10^30 | random, b~10^30 |
|---|---|---|---|---|---|---|---|---|---|---|
| `default` | 3 / 8 | 18 / 21 | 11 / 19 | 8.5 (14) / 18 | 30 / 39 | 16.5 (32) / 40 | 60 / 107 | 35 (48) / 106 | 139 / 349 | 85.5 (106) / 322 |
| `--merge` | 3 / 8 | 17 / 21 | 7 / 9 | 5 (10) / 18 | 22 / 39 | 8.5 (21) / 39.5 | 56 / 107 | 21 (33) / 101.5 | 130 / 349 | 56.5 (77) / 318.5 |
| `--reverse --merge` | 3 / 8 | 17 / 21 | 7 / 9 | 5.5 (10) / 18 | 22 / 39 | 9 (21) / 39.5 | 56 / 107 | 20.5 (32) / 101.5 | 130 / 349 | 56.5 (77) / 318.5 |
| `--greedy` | 3 / 8 | 11 / 23 | 9 / 19 | 8 (14) / 18 | 19 / 39 | 15.5 (27) / 40 | 30 / 100 | 32.5 (45) / 123.5 | 89 / 921 | 72.5 (95) / 640.5 |
| `--limit 2` | 3 / 8 | 12 / 22 | 9 / 19 | 7 (10) / 18.5 | 19 / 39 | 12 (24) / 59 | 52 / 107 | 24.5 (36) / 125 | 117 / 383 | 66.5 (78) / 385.5 |
| `--descent` | 3 / 10 | 6 / 27 | 6 / 15 | 4.5 (6) / 19.5 | 8 / 43 | 6 (8) / 58.5 | 11 / 242 | 9 (10) / 234 | 14 / 1598 | 13 (16) / 2304 |

Larger inputs (one random a/b each): 10⁴⁰ gives 16 terms and 8059 bits in 0.14 s, 10⁶⁰ gives 17
terms and 51887 bits, and 10¹⁰⁰ gives 25 terms and 15.1M bits in 0.6 s. At 10¹⁵⁰ the run exits
with status 2.

Structured inputs:

- (2³¹−1)/(2⁶¹−1). The CF is [0; 2³⁰, 2, 2³⁰−1], and the default gives 3 terms,
  1/(n+1) + 1/((n+1)(2n+1)) + 1/((2n+1)(2⁶¹−1)) with n = 2³⁰ (92 bits). That is optimal. 2⁶¹−1 is
  prime, so a/b = 1/x + 1/y would need a divisor d ∈ {1, b, b²} of b² with a | b + d (from
  (ax − b)(ay − b) = b²), and none exists. `--descent` gives 7 terms with up to 456 bits.
- M107/M127. The default gives 169 terms (540 bits), `-m` 152, `-l 2` 130, `-g` 68; `--descent`
  gives 12 terms with up to 15089 bits. In round 1 there were 13 terms at K = 12; the match with
  169 = 13² was a coincidence, since other K and other `--limit` values break it. Those 13 terms also
  could not be a regrouping of the 169: a group summing to 1/d has d below each member's
  denominator, but 7 of the 13 denominators exceeded the default maximum.

## Primorial denominators

The alternating factorial sums S_k = ½·Σ_{j≤k} (−1)^j j!/(2j+1) of the `orbit` primorial paper
have denominator m# with m = 2k+1, the case the divisor method suits best. For a/m# = frac(S_k):

| m | default | `-m` | round 1, K = 12 | `--descent` | divisor greedy (d \| m#) |
|---|---|---|---|---|---|
| 37 | 33 / 85 | 17 / 85 | 7 / 71 | 8 / 92 | 7 / 40 |
| 53 | 69 / 200 | 44 / 200 | 7 / 174 | 6 / 120 | 7 / 62 |
| 61 | 62 / 153 | 29 / 135 | 8 / 322 | 7 / 97 | 7 / 71 |
| 71 | 49 / 316 | 37 / 316 | 9 / 688 | 7 / 181 | 8 / 87 |

The divisor greedy (Erdős, Tenenbaum–Yokota) takes the largest divisor of m# not above the
remaining numerator, so every denominator divides m#. Here it is nearly as short as the descent:
with Q = 1, a z = 1 step is always available.

The predictable structure is in signed partial fractions, not in a positive expansion. This
identity was verified for every prime m in [11, 400]:

    S_k ≡ 1/6 + Σ_{5≤p≤m} c_p/p  (mod 1),   c_p = centered residue of (−1)^((p−1)/2)·((p−1)/2)!/2 mod p

For p ≡ 3 (mod 4), c_p = ε_p·(p−1)/2, so the contribution is ε_p·(1/2 − 1/(2p)): a signed unit
fraction 1/(2p) whose sign is the class-number bit of Q(√−p) (the paper's Remark on class numbers;
checked by hand for p ≤ 71). For p ≡ 1 (mod 4), c_p is half a square root of −1 mod p. Positive
expansions do not inherit this structure:

- Expansions at consecutive prime entries share 0 or 1 terms (usually 1/2). Each new prime moves
  frac(S) by c_p/p, which is of order one, so size-driven methods start over.
- An expansion prime by prime with denominators dividing m# is impossible:
  Σ_{d | m#/p} 1/(p·d) = (1/p)·Π_{q≠p}(1 + 1/q) ≈ 1.08·ln m/p, which reaches only c_p ≲ ln m.

## Prime factors of the denominators

These denominators are not prime-factor noise; each method leaves a signature. Measured on random
a/b with b ≈ 10¹² (20 samples), against random integers of the same bit lengths:

| | default | `-m` | `-l 2` | `--descent` | random control |
|---|---|---|---|---|---|
| term divisible by a convergent denominator q_i ≥ 1000 | 39 % | 43 % | 51 % | 3 % | 0 % |
| shares a factor > 1000 with another term | 78 % | 76 % | 80 % | 68 % | 1 % |
| bits of the 997-smooth part | 15.3 | 14.9 | 14.8 | 28.0 | 8.3 |
| divisible by b | 3 % | 5 % | 4 % | 20 % | 0 % |

- **CF methods.** Terms are (u−v+vk)(u+vk) with (u, v) = (q_{2i−2}, q_{2i−1}). Neighbouring terms
  share the factor u+vk (the sum telescopes), and the first and last terms of a tuple contain a
  convergent denominator, since u + v·a_{2i} = q_{2i}. The doubled smooth part only reflects that
  each term is a product of two numbers.
- **Descent.** Terms are (M/t)·z with a large p_K-smooth part. M contains b, and the last term MQ/u
  carries every earlier z, so it shares factors with all of them. The z = ⌈Qt/u⌉ themselves look
  like noise.

## Verdict

- **Length:** a large win from b ≈ 10⁶ upward: 6 vs 8.5 terms at 10⁶, 9 vs 20.5 at 10¹² and 13 vs
  56.5 at 10³⁰ (medians against the best existing strategy). No gain at 10³ or on the README
  examples, and none on inputs with a short CF such as (2³¹−1)/(2⁶¹−1).
- **Denominators:** doubly exponential, about 2^n bits after n terms. Median max bits are about
  1.5× `--merge` at 10⁶, 2.3× at 10¹² and 7× at 10³⁰. `--descent` exits once they would pass 2²⁴
  bits, which happens between 10¹⁰⁰ and 10¹⁵⁰.
- **Cost:** at most 2^16 residues per term, well under a second up to 10¹⁰⁰.

Worth keeping as an opt-in mode for "few terms, size no object". Open ideas: preferring z = 1
while it is available (the Tenenbaum–Yokota regime, which keeps denominators ≤ M·b; with
M = p_K#·rad(b) and K = 12 it gave 11 terms within 76 bits at 10¹², but at K = 8 and 10³⁰ it needed
hundreds of terms); a score trading log h against
log z; and cheap trial division of b to add its small primes to M without full factoring.
