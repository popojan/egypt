mod rpn;

use crate::rpn::{_parse_rpn, parse_rpn};

use std::io::{self, BufWriter, StdoutLock, Write};
use std::ops::{Add, Sub, Div, Mul, SubAssign, Neg};
use clap::Parser;
use rug::{Assign, Complete, Integer, Rational};

/// Egyptian Fractions

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
#[command(arg_required_else_help(true))]
struct Args {
    /// Reverse merge strategy
    #[clap(short, long, value_parser, default_value_t = false)]
    reverse: bool,

    /// Extra O(n^2) merge step possibly reducing number of terms
    #[clap(short, long, value_parser, default_value_t = false)]
    merge: bool,

    /// Output minimal number of raw quadruplets (aka symbolic sums)
    #[clap(long, value_parser, default_value_t = false)]
    raw: bool,

    /// Output raw quadruplets bisected according to --limit
    #[clap(long, value_parser, default_value_t = false)]
    bisect: bool,

    /// Split long symbolic sums greedily (largest unit fraction first) instead of symmetrically
    #[clap(short, long, value_parser, default_value_t = false)]
    greedy: bool,

    /// No output
    #[clap(short, long, value_parser, default_value_t = false)]
    silent: bool,

    /// Print integers in hexadecimal (decimal conversion dominates large outputs)
    #[clap(short = 'x', long, value_parser, default_value_t = false)]
    hex: bool,

    /// Batch mode (expects numerator and denominator on each line of stdin)
    #[clap(long, value_parser, default_value_t = false)]
    batch: bool,

    #[clap(value_parser, default_value_t = String::from("1"))]
    numerator: String,

    #[clap(value_parser, default_value_t = String::from("1"))]
    denominator: String,

    /// Maximum number of terms for breaking large symbolic sums
    #[clap(short, long, value_parser, default_value_t = 8)]
    limit: usize,

    /// Precision in bits for irrational/transcendental constants (pi, e, phi, etc.)
    /// Higher precision = more CF terms = more Egypt tuples
    #[clap(short, long, value_parser, default_value_t = 256)]
    precision: u32,

    /// Pell equation mode: output (q, p, norm) for sqrt(D)/1 input
    /// Finds solutions to p² - D·q² = ±1
    #[clap(long, value_parser, default_value_t = false)]
    pell: bool,

    /// Divisor-residue descent: fewer terms, much larger denominators (experimental, doc/divisor-descent.md)
    #[clap(long, value_parser, default_value_t = false)]
    descent: bool,
}

fn merge(eg: &Vec<(Integer, Integer, Integer, Integer)>) -> Vec<(Integer, Integer, Integer, Integer)> {
    let mut i = 0_usize;
    let mut ret = vec![];
    while i < eg.len() {
        let mut q = Rational::from((eg[i].0.clone(), eg[i].1.clone()));
        let (mut ones_i, mut ones_q) = (i, q.clone());
        for j in (i + 1)..eg.len() {
            let p = &eg[j];
            q += Rational::from((p.0.clone(), p.1.clone()));
            if q.numer() == &Integer::from(1) {
                (ones_i, ones_q) = (j, q.clone());
            }
        }
        let (j, q) = (ones_i, ones_q);
        //println!("{:?} {:?} {:?}", ones, x, y);
        let (x, y) = q.into_numer_denom();
        ret.push((x, y, Integer::from(0), Integer::from(0)));
        i = j + 1;
    }
    ret
}

/// Compute CF quotients via single GCD pass: O(log p) instead of O(log² p)
fn cf_quotients(mut a: Integer, mut b: Integer) -> Vec<Integer> {
    let mut quotients = vec![];
    while !b.is_zero() {
        let (q, r) = a.div_rem(b.clone());
        quotients.push(q);
        a = b;
        b = r;
    }
    quotients
}

/// Build convergent denominators from CF quotients
fn convergent_denominators(cf: &[Integer]) -> Vec<Integer> {
    if cf.is_empty() {
        return vec![Integer::from(1)];
    }
    let mut qs = vec![Integer::from(1)];
    if cf.len() > 1 {
        qs.push(cf[1].clone());
    }
    for i in 2..cf.len() {
        let next = cf[i].clone() * &qs[qs.len() - 1] + &qs[qs.len() - 2];
        qs.push(next);
    }
    qs
}

/// ModInv-based computation, kept as the reference form of the paper; the CF pass below is 3-8x faster
#[allow(dead_code)]
fn as_egyptian_fraction_symbolic_modinv(x0: &Integer, y0: &Integer, _expand: bool, ret: &mut Vec<(Integer, Integer, Integer, Integer)>) {
    let gcd = x0.clone().gcd(&y0);
    let mut x = x0.clone().div(&gcd);
    let mut y = y0.clone().div(&gcd);
    if x.ge(&y) {
        ret.push((x.clone().div(&y), 0.into(), 0.into(), 0.into()));
        x.sub_assign(x.clone().div(&y).mul(&y));
    }
    while x.gt(&Integer::from(0)) && y.gt(&Integer::from(1)) {
        let v = x.clone().neg().invert(&y).unwrap();
        let t;
        (t, x) = x.clone().div_rem((x.clone() * &v + 1) / &y);
        y -= v.clone() * &t;
        ret.push((y.clone(), v, 1.into(), t));
    }
    if !x.is_zero() {
        ret.push((y, Integer::from(1), Integer::from(0), Integer::from(0)));
    }
}

/// Single Euclid pass plus the CF-Egypt bijection: one division per quotient instead of one
/// modular inverse per tuple. Tuples come out in CF order (smallest denominators first).
fn as_egyptian_fraction_symbolic_cf(x0: &Integer, y0: &Integer, _expand: bool, ret: &mut Vec<(Integer, Integer, Integer, Integer)>) {
    let gcd = x0.clone().gcd(&y0);
    let mut x = x0.clone().div(&gcd);
    let y = y0.clone().div(&gcd);

    // Handle integer part
    if x.ge(&y) {
        ret.push((x.clone().div(&y), 0.into(), 0.into(), 0.into()));
        x.sub_assign(x.clone().div(&y).mul(&y));
    }

    if x.is_zero() {
        return;
    }

    // Get CF quotients in single pass
    let cf = cf_quotients(x.clone(), y.clone());
    // n = number of partial quotients (CF length excluding leading 0)
    let n = cf.len() - 1;

    if n == 0 {
        // CF = [0], meaning x = 0 (should have returned earlier)
        return;
    }

    if n == 1 {
        // CF = [0; a_1], meaning x/y = 1/y → unit fraction
        // Raw tuple (1, y-1, 1, 1): T = 1/((1)(1+y-1)) = 1/y
        ret.push((Integer::from(1), y.clone() - 1, Integer::from(1), Integer::from(1)));
        return;
    }

    // Build convergent denominators
    let qs = convergent_denominators(&cf);
    let num_tuples = (n + 1) / 2;

    // Apply CF-Egypt bijection formula (0-based indexing)
    // Formula: u_k = q_{2k-2}, v_k = q_{2k-1}, j_k = a_{2k}
    // In 0-based: u = qs[2k-2], v = qs[2k-1], j = cf[2k] (since cf[i] = a_i)
    for k in 1..=num_tuples {
        if k < num_tuples || n % 2 == 0 {
            // Regular case: (q_{2k-2}, q_{2k-1}, 1, a_{2k})
            let u = qs[2 * k - 2].clone();
            let v = qs[2 * k - 1].clone();
            let j = cf[2 * k].clone();  // a_{2k} = cf[2k] in 0-based
            ret.push((u, v, Integer::from(1), j));
        } else {
            // Last tuple, odd CF: (q_{n-1}, q_n - q_{n-1}, 1, 1)
            let u = qs[n - 1].clone();
            let v = qs[n].clone() - &qs[n - 1];
            ret.push((u, v, Integer::from(1), Integer::from(1)));
        }
    }
}

/// Rational entry point: CF pass, then the fractional tuples reversed to keep the
/// largest-first order of the modinv version (raw output and split order depend on it)
fn as_egyptian_fraction_symbolic(x0: &Integer, y0: &Integer, expand: bool, ret: &mut Vec<(Integer, Integer, Integer, Integer)>) {
    let start = ret.len();
    as_egyptian_fraction_symbolic_cf(x0, y0, expand, ret);
    let is_integer_part = ret.len() > start && ret[start].1.is_zero() && ret[start].3.is_zero();
    ret[start + usize::from(is_integer_part)..].reverse();
}

/// CF version - use when stability analysis needed (irrationals)
fn as_egyptian_fraction_symbolic_for_irrational(x0: &Integer, y0: &Integer, expand: bool, ret: &mut Vec<(Integer, Integer, Integer, Integer)>) {
    as_egyptian_fraction_symbolic_cf(x0, y0, expand, ret)
}

fn expand(eg: &Vec<(Integer, Integer, Integer, Integer)>) -> Vec<(Integer, Integer, Integer, Integer)> {
    let mut ret = vec![];
    for (b,v,i,j) in eg.iter() {
        if v.is_zero() && i.is_zero() && j.is_zero() {
            ret.push((b.clone(), Integer::from(1), Integer::from(0), Integer::from(0)));
        } else {
            for k in i.to_usize().unwrap()..j.to_usize().unwrap() + 1 {
                ret.push((
                    Integer::from(1),
                    b.clone().sub(v).add(v.clone().mul(&k))
                        .mul(b.add(v.clone().mul(k.clone()))),
                    Integer::from(0),
                    Integer::from(0)
                ));
            }
        }
    }
    ret
}

fn as_egyptian_fraction(a:&Integer, b:&Integer, args: &Args)->Vec<(Integer, Integer,Integer,Integer)> {
    if args.descent {
        return as_egyptian_fraction_descent(a, b, args);
    }
    let mut res = vec![];
    as_egyptian_fraction_symbolic(
        &a,
        &b,
        args.reverse, &mut res);
    let limit = args.limit.max(2);
    if !args.raw {
        res = halve_symbolic_sums(&res, limit, args.greedy);
        res = expand(&res);
        res.sort_by(|x, y| { x.1.cmp(&y.1)});
        if args.merge {
            if !args.reverse {
                res.reverse();
            }
            res = merge(&res);
        }
        res = fix_duplicates(&res);
    } else {
        if args.bisect {
            res = halve_symbolic_sums(&res, limit, args.greedy);
        }
    }
    res
}

/// Divisors t = d·e of M = p_k#·c with d | p_k# and e in {1, c}, ascending
fn descent_divisors(k: usize, c: &Integer) -> (Integer, Vec<Integer>) {
    let mut gens: Vec<Integer> = (2u32..).filter(|n| (2..*n).take_while(|d| d * d <= *n).all(|d| n % d != 0))
        .take(k).map(Integer::from).collect();
    if *c > 1 {
        gens.push(c.clone());
    }
    let mut divisors = vec![Integer::from(1)];
    for g in &gens {
        let more: Vec<Integer> = divisors.iter().map(|d| Integer::from(d * g)).collect();
        divisors.extend(more);
    }
    divisors.sort();
    (divisors.last().unwrap().clone(), divisors)
}

/// Remainder denominators double in bit length per term; past this the input is too large for the descent
const DESCENT_MAX_BITS: u64 = 1 << 24;

/// Divisor-residue descent after "Short Egyptian fractions" (OpenAI, 2026,
/// https://github.com/openai/math/blob/main/preprints/Short-Egyptian-fractions-September-25-2026/Short-Egyptian-fractions-September-25-2026.pdf):
/// for t | M, u/(MQ) = 1/((M/t)z) + h/(zMQ), z = ceil(Qt/u), h = (-Qt) mod u; each step takes the t
/// with the smallest h (ties: smaller t, i.e. smaller remainder denominator) whose term is unused.
/// M = lcm(p_K#, b) gives Q = 1 without factoring b; K = bits(a)/3 within [4, 16] (doc/divisor-descent.md)
fn as_egyptian_fraction_descent(a: &Integer, b: &Integer, args: &Args) -> Vec<(Integer, Integer, Integer, Integer)> {
    let g = a.clone().gcd(b);
    let (mut a, b) = (a.clone() / &g, b.clone() / &g);
    let mut res = vec![];
    if a >= b {
        res.push((a.clone() / &b, Integer::from(1), Integer::from(0), Integer::from(0)));
        a %= &b;
    }
    let k = (a.significant_bits() as usize).div_ceil(3).clamp(4, 16);
    let (pk, _) = descent_divisors(k, &Integer::from(1));
    let c = b.clone() / b.clone().gcd(&pk);
    let (m, divisors) = descent_divisors(k, &c);
    let (mut u, mut q) = (a * (m.clone() / &b), Integer::from(1));
    let mut seen = std::collections::HashSet::new();
    let mut tmp = Integer::new();
    while !u.is_zero() {
        let qm = q.clone() % &u;
        let h_of = |t: &Integer, r: &mut Integer| {
            r.assign(&qm * t);
            *r %= &u;
            if !r.is_zero() {
                *r = u.clone() - &*r;
            }
        };
        // allocation-free scan for the smallest h: r = Qt mod u is 0 there, or else the largest r
        let (mut best_r, mut order) = (Integer::from(-1), vec![]);
        for (i, t) in divisors.iter().enumerate() {
            tmp.assign(&qm * t);
            tmp %= &u;
            let better = if best_r.is_zero() { false } else { tmp.is_zero() || tmp > best_r };
            if better {
                best_r.assign(&tmp);
                order.clear();
            }
            if better || tmp == best_r {
                order.push(i);
            }
        }
        let try_order = |order: &[usize]| order.iter().find_map(|&i| {
            let t = &divisors[i];
            let mut h = Integer::new();
            h_of(t, &mut h);
            let z = (q.clone() * t + &h) / &u;
            let d = m.clone() / t * &z;
            (!seen.contains(&d)).then(|| (d, z, h))
        });
        // collisions are rare: sort all residues only if every smallest-h term is already used
        let step = try_order(&order).or_else(|| {
            let hs: Vec<Integer> = divisors.iter().map(|t| { let mut h = Integer::new(); h_of(t, &mut h); h }).collect();
            let mut rest: Vec<usize> = (0..divisors.len()).collect();
            rest.sort_by(|&i, &j| hs[i].cmp(&hs[j]).then(i.cmp(&j)));
            try_order(&rest)
        });
        let Some((d, z, h)) = step else {
            eprintln!("Error: --descent found no unused divisor term");
            std::process::exit(2);
        };
        // bit-length sum bounds the product within 2 bits; each factor stays below 2^32 bits
        let bits = [&z, &m, &q].iter().map(|x| u64::from(x.significant_bits())).sum::<u64>();
        if !h.is_zero() && bits > DESCENT_MAX_BITS {
            eprintln!("Error: --descent denominators exceed 2^24 bits, input too large");
            std::process::exit(2);
        }
        res.push((Integer::from(1), d.clone(), Integer::from(0), Integer::from(0)));
        seen.insert(d);
        (u, q) = (h, q * z);
    }
    res.sort_by(|x, y| x.1.cmp(&y.1));
    if args.raw {
        // a unit fraction 1/d is the one-term tuple (1, d-1, 1, 1), an integer n is (n, 0, 0, 0)
        res = res.into_iter().map(|(n, d, _, _)| if d == 1 {
            (n, Integer::new(), Integer::new(), Integer::new())
        } else {
            (n, d - 1u32, Integer::from(1), Integer::from(1))
        }).collect();
    }
    res
}

/// Egyptian fraction for irrationals: uses CF, sorts by fraction size
/// No truncation - user controls precision, we output full CF expansion
fn as_egyptian_fraction_irrational(a: &Integer, b: &Integer, args: &Args) -> Vec<(Integer, Integer, Integer, Integer)> {
    if args.descent {
        return as_egyptian_fraction_descent(a, b, args);
    }
    let mut res = vec![];
    as_egyptian_fraction_symbolic_for_irrational(&a, &b, args.reverse, &mut res);

    // Sort by unit fraction size (large → small)
    sort_by_fraction_size(&mut res);

    let limit = args.limit.max(2);
    if !args.raw {
        res = halve_symbolic_sums(&res, limit, args.greedy);
        res = expand(&res);
        res.sort_by(|x, y| { x.1.cmp(&y.1)});  // expanded: sort by denominator
        if args.merge {
            if !args.reverse {
                res.reverse();
            }
            res = merge(&res);
        }
        res = fix_duplicates(&res);
    } else {
        if args.bisect {
            res = halve_symbolic_sums(&res, limit, args.greedy);
        }
    }
    res
}

/// Replace every run of n equal terms by the expansion of n times that term until no run is left
fn fix_duplicates(eg: &Vec<(Integer, Integer, Integer, Integer)>)
    -> Vec<(Integer, Integer, Integer, Integer)> {
    let mut eg = eg.clone();
    loop {
        eg.sort_by(|x, y| { y.1.cmp(&x.1)});
        let mut run = None;
        let mut i = 0;
        while i < eg.len() {
            let mut j = i + 1;
            while j < eg.len() && eg[j] == eg[i] {
                j += 1;
            }
            if j - i > 1 {
                run = Some((i, j));
                break;
            }
            i = j;
        }
        let Some((i, j)) = run else { break };
        let mut new = vec![];
        as_egyptian_fraction_symbolic(&Integer::from(j - i), &eg[i].1, false, &mut new);
        eg.splice(i..j, expand(&new));
    }
    eg.sort_by(|x, y| { x.1.cmp(&y.1)});
    eg
}

fn calculate_raw_sum(u:&Integer, v:&Integer, i:&Integer, j:&Integer) -> (Integer, Integer) {
    let num = Integer::from(1).sub(i).add(j);
    let den = u.clone().sub(v).add(v.clone().mul(i)).mul(u.add(v.mul(j)).complete());
    let gcd = num.clone().gcd(&den);
    (num.div(&gcd), den.div(&gcd))
}

fn halve_symbolic_sums(a: &Vec<(Integer, Integer, Integer, Integer)>, limit: usize, greedy: bool)
    -> Vec<(Integer, Integer, Integer, Integer)>
{
    let mut stack = a.clone();
    let mut ret = vec![] ;
    let limit = Integer::from(limit);
    let two = Integer::from(2);
    while !stack.is_empty() {
        let (u, v, i, j) = stack.pop().unwrap();
        let term_count = j.clone().sub(&i).add(&Integer::from(1));
        if term_count.le(&limit) {
            ret.push((u, v, i, j));
        } else {
            let (a, b) = calculate_raw_sum(&u, &v, &i, &j);
            if greedy {
                // largest unit fraction 1/d <= a/b; remainder numerator is < a, so this terminates
                let d = (b.clone() + &a - 1u32) / &a;
                ret.push((Integer::from(1), d.clone() - 1u32, Integer::from(1), Integer::from(1)));
                let num = a * &d - &b;
                if !num.is_zero() {
                    as_egyptian_fraction_symbolic(&num, &(b * &d), false, &mut stack);
                }
            } else if a.is_odd() {
                let a1 = a.sub(&Integer::from(1)).div(&two);
                let a2 = a1.clone().add(&Integer::from(1));
                as_egyptian_fraction_symbolic(&a1, &b, false, &mut stack);
                as_egyptian_fraction_symbolic(&a2, &b, false, &mut stack);
            } else {
                let a1 = a.div(&two).sub(&Integer::from(1));
                let a2 = a1.clone().add(&two);
                as_egyptian_fraction_symbolic(&a1, &b, false, &mut stack);
                as_egyptian_fraction_symbolic(&a2, &b, false, &mut stack);
            }
        }
    }
    ret
}
/// Pell equation solver: extract (q, p, norm) from Egypt tuples
/// Tuples encode q_{2k-2}, q_{2k-1} as u, v parameters
/// Returns: Vec<(q, p, norm)> where norm = p² - D·q²
fn solve_pell(tuples: &[(Integer, Integer, Integer, Integer)], d: &Integer) -> Vec<(Integer, Integer, Integer)> {
    let mut results = vec![];

    // Extract integer part (a0) and q values from tuples
    let mut a0 = Integer::from(0);
    let mut qs: Vec<Integer> = vec![];

    for (u, v, i, j) in tuples {
        if v.is_zero() && i.is_zero() && j.is_zero() {
            // Integer part tuple
            a0 = u.clone();
        } else {
            // Fractional tuple: u = q_{2k-2}, v = q_{2k-1}
            if qs.is_empty() || qs.last() != Some(u) {
                qs.push(u.clone());
            }
            qs.push(v.clone());
        }
    }

    if a0.is_zero() {
        return results;
    }

    // Compute p values using CF recurrence:
    // p_{-1} = 1, p_0 = a0
    // p_n = a_n * p_{n-1} + p_{n-2}
    // where a_n = (q_n - q_{n-2}) / q_{n-1}
    let mut ps = vec![Integer::from(1), a0.clone()];

    // First: q_0 = 1, p_0 = a0
    let norm0 = a0.clone().square() - d;
    results.push((Integer::from(1), a0.clone(), norm0));

    for n in 1..qs.len() {
        let a_n = if n == 1 {
            qs[1].clone()
        } else {
            (qs[n].clone() - &qs[n - 2]) / &qs[n - 1]
        };

        let p_n = a_n * &ps[ps.len() - 1] + &ps[ps.len() - 2];
        ps.push(p_n.clone());

        let q = qs[n].clone();
        let norm = p_n.clone().square() - d * q.clone().square();
        results.push((q, p_n, norm));
    }

    results
}

/// Check if RPN expression contains irrational/transcendental constants
fn contains_irrational(s: &str) -> bool {
    let lower = s.to_lowercase();
    let trimmed = lower.trim();
    lower.contains("pi") ||
    lower.contains(" e ") || lower.starts_with("e ") || lower.ends_with(" e") || trimmed == "e" ||
    lower.contains("phi") || lower.contains("sqrt2") || lower.contains("gamma") ||
    lower.contains("sqrt")  // sqrt of non-perfect square
}

/// Evaluate both RPN expressions as rationals and divide: (a/b) / (c/d) = ad / bc
/// Returns (numerator, denominator, is_irrational)
fn parse_rpn_auto(num_str: &str, den_str: &str, precision: u32) -> (Integer, Integer, bool) {
    let is_irrational = contains_irrational(num_str) || contains_irrational(den_str);
    let (num_a, den_a) = parse_rpn(num_str, precision);
    let (num_b, den_b) = parse_rpn(den_str, precision);
    let final_num = (num_a * &den_b).abs();
    let final_den = (den_a * &num_b).abs();
    if final_den.is_zero() {
        eprintln!("Error: denominator evaluates to zero");
        std::process::exit(2);
    }
    let gcd = final_num.clone().gcd(&final_den);
    (final_num / &gcd, final_den / &gcd, is_irrational)
}

/// Calculate the smallest denominator a raw tuple produces (for sorting)
/// For (u, v, i, j): smallest denom = (u-v+v*i)*(u+v*i) when i=1
fn tuple_min_denominator(u: &Integer, v: &Integer) -> Integer {
    // First unit fraction: 1/(u*(u+v))
    u.clone() * (u.clone() + v)
}

/// Sort tuples by unit fraction size (large → small = small denominator first)
fn sort_by_fraction_size(tuples: &mut Vec<(Integer, Integer, Integer, Integer)>) {
    tuples.sort_by(|a, b| {
        let denom_a = tuple_min_denominator(&a.0, &a.1);
        let denom_b = tuple_min_denominator(&b.0, &b.1);
        denom_a.cmp(&denom_b)  // smaller denominator = larger fraction = comes first
    });
}

/// Decimal conversion is superlinear and dominates large outputs; hex is linear
fn fmt_int(i: &Integer, hex: bool) -> String {
    if hex { format!("{:x}", i) } else { i.to_string() }
}

fn stdout_writer() -> BufWriter<StdoutLock<'static>> {
    BufWriter::with_capacity(1 << 20, io::stdout().lock())
}

fn main() {
    let args = Args::parse();

    if args.batch {
        let mut out = stdout_writer();
        for line in io::stdin().lines() {
            if let Ok(line) = line {
                let num_den = line.split("\t").take(2).collect::<Vec<&str>>();
                if num_den.len() < 2 {
                    writeln!(out, "expecting tab delimited numerator and denominator").unwrap();
                    continue;
                }

                let (num, den, is_irrational) = parse_rpn_auto(num_den[0], num_den[1], args.precision);
                let fractions = if is_irrational {
                    as_egyptian_fraction_irrational(&num, &den, &args)
                } else {
                    as_egyptian_fraction(&num, &den, &args)
                };
                if !args.silent {
                    let hex = args.hex;
                    let mut gt0 = false;
                    write!(out, "{}\t{}\t", fmt_int(&num, hex), fmt_int(&den, hex)).unwrap();
                    for (i, (a, b, c, d)) in fractions.iter().enumerate() {
                        let is_natural = (args.raw && b.is_zero() && c.is_zero() && d.is_zero())
                            || (!args.raw && *b == Integer::from(1));
                        let term = if args.raw {
                            format!("{},{},{},{}", fmt_int(a, hex), fmt_int(b, hex), fmt_int(c, hex), fmt_int(d, hex))
                        } else {
                            fmt_int(b, hex)
                        };
                        if i == 0 && is_natural {
                            write!(out, "{}\t", fmt_int(a, hex)).unwrap();
                            gt0 = true;
                        } else if i == 0 {
                            write!(out, "0\t{}", term).unwrap();
                        } else if i == 1 && gt0 {
                            write!(out, "{}", term).unwrap();
                        } else {
                            write!(out, " {}", term).unwrap();
                        }
                    }
                    writeln!(out).unwrap();
                }
                // one line per input, so a coprocess feeding stdin sees each answer at once
                out.flush().unwrap();
            }
        }
    } else {
        let (num, den, is_irrational) = parse_rpn_auto(&args.numerator, &args.denominator, args.precision);

        // Pell mode: extract D from "D sqrt" pattern
        if args.pell {
            // Parse D from numerator (expecting "D sqrt" RPN)
            let d = extract_pell_d(&args.numerator);
            if d.is_none() {
                eprintln!("Error: --pell requires input like 'egypt \"D sqrt\" 1 --pell'");
                std::process::exit(2);
            }
            let d = d.unwrap();

            // Get raw tuples (always use raw mode internally for Pell)
            let mut res = vec![];
            as_egyptian_fraction_symbolic_for_irrational(&num, &den, args.reverse, &mut res);
            sort_by_fraction_size(&mut res);

            // Solve Pell equation
            let pell_results = solve_pell(&res, &d);

            if !args.silent {
                println!("q\tp\tnorm");
                let mut found_fundamental = false;
                let mut found_quasi = false;

                for (q, p, norm) in pell_results {
                    println!("{}\t{}\t{}", fmt_int(&q, args.hex), fmt_int(&p, args.hex), fmt_int(&norm, args.hex));

                    if !found_quasi && norm == Integer::from(-1) {
                        eprintln!("# Quasi-solution (norm=-1): p={}, q={}", p, q);
                        found_quasi = true;
                    }
                    if !found_fundamental && norm == Integer::from(1) {
                        eprintln!("# Fundamental solution (norm=1): p={}, q={}", p, q);
                        found_fundamental = true;
                        break;
                    }
                }

                if !found_fundamental {
                    if found_quasi {
                        eprintln!("# No fundamental solution, but quasi-solution exists");
                    } else {
                        eprintln!("# No Pell solution found (increase -p precision)");
                    }
                    std::process::exit(1);
                }
            }
        } else {
            let fractions = if is_irrational {
                as_egyptian_fraction_irrational(&num, &den, &args)
            } else {
                as_egyptian_fraction(&num, &den, &args)
            };
            if !args.silent {
                let (hex, mut out) = (args.hex, stdout_writer());
                for (a, b, c, d) in fractions.iter() {
                    if args.raw {
                        writeln!(out, "{}\t{}\t{}\t{}", fmt_int(a, hex), fmt_int(b, hex), fmt_int(c, hex), fmt_int(d, hex)).unwrap();
                    } else {
                        writeln!(out, "{}\t{}", fmt_int(a, hex), fmt_int(b, hex)).unwrap();
                    }
                }
                out.flush().unwrap();
            }
        }
    }
}

/// Extract D from "D sqrt" RPN pattern for Pell mode
fn extract_pell_d(s: &str) -> Option<Integer> {
    let parts: Vec<&str> = s.trim().split_whitespace().collect();
    if parts.len() >= 2 && parts.last() == Some(&"sqrt") {
        // Try to parse the number before "sqrt"
        let num_str = parts[..parts.len()-1].join(" ");
        Some(_parse_rpn(&num_str))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::ops::Pow;
    use std::collections::HashSet;

    /// Expansion must sum exactly to the input and use each unit fraction once
    fn check(cli: &[&str]) {
        let args = Args::parse_from(std::iter::once("egypt").chain(cli.iter().copied()));
        let (num, den, irrational) = parse_rpn_auto(&args.numerator, &args.denominator, args.precision);
        let fractions = if irrational {
            as_egyptian_fraction_irrational(&num, &den, &args)
        } else {
            as_egyptian_fraction(&num, &den, &args)
        };
        let mut sum = Rational::new();
        let mut seen = HashSet::new();
        for (a, b, _, _) in &fractions {
            sum += Rational::from((a.clone(), b.clone()));
            if *a == 1 {
                assert!(seen.insert(b.clone()), "duplicate 1/{} in {:?}", b, cli);
            }
        }
        assert_eq!(sum, Rational::from((num, den)), "{:?}", cli);
    }

    #[test]
    fn expansions_sum_to_input_with_distinct_terms() {
        let inputs = [
            ("7", "19"), ("2023", "2024"), ("999999", "1000000"), ("58", "3511471"),
            ("5", "3"), ("2 107 ^ 1 -", "2 127 ^ 1 -"), ("1 3 /", "0.7"), ("49 sqrt", "11"),
            // at -l 2 the duplicate 1/37 expands to 1/19 + 1/703 and 1/19 is itself a repeat
            ("24285564808518", "223702468834263"),
        ];
        let modes: [&[&str]; 5] = [&[], &["-g"], &["-m"], &["-m", "-r"], &["--descent"]];
        for (a, b) in inputs {
            for limit in ["2", "8"] {
                for mode in modes {
                    let mut cli = vec!["-l", limit, a, b];
                    cli.extend_from_slice(mode);
                    check(&cli);
                }
            }
        }
    }

    #[test]
    fn big_literals_parse_exactly() {
        let (n, d, irrational) = parse_rpn_auto(
            "162259276829213363391578010288127", "170141183460469231731687303715884105727", 64);
        assert_eq!(n, Integer::from(2u32).pow(107u32) - 1u32);
        assert_eq!(d, Integer::from(2u32).pow(127u32) - 1u32);
        assert!(!irrational);
    }
}