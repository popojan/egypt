use std::cmp::Ordering;
use num_prime::nt_funcs::nth_prime;
use rug::{Integer, Float, Rational};
use rug::float::{Constant, Round};
use rug::ops::{AssignRound, Pow};

/// Value of an RPN expression: exact, or an interval that rigorously encloses it
#[derive(Clone, Debug, PartialEq)]
pub enum Enclosure {
    Exact(Rational),
    /// lo < hi; the true value lies in [lo, hi]
    Interval(Rational, Rational),
    /// nothing is known at this precision (e.g. a division by an interval containing 0)
    Unknown,
}

#[derive(Clone)]
enum Value {
    Exact(Rational),
    Interval(Float, Float),
    Unknown,
}

fn dn<T>(prec: u32, v: T) -> Float where Float: AssignRound<T, Round = Round, Ordering = Ordering> {
    Float::with_val_round(prec, v, Round::Down).0
}

fn up<T>(prec: u32, v: T) -> Float where Float: AssignRound<T, Round = Round, Ordering = Ordering> {
    Float::with_val_round(prec, v, Round::Up).0
}

fn fail(msg: String) -> ! {
    eprintln!("Error: {}", msg);
    std::process::exit(1);
}

impl Value {
    fn bounds(&self, prec: u32) -> Option<(Float, Float)> {
        match self {
            Value::Exact(r) => Some((dn(prec, r), up(prec, r))),
            Value::Interval(lo, hi) => Some((lo.clone(), hi.clone())),
            Value::Unknown => None,
        }
    }

    fn to_integer(&self) -> Integer {
        match self {
            Value::Exact(r) => r.numer().clone(),
            Value::Interval(lo, hi) => match (lo.to_integer(), hi.to_integer()) {
                (Some(a), Some(b)) if a == b => a,
                _ => fail("integer argument is not determined at this precision".into()),
            },
            Value::Unknown => fail("integer argument is not determined at this precision".into()),
        }
    }

    fn to_u64(&self) -> u64 {
        self.to_integer().to_u64().unwrap_or(0)
    }
}

/// Interval hull of values computed with downward and upward rounding
fn hull(lows: Vec<Float>, highs: Vec<Float>) -> Value {
    let lo = lows.into_iter().reduce(|a, b| if b < a { b } else { a }).unwrap();
    let hi = highs.into_iter().reduce(|a, b| if b > a { b } else { a }).unwrap();
    if lo.is_nan() || hi.is_nan() || lo.is_infinite() || hi.is_infinite() {
        Value::Unknown
    } else {
        Value::Interval(lo, hi)
    }
}

fn contains_zero(lo: &Float, hi: &Float) -> bool {
    *lo <= 0 && *hi >= 0
}

fn binop(op: &str, a: &Value, b: &Value, prec: u32) -> Value {
    if let (Value::Exact(ra), Value::Exact(rb)) = (a, b) {
        // exact unless the exponent is fractional
        if op != "^" || *rb.denom() == 1 {
            return Value::Exact(match op {
                "+" => Rational::from(ra + rb),
                "-" => Rational::from(ra - rb),
                "*" => Rational::from(ra * rb),
                "/" => {
                    if rb.numer().is_zero() {
                        fail("division by zero".into());
                    }
                    Rational::from(ra / rb)
                }
                _ => match rb.numer().to_i32() {
                    Some(exp) if exp >= 0 || !ra.numer().is_zero() => Rational::from(ra.pow(exp)),
                    _ => fail(format!("exponent {} is out of range for base {}", rb, ra)),
                },
            });
        }
    }
    if let (Value::Exact(rb), "^") = (b, op) {
        if *rb.denom() == 1 {
            return int_pow(a, rb.numer(), prec);
        }
    }
    let (Some((al, ah)), Some((bl, bh))) = (a.bounds(prec), b.bounds(prec)) else {
        return Value::Unknown;
    };
    match op {
        "+" => Value::Interval(dn(prec, &al + &bl), up(prec, &ah + &bh)),
        "-" => Value::Interval(dn(prec, &al - &bh), up(prec, &ah - &bl)),
        "*" => {
            let c = [(&al, &bl), (&al, &bh), (&ah, &bl), (&ah, &bh)];
            hull(c.iter().map(|&(x, y)| dn(prec, x * y)).collect(), c.iter().map(|&(x, y)| up(prec, x * y)).collect())
        }
        "/" => {
            if contains_zero(&bl, &bh) {
                if let Value::Exact(rb) = b { if rb.numer().is_zero() { fail("division by zero".into()); } }
                return Value::Unknown;
            }
            let c = [(&al, &bl), (&al, &bh), (&ah, &bl), (&ah, &bh)];
            hull(c.iter().map(|&(x, y)| dn(prec, x / y)).collect(), c.iter().map(|&(x, y)| up(prec, x / y)).collect())
        }
        _ => {
            // x^y is monotone in each argument for x > 0, so the extremes sit at the corners
            if ah <= 0 {
                fail("fractional power of a non-positive base".into());
            }
            if al <= 0 {
                return Value::Unknown;
            }
            let c = [(&al, &bl), (&al, &bh), (&ah, &bl), (&ah, &bh)];
            hull(c.iter().map(|&(x, y)| dn(prec, x.pow(y))).collect(), c.iter().map(|&(x, y)| up(prec, x.pow(y))).collect())
        }
    }
}

/// Interval to an integer power: monotone for odd n, monotone in |x| for even n
fn int_pow(a: &Value, n: &Integer, prec: u32) -> Value {
    let Some(e) = n.to_i32() else { fail(format!("exponent {} is out of range", n)) };
    let Some((lo, hi)) = a.bounds(prec) else { return Value::Unknown };
    if e == 0 {
        return Value::Exact(Rational::from(1));
    }
    let m = e.unsigned_abs();
    // borrowed receivers give incomplete values, rounded once in the requested direction
    let (pl, ph) = if m % 2 == 1 || lo >= 0 {
        (dn(prec, (&lo).pow(m)), up(prec, (&hi).pow(m)))
    } else if hi <= 0 {
        (dn(prec, (&hi).pow(m)), up(prec, (&lo).pow(m)))
    } else {
        let top = if Float::with_val(prec, -&lo) > hi { Float::with_val(prec, -&lo) } else { hi };
        (Float::with_val(prec, 0), up(prec, (&top).pow(m)))
    };
    let p = Value::Interval(pl, ph);
    if e > 0 { p } else { binop("/", &Value::Exact(Rational::from(1)), &p, prec) }
}

/// Unified RPN calculator supporting both exact rational and irrational operations.
///
/// Supported tokens:
/// - Numbers: integers of any size, decimals and scientific notation (all exact)
/// - Constants: pi, e, phi (golden ratio), sqrt2, gamma (Euler-Mascheroni)
/// - Binary ops: +, -, *, /, ^
/// - Unary ops: sqrt, inv, !, # (primorial: product of primes <= n), p (nth prime), np (next prime), pp (prev prime), fib
/// - Aggregators: seq (range), sum, prod, lcm
///
/// Irrational values are carried as intervals with outward rounding at the given precision (in bits).
pub fn eval_rpn(s: &str, precision: u32) -> Enclosure {
    let prec = precision;
    let mut stack: Vec<Value> = Vec::new();
    let pop = |stack: &mut Vec<Value>| stack.pop().unwrap_or_else(|| fail("RPN stack underflow".into()));

    for el in s.split_whitespace() {
        let el_lower = el.to_lowercase();
        let one = Float::with_val(prec, 1);

        // Constants (irrational)
        if el_lower == "pi" {
            stack.push(Value::Interval(dn(prec, Constant::Pi), up(prec, Constant::Pi)));
        } else if el_lower == "e" {
            stack.push(Value::Interval(dn(prec, one.exp_ref()), up(prec, one.exp_ref())));
        } else if el_lower == "phi" {
            let five = Float::with_val(prec, 5);
            let (sl, sh) = (dn(prec, five.sqrt_ref()), up(prec, five.sqrt_ref()));
            stack.push(Value::Interval(dn(prec, &sl + 1u32) / 2u32, up(prec, &sh + 1u32) / 2u32));
        } else if el_lower == "sqrt2" {
            let two = Float::with_val(prec, 2);
            stack.push(Value::Interval(dn(prec, two.sqrt_ref()), up(prec, two.sqrt_ref())));
        } else if el_lower == "gamma" {
            stack.push(Value::Interval(dn(prec, Constant::Euler), up(prec, Constant::Euler)));
        }
        // Binary operators
        else if el == "^" || el == "-" || el == "+" || el == "*" || el == "/" {
            let b = pop(&mut stack);
            let a = pop(&mut stack);
            stack.push(binop(el, &a, &b, prec));
        }
        // Sequence generator: a b seq -> pushes a, a+1, ..., b
        else if el == "seq" {
            let b = pop(&mut stack).to_u64();
            let a = pop(&mut stack).to_u64();
            for i in a..=b {
                stack.push(Value::Exact(Rational::from(i)));
            }
        }
        // Aggregators: consume entire stack
        else if el == "sum" || el == "prod" {
            let (op, unit) = if el == "sum" { ("+", 0) } else { ("*", 1) };
            let mut c = Value::Exact(Rational::from(unit));
            while let Some(v) = stack.pop() {
                c = binop(op, &c, &v, prec);
            }
            stack.push(c);
        } else if el == "lcm" {
            let mut c = Integer::from(1);
            while let Some(v) = stack.pop() {
                c.lcm_mut(&v.to_integer());
            }
            stack.push(Value::Exact(Rational::from(c)));
        }
        // Factorial
        else if el == "!" {
            let n = pop(&mut stack).to_u64();
            let mut result = Integer::from(1);
            for i in 2..=n {
                result *= i;
            }
            stack.push(Value::Exact(Rational::from(result)));
        }
        // Primorial: n# = product of all primes <= n
        else if el == "#" {
            let n = pop(&mut stack).to_integer();
            let mut result = Integer::from(1);
            let mut p = Integer::from(2);
            while p <= n {
                result *= &p;
                p.next_prime_mut();
            }
            stack.push(Value::Exact(Rational::from(result)));
        }
        // Fibonacci (fast doubling)
        else if el == "fib" {
            let n = pop(&mut stack).to_u64();

            fn fib_pair(n: u64) -> (Integer, Integer) {
                if n == 0 {
                    return (Integer::from(0), Integer::from(1));
                }
                let (a, b) = fib_pair(n / 2);
                let c = &a * (Integer::from(2) * &b - &a);
                let d = a.clone() * &a + b.clone() * &b;
                if n % 2 == 0 {
                    (c, d)
                } else {
                    (d.clone(), c + d)
                }
            }

            stack.push(Value::Exact(Rational::from(fib_pair(n).0)));
        }
        // Prime functions
        else if el == "p" {
            let result = nth_prime(pop(&mut stack).to_u64());
            stack.push(Value::Exact(Rational::from(result)));
        } else if el == "np" {
            let mut n = pop(&mut stack).to_integer();
            n.next_prime_mut();
            stack.push(Value::Exact(Rational::from(n)));
        } else if el == "pp" {
            let n = pop(&mut stack).to_integer();
            // prev_prime: find largest prime < n
            let mut candidate = if n > 2 { Integer::from(&n - 1) } else { Integer::from(2) };
            while candidate > 1 {
                match candidate.is_probably_prime(25) {
                    rug::integer::IsPrime::Yes | rug::integer::IsPrime::Probably => break,
                    rug::integer::IsPrime::No => candidate -= 1,
                }
            }
            stack.push(Value::Exact(Rational::from(candidate)));
        }
        // Square root - exact for perfect squares of rationals
        else if el == "sqrt" {
            let a = pop(&mut stack);
            if let Value::Exact(r) = &a {
                if *r < 0 {
                    fail(format!("square root of negative {}", r));
                }
                let (n, nr) = r.numer().clone().sqrt_rem(Integer::new());
                let (d, dr) = r.denom().clone().sqrt_rem(Integer::new());
                if nr == 0 && dr == 0 {
                    stack.push(Value::Exact(Rational::from((n, d))));
                    continue;
                }
            }
            stack.push(match a.bounds(prec) {
                None => Value::Unknown,
                Some((lo, hi)) => {
                    if hi < 0 {
                        fail("square root of a negative value".into());
                    }
                    // an enclosure reaching below 0 still has a non-negative true value
                    let l = if lo <= 0 { Float::with_val(prec, 0) } else { dn(prec, lo.sqrt_ref()) };
                    Value::Interval(l, up(prec, hi.sqrt_ref()))
                }
            });
        }
        // Inverse
        else if el == "inv" {
            let a = pop(&mut stack);
            stack.push(binop("/", &Value::Exact(Rational::from(1)), &a, prec));
        }
        // Number literal: arbitrary-size integer, exact decimal or scientific notation
        else if let Ok(n) = el.parse::<Integer>() {
            stack.push(Value::Exact(Rational::from(n)));
        } else if let Some(r) = parse_decimal(el) {
            stack.push(Value::Exact(r));
        } else {
            eprintln!("Warning: unknown token '{}', using 0", el);
            stack.push(Value::Exact(Rational::from(0)));
        }
    }

    match stack.pop().unwrap_or(Value::Exact(Rational::from(0))) {
        Value::Exact(r) => Enclosure::Exact(r),
        Value::Unknown => Enclosure::Unknown,
        Value::Interval(lo, hi) => {
            let (lo, hi) = (lo.to_rational().unwrap(), hi.to_rational().unwrap());
            // zero width only when every rounding was exact
            if lo == hi { Enclosure::Exact(lo) } else { Enclosure::Interval(lo, hi) }
        }
    }
}

/// Point value of an RPN expression (exact, or the lower end of its enclosure)
pub fn parse_rpn(s: &str, precision: u32) -> (Integer, Integer) {
    let r = match eval_rpn(s, precision) {
        Enclosure::Exact(r) | Enclosure::Interval(r, _) => r,
        Enclosure::Unknown => Rational::from(0),
    };
    r.into_numer_denom()
}

/// Exact value of a decimal literal such as "2.50", "-.5" or "1.5e-3"; None for anything else
fn parse_decimal(s: &str) -> Option<Rational> {
    if let Some((mantissa, exp)) = s.split_once(['e', 'E']) {
        let m = mantissa.parse::<Integer>().ok().map(Rational::from).or_else(|| parse_plain_decimal(mantissa))?;
        let exp: i32 = exp.parse().ok()?;
        let scale = Rational::from(Integer::from(10).pow(exp.unsigned_abs()));
        return Some(if exp >= 0 { m * scale } else { m / scale });
    }
    parse_plain_decimal(s)
}

fn parse_plain_decimal(s: &str) -> Option<Rational> {
    let (negative, body) = match s.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, s.strip_prefix('+').unwrap_or(s)),
    };
    let (int_part, frac_part) = body.split_once('.')?;
    let digits = format!("{}{}", int_part, frac_part);
    if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
        return None;
    }
    let mut num = digits.parse::<Integer>().ok()?;
    if negative {
        num = -num;
    }
    let den = Integer::from(10).pow(frac_part.len() as u32);
    Some(Rational::from((num, den)))
}

/// Integer-valued RPN (used for Pell's D); a fractional result is truncated
pub fn _parse_rpn(s: &str) -> Integer {
    let (num, den) = parse_rpn(s, 256);
    if den == 1 {
        num
    } else {
        num / den
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn exact(expr: &str, num: &str, den: &str) {
        let (n, d) = parse_rpn(expr, 64);
        assert_eq!((n.to_string(), d.to_string()), (num.to_string(), den.to_string()), "{}", expr);
    }

    fn close_to_sqrt2(expr: &str) {
        let (n, d) = parse_rpn(expr, 64);
        let x = Rational::from((n, d));
        let err = (Rational::from(&x * &x) - 2u32).abs();
        assert!(err < Rational::from((1u32, Integer::from(1) << 60u32)), "{}", expr);
    }

    #[test]
    fn integer_literals_of_any_size() {
        exact("7", "7", "1");
        exact("-3", "-3", "1");
        exact("162259276829213363391578010288127", "162259276829213363391578010288127", "1");
        exact("162259276829213363391578010288127 170141183460469231731687303715884105727 /",
              "162259276829213363391578010288127", "170141183460469231731687303715884105727");
    }

    #[test]
    fn decimals_are_exact() {
        exact("0.1", "1", "10");
        exact("2.50", "5", "2");
        exact("-0.25", "-1", "4");
        exact(".5 2 *", "1", "1");
    }

    #[test]
    fn arithmetic() {
        exact("2 3 +", "5", "1");
        exact("2 3 -", "-1", "1");
        exact("2 3 *", "6", "1");
        exact("1 3 /", "1", "3");
        exact("2 10 ^", "1024", "1");
        exact("2 -2 ^", "1", "4");
        exact("2 107 ^ 1 -", "162259276829213363391578010288127", "1");
        exact("4 inv", "1", "4");
    }

    #[test]
    fn aggregators_and_sequences() {
        exact("1 10 seq sum", "55", "1");
        exact("1 5 seq prod", "120", "1");
        exact("1 10 seq lcm", "2520", "1");
    }

    #[test]
    fn number_theory_tokens() {
        exact("20 !", "2432902008176640000", "1");
        exact("10 #", "210", "1");
        exact("5 p", "11", "1");
        exact("10 np", "11", "1");
        exact("10 pp", "7", "1");
        exact("100 fib", "354224848179261915075", "1");
        exact("49 sqrt", "7", "1");
    }

    #[test]
    fn irrationals_are_approximated_to_precision() {
        close_to_sqrt2("2 sqrt");
        close_to_sqrt2("sqrt2");
        close_to_sqrt2("2 0.5 ^");
        let (n, d) = parse_rpn("pi", 64);
        let pi = Rational::from((n, d));
        assert!(pi > Rational::from((314159u32, 100000u32)) && pi < Rational::from((314160u32, 100000u32)));
    }

    #[test]
    fn scientific_notation_is_exact() {
        exact("1e-3", "1", "1000");
        exact("1.5e3", "1500", "1");
        exact("-2.5E-1", "-1", "4");
    }

    #[test]
    fn intervals_enclose_the_value() {
        let Enclosure::Interval(lo, hi) = eval_rpn("pi", 64) else { panic!() };
        let pi = Rational::from((Integer::from_str_radix("3243f6a8885a308d313198a2e037", 16).unwrap(), Integer::from(1) << 108u32));
        assert!(lo < pi && pi < hi && Rational::from(&hi - &lo) < Rational::from((1u32, Integer::from(1) << 60u32)));
        for expr in ["2 sqrt 2 ^", "2 sqrt 2 sqrt *", "2 0.5 ^ 2 0.5 ^ *"] {
            let Enclosure::Interval(lo, hi) = eval_rpn(expr, 64) else { panic!("{}", expr) };
            assert!(lo < 2 && hi > 2, "{}", expr);
        }
        let Enclosure::Interval(lo, hi) = eval_rpn("pi pi -", 64) else { panic!() };
        assert!(lo < 0 && hi > 0);
        assert_eq!(eval_rpn("1 pi pi - /", 64), Enclosure::Unknown);
        assert_eq!(eval_rpn("pi 0 *", 64), Enclosure::Exact(Rational::new()));
    }

    #[test]
    fn legacy_integer_entry_point() {
        assert_eq!(_parse_rpn("4729494"), 4729494);
        assert_eq!(_parse_rpn("2 61 ^ 1 -"), Integer::from(2305843009213693951u64));
    }
}
