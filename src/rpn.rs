use std::convert::TryFrom;
use num_prime::nt_funcs::nth_prime;
use rug::{Integer, Float, Rational};
use rug::float::Constant;
use rug::ops::Pow;

/// Unified RPN calculator supporting both exact rational and irrational operations.
/// Returns (numerator, denominator) pair.
///
/// Supported tokens:
/// - Numbers: integers of any size and decimals (both exact), scientific notation (approximate)
/// - Constants: pi, e, phi (golden ratio), sqrt2, gamma (Euler-Mascheroni)
/// - Binary ops: +, -, *, /, ^
/// - Unary ops: sqrt, inv, !, # (primorial: product of primes <= n), p (nth prime), np (next prime), pp (prev prime), fib
/// - Aggregators: seq (range), sum, prod, lcm
///
/// When irrational operations are used (sqrt on non-perfect-square, pi, e, etc.),
/// the result is approximated using the given precision (in bits).
pub fn parse_rpn(s: &str, precision: u32) -> (Integer, Integer) {
    let parts: Vec<&str> = s.split_whitespace().collect();

    // Use Rational for exact arithmetic, convert to Float only when needed
    enum Value {
        Exact(Rational),
        Approx(Float),
    }

    impl Value {
        fn to_float(&self, precision: u32) -> Float {
            match self {
                Value::Exact(r) => Float::with_val(precision, r),
                Value::Approx(f) => f.clone(),
            }
        }

        fn to_rational(&self) -> Rational {
            match self {
                Value::Exact(r) => r.clone(),
                Value::Approx(f) => Rational::try_from(f).unwrap_or(Rational::from(0)),
            }
        }

        fn to_integer(&self) -> Integer {
            match self {
                Value::Exact(r) => r.numer().clone(),
                Value::Approx(f) => f.to_integer().unwrap_or(Integer::from(0)),
            }
        }

        fn to_u64(&self) -> u64 {
            self.to_integer().to_u64().unwrap_or(0)
        }
    }

    let mut stack: Vec<Value> = Vec::new();

    for el in parts.iter() {
        let el_lower = el.to_lowercase();

        // Constants (irrational)
        if el_lower == "pi" {
            stack.push(Value::Approx(Float::with_val(precision, Constant::Pi)));
        } else if el_lower == "e" {
            stack.push(Value::Approx(Float::with_val(precision, 1).exp()));
        } else if el_lower == "phi" {
            let sqrt5 = Float::with_val(precision, 5).sqrt();
            stack.push(Value::Approx((Float::with_val(precision, 1) + sqrt5) / 2));
        } else if el_lower == "sqrt2" {
            stack.push(Value::Approx(Float::with_val(precision, 2).sqrt()));
        } else if el_lower == "gamma" {
            stack.push(Value::Approx(Float::with_val(precision, Constant::Euler)));
        }
        // Binary operators
        else if el == &"^" || el == &"-" || el == &"+" || el == &"*" || el == &"/" {
            let b = stack.pop().unwrap();
            let a = stack.pop().unwrap();

            let result = match (&a, &b) {
                // exact unless the exponent is fractional
                (Value::Exact(ra), Value::Exact(rb)) if *el != "^" || *rb.denom() == 1 => {
                    let c = if *el == "+" {
                        Rational::from(ra + rb)
                    } else if *el == "-" {
                        Rational::from(ra - rb)
                    } else if *el == "*" {
                        Rational::from(ra * rb)
                    } else if *el == "/" {
                        if rb.numer().is_zero() {
                            eprintln!("Error: division by zero");
                            std::process::exit(1);
                        }
                        Rational::from(ra / rb)
                    } else {
                        match rb.numer().to_i32() {
                            Some(exp) if exp >= 0 || !ra.numer().is_zero() => Rational::from(ra.pow(exp)),
                            _ => {
                                eprintln!("Error: exponent {} is out of range for base {}", rb, ra);
                                std::process::exit(1);
                            }
                        }
                    };
                    Value::Exact(c)
                }
                _ => {
                    let fa = a.to_float(precision);
                    let fb = b.to_float(precision);
                    let c = if *el == "+" {
                        fa + fb
                    } else if *el == "-" {
                        fa - fb
                    } else if *el == "*" {
                        fa * fb
                    } else if *el == "/" {
                        fa / fb
                    } else {
                        fa.pow(&fb)
                    };
                    Value::Approx(c)
                }
            };
            stack.push(result);
        }
        // Sequence generator: a b seq -> pushes a, a+1, ..., b
        else if *el == "seq" {
            let b = stack.pop().unwrap().to_u64();
            let a = stack.pop().unwrap().to_u64();
            for i in a..=b {
                stack.push(Value::Exact(Rational::from(i)));
            }
        }
        // Aggregators: consume entire stack
        else if *el == "sum" {
            let mut c = Rational::from(0);
            while !stack.is_empty() {
                c += stack.pop().unwrap().to_rational();
            }
            stack.push(Value::Exact(c));
        } else if *el == "prod" {
            let mut c = Rational::from(1);
            while !stack.is_empty() {
                c *= stack.pop().unwrap().to_rational();
            }
            stack.push(Value::Exact(c));
        } else if *el == "lcm" {
            let mut c = Integer::from(1);
            while !stack.is_empty() {
                let val = stack.pop().unwrap().to_integer();
                c.lcm_mut(&val);
            }
            stack.push(Value::Exact(Rational::from(c)));
        }
        // Factorial
        else if *el == "!" {
            let a = stack.pop().unwrap();
            let n = a.to_u64();
            let mut result = Integer::from(1);
            for i in 2..=n {
                result *= i;
            }
            stack.push(Value::Exact(Rational::from(result)));
        }
        // Primorial: n# = product of all primes <= n
        else if *el == "#" {
            let a = stack.pop().unwrap();
            let n = a.to_integer();
            let mut result = Integer::from(1);
            let mut p = Integer::from(2);
            while p <= n {
                result *= &p;
                p.next_prime_mut();
            }
            stack.push(Value::Exact(Rational::from(result)));
        }
        // Fibonacci (fast doubling)
        else if *el == "fib" {
            let a = stack.pop().unwrap();
            let n = a.to_u64();

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

            let result = fib_pair(n).0;
            stack.push(Value::Exact(Rational::from(result)));
        }
        // Prime functions
        else if *el == "p" {
            let a = stack.pop().unwrap();
            let result = nth_prime(a.to_u64());
            stack.push(Value::Exact(Rational::from(result)));
        } else if *el == "np" {
            let a = stack.pop().unwrap();
            let mut n = a.to_integer();
            n.next_prime_mut();
            stack.push(Value::Exact(Rational::from(n)));
        } else if *el == "pp" {
            let a = stack.pop().unwrap();
            let n = a.to_integer();
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
        // Square root - check if exact
        else if *el == "sqrt" {
            let a = stack.pop().unwrap();
            match &a {
                Value::Exact(r) if r.is_integer() => {
                    let n = r.numer();
                    let (root, rem) = n.clone().sqrt_rem(Integer::new());
                    if rem == 0 {
                        // Perfect square
                        stack.push(Value::Exact(Rational::from(root)));
                    } else {
                        // Not perfect square, use float
                        let f = Float::with_val(precision, n).sqrt();
                        stack.push(Value::Approx(f));
                    }
                }
                _ => {
                    let f = a.to_float(precision).sqrt();
                    stack.push(Value::Approx(f));
                }
            }
        }
        // Inverse
        else if *el == "inv" {
            let a = stack.pop().unwrap();
            match a {
                Value::Exact(r) => {
                    stack.push(Value::Exact(Rational::from(r.recip())));
                }
                Value::Approx(f) => {
                    stack.push(Value::Approx(Float::with_val(precision, 1) / f));
                }
            }
        }
        // Number literal: arbitrary-size integer, exact decimal, then f64 as last resort
        else if let Ok(n) = el.parse::<Integer>() {
            stack.push(Value::Exact(Rational::from(n)));
        } else if let Some(r) = parse_decimal(el) {
            stack.push(Value::Exact(r));
        } else if let Ok(f) = el.parse::<f64>() {
            stack.push(Value::Approx(Float::with_val(precision, f)));
        } else {
            eprintln!("Warning: unknown token '{}', using 0", el);
            stack.push(Value::Exact(Rational::from(0)));
        }
    }

    let result = stack.pop().unwrap_or(Value::Exact(Rational::from(0)));
    let rational = result.to_rational();
    let (num, den) = rational.into_numer_denom();

    (num, den)
}

/// Exact value of a plain decimal literal such as "2.50" or "-.5"; None for anything else
fn parse_decimal(s: &str) -> Option<Rational> {
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
    fn legacy_integer_entry_point() {
        assert_eq!(_parse_rpn("4729494"), 4729494);
        assert_eq!(_parse_rpn("2 61 ^ 1 -"), Integer::from(2305843009213693951u64));
    }
}
