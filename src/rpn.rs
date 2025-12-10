use std::convert::TryFrom;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};
use std::str::FromStr;
use num_prime::nt_funcs::nth_prime;
use rug::{Integer, Float, Rational};
use rug::float::Constant;
use rug::ops::Pow;

pub fn _parse_rpn(s: &str) -> Integer { // TODO error handling
    let parts = s.split(" ")
        .into_iter().map(|x| x.to_string()).collect::<Vec::<String>>();
    let mut stack = Vec::<String>::new();
    for el in parts.iter() {
        if el == "^" || el == "-" || el == "+" || el == "*" {
            let b = Integer::from_str(&stack.pop().unwrap()).unwrap();
            let a = Integer::from_str(&stack.pop().unwrap()).unwrap();
            let c = if el == "^" {
                a.pow(&b.to_u32().unwrap())
            } else if el == "+" {
                a.add(&b)
            } else if el == "-" {
                a.sub(&b)
            } else if el == "*" {
                a.mul(&b)
            } else {
                Integer::from(0)
            };
            stack.push(c.to_string());
        } else if el == "seq" {
            let b = Integer::from_str(&stack.pop().unwrap()).unwrap().to_u64().unwrap();
            let mut a = Integer::from_str(&stack.pop().unwrap()).unwrap().to_u64().unwrap();
            while a <= b {
                stack.push(a.clone().to_string());
                a.add_assign(1);
            }
        } else if el == "sum" || el == "prod" || el == "lcm" {
            let mut c = Integer::from(0);
            if el == "sum" {
                while !stack.is_empty() {
                    let a = Integer::from_str(&stack.pop().unwrap()).unwrap();
                    c.add_assign(&a);
                }
            } else if el == "prod" {
                c.add_assign(1);
                while !stack.is_empty() {
                    let a = Integer::from_str(&stack.pop().unwrap()).unwrap();
                    c.mul_assign(&a);
                }
            } else if el == "lcm" {
                c.add_assign(1);
                while !stack.is_empty() {
                    let a = Integer::from_str(&stack.pop().unwrap()).unwrap();
                    c.lcm_mut(&a);
                }
            }
            stack.push(c.to_string());

        } else if el == "!"  || el == "p" || el == "sqrt" || el == "np" || el == "pp" || el == "fib" {
            let mut a = Integer::from_str(&stack.pop().unwrap()).unwrap();
            let c = if el == "fib" {
                // Fibonacci using fast doubling
                fn fib(n: u64) -> (Integer, Integer) {
                    if n == 0 {
                        return (Integer::from(0), Integer::from(1));
                    }
                    let (a, b) = fib(n / 2);
                    let c = &a * (Integer::from(2) * &b - &a);
                    let d = a.clone() * &a + b.clone() * &b;
                    if n % 2 == 0 {
                        (c, d)
                    } else {
                        (d.clone(), c + d)
                    }
                }
                fib(a.to_u64().unwrap()).0.to_string()
            } else if el == "!" {
                let mut c = Integer::from(1);
                while a > 1 {
                    c.mul_assign(&a);
                    a.sub_assign(Integer::from(1));
                }
                c.to_string()
            } else if el == "p" {
                nth_prime(a.to_u64().unwrap()).to_string()
            } else if el == "np" {
                a.next_prime().to_string()
            } else if el == "pp" {
                a.prev_prime().to_string()
            } else if el == "sqrt" {
                a.sqrt().to_string()
            } else {
                "0".to_string()
            };
            stack.push(c);
        }  else {
            stack.push(el.to_string());
        }
    };
    Integer::from_str(&stack.pop().unwrap()).unwrap()
}

/// Parse RPN expression with irrational constants, returning (numerator, denominator)
/// Supports: pi, e, phi (golden ratio), sqrt2, gamma (Euler-Mascheroni)
/// precision: number of bits for Float computation
pub fn _parse_rpn_irrational(s: &str, precision: u32) -> (Integer, Integer) {
    let parts = s.split(" ")
        .into_iter().map(|x| x.to_string()).collect::<Vec<String>>();
    let mut stack = Vec::<Float>::new();

    for el in parts.iter() {
        let el_lower = el.to_lowercase();

        // Constants
        if el_lower == "pi" {
            stack.push(Float::with_val(precision, Constant::Pi));
        } else if el_lower == "e" {
            stack.push(Float::with_val(precision, 1).exp());
        } else if el_lower == "phi" {
            // Golden ratio: (1 + sqrt(5)) / 2
            let sqrt5 = Float::with_val(precision, 5).sqrt();
            stack.push((Float::with_val(precision, 1) + sqrt5) / 2);
        } else if el_lower == "sqrt2" {
            stack.push(Float::with_val(precision, 2).sqrt());
        } else if el_lower == "gamma" {
            // Euler-Mascheroni constant
            stack.push(Float::with_val(precision, Constant::Euler));
        } else if el == "^" || el == "-" || el == "+" || el == "*" || el == "/" {
            let b = stack.pop().unwrap();
            let a = stack.pop().unwrap();
            let c = if el == "^" {
                a.pow(b.to_u32_saturating().unwrap_or(1))
            } else if el == "+" {
                a + b
            } else if el == "-" {
                a - b
            } else if el == "*" {
                a * b
            } else if el == "/" {
                a / b
            } else {
                Float::with_val(precision, 0)
            };
            stack.push(c);
        } else if el == "sqrt" {
            let a = stack.pop().unwrap();
            stack.push(a.sqrt());
        } else if el == "inv" {
            let a = stack.pop().unwrap();
            stack.push(Float::with_val(precision, 1) / a);
        } else {
            // Try parsing as number
            if let Ok(n) = el.parse::<i64>() {
                stack.push(Float::with_val(precision, n));
            } else if let Ok(f) = el.parse::<f64>() {
                stack.push(Float::with_val(precision, f));
            } else {
                eprintln!("Warning: unknown token '{}', using 0", el);
                stack.push(Float::with_val(precision, 0));
            }
        }
    }

    let result = stack.pop().unwrap_or(Float::with_val(precision, 0));

    // Convert to Rational preserving full precision (not via f64!)
    let rational = Rational::try_from(&result).unwrap_or(Rational::from(0));
    let (num, den) = rational.into_numer_denom();

    (num, den)
}