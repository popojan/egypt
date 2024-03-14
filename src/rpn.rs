use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};
use std::str::FromStr;
use num_prime::nt_funcs::nth_prime;
use rug::Integer;
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

        } else if el == "!"  || el == "p" || el == "sqrt" {
            let mut a = Integer::from_str(&stack.pop().unwrap()).unwrap();
            let c = if el == "!" {
                let mut c = Integer::from(1);
                while a > 1 {
                    c.mul_assign(&a);
                    a.sub_assign(Integer::from(1));
                }
                c.to_string()
            } else if el == "p" {
                nth_prime(a.to_u64().unwrap()).to_string()
            } else if el == "sqrt" {
                a.sqrt().to_string()
            } else {
                "0".to_string()
            };
            stack.push(c);
        } else {
            stack.push(el.to_string());
        }
    };
    Integer::from_str(&stack.pop().unwrap()).unwrap()
}