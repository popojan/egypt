use num_bigint::BigInt;
use num_traits::identities::One;
use std::ops::{Add, Sub, Div, Mul, AddAssign};
use std::env;
use std::str::FromStr;
use num_integer::Integer;

fn low(x: &BigInt, y: &BigInt) -> (BigInt, BigInt) {
    let b = if y.sub(x).eq(&BigInt::one()) {
        x.clone()
    } else {
        x.extended_gcd(&y).x.add(y).mod_floor(&y)
    };
    let m = x.clone().mul(&b);
    let mut a = m.clone().div(y);
    if m.sub(a.clone().mul(y)).gt(&y.clone().div(&BigInt::from(2_u32))) {
        a.add_assign(&BigInt::one());
    }
    let gcd = a.gcd(&b);
    return (a.div(&gcd), b.div(&gcd))
}

fn as_egyptian_fraction_raw(x0: &BigInt, y0: &BigInt) -> Vec<(BigInt, BigInt)> {
    let mut whole = Vec::<(BigInt, BigInt)>::new();
    let mut ret = Vec::<(BigInt, BigInt)>::new();
    let mut x = x0.clone();
    let mut y = y0.clone();
     if x.ge(&y) {
        whole.push((x.clone().div(&y), BigInt::one()));
        x = x.clone().sub(x.div(&y).mul(y.clone()));
    }
    while x.gt(&BigInt::one()) {
        let (x2, y2) = low(&x, &y);
        let (a, b) = (x.clone().mul(&y2).sub(x2.clone().mul(&y)), y.clone().mul(&y2));
        let gcd = a.gcd(&b);
        ret.push((a.div(&gcd), b.div(&gcd)));
        (x, y) = (x2, y2);
    }
    ret.extend(vec![(x.clone(), y.clone())]);
    ret.extend(whole);
    ret.reverse();
    ret
}

fn merge(eg: &Vec<(BigInt, BigInt)>) -> Vec<(BigInt, BigInt)> {
    let mut i = 0_usize;
    let mut ret = Vec::<(BigInt, BigInt)>::new();
    while i < eg.len() {
        let (mut x, mut y) = eg[i].clone();
        let mut ones = vec![(i, x.clone(), y.clone())];
        for j in (i + 1)..eg.len() {
            let (x2, y2) = &eg[j];
            (x, y) = (x.mul(y2).add(x2.mul(y.clone())), y.clone().mul(y2));
            let gcd = x.gcd(&y);
            (x, y) = (x.div(&gcd), y.div(&gcd));
            if x.eq(&BigInt::one()) {
                ones.push((j, x.clone(), y.clone()));
            }
        }
        let (j, x, y) = ones.last().unwrap();
        ret.push((x.clone(), y.clone()));
        i = j + 1;
    }
    return ret
}

fn as_egyptian_fraction(x0: &BigInt, y0: &BigInt) -> Vec<(BigInt, BigInt)> {
    merge(&as_egyptian_fraction_raw(x0, y0))
}

/// Converts rational numbers to egyptian fractions
/// featuring quite small denominators compared to competing algorithms
/// Largest factor of any denominator is not greater than the original denominator

fn main() {
    let mut num = BigInt::from(7_u32);
    let mut den = BigInt::from(11_u32);

    if env::args().len() < 3 {
        println!("Usage: ./egypt <numerator> <denominator>");
        return;
    }
    for (i, arg) in env::args().enumerate() {
        if i == 1 {
            num = BigInt::from_str(&arg).unwrap();
        }
        if i == 2 {
            den = BigInt::from_str(&arg).unwrap();
        }
    }
    for (a, b) in as_egyptian_fraction(&num, &den).iter() {
        println!("{:?}\t{:?}", a, b);
    }
    return;
}