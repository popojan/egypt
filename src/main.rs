use std::collections::HashSet;
use std::io;
use num_bigint::BigInt;
use num_traits::identities::One;
use std::ops::{Add, Sub, Div, Mul, SubAssign, MulAssign, DivAssign};
use num_integer::Integer;
use clap::Parser;
use std::str::FromStr;
use num_traits::{ToPrimitive, Zero};

/// Egyptian Fractions

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
#[command(arg_required_else_help(true))]
struct Args {
    /// Reverse merge strategy
    #[clap(short, long, value_parser, default_value_t = false)]
    reverse: bool,

    /// Do not merge
    #[clap(long, value_parser, default_value_t = false)]
    raw: bool,

    /// No output
    #[clap(short, long, value_parser, default_value_t = false)]
    silent: bool,

    /// Batch mode (expects numerator and denominator on each line of standard input)
    #[clap(long, value_parser, default_value_t = false)]
    batch: bool,

    #[clap(value_parser, default_value_t = BigInt::from(1))]
    numerator: BigInt,

    #[clap(value_parser, default_value_t = BigInt::from(1))]
    denominator: BigInt,
}

fn low(x: &BigInt, y: &BigInt) -> (BigInt, BigInt) {
    let b = if y.sub(x).eq(&BigInt::one()) {
        x.clone()
    } else {
        x.extended_gcd(&y).x.add(y).mod_floor(&y)
    };
    let m = x.clone().mul(&b);
    let a = m.clone().div(y);
    let gcd = a.gcd(&b);
    return (a.div(&gcd), b.div(&gcd))
}

fn _as_egyptian_fraction_raw(x0: &BigInt, y0: &BigInt, _reverse: bool) -> Vec<(BigInt, BigInt)> {
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
    if !x.is_zero() {
        ret.extend(vec![(x.clone(), y.clone())]);
    }
    ret.extend(whole);
    ret.sort_by(|x, y| { x.1.cmp(&y.1)});
    ret
}

fn merge(eg: &Vec<(BigInt, BigInt, BigInt, BigInt)>) -> Vec<(BigInt, BigInt, BigInt, BigInt)> {
    let mut i = 0_usize;
    let mut ret = HashSet::<(BigInt, BigInt, BigInt, BigInt)>::new();
    let egs = HashSet::<(BigInt, BigInt, BigInt, BigInt)>::from_iter(
        eg.clone().into_iter());
    while i < eg.len() {
        let (mut x, mut y, _, _) = eg[i].clone();
        let mut ones = vec![(i, x.clone(), y.clone())];
        for j in (i + 1)..eg.len() {
            let (x2, y2, _, _) = &eg[j];
            (x, y) = (x.mul(y2).add(x2.mul(y.clone())), y.clone().mul(y2));
            let gcd = x.gcd(&y);
            (x, y) = (x.div(&gcd), y.div(&gcd));
            let tup = (x.clone(), y.clone(), BigInt::zero(), BigInt::zero());
            //println!("{:?}/{:?}: {} / {}, {}", x, y, j, i, eg.len());
            if x.eq(&BigInt::one()) && !ret.contains(&tup)
                && !egs.contains(&tup)
            {
                ones.push((j, x.clone(), y.clone()));
            }
        }
        let (j, x, y) = ones.last().unwrap();
        //println!("{:?} {:?} {:?}", ones, x, y);
        ret.insert((x.clone(), y.clone(), BigInt::zero(), BigInt::zero()));
        i = j + 1;
    }
    let mut ret = ret
        .into_iter().collect::<Vec<(BigInt, BigInt, BigInt, BigInt)>>();
    ret.sort_by(|x, y| { x.1.cmp(&y.1)});
    ret
}

fn as_egyptian_fraction_symbolic(x0: &BigInt, y0: &BigInt, _expand: bool) -> Vec<(BigInt, BigInt, BigInt, BigInt)> {
    let mut whole = Vec::<(BigInt, BigInt, BigInt, BigInt)>::new();
    let mut ret = Vec::<(BigInt, BigInt, BigInt, BigInt)>::new();
    let mut x = x0.clone();
    let mut y = y0.clone();
    if x.ge(&y) {
        whole.push((BigInt::one(), BigInt::zero(), BigInt::one(), x.clone().div(&y)));
        x.sub_assign(x.clone().div(&y).mul(y.clone()));
    }
    while x.gt(&BigInt::zero()) && y.gt(&BigInt::one()) {
        let v = y.clone() - low(&x, &y).1;
        let t = x.clone().mul(&y).div(x.clone().mul(&v).add(&BigInt::one()));
        ret.push((y.clone() - t.clone().mul(&v), v.clone(), BigInt::one(), t.clone()));
        let mul = y.clone() - t.clone().mul(&v);
        y.mul_assign(&mul);
        x.mul_assign(&mul);
        x.sub_assign(&t);
        let gcd = x.clone().gcd(&y);
        x.div_assign(&gcd);
        y.div_assign(&gcd);
    }
    if !x.is_zero() {
        ret.push((y.clone(), BigInt::one(), BigInt::one(), BigInt::one()));
    }
    ret.extend(whole);
    ret.reverse();//sort_by(|x, y| { x.1.cmp(&y.1)});
    ret
}

fn expand(eg: &Vec<(BigInt, BigInt, BigInt, BigInt)>) -> Vec<(BigInt, BigInt, BigInt, BigInt)> {
    let mut ret = vec![];
    for (b,v,i,j) in eg.iter() {
        for k in i.to_usize().unwrap() .. j.to_usize().unwrap() + 1 {
            let num = BigInt::one();
            let den = b.sub(v).add(v.clone().mul(&k))
                .mul(b.add(v.clone().mul(k.clone())));
            ret.push((num, den, BigInt::zero(), BigInt::zero()));
        }
    }
    ret
}
/// Converts rational numbers to egyptian fractions

fn as_egyptian_fraction(a:&BigInt, b:&BigInt, raw:bool, reverse:bool)->Vec<(BigInt, BigInt,BigInt,BigInt)> {
    let mut res = as_egyptian_fraction_symbolic(
        &a,
        &b,
        reverse);
    if !raw {
        res = expand(&res);
        if !reverse {
            res.reverse();
        }
        res = merge(&res);
    }
    res
}

fn main() {
    let args = Args::parse();

    if args.batch {
        for line in io::stdin().lines() {
            if let Ok(line) = line {
                let num_den = line.split_ascii_whitespace().take(2).collect::<Vec<&str>>();
                let num = BigInt::from_str(num_den[0]).unwrap();
                let den = BigInt::from_str(num_den[1]).unwrap();
                if !args.silent {
                    let mut gt0 = false;
                    print!("{}\t{}\t", num.to_str_radix(10), den.to_str_radix(10));
                    for (i, (a, b, c, d)) in as_egyptian_fraction(&num, &den, args.raw, args.reverse).iter().enumerate() {
                        if i == 0 && b.is_one() {
                            print!("{}\t", a);
                            gt0 = true;
                        } else if i == 0 {
                            print!("0\t{}", b);
                        } else if i == 1 && gt0 {
                            if args.raw {
                                print!("{},{},{},{}", a, b, c, d);
                            } else {
                                print!("{}", b);
                            }
                        } else {
                            if args.raw {
                                print!(" {},{},{},{}", a, b, c, d);
                            } else {
                                print!(" {}", b);
                            }
                        }
                    }
                    print!("\n");
                }
            }
        }
    } else {
        for (a, b, c, d) in as_egyptian_fraction(
            &args.numerator, &args.denominator, args.raw, args.reverse).iter() {
            if !args.silent {
                if c.is_zero() && d.is_zero() {
                    println!("{:?}\t{:?}", a, b);
                } else {
                    println!("{:?}\t{:?}\t{:?}\t{:?}", a, b, c, d);
                }
            }
        }
    }
    return;
}