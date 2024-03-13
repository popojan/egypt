use std::collections::HashSet;
use std::io;
use num_bigint::BigInt;
use num_traits::identities::One;
use std::ops::{Add, Sub, Div, Mul, SubAssign, Neg};
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

    /// Extra merge step with quadratic complexity possibly reducing number of terms
    #[clap(short, long, value_parser, default_value_t = false)]
    merge: bool,

    /// Output minimal number of raw quadruplets (aka symbolic sums)
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

    /// Maximum number of terms for breaking large symbolic sums
    #[clap(short, long, value_parser, default_value_t = 2)]
    limit: usize,

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

fn add_fractions(a1:&BigInt, b1:&BigInt, a2:&BigInt, b2:&BigInt) -> (BigInt, BigInt) {
    let den = b1.lcm(b2);
    let num = a1.mul(den.clone().div(b1))
                + a2.mul(den.clone().div(b2));
    let gcd = num.gcd(&den);
    (num.div(&gcd), den.div(&gcd))
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
            (x, y) = add_fractions(&x, &y, &x2, &y2);
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
    let gcd = x0.gcd(&y0);
    let mut x = x0.div(&gcd).clone();
    let mut y = y0.div(&gcd).clone();
    if x.ge(&y) {
        whole.push((x.clone().div(&y), BigInt::zero(), BigInt::zero(), BigInt::zero()));
        x.sub_assign(x.clone().div(&y).mul(y.clone()));
    }
    while x.gt(&BigInt::zero()) && y.gt(&BigInt::one()) {
        let v = y.clone() - low(&x, &y).1;
        let t = x.clone().mul(&y).div(x.clone().mul(&v).add(&BigInt::one()));
        ret.push((y.clone() - t.clone().mul(&v), v.clone(), BigInt::one(), t.clone()));
        let den = y.clone().mul(y.clone() - t.clone().mul(&v));
        (x, y) = add_fractions(&x, &y, &t.clone().neg(), &den);
    }
    if !x.is_zero() {
        ret.push((y.clone(), BigInt::one(), BigInt::zero(), BigInt::zero()));
    }
    ret.extend(whole);
    ret.reverse();
    ret
}

fn expand(eg: &Vec<(BigInt, BigInt, BigInt, BigInt)>) -> Vec<(BigInt, BigInt, BigInt, BigInt)> {
    let mut ret = vec![];
    for (b,v,i,j) in eg.iter() {
        if v.is_zero() && i.is_zero() && j.is_zero() {
            ret.push((b.clone(), BigInt::one(), BigInt::zero(), BigInt::zero()));
            continue;
        }
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

fn as_egyptian_fraction(a:&BigInt, b:&BigInt, args: &Args)->Vec<(BigInt, BigInt,BigInt,BigInt)> {
    let mut res = as_egyptian_fraction_symbolic(
        &a,
        &b,
        args.reverse);
    if !args.raw {
        res = halve_symbolic_sums(&res, args.limit);
        res = expand(&res);
        res = fix_duplicates(&res);
        if args.merge {
            if !args.reverse {
                res.reverse();
            }
            res = merge(&res);
        }
    }
    res
}

fn fix_duplicates(eg: &Vec<(BigInt, BigInt, BigInt, BigInt)>)
    -> Vec<(BigInt, BigInt, BigInt, BigInt)> {
      if eg.is_empty() {
          return eg.clone();
      }
    let mut last_i = 0;
    let mut eg = eg.clone();
    while last_i < eg.len() {
        eg.sort_by(|x, y| { x.1.cmp(&y.1)});
        let mut ret = vec![];
        let mut cnt = 1;
        let mut prev = eg.first().unwrap();
        last_i = eg.len();
        for (i, current) in eg.iter().enumerate().skip(1) {
            if current == prev {
                cnt += 1;
            } else if cnt > 1 {
                last_i = i;
                break;
            } else {
                cnt = 1;
                ret.push(prev.clone());
                prev = current;
            }
        }
        if last_i < eg.len() {
            let a = BigInt::from(cnt);
            let b = prev.clone();
            let gcd = a.gcd(&b.1);
            let new = as_egyptian_fraction_symbolic(&a.div(&gcd), &b.1.div(&gcd), false);
            ret.extend(expand(&new));
            ret.extend(eg[last_i..eg.len()].to_vec());
        } else {
            ret.push(prev.clone());
        }
        if eg == ret {
            break;
        }
        eg.clear();
        eg.extend(ret);
    }
    eg
}

fn calculate_raw_sum(u:&BigInt, v:&BigInt, i:&BigInt, j:&BigInt) -> (BigInt, BigInt) {
    let num = BigInt::one().sub(i).add(j);
    let den = u.clone().sub(v).add(v.clone().mul(i)).mul(u.add(v.mul(j)));
    let gcd = num.gcd(&den);
    (num.div(&gcd), den.div(&gcd))
}

fn halve_symbolic_sums(a: &Vec<(BigInt, BigInt, BigInt, BigInt)>, limit: usize)
    -> Vec<(BigInt, BigInt, BigInt, BigInt)>
{
    let mut stack = a.clone();
    let mut ret = vec![] ;
    let limit = BigInt::from(limit);
    let two = BigInt::from(2);
    while !stack.is_empty() {
        let (u, v, i, j) = stack.pop().unwrap();
        let term_count = j.clone().sub(&i).add(&BigInt::one());
        if term_count.le(&limit) {
            ret.push((u, v, i, j));
        } else {
            let (a, b) = calculate_raw_sum(&u, &v, &i, &j);
            if a.is_odd() {
                let a1 = a.sub(&BigInt::one()).div(&two);
                let a2 = a1.clone().add(&BigInt::one());
                let gcd1 = a1.gcd(&b);
                let gcd2 = a2.gcd(&b);
                stack.extend(as_egyptian_fraction_symbolic(&a1.div(&gcd1), &b.clone().div(&gcd1), false));
                stack.extend(as_egyptian_fraction_symbolic(&a2.div(&gcd2), &b.div(&gcd2), false));
            } else {
                let a1 = a.div(&two).sub(&BigInt::one());
                let a2 = a1.clone().add(&two);
                let gcd1 = a1.gcd(&b);
                let gcd2 = a2.gcd(&b);
                stack.extend(as_egyptian_fraction_symbolic(&a1.div(&gcd1), &b.clone().div(&gcd1), false));
                stack.extend(as_egyptian_fraction_symbolic(&a2.div(&gcd2), &b.div(&gcd2), false));
            }
        }
    }
    ret
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
                    for (i, (a, b, c, d))
                        in as_egyptian_fraction(&num, &den, &args).iter().enumerate() {
                        let is_natural = (args.raw && b.is_zero() && c.is_zero() && d.is_zero())
                            || (!args.raw && b.is_one());
                        if i == 0 && is_natural {
                            print!("{}\t", a);
                            gt0 = true;
                        } else if i == 0 {
                            if !args.raw {
                                print!("0\t{}", b);
                            } else {
                                print!("0\t{},{},{},{}", a, b, c, d);
                            }
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
            &args.numerator, &args.denominator, &args).iter() {
            if !args.silent {
                if !args.raw {
                    println!("{:?}\t{:?}", a, b);
                } else {
                    println!("{:?}\t{:?}\t{:?}\t{:?}", a, b, c, d);
                }
            }
        }
    }
    return;
}
