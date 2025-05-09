mod rpn;

use crate::rpn::_parse_rpn;

use std::io;
use std::ops::{Add, Sub, Div, Mul, SubAssign, Neg};
use clap::Parser;
use rug::{Complete, Integer, Rational};

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

    /// No output
    #[clap(short, long, value_parser, default_value_t = false)]
    silent: bool,

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

fn as_egyptian_fraction_symbolic(x0: &Integer, y0: &Integer, _expand: bool, ret: &mut Vec<(Integer, Integer, Integer, Integer)>) {
    let gcd = x0.clone().gcd(&y0);
    let mut x = x0.clone().div(&gcd);
    let mut y = y0.clone().div(&gcd);
    if x.ge(&y) {
        ret.push((x.clone().div(&y),0.into(), 0.into(), 0.into()));
        x.sub_assign(x.clone().div(&y).mul(&y));
    }
    while x.gt(&Integer::from(0)) && y.gt((&1).into()) {
        let v = x.clone().neg().invert(&y).unwrap();
        let t;
        (t, x) = x.clone().div_rem( (x.clone() * &v + 1) / &y);
        y -= v.clone() * &t;
        ret.push((y.clone(), v, 1.into(), t));
    }
    if !x.is_zero() {
        ret.push((y, Integer::from(1), Integer::from(0), Integer::from(0)));
    }
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
    let mut res = vec![];
    as_egyptian_fraction_symbolic(
        &a,
        &b,
        args.reverse, &mut res);
    let limit = args.limit.max(2);
    if !args.raw {
        res = halve_symbolic_sums(&res, limit);
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
            res = halve_symbolic_sums(&res, limit);
        }
    }
    res
}

fn fix_duplicates(eg: &Vec<(Integer, Integer, Integer, Integer)>)
    -> Vec<(Integer, Integer, Integer, Integer)> {
      if eg.is_empty() {
          return eg.clone();
      }
    let mut last_i = 0;
    let mut eg = eg.clone();
    while last_i < eg.len() {
        eg.sort_by(|x, y| { y.1.cmp(&x.1)});
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
            let a = Integer::from(cnt);
            let b = prev.clone();
            let gcd = a.clone().gcd(&b.1);
            let mut new = vec![];
            as_egyptian_fraction_symbolic(&a.div(&gcd), &b.1.div(&gcd), false, &mut new);
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
    eg.sort_by(|x, y| { x.1.cmp(&y.1)});
    eg
}

fn calculate_raw_sum(u:&Integer, v:&Integer, i:&Integer, j:&Integer) -> (Integer, Integer) {
    let num = Integer::from(1).sub(i).add(j);
    let den = u.clone().sub(v).add(v.clone().mul(i)).mul(u.add(v.mul(j)).complete());
    let gcd = num.clone().gcd(&den);
    (num.div(&gcd), den.div(&gcd))
}

fn halve_symbolic_sums(a: &Vec<(Integer, Integer, Integer, Integer)>, limit: usize)
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
            if a.is_odd() {
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
fn main() {
    let args = Args::parse();

    if args.batch {
        for line in io::stdin().lines() {
            if let Ok(line) = line {
                let num_den = line.split("\t").take(2).collect::<Vec<&str>>();
                if num_den.len() < 2 {
                    println!("expecting tab delimited numerator and denominator");
                    continue;
                }

                let num = _parse_rpn(num_den[0]).abs();
                let den = _parse_rpn(num_den[1]).abs();
                if !args.silent {
                    let mut gt0 = false;
                    print!("{}\t{}\t", num.to_string(), den.to_string());
                    for (i, (a, b, c, d))
                        in as_egyptian_fraction(&num, &den, &args).iter().enumerate() {
                        let is_natural = (args.raw && b.is_zero() && c.is_zero() && d.is_zero())
                            || (!args.raw && *b == Integer::from(1));
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
            &_parse_rpn(&args.numerator).abs(), &_parse_rpn(&args.denominator).abs(), &args).iter() {
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