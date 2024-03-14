use std::collections::HashSet;
use std::io;
use std::ops::{Add, Sub, Div, Mul, SubAssign, Neg};
use clap::Parser;
use std::str::FromStr;
use rug::{Complete, Integer};
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

    /// Output raw quadruplets bisected according to --limit
    #[clap(long, value_parser, default_value_t = false)]
    bisect: bool,

    /// No output
    #[clap(short, long, value_parser, default_value_t = false)]
    silent: bool,

    /// Batch mode (expects numerator and denominator on each line of standard input)
    #[clap(long, value_parser, default_value_t = false)]
    batch: bool,

    #[clap(value_parser, default_value_t = Integer::from(1))]
    numerator: Integer,

    #[clap(value_parser, default_value_t = Integer::from(1))]
    denominator: Integer,

    /// Maximum number of terms for breaking large symbolic sums
    #[clap(short, long, value_parser, default_value_t = 8)]
    limit: usize,

}

fn low(x: &Integer, y: &Integer) -> (Integer, Integer) {
    let b = if (y - x).complete() == Integer::from(1) {
        x.clone()
    } else {
        x.clone().extended_gcd(y.clone(), Integer::new()).1.add(y).modulo(y)
    };
    let m = x.clone().mul(&b);
    let a = m.clone().div(y);
    let gcd = a.clone().gcd(&b);
    return (a.div(&gcd), b.clone().div(&gcd))
}

fn add_fractions(a1:&Integer, b1:&Integer, a2:&Integer, b2:&Integer) -> (Integer, Integer) {
    let den = b1.clone().lcm(b2);
    let num = a1.mul(den.clone().div(b1))
                + a2.mul(den.clone().div(b2));
    let gcd = num.clone().gcd(&den);
    (num.div(&gcd), den.div(&gcd))
}

fn _as_egyptian_fraction_raw(x0: &Integer, y0: &Integer, _reverse: bool) -> Vec<(Integer, Integer)> {
    let mut whole = Vec::<(Integer, Integer)>::new();
    let mut ret = Vec::<(Integer, Integer)>::new();
    let mut x = x0.clone();
    let mut y = y0.clone();
     if x.ge(&y) {
        whole.push((x.clone().div(&y), Integer::from(1)));
        x = x.clone().sub(x.div(&y).mul(y.clone()));
    }
    while x.gt(&Integer::from(1)) {
        let (x2, y2) = low(&x, &y);
        let (a, b) = (x.clone().mul(&y2).sub(x2.clone().mul(&y)), y.clone().mul(&y2));
        let gcd = a.clone().gcd(&b);
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

fn merge(eg: &Vec<(Integer, Integer, Integer, Integer)>) -> Vec<(Integer, Integer, Integer, Integer)> {
    let mut i = 0_usize;
    let mut ret = HashSet::<(Integer, Integer, Integer, Integer)>::new();
    while i < eg.len() {
        let (mut x, mut y, _, _) = eg[i].clone();
        let (mut ones_i, mut ones_x, mut ones_y) = (i, x.clone(), y.clone());
        for j in (i + 1)..eg.len() {
            let (x2, y2, _, _) = &eg[j];
            (x, y) = add_fractions(&x, &y, &x2, &y2);
            if x == Integer::from(1) {
                (ones_i, ones_x, ones_y) = (j, x.clone(), y.clone());
            }
        }
        let (j, x, y) = (ones_i, ones_x, ones_y);
        //println!("{:?} {:?} {:?}", ones, x, y);
        ret.insert((x.clone(), y.clone(), Integer::from(0), Integer::from(0)));
        i = j + 1;
    }
    let mut ret = ret
        .into_iter().collect::<Vec<(Integer, Integer, Integer, Integer)>>();
    ret.sort_by(|x, y| { x.1.cmp(&y.1)});
    ret
}

fn as_egyptian_fraction_symbolic(x0: &Integer, y0: &Integer, _expand: bool, ret: &mut Vec<(Integer, Integer, Integer, Integer)>) {
    let mut whole = Vec::<(Integer, Integer, Integer, Integer)>::new();
    let gcd = x0.clone().gcd(&y0);
    let mut x = x0.div(&gcd).complete();
    let mut y = y0.div(&gcd).complete();
    if x.ge(&y) {
        whole.push((x.clone().div(&y), Integer::from(0), Integer::from(0), Integer::from(0)));
        x.sub_assign(x.clone().div(&y).mul(y.clone()));
    }
    while x.gt(&Integer::from(0)) && y.gt(&Integer::from(1)) {
        let v = y.clone() - low(&x, &y).1;
        let t = x.clone().mul(&y).div(x.clone().mul(&v).add(&Integer::from(1)));
        ret.push((y.clone() - t.clone().mul(&v), v.clone(), Integer::from(1), t.clone()));
        let den = y.clone().mul(y.clone() - t.clone().mul(&v));
        (x, y) = add_fractions(&x, &y, &t.clone().neg(), &den);
    }
    if !x.is_zero() {
        ret.push((y.clone(), Integer::from(1), Integer::from(0), Integer::from(0)));
    }
    ret.extend(whole);
    ret.reverse();
}

fn expand(eg: &Vec<(Integer, Integer, Integer, Integer)>) -> Vec<(Integer, Integer, Integer, Integer)> {
    let mut ret = vec![];
    for (b,v,i,j) in eg.iter() {
        if v.is_zero() && i.is_zero() && j.is_zero() {
            ret.push((b.clone(), Integer::from(1), Integer::from(0), Integer::from(0)));
            continue;
        }
        for k in i.to_usize().unwrap() .. j.to_usize().unwrap() + 1 {
            let num = Integer::from(1);
            let den = b.sub(v).complete().add(v.clone().mul(&k))
                .mul(b.add(v.clone().mul(k.clone())));
            ret.push((num, den, Integer::from(0), Integer::from(0)));
        }
    }
    ret
}
/// Converts rational numbers to egyptian fractions

fn as_egyptian_fraction(a:&Integer, b:&Integer, args: &Args)->Vec<(Integer, Integer,Integer,Integer)> {
    let mut res = vec![];
    as_egyptian_fraction_symbolic(
        &a,
        &b,
        args.reverse, &mut res);
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
    } else {
        if args.bisect {
            res = halve_symbolic_sums(&res, args.limit);
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
                let num_den = line.split_ascii_whitespace().take(2).collect::<Vec<&str>>();
                let num = Integer::from_str(num_den[0]).unwrap();
                let den = Integer::from_str(num_den[1]).unwrap();
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
