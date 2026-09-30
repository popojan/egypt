//! CLI regression tests: golden outputs per switch, plus the invariants of
//! wl/test.wls (unit numerators, no duplicate denominators, exact sum).
//! Regenerate goldens with `UPDATE_GOLDEN=1 cargo test --release --test cli`.

use rug::{Integer, Rational};
use std::collections::HashSet;
use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

struct Run {
    stdout: String,
    stderr: String,
    code: i32,
}

fn egypt_stdin(args: &[&str], stdin: &str) -> Run {
    let mut child = Command::new(env!("CARGO_BIN_EXE_egypt"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn egypt");
    child.stdin.take().unwrap().write_all(stdin.as_bytes()).unwrap();
    let out = child.wait_with_output().unwrap();
    Run {
        stdout: String::from_utf8(out.stdout).unwrap(),
        stderr: String::from_utf8(out.stderr).unwrap(),
        code: out.status.code().unwrap_or(-1),
    }
}

fn egypt(args: &[&str]) -> Run {
    egypt_stdin(args, "")
}

fn golden_stdin(name: &str, args: &[&str], stdin: &str) -> Run {
    let run = egypt_stdin(args, stdin);
    assert_eq!(run.code, 0, "{:?}: {}", args, run.stderr);
    let path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "tests", "golden", &format!("{}.out", name)]
        .iter()
        .collect();
    if std::env::var_os("UPDATE_GOLDEN").is_some() {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, &run.stdout).unwrap();
    }
    let expected = std::fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("missing {}; run with UPDATE_GOLDEN=1", path.display()));
    assert_eq!(run.stdout, expected, "{:?}", args);
    run
}

fn golden(name: &str, args: &[&str]) -> Run {
    golden_stdin(name, args, "")
}

fn ratio(a: &str, b: &str) -> Rational {
    Rational::from((a.parse::<Integer>().unwrap(), b.parse::<Integer>().unwrap()))
}

/// "a\tb" lines of a non-raw run
fn fractions(stdout: &str) -> Vec<(Integer, Integer)> {
    stdout
        .lines()
        .map(|line| {
            let mut cols = line.split('\t').map(|c| c.parse::<Integer>().unwrap());
            (cols.next().unwrap(), cols.next().unwrap())
        })
        .collect()
}

/// Same checks as wl/test.wls, plus ascending denominators
fn assert_expansion(args: &[&str], value: &Rational) {
    let run = egypt(args);
    assert_eq!(run.code, 0, "{:?}: {}", args, run.stderr);
    let mut sum = Rational::new();
    let mut seen = HashSet::new();
    let mut prev = Integer::from(0);
    for (i, (a, b)) in fractions(&run.stdout).iter().enumerate() {
        sum += Rational::from((a.clone(), b.clone()));
        if i == 0 && *b == 1 {
            continue; // integer part
        }
        assert_eq!(*a, 1, "non-unit term {}/{} for {:?}", a, b, args);
        assert!(*b > prev, "denominators not ascending at {} for {:?}", b, args);
        assert!(seen.insert(b.clone()), "duplicate 1/{} for {:?}", b, args);
        prev = b.clone();
    }
    assert_eq!(&sum, value, "{:?}", args);
}

#[test]
fn readme_examples() {
    golden("readme_7_19_merge", &["--merge", "--limit", "19", "7", "19"]);
    golden("readme_2023_2024_merge", &["--merge", "--limit", "2023", "2023", "2024"]);
    golden("readme_2023_2024_reverse_merge", &["--reverse", "--merge", "--limit", "2023", "2023", "2024"]);
    golden("readme_2023_2024_limit2", &["--limit", "2", "2023", "2024"]);
    golden("readme_999999_limit2", &["--limit", "2", "999999", "1000000"]);
}

#[test]
fn raw_bisect_and_limit_switches() {
    golden("raw_7_19", &["--raw", "7", "19"]);
    golden("raw_bisect_limit3", &["--raw", "--bisect", "-l", "3", "58", "3511471"]);
    golden("raw_bisect_greedy", &["--raw", "--bisect", "-l", "8", "-g", "58", "3511471"]);
    golden("limit_1", &["-l", "1", "58", "3511471"]);
    golden("integer_part", &["22", "7"]);
    golden("exact_expression", &["1 3 /", "0.7"]);
}

#[test]
fn greedy_switch() {
    golden("greedy_limit8", &["-l", "8", "-g", "58", "3511471"]);
    golden("greedy_limit2_2023", &["-l", "2", "-g", "2023", "2024"]);
}

#[test]
fn hex_switch() {
    let args = ["-l", "8", "2 107 ^ 1 -", "2 127 ^ 1 -"];
    let dec = egypt(&args);
    let hex_args: Vec<&str> = [&["-x"][..], &args[..]].concat();
    let hex = egypt(&hex_args);
    let decoded: String = hex
        .stdout
        .lines()
        .map(|line| {
            let cols: Vec<String> = line
                .split('\t')
                .map(|c| Integer::from_str_radix(c, 16).unwrap().to_string())
                .collect();
            cols.join("\t") + "\n"
        })
        .collect();
    assert_eq!(decoded, dec.stdout);
    golden("hex_raw_58", &["--hex", "--raw", "58", "3511471"]);
    golden_stdin("batch_hex", &["--batch", "--hex", "-l", "2"], "22\t7\n2023\t2024\n");
}

#[test]
fn unbisected_output_matches_raw_tuples() {
    // with the limit above every tuple length, the terms are exactly (u+v(k-1))(u+vk)
    let raw = egypt(&["--raw", "58", "3511471"]);
    let mut expected: Vec<Integer> = vec![];
    for line in raw.stdout.lines() {
        let t: Vec<Integer> = line.split('\t').map(|c| c.parse().unwrap()).collect();
        let (u, v) = (&t[0], &t[1]);
        for k in t[2].to_u64().unwrap()..=t[3].to_u64().unwrap() {
            let s_prev = Integer::from(v * (k - 1)) + u;
            let s_k = Integer::from(v * k) + u;
            expected.push(s_prev * s_k);
        }
    }
    expected.sort();
    let full = egypt(&["-l", "100000", "58", "3511471"]);
    let got: Vec<Integer> = fractions(&full.stdout).into_iter().map(|(_, b)| b).collect();
    assert_eq!(got, expected);
}

#[test]
fn precision_switch_and_irrational_goldens() {
    let p64 = egypt(&["pi", "1", "--raw", "-p", "64"]);
    let p128 = egypt(&["pi", "1", "--raw", "-p", "128"]);
    let (a, b): (Vec<&str>, Vec<&str>) = (p64.stdout.lines().collect(), p128.stdout.lines().collect());
    assert!(b.len() > a.len(), "more precision must give more tuples");
    assert_eq!(a[..5], b[..5], "leading tuples are stable under a precision increase");
    golden("pi_4_raw_p64", &["pi", "4", "--raw", "-p", "64"]);
    golden("phi_raw_p64", &["phi", "1", "--raw", "-p", "64"]);
    golden("sqrt7_11_p32", &["-p", "32", "7 sqrt", "11"]);
    golden("sqrt7_11_p32_raw", &["-p", "32", "--raw", "7 sqrt", "11"]);
}

#[test]
fn pell_switch() {
    let run = golden("pell_13", &["13 sqrt", "1", "--pell", "-p", "64"]);
    assert!(run.stderr.contains("Fundamental solution (norm=1): p=649, q=180"), "{}", run.stderr);
    let run = golden("pell_cattle", &["4729494 sqrt", "1", "--pell", "-p", "512"]);
    assert!(run.stdout.lines().last().unwrap().ends_with("\t1"));
}

#[test]
fn batch_switch() {
    let input = "7\t19\n22\t7\n2023\t2024\npi\t4\n2 107 ^ 1 -\t2 127 ^ 1 -\n7 19\n";
    golden_stdin("batch_limit8", &["--batch", "-l", "8"], input);
    golden_stdin("batch_raw", &["--batch", "--raw"], input);
    golden_stdin("batch_greedy", &["--batch", "-l", "2", "-g"], input);
}

#[test]
fn silent_and_errors() {
    let run = egypt(&["-s", "7", "19"]);
    assert_eq!((run.code, run.stdout.as_str()), (0, ""));
    let run = egypt(&["1", "0"]);
    assert_eq!(run.code, 2);
    assert!(run.stderr.contains("zero"), "{}", run.stderr);
    let run = egypt(&["7", "1", "--pell"]);
    assert_eq!(run.code, 2);
    let run = egypt(&[]);
    assert_eq!(run.code, 2, "no arguments prints help");
    assert!(run.stdout.contains("Usage") || run.stderr.contains("Usage"));
}

#[test]
fn big_literal_equals_rpn_form() {
    let m107 = "162259276829213363391578010288127";
    let m127 = "170141183460469231731687303715884105727";
    let literal = egypt(&["-l", "8", m107, m127]);
    let rpn = egypt(&["-l", "8", "2 107 ^ 1 -", "2 127 ^ 1 -"]);
    assert_eq!(literal.stdout, rpn.stdout);
    assert_eq!(literal.stdout.lines().count(), 169);
    assert_expansion(&["-l", "8", m107, m127], &ratio(m107, m127));
}

#[test]
fn switch_matrix_keeps_invariants() {
    let inputs = [
        ("7", "19", ratio("7", "19")),
        ("22", "7", ratio("22", "7")),
        ("1", "1", ratio("1", "1")),
        ("2023", "2024", ratio("2023", "2024")),
        ("999999", "1000000", ratio("999999", "1000000")),
        ("58", "3511471", ratio("58", "3511471")),
        ("49 sqrt", "11", ratio("7", "11")),
        ("1 3 /", "0.7", ratio("10", "21")),
    ];
    let modes: [&[&str]; 8] = [
        &[], &["-l", "2"], &["-g"], &["-l", "2", "-g"],
        &["-m"], &["-r", "-m"], &["-m", "-l", "2"], &["-m", "-g", "-l", "3"],
    ];
    for (a, b, value) in &inputs {
        for mode in modes {
            let mut args = mode.to_vec();
            args.extend([*a, *b]);
            assert_expansion(&args, value);
        }
    }
}

/// Deterministic stand-in for wl/test.wls: random rationals with 48-bit parts
#[test]
fn random_rationals_keep_invariants() {
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    for _ in 0..40 {
        let a = ((next() >> 16) + 1).to_string();
        let b = ((next() >> 16) + 1).to_string();
        let value = ratio(&a, &b);
        for mode in [&[][..], &["-l", "2"][..], &["-g"][..], &["-m"][..]] {
            let mut args = mode.to_vec();
            args.extend([a.as_str(), b.as_str()]);
            assert_expansion(&args, &value);
        }
    }
}
