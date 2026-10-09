#!/usr/bin/env python3
"""Terms and max denominator bits per strategy (doc/divisor-descent.md).
Usage: scripts/bench_descent.py [path/to/egypt] [samples]"""
import random, statistics, subprocess, sys, time
from fractions import Fraction

EGYPT = sys.argv[1] if len(sys.argv) > 1 else "target/release/egypt"
SAMPLES = int(sys.argv[2]) if len(sys.argv) > 2 else 20

STRATEGIES = [
    ("default", []), ("--merge", ["-m"]), ("--reverse --merge", ["-r", "-m"]),
    ("--greedy", ["-g"]), ("--limit 2", ["-l", "2"]),
    ("--descent", ["--descent"]),
]


def run(args, a, b):
    t0 = time.perf_counter()
    out = subprocess.run([EGYPT, "-x", *args, str(a), str(b)], capture_output=True, text=True, check=True).stdout
    dt = time.perf_counter() - t0
    terms = [tuple(int(c, 16) for c in line.split("\t")) for line in out.splitlines()]
    assert sum(Fraction(n, d) for n, d in terms) == Fraction(a, b), (args, a, b)
    dens = [d for n, d in terms]
    assert len(set(dens)) == len(dens) and all(n == 1 for n, _ in terms), (args, a, b)
    return len(terms), max(dens).bit_length(), dt


def row(label, cases, args):
    res = [run(args, a, b) for a, b in cases]
    terms = [r[0] for r in res]
    bits = [r[1] for r in res]
    secs = max(r[2] for r in res)
    if len(cases) == 1:
        return f"{terms[0]} / {bits[0]}"
    return f"{statistics.median(terms):g} ({max(terms)}) / {statistics.median(bits):g}" + (f" [{secs:.1f}s]" if secs > 1 else "")


def main():
    readme = [("7/19", ["-m", "-l", "19"], 7, 19), ("2023/2024", ["-m", "-l", "2023"], 2023, 2024),
              ("2023/2024", ["-r", "-m", "-l", "2023"], 2023, 2024)]
    for name, args, a, b in readme:
        print(f"README {name} `{' '.join(args)}`: {row('', [(a, b)], args)}")
    rng = random.Random(304)
    groups = [("7/19", [(7, 19)]), ("2023/2024", [(2023, 2024)])]
    for e in (3, 6, 12, 30):
        b = 10 ** e
        groups.append((f"(10^{e}-1)/10^{e}", [(b - 1, b)]))
        cases = []
        while len(cases) < SAMPLES:
            bb = rng.randrange(b, 2 * b)
            a = rng.randrange(1, bb)
            cases.append((a, bb))
        groups.append((f"random, b~10^{e}", cases))
    names = [g[0] for g in groups]
    print("| strategy | " + " | ".join(names) + " |")
    print("|---" * (len(names) + 1) + "|")
    for label, args in STRATEGIES:
        cells = [row(label, cases, args) for _, cases in groups]
        print(f"| `{label}` | " + " | ".join(cells) + " |", flush=True)


if __name__ == "__main__":
    main()
