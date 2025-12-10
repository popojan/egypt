#!/usr/bin/env python3
"""
Pell equation solver using Egypt output.
Usage: egypt "D sqrt" 1 --raw -p PREC | python3 scripts/pell.py D

Input: Egypt raw tuples (u v i j) from stdin
Output: q  p  norm (tab-separated)
Exit 0 if Pell solution found, 1 otherwise.
"""

import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 pell.py D < egypt_output", file=sys.stderr)
        sys.exit(2)

    D = int(sys.argv[1])

    # Parse Egypt output, extract q values
    # Tuple format: u v i j (tab-separated)
    # u = q_{2k-2}, v = q_{2k-1}
    # Integer part: (n, 0, 0, 0)

    a0 = 0  # floor(sqrt(D))
    qs = []  # q values: q_0, q_1, q_2, ...

    for line in sys.stdin:
        parts = line.strip().split('\t')
        if len(parts) != 4:
            continue
        u, v, i, j = map(int, parts)

        if v == 0 and i == 0 and j == 0:
            a0 = u  # integer part
        else:
            # Fractional tuple: u = q_{2k-2}, v = q_{2k-1}
            if not qs or qs[-1] != u:
                qs.append(u)
            qs.append(v)

    if a0 == 0:
        print("No integer part found", file=sys.stderr)
        sys.exit(2)

    # Compute p values using CF recurrence:
    # p_{-1} = 1, p_0 = a0
    # p_n = a_n * p_{n-1} + p_{n-2}
    # where a_n = (q_n - q_{n-2}) / q_{n-1}

    # For q indices: q[0] = q_0, q[1] = q_1, etc. (0-based in our list)
    # We need q_{-1} = 0 implicitly

    ps = [1, a0]  # p_{-1}, p_0

    for n in range(len(qs)):
        q_n = qs[n]
        if n == 0:
            # a_1 = q_1 / q_0 = q_1 (since q_{-1}=0: a_1 = (q_1 - 0)/q_0 = q_1/1)
            # But wait, q[0] is q_0, not q_1. Let me reconsider.
            # Actually qs[0] could be q_0=1 or q_1 depending on tuples
            # From Egypt: first tuple has u=q_0=1, v=q_1
            # So qs = [1, q_1, q_2, q_3, ...]
            # For n=0: q_n = qs[0] = q_0 = 1
            # a_0 is already handled (integer part)
            # For n=1 (computing p_1):
            #   a_1 = (q_1 - q_{-1}) / q_0 = q_1 / 1 = q_1
            continue  # q_0 already covered by a0
        elif n == 1:
            # a_1 = q_1 (since q_{-1}=0)
            a_n = qs[1]
        else:
            # a_n = (q_n - q_{n-2}) / q_{n-1}
            a_n = (qs[n] - qs[n-2]) // qs[n-1]

        p_n = a_n * ps[-1] + ps[-2]
        ps.append(p_n)

    # Output
    print("q\tp\tnorm")

    found_solution = False

    # First: q_0 = 1, p_0 = a0
    norm = a0*a0 - D
    print(f"1\t{a0}\t{norm}")
    if abs(norm) == 1:
        found_solution = True

    # Rest: q_n, p_n for n >= 1
    quasi_found = False
    for n in range(1, len(qs)):
        q = qs[n]
        p = ps[n+1]  # ps[0]=p_{-1}, ps[1]=p_0, ps[2]=p_1, so p_n = ps[n+1]
        norm = p*p - D*q*q
        print(f"{q}\t{p}\t{norm}")
        if norm == -1 and not quasi_found:
            print(f"# Quasi-solution (norm=-1): p={p}, q={q}", file=sys.stderr)
            quasi_found = True
        if norm == 1:
            print(f"# Fundamental solution (norm=1): p={p}, q={q}", file=sys.stderr)
            found_solution = True
            break

    if found_solution:
        sys.exit(0)
    elif quasi_found:
        print("# No fundamental solution, but quasi-solution exists", file=sys.stderr)
        sys.exit(1)
    else:
        print("# No Pell solution found", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
