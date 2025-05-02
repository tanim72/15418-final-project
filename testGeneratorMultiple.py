
import random, pathlib

M        = 1                # number of independent systems
N        = 131072                # size of each tridiagonal system
OUTFILE  = f"inputs/test{M}x{N}_input.txt"   # output filename

SEED     = None # random seed


def gen_system(n: int, rng: random.Random):
    b, a, c, d = [], [], [], []

    for i in range(n):
        if i > 0: # sub-diag 
            a.append(float(rng.choice([-2, -1, 1, 2])))
        if i < n - 1:  # super-diag 
            c.append(float(rng.choice([-2, -1, 1, 2])))

    for i in range(n):
        dom = (abs(a[i-1]) if i else 0) + (abs(c[i]) if i < n-1 else 0) + 1.0
        b.append(float(dom))
        d.append(float(rng.uniform(-10.0, 10.0)))

    return b, a, c, d


def write_case(path: pathlib.Path, m: int, n: int, seed):
    rng = random.Random(seed)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{m} {n}\n")
        for _ in range(m):
            b, a, c, d = gen_system(n, rng)
            f.write(" ".join(f"{x:.6g}" for x in b) + "\n")
            f.write(" ".join(f"{x:.6g}" for x in a) + "\n")
            f.write(" ".join(f"{x:.6g}" for x in c) + "\n")
            f.write(" ".join(f"{x:.6g}" for x in d) + "\n")
    print(f"Wrote {m}x{n} test case → {path}")


def main():
    if N < 2 or M < 1:
        raise ValueError("Require N ≥ 2 and M ≥ 1.")
    write_case(pathlib.Path(OUTFILE), M, N, SEED)

if __name__ == "__main__":
    main()
