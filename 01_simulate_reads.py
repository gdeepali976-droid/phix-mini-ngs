"""Step 1: simulate Illumina-style single-end reads from phiX174.

Reads are 100 bp, drawn from random positions on both strands (wrapping
round the circular genome). Quality falls along the read, a fraction of
reads have a bad 3' tail, and every base is miscalled with exactly the
probability its Phred score claims. Because the true position of each read
is stored in its name, later steps can grade QC and assembly against truth.

Usage: python 01_simulate_reads.py [coverage]    (default 30)
"""
import random
import sys
from ngstools import path, read_fasta, revcomp, write_fastq, qstring

READ_LEN = 100


def simulate(ref, coverage, seed=42, read_len=READ_LEN):
    rng = random.Random(seed)
    L = len(ref)
    ref2 = ref + ref[:read_len]                 # circular genome
    n = round(coverage * L / read_len)
    reads = []
    for i in range(n):
        pos = rng.randrange(L)
        frag = ref2[pos:pos + read_len]
        fwd = rng.random() < 0.5
        true = frag if fwd else revcomp(frag)
        bad_from = rng.randint(55, 90) if rng.random() < 0.15 else None
        quals, seq = [], []
        for j, base in enumerate(true):
            q = 38 - 12 * (j / read_len) ** 2 + rng.gauss(0, 2)
            if bad_from is not None and j >= bad_from:
                q -= 16
            q = max(2, min(41, round(q)))
            quals.append(q)
            if rng.random() < 10 ** (-q / 10):  # the base really is wrong
                base = rng.choice([b for b in "ACGT" if b != base])
            seq.append(base)
        reads.append((f"r{i}_{pos}_{'+' if fwd else '-'}", "".join(seq), qstring(quals)))
    return reads


if __name__ == "__main__":
    cov = float(sys.argv[1]) if len(sys.argv) > 1 else 30
    ref = read_fasta(path("phix174.fasta"))[0][1]
    reads = simulate(ref, cov)
    write_fastq(path("reads.fastq"), reads)
    print(f"{len(reads)} reads x {READ_LEN} bp = {len(reads)*READ_LEN} bases "
          f"(~{cov:.0f}x of a {len(ref)} bp genome)")
