"""Sanity checks. Run with:  python test_known_values.py

Most of these exist because the bug they guard against is silent: a wrong
strand, an off-by-one in a k-mer, or quality scores that lie all give
plausible-looking output.
"""
import importlib
import random
from assembler import assemble
from ngstools import (revcomp, canonical, n50, phred, qstring, kmer_set,
                      path, read_gbk, read_fasta, trim_3prime)

sim = importlib.import_module("01_simulate_reads")


def test_revcomp_and_canonical():
    assert revcomp("AACG") == "CGTT"
    assert revcomp(revcomp("ACGTTGCA")) == "ACGTTGCA"
    assert canonical("TTT") == canonical("AAA") == "AAA"


def test_phred_roundtrip():
    q = [2, 20, 30, 41]
    assert phred(qstring(q)) == q


def test_n50():
    assert n50([100]) == 100
    assert n50([80, 70, 50, 30, 20]) == 70      # 250 total; 80+70 = 150 >= 125


def test_trim_cuts_bad_tail():
    seq, qual = trim_3prime("A" * 10, qstring([35] * 6 + [5] * 4))
    assert len(seq) < 10 and len(qual) == len(seq)


def test_reference_parses():
    seq, genes = read_gbk(path("phix174.gbk"))
    assert len(seq) == 5386 and len(genes) == 11
    assert all(len(g["seq"]) % 3 == 0 for g in genes)        # whole codons
    assert all(g["seq"][:3] == "ATG" for g in genes)         # all start with ATG


def test_assembler_recovers_random_sequence():
    rng = random.Random(1)
    s = "".join(rng.choice("ACGT") for _ in range(600))
    reads = [s[i:i + 60] for i in range(0, len(s) - 60 + 1, 5)]
    contigs, circ, _, _ = assemble(reads, k=25, min_count=1)
    assert len(contigs) == 1 and contigs[0] in (s, revcomp(s)) and not circ[0]


def test_simulation_is_reproducible_and_honest():
    ref = read_fasta(path("phix174.fasta"))[0][1]
    a, b = sim.simulate(ref, 10), sim.simulate(ref, 10)
    assert a == b
    ref2 = ref + ref[:200]
    claimed = observed = n = 0
    for name, s, q in a:
        _, pos, st = name.split("_")
        frag = ref2[int(pos):int(pos) + len(s)]
        true = frag if st == "+" else revcomp(frag)
        for x, y, qq in zip(s, true, phred(q)):
            claimed += 10 ** (-qq / 10)
            observed += x != y
            n += 1
    assert abs(observed - claimed) / claimed < 0.2     # Phred scores tell the truth


def test_final_assembly_is_the_genome():
    ref = read_fasta(path("phix174.fasta"))[0][1]
    (name, contig), = read_fasta(path("contigs.fasta"))
    assert len(contig) == len(ref)
    assert contig in ref + ref or revcomp(contig) in ref + ref   # same circle


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok  ", name)
