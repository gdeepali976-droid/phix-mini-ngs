"""Helpers shared by every script: sequences, GenBank parsing, FASTQ, k-mers.

Standard library only. The point of the repo is to see what FastQC and an
assembler actually compute, so nothing here is hidden inside a package.
"""
import os
import re
from collections import Counter

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
COMP = str.maketrans("ACGT", "TGCA")


def path(name):
    return os.path.join(DATA, name)


def revcomp(s):
    return s.translate(COMP)[::-1]


def canonical(s):
    r = revcomp(s)
    return s if s < r else r


def gc(s):
    return (s.count("G") + s.count("C")) / len(s)


# ---------------------------------------------------------------- GenBank
def read_gbk(fname):
    """Return (sequence, genes). genes = list of dicts with name, seq, strand,
    parts. Handles complement(...) and join(...), including genes that wrap
    around the end of the circular genome (phiX174 has one)."""
    text = open(fname).read()
    seq = "".join(re.findall(
        r"[acgt]+", text.split("\nORIGIN")[1].split("//")[0])).upper()
    genes = []
    feats = text.split("\nFEATURES")[1].split("\nORIGIN")[0]
    for block in re.split(r"\n {5}(?=\S)", feats):
        if not block.startswith("CDS"):
            continue
        loc_txt = re.split(r"\n {21}/", block)[0]
        loc = "".join(loc_txt.split())[3:]
        strand = "-" if loc.startswith("complement") else "+"
        parts = [(int(a), int(b)) for a, b in
                 re.findall(r"(\d+)\.\.>?(\d+)", loc)]
        s = "".join(seq[a - 1:b] for a, b in parts)
        if strand == "-":
            s = revcomp(s)
        m = re.search(r'/product="([^"]+)"', block)
        genes.append({"name": m.group(1) if m else f"cds{len(genes)+1}",
                      "seq": s, "strand": strand, "parts": parts})
    return seq, genes


def write_fasta(fname, records):
    with open(fname, "w") as fh:
        for name, s in records:
            fh.write(f">{name}\n")
            for i in range(0, len(s), 70):
                fh.write(s[i:i + 70] + "\n")


def read_fasta(fname):
    out, name, buf = [], None, []
    for ln in open(fname):
        ln = ln.strip()
        if ln.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = ln[1:], []
        elif ln:
            buf.append(ln)
    if name is not None:
        out.append((name, "".join(buf)))
    return out


# ------------------------------------------------------------------ FASTQ
def phred(qstr):
    """Phred+33 string -> list of ints."""
    return [ord(c) - 33 for c in qstr]


def qstring(qs):
    return "".join(chr(q + 33) for q in qs)


def write_fastq(fname, reads):
    with open(fname, "w") as fh:
        for name, s, q in reads:
            fh.write(f"@{name}\n{s}\n+\n{q}\n")


def read_fastq(fname):
    with open(fname) as fh:
        while True:
            h = fh.readline().rstrip()
            if not h:
                return
            s = fh.readline().rstrip()
            fh.readline()
            q = fh.readline().rstrip()
            yield h[1:], s, q


def trim_3prime(seq, qual, window=4, min_q=20):
    """Sliding-window trim from the 5' end: cut at the first window whose mean
    quality falls below min_q (the idea behind Trimmomatic SLIDINGWINDOW)."""
    qs = phred(qual)
    for i in range(len(qs) - window + 1):
        if sum(qs[i:i + window]) / window < min_q:
            return seq[:i], qual[:i]
    return seq, qual


# --------------------------------------------------------------- k-mers
def kmer_set(s, k):
    return {canonical(s[i:i + k]) for i in range(len(s) - k + 1)}


def n50(lengths):
    """Length L such that contigs of length >= L hold half the assembly."""
    total, run = sum(lengths), 0
    for L in sorted(lengths, reverse=True):
        run += L
        if run * 2 >= total:
            return L
    return 0
