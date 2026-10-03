"""A small de Bruijn graph assembler (the idea behind ABySS, Velvet, SPAdes).

1. Count every k-mer in every read, on both strands.
2. Throw away k-mers seen fewer than `min_count` times. A sequencing error
   creates up to k brand-new k-mers that almost nobody else shares.
3. Link k-mers that overlap by k-1 bases. Read off the unbranched paths
   ("unitigs"). Those are the contigs.

No bubble popping, no tip clipping, no paired-end scaffolding. That is why
real assemblers exist, and it is also why this one is readable.
"""
from collections import Counter
from ngstools import revcomp, canonical

BASES = "ACGT"


def count_kmers(seqs, k):
    c = Counter()
    for s in seqs:
        for t in (s, revcomp(s)):
            for i in range(len(t) - k + 1):
                c[t[i:i + k]] += 1
    return c


def assemble(seqs, k=31, min_count=2):
    counts = count_kmers(seqs, k)
    solid = {m for m, n in counts.items() if n >= min_count}

    def succ(m):
        return [m[1:] + b for b in BASES if m[1:] + b in solid]

    def pred(m):
        return [b + m[:-1] for b in BASES if b + m[:-1] in solid]

    seen, contigs, circular = set(), [], []

    def walk(start):
        path, m = [start], start
        seen.add(start)
        while True:
            nxt = succ(m)
            if len(nxt) != 1:
                break
            n = nxt[0]
            if n in seen or len(pred(n)) != 1:
                break
            path.append(n)
            seen.add(n)
            m = n
        return path[0] + "".join(p[-1] for p in path[1:])

    # unitig starts: no predecessor, a branch behind, or a fork ahead of it
    for m in solid:
        if m in seen:
            continue
        p = pred(m)
        if len(p) == 1 and len(succ(p[0])) == 1:
            continue                     # middle of a path, reached from its start
        contigs.append(walk(m))
        circular.append(False)
    for m in solid:                      # anything left is a closed circle
        if m not in seen:
            c = walk(m)
            contigs.append(c[:len(c) - (k - 1)])   # drop the k-1 overlap that closes the loop
            circular.append(True)

    # every contig exists once per strand; keep one
    uniq = {}
    for c, circ in zip(contigs, circular):
        uniq[canon_circular(c) if circ else canonical(c)] = circ
    items = sorted(uniq.items(), key=lambda kv: len(kv[0]), reverse=True)
    return [c for c, _ in items], [f for _, f in items], len(solid), len(counts)


def canon_circular(s):
    """Same circle read from any start or either strand -> one fixed string."""
    r = revcomp(s)
    return min(min(s[i:] + s[:i] for i in range(len(s))),
               min(r[i:] + r[:i] for i in range(len(r))))
