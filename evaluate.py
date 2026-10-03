"""Grade an assembly against the known reference, using k-mers (no aligner).

genome fraction : share of the reference's k-mers found in the assembly
false k-mers    : share of the assembly's k-mers that are NOT in the reference
                  (errors and misjoins show up here)
genes           : a BUSCO-style Complete / Fragmented / Missing call for each
                  annotated gene. Complete = at least 99% of the gene's k-mers
                  lie in ONE contig. Fragmented = at least 50% of them are in
                  the assembly, but spread over several contigs or only partly
                  covered. Missing = anything less.
"""
from ngstools import kmer_set, n50

K = 31


def evaluate(contigs, ref, genes, k=K, circular=None):
    circular = circular or [False] * len(contigs)
    sets = [kmer_set(c + c[:k - 1] if circ else c, k)
            for c, circ in zip(contigs, circular)]
    allk = set().union(*sets) if sets else set()
    refk = kmer_set(ref + ref[:k - 1], k)       # circular genome
    lens = [len(c) for c in contigs]
    res = {
        "contigs": len(contigs),
        "total_len": sum(lens),
        "largest": max(lens) if lens else 0,
        "n50": n50(lens) if lens else 0,
        "genome_fraction": len(refk & allk) / len(refk),
        "false_kmer_fraction": len(allk - refk) / len(allk) if allk else 0.0,
        "genes": [],
    }
    for g in genes:
        gk = kmer_set(g["seq"], k)
        best = max((len(gk & s) / len(gk) for s in sets), default=0.0)
        total = len(gk & allk) / len(gk)
        status = ("Complete" if best >= 0.99
                  else "Fragmented" if total >= 0.5 else "Missing")
        res["genes"].append((g["name"], status, best, total))
    res["complete"] = sum(1 for g in res["genes"] if g[1] == "Complete")
    res["fragmented"] = sum(1 for g in res["genes"] if g[1] == "Fragmented")
    res["missing"] = sum(1 for g in res["genes"] if g[1] == "Missing")
    return res
