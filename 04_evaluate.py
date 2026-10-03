"""Step 4: grade the assembly against the true phiX174 genome.

Usage: python 04_evaluate.py [contigs.fasta]
"""
import sys
from ngstools import path, read_fasta, read_gbk
from evaluate import evaluate, K

fa = sys.argv[1] if len(sys.argv) > 1 else path("contigs.fasta")
ref, genes = read_gbk(path("phix174.gbk"))
records = read_fasta(fa)
contigs = [s for _, s in records]
circ = ["circular=true" in n for n, _ in records]
r = evaluate(contigs, ref, genes, circular=circ)

lines = [
    f"{r['contigs']} contigs, {r['total_len']} bp (reference is {len(ref)} bp), "
    f"largest {r['largest']} bp, N50 {r['n50']} bp",
    f"genome fraction (share of reference {K}-mers recovered): {100*r['genome_fraction']:.2f}%",
    f"false {K}-mers (in assembly, not in reference):        {100*r['false_kmer_fraction']:.2f}%",
    "",
    f"Gene completeness, BUSCO-style: {r['complete']} complete, "
    f"{r['fragmented']} fragmented, {r['missing']} missing (of {len(genes)})",
] + [f"  {n[:40]:<40} {s:<11} best contig {100*b:5.1f}%   assembly {100*t:5.1f}%"
     for n, s, b, t in r["genes"]]
print("\n".join(lines))
with open("results/04_evaluate.txt", "w") as fh:
    fh.write("\n".join(lines) + "\n")
