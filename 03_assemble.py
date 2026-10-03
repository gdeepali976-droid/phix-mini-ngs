"""Step 3: assemble the reads with the de Bruijn graph assembler.

Usage: python 03_assemble.py [reads.fastq] [k] [min_count]
Defaults: the trimmed reads, k = 31, min_count = 2.
"""
import sys
from ngstools import path, read_fastq, write_fasta, n50
from assembler import assemble

fq = sys.argv[1] if len(sys.argv) > 1 else path("reads_trimmed.fastq")
k = int(sys.argv[2]) if len(sys.argv) > 2 else 31
min_count = int(sys.argv[3]) if len(sys.argv) > 3 else 2

seqs = [s for _, s, _ in read_fastq(fq)]
contigs, circ, solid, distinct = assemble(seqs, k=k, min_count=min_count)
lens = [len(c) for c in contigs]
lines = [
    f"reads: {len(seqs)}   k = {k}   min k-mer count = {min_count}",
    f"distinct k-mers: {distinct}   kept (seen >= {min_count} times): {solid}",
    f"contigs: {len(contigs)}   total length: {sum(lens)} bp   "
    f"largest: {max(lens)} bp   N50: {n50(lens)} bp   circular contigs: {sum(circ)}",
    "five longest contigs (bp): " + ", ".join(str(x) for x in lens[:5]),
]
print("\n".join(lines))
write_fasta(path("contigs.fasta"),
            [(f"contig_{i+1}_len{len(c)}" + (" circular=true" if f else ""), c)
             for i, (c, f) in enumerate(zip(contigs, circ))])
with open("results/03_assemble.txt", "w") as fh:
    fh.write("\n".join(lines) + "\n")
