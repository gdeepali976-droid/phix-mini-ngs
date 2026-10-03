"""Step 0: download the phiX174 reference (5,386 bp) and its gene annotation.

phiX174 is the control genome Illumina spikes into sequencing runs, and it
was the first DNA genome ever sequenced (Sanger, 1977). Small enough to
assemble in seconds, real enough to be honest.
"""
import os
import urllib.request as U
from ngstools import path, read_gbk, write_fasta

URL = ("https://raw.githubusercontent.com/biopython/biopython/master/"
       "Tests/GFF/NC_001422.gbk")
os.makedirs(os.path.dirname(path("x")), exist_ok=True)
if not os.path.exists(path("phix174.gbk")):
    U.urlretrieve(URL, path("phix174.gbk"))

seq, genes = read_gbk(path("phix174.gbk"))
write_fasta(path("phix174.fasta"), [("NC_001422.1 phiX174", seq)])
print(f"reference: {len(seq)} bp, GC {100*(seq.count('G')+seq.count('C'))/len(seq):.1f}%")
print(f"annotated genes: {len(genes)}")
for g in genes:
    print(f"  {g['name'][:42]:<42} {g['strand']} {len(g['seq']):5d} bp  {g['parts']}")
