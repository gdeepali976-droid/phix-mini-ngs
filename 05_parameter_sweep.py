"""Step 5: what actually decides whether the assembly works?

Three knobs, one at a time, everything else held at the defaults
(30x coverage, trimmed reads, k = 31, min k-mer count = 2):

  A. sequencing depth      B. k-mer size      C. error filtering (raw vs
                                              trimmed reads, min count 1/2/3)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ngstools import path, read_gbk, trim_3prime
from assembler import assemble
from evaluate import evaluate
import importlib
sim = importlib.import_module("01_simulate_reads")

ref, genes = read_gbk(path("phix174.gbk"))
lines = []


def out(s=""):
    print(s)
    lines.append(s)


def reads_at(cov, trim=True):
    rs = sim.simulate(ref, cov)
    if not trim:
        return [s for _, s, _ in rs]
    kept = []
    for _, s, q in rs:
        s2, _ = trim_3prime(s, q)
        if len(s2) >= 50:
            kept.append(s2)
    return kept


def run(seqs, k=31, min_count=2):
    contigs, circ, solid, distinct = assemble(seqs, k=k, min_count=min_count)
    r = evaluate(contigs, ref, genes, circular=circ)
    r["solid"], r["distinct"] = solid, distinct
    return r


def row(label, r):
    out(f"{label:>10}  contigs {r['contigs']:>4}  N50 {r['n50']:>5}  "
        f"genome {100*r['genome_fraction']:6.2f}%  false k-mers {100*r['false_kmer_fraction']:5.2f}%  "
        f"genes C/F/M {r['complete']}/{r['fragmented']}/{r['missing']}")


out("A. Sequencing depth (trimmed reads, k = 31, min count 2)")
covs = [3, 5, 8, 10, 15, 20, 30, 50, 100]
A = []
for c in covs:
    r = run(reads_at(c)); A.append(r); row(f"{c}x", r)

out()
out("B. k-mer size (30x, trimmed reads, min count 2)")
ks = [15, 21, 31, 41, 51, 71]
seqs30 = reads_at(30)
B = []
for k in ks:
    r = run(seqs30, k=k); B.append(r); row(f"k={k}", r)

out()
out("C. Error filtering (30x, k = 31)")
C = []
for trim in (False, True):
    sq = reads_at(30, trim=trim)
    for m in (1, 2, 3):
        r = run(sq, min_count=m); C.append((trim, m, r))
        row(f"{'trim' if trim else 'raw'} m={m}", r)
        out(f"{'':>10}  k-mers seen: {r['distinct']}, kept: {r['solid']}")

fig, ax = plt.subplots(1, 3, figsize=(15, 4))
ax[0].plot(covs, [100 * r["genome_fraction"] for r in A], "o-", label="genome recovered %")
ax[0].plot(covs, [100 * r["complete"] / len(genes) for r in A], "s--", label="genes complete %")
ax[0].set_xscale("log"); ax[0].set_xlabel("coverage (x)"); ax[0].set_ylim(0, 105)
ax[0].set_title("A. Depth"); ax[0].legend(fontsize=8)
ax[1].plot(ks, [100 * r["genome_fraction"] for r in B], "o-", label="genome recovered %")
ax[1].plot(ks, [r["n50"] / len(ref) * 100 for r in B], "s--", label="N50 as % of genome")
ax[1].set_xlabel("k"); ax[1].set_ylim(0, 105); ax[1].set_title("B. k-mer size"); ax[1].legend(fontsize=8)
labels = [f"{'trim' if t else 'raw'}\nm={m}" for t, m, _ in C]
ax[2].bar(labels, [r["contigs"] for _, _, r in C], color=["#d9694a"] * 3 + ["#4a90d9"] * 3)
ax[2].set_ylabel("number of contigs (fewer is better)")
ax[2].set_title("C. Error filtering")
fig.tight_layout(); fig.savefig("figures/05_parameter_sweep.png", dpi=150)

with open("results/05_parameter_sweep.txt", "w") as fh:
    fh.write("\n".join(lines) + "\n")
