"""Step 2: FastQC-style quality control, then trimming.

Computes the same headline numbers FastQC reports (per-base quality, Q20/Q30,
GC content, duplicates) and then does something FastQC cannot: because we
know where every read really came from, it measures the true error rate and
checks that the quality scores told the truth.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ngstools import (path, read_fasta, read_fastq, revcomp, phred, gc,
                      trim_3prime, write_fastq)

ref = read_fasta(path("phix174.fasta"))[0][1]
ref2 = ref + ref[:200]
reads = list(read_fastq(path("reads.fastq")))
RL = max(len(s) for _, s, _ in reads)
lines = []


def out(s=""):
    print(s)
    lines.append(s)


def truth(name, n):
    """True sequence of the first n bases of the read as it was sequenced."""
    _, pos, strand = name.split("_")
    frag = ref2[int(pos):int(pos) + RL]
    return (frag if strand == "+" else revcomp(frag))[:n]


def profile(rs):
    """per-position mean Q, observed error rate, expected error rate."""
    sumq, err, exp, cnt = [0] * RL, [0] * RL, [0.0] * RL, [0] * RL
    for name, s, q in rs:
        t = truth(name, len(s))
        for j, qq in enumerate(phred(q)):
            sumq[j] += qq
            cnt[j] += 1
            exp[j] += 10 ** (-qq / 10)
            err[j] += s[j] != t[j]
    ok = [j for j in range(RL) if cnt[j]]
    return ([sumq[j] / cnt[j] for j in ok], [err[j] / cnt[j] for j in ok],
            [exp[j] / cnt[j] for j in ok])


def summary(label, rs):
    bases = sum(len(s) for _, s, _ in rs)
    qs = [x for _, _, q in rs for x in phred(q)]
    q20 = 100 * sum(x >= 20 for x in qs) / len(qs)
    q30 = 100 * sum(x >= 30 for x in qs) / len(qs)
    errs = sum(a != b for n_, s, _ in rs for a, b in zip(s, truth(n_, len(s))))
    dup = 100 * (1 - len({s for _, s, _ in rs}) / len(rs))
    out(f"{label:<10}{len(rs):>7} reads {bases:>9} bases   Q20 {q20:5.1f}%  Q30 {q30:5.1f}%   "
        f"true error {100*errs/bases:.3f}%   duplicates {dup:.1f}%")


out("QC summary")
summary("raw", reads)

trimmed = []
for name, s, q in reads:
    s2, q2 = trim_3prime(s, q, window=4, min_q=20)
    if len(s2) >= 50:
        trimmed.append((name, s2, q2))
summary("trimmed", trimmed)
out(f"kept {100*len(trimmed)/len(reads):.1f}% of reads, "
    f"{100*sum(len(s) for _,s,_ in trimmed)/sum(len(s) for _,s,_ in reads):.1f}% of bases")
write_fastq(path("reads_trimmed.fastq"), trimmed)

out()
out(f"GC content: genome {100*gc(ref):.1f}%, reads (mean) "
    f"{100*sum(gc(s) for _,s,_ in reads)/len(reads):.1f}%")

mq, oe, ee = profile(reads)
mq2, oe2, ee2 = profile(trimmed)
out()
out("Did the quality scores tell the truth? (error rate, raw reads)")
out(f"  claimed by Phred scores: {100*sum(ee)/len(ee):.3f}%    observed against reference: {100*sum(oe)/len(oe):.3f}%")
out(f"  last 10 positions  claimed {100*sum(ee[-10:])/10:.2f}%  observed {100*sum(oe[-10:])/10:.2f}%")

fig, ax = plt.subplots(1, 2, figsize=(12, 4))
ax[0].plot(mq, label="raw", color="tab:red")
ax[0].plot(mq2, label="trimmed", color="tab:green")
ax[0].axhline(30, ls=":", color="grey"); ax[0].axhline(20, ls=":", color="grey")
ax[0].set_xlabel("position in read"); ax[0].set_ylabel("mean Phred quality")
ax[0].set_title("Per-base quality"); ax[0].legend()
ax[1].semilogy(oe, label="observed (vs true sequence)", color="k")
ax[1].semilogy(ee, label="claimed by Phred score", color="tab:orange", ls="--")
ax[1].set_xlabel("position in read"); ax[1].set_ylabel("error rate")
ax[1].set_title("Are the quality scores honest?"); ax[1].legend()
fig.tight_layout(); fig.savefig("figures/02_read_qc.png", dpi=150)

with open("results/02_read_qc.txt", "w") as fh:
    fh.write("\n".join(lines) + "\n")
