# A mini NGS pipeline from scratch: reads, QC, assembly, evaluation

A small end-to-end sequencing workflow written in plain Python (standard library, plus matplotlib for figures). It simulates Illumina-style reads from the phiX174 genome, runs FastQC-style quality control and trimming, assembles the reads with a de Bruijn graph assembler I wrote, and grades the result against the known genome.

I built it to understand what FastQC, ABySS and BUSCO actually compute, since those were the tools I used during my internship at CDFD. phiX174 (5,386 bp, 11 genes) is the control genome Illumina spikes into sequencing runs. It was also the first DNA genome ever sequenced.

## Pipeline

| Step | Script | What it does |
|------|--------|--------------|
| 0 | `00_fetch_reference.py` | Downloads the phiX174 sequence and its 11 annotated genes |
| 1 | `01_simulate_reads.py` | Simulates 100 bp reads from both strands of the circular genome, with a Phred quality profile and errors that match it |
| 2 | `02_read_qc.py` | FastQC-style metrics, sliding-window trimming, and a check that the quality scores tell the truth |
| 3 | `03_assemble.py` | De Bruijn graph assembly (`assembler.py`) |
| 4 | `04_evaluate.py` | Genome recovery, false k-mers, and a BUSCO-style Complete/Fragmented/Missing call per gene (`evaluate.py`) |
| 5 | `05_parameter_sweep.py` | Which settings decide whether the assembly works |

## Run it

```
pip install matplotlib
python run_all.py
```

It downloads about 23 KB, runs all six steps and the tests, and takes under a minute. Output goes to `results/` and `figures/`.

## Results (30x coverage, seed 42)

**QC.** 95.9% of bases are Q20 and 81.9% are Q30. Because I know where each read came from, I can measure the real error rate: 0.300%, against 0.289% claimed by the quality scores, so the scores are honest. Sliding-window trimming (window 4, mean Q20) drops the real error rate to 0.065% while keeping 95.8% of the bases.

**Assembly.** From the trimmed reads the assembler produces one circular contig of 5,386 bp. It recovers 100% of the reference 31-mers with 0 false k-mers, and all 11 genes are complete.

**What decides whether it works** (`figures/05_parameter_sweep.png`):

- **Depth.** 15x is enough for one complete contig. At 10x there are 7 contigs and only 8 of 11 genes are complete, and at 5x no gene is complete. At 100x it gets worse again, with 5 contigs and 9 of 11 genes complete, because sequencing errors start repeating often enough to pass the min-count filter. A fixed filter does not scale with depth.
- **k-mer size.** k = 15 to 51 all give the same perfect result. At k = 71 the genome splits in two, because short trimmed reads contribute too few k-mers that long.
- **Error filtering.** Keeping every k-mer (min count 1) on raw reads gives 588 contigs and 46% false k-mers. Raw reads need min count 3 to assemble cleanly, but trimmed reads only need 2. Trimming and filtering do the same job from two sides.

## How the assembler works

1. Count every k-mer in every read, on both strands.
2. Drop k-mers seen fewer than `min_count` times. One sequencing error creates up to k new k-mers that nobody else shares.
3. Link k-mers that overlap by k−1 bases and read off the unbranched paths (unitigs). Closed loops are reported as circular contigs.

## Limitations

- **The reads are simulated.** The error model is substitutions only. There are no indels, no GC bias, no adapters, no paired ends and no contamination. Real data is messier, and this project does not claim otherwise.
- **phiX174 is an easy genome.** It is small and has few repeats. The assembler has no bubble popping, tip clipping or scaffolding, so it would fail on a repeat-rich genome.
- **N50 is trivial here.** One circular contig means N50 equals the genome length. It matters in the depth sweep, where the assembly is fragmented.
- **Evaluation uses k-mers, not alignment.** That is enough to say what was recovered and what was invented, but it cannot locate individual mismatches. The circular-contig deduplication is O(n²), which is fine for 5 kb and not for real genomes.

## Next steps

- Swap the simulated reads for real phiX174 reads from a public run and compare the QC numbers with FastQC.
- Add paired-end reads and indel errors to the simulator.
- Assemble the same reads with ABySS or SPAdes and compare against my assembler.
- Test on a genome with repeats.
