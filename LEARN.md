# LEARN — Bioinformatics for people who already know Python

You know Python. You do not know biology. This file is the shortest path
from there to being useful.

Twenty points. Each one tells you **what it is**, **why you will hit it**,
and gives a snippet you can paste and run.

Companion files: [GLOSSARY.md](GLOSSARY.md) for the vocabulary that blocks
people, and [PROTEINS.md](PROTEINS.md) for protein and structure work.

## How to run the snippets

```bash
git clone https://github.com/hadeelsharaf/Bioinformatic-docs.git
cd Bioinformatic-docs
python -m venv .venv
.venv/Scripts/activate        # Windows
# source .venv/bin/activate   # macOS / Linux
pip install -r requirements.txt
```

Run every snippet from the repo root. They read the tiny files in `data/`.

Every snippet below was run on **Python 3.14.3 / Biopython 1.88 / pandas 3.0.5**.
The outputs shown are the real outputs. Two snippets are marked *not run* and
say why.

## Contents

**Part 1 — the biology you must hold in your head**
1. [The central dogma in 60 seconds](#1-the-central-dogma-in-60-seconds)
2. [The alphabet is not just ACGT](#2-the-alphabet-is-not-just-acgt)
3. [Strands: reverse complement is not `[::-1]`](#3-strands-reverse-complement-is-not--1)
4. [Coordinates: the off-by-one that eats everyone](#4-coordinates-the-off-by-one-that-eats-everyone)

**Part 2 — the file formats**
5. [FASTA](#5-fasta)
6. [FASTQ and Phred quality scores](#6-fastq-and-phred-quality-scores)
7. [Everything is gzipped: stream, never `read()`](#7-everything-is-gzipped-stream-never-read)
8. [SAM / BAM: reads glued to a genome](#8-sam--bam-reads-glued-to-a-genome)
9. [VCF, BED and GFF: variants and annotations](#9-vcf-bed-and-gff-variants-and-annotations)

**Part 3 — the core operations**
10. [Translation and codon tables](#10-translation-and-codon-tables)
11. [GC content (and a removed function that breaks old tutorials)](#11-gc-content-and-a-removed-function-that-breaks-old-tutorials)
12. [Open reading frames in all six frames](#12-open-reading-frames-in-all-six-frames)
13. [Pairwise alignment](#13-pairwise-alignment)
14. [k-mers](#14-k-mers)

**Part 4 — doing real work**
15. [Get real data: NCBI Entrez](#15-get-real-data-ncbi-entrez)
16. [BLAST: search by similarity](#16-blast-search-by-similarity)
17. [pandas and matplotlib on sequences](#17-pandas-and-matplotlib-on-sequences)
18. [Speed: pure Python loops die on a genome](#18-speed-pure-python-loops-die-on-a-genome)
19. [Reproducibility: bioconda and workflow managers](#19-reproducibility-bioconda-and-workflow-managers)
20. [Ethics and licensing](#20-ethics-and-licensing)

**Part 5 — the NGS pipeline, stage by stage**
21. [Quality control](#21-quality-control)
22. [Trimming](#22-trimming)
23. [Mapping](#23-mapping)
24. [Duplicate marking](#24-duplicate-marking)
25. [Variant calling](#25-variant-calling)
26. [Coverage](#26-coverage)
27. [RNA-seq counting](#27-rna-seq-counting)
28. [Differential expression](#28-differential-expression)

**Also here**
- [Beyond Biopython: scikit-bio](#beyond-biopython-scikit-bio)

---

# Part 1 — the biology you must hold in your head

## 1. The central dogma in 60 seconds

DNA is a long text. A cell copies part of it into RNA (**transcription**),
then reads the RNA three letters at a time to build a protein
(**translation**). Three letters = one **codon** = one amino acid.

Why you care: almost every function in this field is one of those two steps,
or a search over the text they act on.

```python
from Bio.Seq import Seq

dna = Seq("ATGGCCATTGTAATGGGCCGCTGA")
rna = dna.transcribe()      # T becomes U
protein = rna.translate()   # 3 bases become 1 amino acid
print(dna, rna, protein, sep="\n")
```

Out:
```
ATGGCCATTGTAATGGGCCGCTGA
AUGGCCAUUGUAAUGGGCCGCUGA
MAIVMGR*
```

`*` is a stop codon. It means "protein ends here".

Link: [Khan Academy AP Biology](https://www.khanacademy.org/science/ap-biology) (units 6 and 7)

## 2. The alphabet is not just ACGT

Real sequence files contain **IUPAC ambiguity codes**. `N` means "any base".
`R` means "A or G". `Y` means "C or T". They appear when the machine was not
sure.

Why you care: your `if base in "ACGT"` check will silently drop real data, and
your GC calculation will quietly change meaning depending on how you treat them.

```python
from Bio.SeqUtils import gc_fraction

seq = "ACGTNNNRYACGT"
print(sorted(set(seq) - set("ACGT")))        # what is actually in there
print(gc_fraction(seq))                      # default: ambiguous bases removed
print(gc_fraction(seq, ambiguous="ignore"))  # ambiguous bases kept in the total
```

Out:
```
['N', 'R', 'Y']
0.5
0.3076923076923077
```

Same sequence, two different answers. Know which one you meant.

Link: [IUPAC codes](https://www.bioinformatics.org/sms/iupac.html)

## 3. Strands: reverse complement is not `[::-1]`

DNA is two paired strands running in opposite directions. Sequence is always
written 5'→3' (say "five prime to three prime"). To read the other strand you
must reverse **and** swap each base for its pair: A↔T, C↔G.

Why you care: a gene can sit on either strand. Reverse only, and you get
nonsense.

```python
from Bio.Seq import Seq

dna = Seq("ATGGCC")
print(dna[::-1])                  # reversed only  -> WRONG
print(dna.reverse_complement())   # reversed and complemented -> right
```

Out:
```
CCGGTA
GGCCAT
```

Link: [Biopython Seq objects](https://biopython.org/docs/latest/Tutorial/chapter_seq_objects.html)

## 4. Coordinates: the off-by-one that eats everyone

Two coordinate systems are in daily use and they disagree:

| System | Used by | Start counts from | End |
|---|---|---|---|
| 0-based, half-open | BED, Python slices, pysam | 0 | excluded |
| 1-based, inclusive | GFF, GTF, VCF, SAM text, Ensembl | 1 | included |

Why you care: this is the single most common bug in the field. The same three
bases have two different addresses.

```python
seq = "ACGTACGTAC"

# BED / Python: 0-based, end excluded
print(seq[2:5])            # BED line would read: chr1  2  5

# GFF / VCF: 1-based, end included -> the SAME bases are 3..5
start, end = 3, 5
print(seq[start - 1:end])
```

Out:
```
GTA
GTA
```

Rule of thumb: when a tool gives you a region, write down which system it uses
before you slice anything.

Link: [UCSC format FAQ](https://genome.ucsc.edu/FAQ/FAQformat.html)

---

# Part 2 — the file formats

## 5. FASTA

A header line starting with `>`, then the sequence on the lines below.
That is the whole format. Extensions: `.fasta`, `.fa`, `.fna` (nucleotides),
`.faa` (amino acids).

Why you care: it is the default way a sequence travels between tools.

```python
from Bio import SeqIO

for rec in SeqIO.parse("data/sample.fasta", "fasta"):
    print(rec.id, len(rec.seq))
```

Out:
```
seq1 39
seq2 36
seq3 23
```

`SeqIO.parse` is a generator, so a 3 GB file costs you one record of memory.
Use `SeqIO.read` only when you know the file holds exactly one record.

Link: [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) · [Bio.SeqIO API](https://biopython.org/docs/latest/api/Bio.SeqIO.html)

## 6. FASTQ and Phred quality scores

FASTQ is FASTA plus a confidence score for every single base. Four lines per
record: id, sequence, `+`, quality. The quality line is the same length as the
sequence line — always check this if you ever write one.

The score is **Phred**: `Q = -10 * log10(P_error)`. So Q20 = 1 wrong base in
100. Q30 = 1 in 1000. Q40 = 1 in 10,000. It is stored as one ASCII character
per base (`chr(Q + 33)`).

Why you care: raw sequencer output is FASTQ, and filtering on quality is the
first step of nearly every pipeline.

```python
from Bio import SeqIO

for rec in SeqIO.parse("data/sample.fastq", "fastq"):
    q = rec.letter_annotations["phred_quality"]   # list of ints
    print(rec.id, "mean Q =", round(sum(q) / len(q), 1), "min Q =", min(q))
```

Out:
```
read1 mean Q = 27.0 min Q = 2
read2 mean Q = 40.0 min Q = 40
```

`read1` ends badly — quality falling off at the end of a read is normal and is
why trimming tools exist.

Link: [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format)

## 7. Everything is gzipped: stream, never `read()`

Real files arrive as `.fastq.gz` and are gigabytes. Open them in text mode
(`"rt"`) and hand the handle straight to the parser.

Why you care: `open(f).read()` on a real FASTQ will kill your process. There is
no polite failure.

```python
import gzip
from Bio import SeqIO

with gzip.open("data/sample.fastq.gz", "rt") as handle:
    for rec in SeqIO.parse(handle, "fastq"):
        print(rec.id, len(rec.seq))
```

Out:
```
read1 36
read2 36
```

Note `"rt"`, not `"r"`. Without it you get bytes and the parser fails.

Link: [Biopython SeqIO tutorial](https://biopython.org/docs/latest/Tutorial/chapter_seq_annot.html)

## 8. SAM / BAM: reads glued to a genome

Once short reads are **aligned** to a reference genome, the result is SAM
(text) or BAM (compressed binary). Each line says: this read sits at this
position, with this confidence (MAPQ), with these differences (CIGAR).

Why you care: this is where the actual biology gets read out. `pysam` is the
standard Python door into it.

```python
# NOT RUN HERE: pysam has no Windows wheel. Use Linux, macOS, or WSL.
# pip install pysam
import pysam

with pysam.AlignmentFile("aln.bam", "rb") as bam:
    for read in bam.fetch("chr1", 1000, 2000):   # 0-based, half-open
        print(read.query_name, read.reference_start, read.mapping_quality)
```

Two traps: `fetch` needs an index file (`aln.bam.bai`) next to the BAM, and
`reference_start` is 0-based even though the SAM text file is 1-based.

Link: [SAM/BAM spec](https://samtools.github.io/hts-specs/) · [pysam API](https://pysam.readthedocs.io/en/latest/api.html)

## 9. VCF, BED and GFF: variants and annotations

- **VCF** — how one sample differs from the reference. One line per variant.
- **BED** — plain regions. 3 required columns: chrom, start, end. 0-based.
- **GFF / GTF** — rich annotations (this region is an exon of that gene). 1-based.

Why you care: BED and GFF are just TSV, so pandas handles them. VCF is not —
it has a header block and typed INFO fields, so use a real parser.

```python
import pandas as pd

bed = pd.read_csv("data/regions.bed", sep="\t", header=None,
                  names=["chrom", "start", "end", "name"])
print(bed)
print(bed["end"].sub(bed["start"]).sum(), "bases covered")
```

Out:
```
  chrom  start  end      name
0  chr1    100  250  promoter
1  chr1    500  900     exon1
2  chr2     50  120     exon2
620 bases covered
```

For VCF use `pysam.VariantFile` or `cyvcf2`. Do not write your own VCF parser —
the spec has more corners than you expect.

Link: [hts-specs (VCF, SAM, BED)](https://samtools.github.io/hts-specs/) · [UCSC format FAQ](https://genome.ucsc.edu/FAQ/FAQformat.html)

---

# Part 3 — the core operations

## 10. Translation and codon tables

64 codons map to 20 amino acids plus "stop". That mapping is a **codon table**,
and there is more than one. Mitochondria use a different one.

Why you care: translate mitochondrial DNA with the standard table and you get
a protein that stops in the middle and is wrong.

```python
from Bio.Seq import Seq
from Bio.Data import CodonTable

t = CodonTable.unambiguous_dna_by_name["Standard"]
print(t.stop_codons)
print(t.forward_table["ATG"])                                  # start codon
print(Seq("ATGGCCTGA").translate())                            # standard
print(Seq("ATGGCCTGA").translate(table="Vertebrate Mitochondrial"))
```

Out:
```
['TAA', 'TAG', 'TGA']
M
MA*
MAW
```

Same DNA. `TGA` is a stop codon in the standard table and the amino acid
tryptophan (`W`) in the mitochondrial one.

Link: [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

## 11. GC content (and a removed function that breaks old tutorials)

GC content is the share of bases that are G or C. High GC means the DNA melts
at a higher temperature, and protein-coding regions tend to sit in GC-rich
areas.

Why you care beyond biology: **`Bio.SeqUtils.GC()` was removed.** Half the
tutorials online still use it. Use `gc_fraction`, and remember it returns a
fraction (0–1), not a percentage.

```python
from Bio.SeqUtils import gc_fraction

print(round(gc_fraction("GCGCGCGCGGCGGCGCGCGCGGCGCGCGCGGCGCGC") * 100, 1))
# from Bio.SeqUtils import GC   # <- AttributeError on modern Biopython
```

Out:
```
100.0
```

Link: [Bio.SeqUtils API](https://biopython.org/docs/latest/api/Bio.SeqUtils.html)

## 12. Open reading frames in all six frames

An **ORF** is a stretch that could code for a protein: starts at `M`, runs to a
stop. But you do not know where the gene starts, so you must try 3 offsets —
and the gene may be on the other strand, so 3 more. Six frames total.

Why you care: the one-frame regex version in most tutorials misses most genes.

```python
from Bio import SeqIO

def find_orfs(seq, min_aa=5):
    """Yield (strand, frame, protein) for ORFs in all six reading frames."""
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            usable = len(nuc) - frame - (len(nuc) - frame) % 3   # keep it a multiple of 3
            protein = nuc[frame:frame + usable].translate()
            for piece in protein.split("*"):
                if "M" in piece:
                    orf = piece[piece.index("M"):]
                    if len(orf) >= min_aa:
                        yield strand, frame, str(orf)

record = next(SeqIO.parse("data/sample.fasta", "fasta"))
for strand, frame, orf in find_orfs(record.seq, min_aa=5):
    print(strand, frame, orf)
```

Out:
```
1 0 MAIVMGR
```

Raise `min_aa` to cut noise: short "ORFs" appear by chance in any random text.

Link: [Biopython tutorial](https://biopython.org/docs/latest/Tutorial/index.html)

## 13. Pairwise alignment

Line two sequences up to see how similar they are, allowing mismatches and gaps.
You set the scores; the algorithm finds the best-scoring arrangement.

Why you care: similarity is how function gets guessed. It is also the engine
under BLAST.

Modern API is `Bio.Align.PairwiseAligner`. **`Bio.pairwise2` is deprecated** —
if a tutorial shows it, the tutorial is old.

```python
from Bio import Align

aligner = Align.PairwiseAligner(match_score=1, mismatch_score=-1,
                                open_gap_score=-2, extend_gap_score=-0.5)
alignment = aligner.align("ACGTACGT", "ACGACGT")[0]
print(alignment)
print("score:", alignment.score)
```

Out:
```
target            0 ACGTACGT 8
                  0 |||-|||| 8
query             0 ACG-ACGT 7

score: 5.0
```

`aligner.mode` is `"global"` by default (align end to end). Set it to
`"local"` to find the best matching stretch inside longer sequences.

Link: [Biopython pairwise alignment](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html)

## 14. k-mers

Chop a sequence into every overlapping substring of length k. That is it.
`ACGT` with k=3 gives `ACG`, `CGT`.

Why you care: genome assembly, read classification, and fast similarity search
are all built on counting k-mers instead of aligning. It turns a hard biology
problem into a dictionary problem, which is exactly what Python is good at.

```python
from collections import Counter

def kmers(seq, k):
    return Counter(seq[i:i + k] for i in range(len(seq) - k + 1))

print(kmers("ACGTACGTAC", 3).most_common(3))
```

Out:
```
[('ACG', 2), ('CGT', 2), ('GTA', 2)]
```

Watch the memory: a human genome at k=31 has billions of distinct k-mers. Real
tools use Bloom filters or disk-backed counts, not a `dict`.

Link: [k-mer](https://en.wikipedia.org/wiki/K-mer)

---

# Part 4 — doing real work

## 15. Get real data: NCBI Entrez

NCBI holds most public sequence data. `Bio.Entrez` is the Python client.
You must set `Entrez.email` — it is not optional politeness, NCBI uses it to
contact you before blocking you.

Why you care: toy strings teach syntax. Real genomes teach the job.

```python
from Bio import Entrez, SeqIO

Entrez.email = "you@example.com"          # put YOUR address here
with Entrez.efetch(db="nucleotide", id="NC_045512.2",
                   rettype="fasta", retmode="text") as handle:
    rec = SeqIO.read(handle, "fasta")

print(rec.id, len(rec.seq), "bp")
print(rec.description)
```

Out:
```
NC_045512.2 29903 bp
NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
```

Rules: 3 requests per second without an API key, 10 with one. Cache what you
download to a local file — do not re-fetch the same record in a loop.

Link: [Biopython Entrez chapter](https://biopython.org/docs/latest/Tutorial/chapter_entrez.html) · [NCBI usage policy](https://www.ncbi.nlm.nih.gov/books/NBK25497/)

## 16. BLAST: search by similarity

BLAST answers "what known sequences look like mine?" It is the most used tool
in biology. You can call the NCBI web service, or install BLAST+ and run it
locally against your own database.

Why you care: web BLAST takes minutes per query and is rate limited. If you are
looping, you want it local.

```python
# NOT RUN HERE: a single web BLAST takes minutes. Run it once, save the result.
from Bio import Blast                      # Biopython >= 1.85

result = Blast.qblast("blastn", "nt", "ACGTACGTACGTACGTACGTACGTACGT")
records = Blast.parse(result)
for hit in next(records)[:5]:
    print(hit.target.id, hit[0].annotations["evalue"])
```

Read the **E-value**, not the percent identity: it is the number of hits this
good you would expect by chance. Lower is better. `1e-50` is strong, `2.0` is
noise.

Link: [Biopython BLAST chapter](https://biopython.org/docs/latest/Tutorial/chapter_blast.html) · [BLAST docs](https://blast.ncbi.nlm.nih.gov/doc/blast-help/)

## 17. pandas and matplotlib on sequences

Once you have parsed sequences, the job becomes ordinary data work: build a
DataFrame, group it, plot it.

Why you care: this is where your existing Python skill pays off immediately,
with no new biology to learn.

```python
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

df = pd.DataFrame(
    [{"id": r.id, "length": len(r.seq), "gc": round(gc_fraction(r.seq), 3)}
     for r in SeqIO.parse("data/sample.fasta", "fasta")]
)
print(df)
# df.plot.bar(x="id", y="gc"); import matplotlib.pyplot as plt; plt.savefig("gc.png")
```

Out:
```
     id  length     gc
0  seq1      39  0.564
1  seq2      36  1.000
2  seq3      23  0.500
```

## 18. Speed: pure Python loops die on a genome

A human genome is ~3.1 billion bases. A per-base Python loop over it takes
hours. Push the loop into C.

Why you care: the difference is not 2x, it is 100x, and it decides whether your
script finishes today.

```python
import numpy as np

big = "ACGT" * 250_000                          # 1 million bases
arr = np.frombuffer(big.encode(), dtype="S1")   # no copy, one byte per base
print(int((arr == b"G").sum() + (arr == b"C").sum()))
```

Out:
```
500000
```

Also worth knowing: `pyfastx` for random access into a huge FASTA without
loading it, and `biotite` when you want numpy-native structures instead of
Biopython objects.

Link: [pyfastx](https://pyfastx.readthedocs.io/en/latest/)

## 19. Reproducibility: bioconda and workflow managers

Most bioinformatics tools are C binaries, not pip packages. `samtools`, `bwa`,
`bedtools` — `pip install` will not get them. **Bioconda** will.

Why you care: "it worked on my laptop" is the default state of this field, and
a result nobody can reproduce is not a result.

```bash
# NOT RUN HERE: shell, and it needs conda installed.
conda create -n bio -c conda-forge -c bioconda python=3.12 biopython samtools
conda activate bio
conda env export --no-builds > environment.yml   # commit this file
```

When your analysis grows past one script, move it into **Snakemake** or
**Nextflow**. They re-run only the steps whose inputs changed, and they record
what ran.

Link: [Bioconda](https://bioconda.github.io/) · [Snakemake](https://snakemake.readthedocs.io/en/stable/)

## 20. Ethics and licensing

This point has a checklist instead of a snippet, because the mistakes here are
not coding mistakes.

- **Human genomic data is identifying.** A genome is not anonymous data, and it
  cannot be de-identified the way a name or an address can. Treat it as personal
  data even when a file has no name in it.
- **Controlled access is real.** Data from dbGaP, EGA and most clinical sources
  needs an approved application, and the agreement usually forbids re-sharing,
  re-uploading to a cloud service, or committing it to git.
- **Check the license on reference data too.** Public does not mean
  unrestricted, and some databases forbid commercial use.
- **Never commit real patient data to a repository.** Git history keeps it after
  you delete the file. Put a `data/` rule in `.gitignore` before your first
  commit, not after.
- **Say what you did.** Record tool versions and parameters. An unreproducible
  clinical claim can do real harm.

Link: [GA4GH](https://www.ga4gh.org/) · [dbGaP](https://www.ncbi.nlm.nih.gov/gap/)

---

---

# Part 5 — the NGS pipeline, stage by stage

These eight stages sit between a sequencing machine and an answer. Knowing the
names tells you what a paper's methods section is talking about.

Be warned: **about half of this work is not Python.** The heavy stages are C
programs you call from a shell. Where that is true below, the block is marked
`SHELL, NOT RUN` and you get the Python that surrounds it instead. Dressing up a
`bwa` command as Python would teach you nothing.

## 21. Quality control

Before anything else, look at your reads. How long are they, how good are the
scores, does quality collapse at the end. The standard tool is **FastQC**, which
produces a report, not data.

Why you care: bad input silently produces confident, wrong output. QC is the
cheapest bug prevention in the field.

```python
from Bio import SeqIO

for rec in SeqIO.parse("data/sample.fastq", "fastq"):
    q = rec.letter_annotations["phred_quality"]
    bad = sum(1 for x in q if x < 20)
    print(f"{rec.id:6} mean {sum(q) / len(q):5.1f}  bases below Q20: {bad}/{len(q)}")
```

Out:
```
read1  mean  27.0  bases below Q20: 6/36
read2  mean  40.0  bases below Q20: 0/36
```

Link: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## 22. Trimming

Cut the bad parts off. Usually that means the end of the read, where quality
falls, plus any leftover adapter — the artificial sequence the machine attaches
to every fragment. Common tools: **fastp**, **Trimmomatic**, **cutadapt**.

Why you care: adapter left in a read will not map, or worse, will map to the
wrong place. Trimming is the difference between usable and misleading data.

```python
from Bio import SeqIO

def trim_tail(rec, cutoff=20):
    """Chop bases off the end while their quality is below the cutoff."""
    q = rec.letter_annotations["phred_quality"]
    end = len(q)
    while end > 0 and q[end - 1] < cutoff:
        end -= 1
    return rec[:end]

for rec in SeqIO.parse("data/sample.fastq", "fastq"):
    print(f"{rec.id:6} {len(rec.seq)} -> {len(trim_tail(rec).seq)} bases")
```

Out:
```
read1  36 -> 30 bases
read2  36 -> 36 bases
```

Slicing a `SeqRecord` slices its quality scores too. That is why `rec[:end]`
just works and you never handle the two lists separately.

Link: [fastp](https://github.com/OpenGene/fastp)

## 23. Mapping

Find where each read came from in the reference genome. This is the expensive
step. Tools: **bwa-mem2** or **minimap2** for DNA, **STAR** or **HISAT2** for
RNA. Output is SAM or BAM.

Why you care: everything downstream is coordinates, and coordinates come from
here. Mapping quality (MAPQ) tells you how sure the aligner is — a read that
fits equally well in three places gets MAPQ 0, and you usually throw it away.

```bash
# SHELL, NOT RUN: these are C programs, not Python packages.
bwa-mem2 index reference.fasta                       # once per reference
bwa-mem2 mem reference.fasta reads.fastq.gz \
  | samtools sort -o aligned.bam
samtools index aligned.bam                           # makes aligned.bam.bai
```

The Python side of this stage is reading the result — see point 8 for `pysam`.

Link: [minimap2](https://github.com/lh3/minimap2) · [samtools](https://www.htslib.org/)

## 24. Duplicate marking

The lab step that copies DNA can copy the same fragment many times. Those copies
are not independent evidence, so they get flagged and ignored. Tools:
**samtools markdup**, **Picard MarkDuplicates**.

Why you care: forget this and one lucky fragment looks like twenty supporting
reads, which turns noise into a confident false variant.

```python
from collections import Counter
from Bio import SeqIO

seqs = [str(r.seq) for r in SeqIO.parse("data/sample.fastq", "fastq")]
seqs.append(seqs[0])                       # pretend read1 got copied

counts = Counter(seqs)
print("total:", len(seqs), " unique:", len(counts))
print("duplicate rate:", f"{1 - len(counts) / len(seqs):.0%}")
```

Out:
```
total: 3  unique: 2
duplicate rate: 33%
```

Real tools compare mapped positions, not sequence text, because two genuinely
different fragments can share a sequence by chance.

Link: [samtools markdup](https://www.htslib.org/doc/samtools-markdup.html)

## 25. Variant calling

Compare the piled-up reads to the reference and decide, position by position,
whether a difference is real. Tools: **bcftools**, **GATK**, **DeepVariant**.
Output is VCF.

Why you care: this is the stage that produces the actual finding. The Python
side is reading and filtering the VCF that comes out.

```bash
# SHELL, NOT RUN
bcftools mpileup -f reference.fasta aligned.bam \
  | bcftools call -mv -Ob -o calls.bcf
```

A VCF data line is tab separated, and the fixed columns are always in this order:

```python
line = "chr1\t150\trs123\tA\tG\t60\tPASS\tDP=32"
chrom, pos, vid, ref, alt, qual, filt, info = line.split("\t")

fields = dict(kv.split("=") for kv in info.split(";"))
print(f"{chrom}:{pos} {ref}>{alt}  qual={qual}  depth={fields['DP']}")
```

Out:
```
chr1:150 A>G  qual=60  depth=32
```

Do that for one line to learn the shape, then use `pysam.VariantFile` for real
files — the header, the INFO types and the sample columns get complicated fast.

Link: [bcftools](https://samtools.github.io/bcftools/)

## 26. Coverage

How many reads sit over each position. Check it before you trust anything: a
region with no reads produces no variants, which looks identical to a region
with nothing wrong.

Why you care: "we found no mutation" and "we could not see that region" are
completely different claims, and only coverage tells them apart.

```python
from collections import Counter

reads = [(100, 136), (110, 146), (120, 156), (500, 536)]   # start, end
depth = Counter(pos for start, end in reads for pos in range(start, end))

print("bases covered:", len(depth))
print("max depth:", max(depth.values()))
print("mean depth:", round(sum(depth.values()) / len(depth), 2))
print("depth at 125:", depth[125], "| at 505:", depth[505], "| at 400:", depth[400])
```

Out:
```
bases covered: 92
max depth: 3
mean depth: 1.57
depth at 125: 3 | at 505: 1 | at 400: 0
```

For a whole genome, use `samtools depth` or `mosdepth` — a Python loop over
3 billion positions is the wrong tool.

Link: [mosdepth](https://github.com/brentp/mosdepth)

## 27. RNA-seq counting

For RNA work the question is not "what changed" but "how much of each gene is
present". So you count how many reads land in each gene. Tools:
**featureCounts**, **HTSeq**, or **salmon**, which skips full mapping.

Why you care: the output is a plain table of genes by samples. Once you have it,
the work becomes pandas.

```python
import pandas as pd

genes = pd.read_csv("data/regions.bed", sep="\t", header=None,
                    names=["chrom", "start", "end", "gene"])
reads = [("chr1", 120), ("chr1", 130), ("chr1", 600), ("chr2", 60), ("chr1", 999)]

counts = {g.gene: 0 for g in genes.itertuples()}
for chrom, pos in reads:
    for g in genes.itertuples():
        if g.chrom == chrom and g.start <= pos < g.end:
            counts[g.gene] += 1

print(counts)
print("reads assigned:", sum(counts.values()), "of", len(reads))
```

Out:
```
{'promoter': 2, 'exon1': 1, 'exon2': 1}
reads assigned: 4 of 5
```

One read landed at chr1:999, outside every gene, so it counts for nothing. That
gap is normal and worth reporting — a low assignment rate means something is
wrong with your annotation or your mapping.

Link: [HTSeq](https://htseq.readthedocs.io/en/master/index.html)

## 28. Differential expression

Given counts for two groups, which genes really changed? Not simply the ones
with the biggest ratio — a gene with 2 reads against 8 is 4x higher and means
nothing. You need the statistics. Tools: **DESeq2** and **edgeR** (both R), or
**PyDESeq2** in Python.

Why you care: this is where a naive fold-change gives confidently wrong answers,
and it is the most common analysis mistake a programmer makes here.

```python
import math

for gene, control, treated in [("geneA", 100, 400), ("geneB", 2, 8)]:
    print(f"{gene}: log2 fold change = {math.log2(treated / control)}")
```

Out:
```
geneA: log2 fold change = 2.0
geneB: log2 fold change = 2.0
```

Identical fold change. Completely different confidence. `geneA` rests on 500
reads, `geneB` on 10. A real method models that spread and returns a p-value —
which is precisely why you should not write this analysis yourself.

Link: [PyDESeq2](https://pydeseq2.readthedocs.io/en/stable/)

---

# Beyond Biopython: scikit-bio

Biopython is not the only library. [**scikit-bio**](https://github.com/scikit-bio/scikit-bio)
is a community-built package that covers the half Biopython does not: statistics.

It is genuinely maintained — version 0.7.3 was released in June 2026, and the
repository had commits within days of this file being written. BSD-3 licensed,
Python 3.10 and newer.

**The split is clean:**

| | Biopython | scikit-bio |
|---|---|---|
| Best at | reading files, fetching records, sequence operations | statistics, distances, diversity, ordination |
| Objects | its own `Seq` and `SeqRecord` | numpy and pandas native |
| Reach for it when | you need to parse or convert something | you need to answer "are these samples different?" |

**Where it earns its place:**

1. **Microbiome and community work.** Alpha diversity (how varied is one sample), beta diversity (how different are two samples), and PERMANOVA to test whether groups truly separate. Biopython has none of this.
2. **Ordination.** PCoA squashes a distance matrix into two dimensions you can plot — the standard way to show sample similarity.
3. **Sequences as data.** Its `DNA` object exposes `gc_content()`, `kmer_frequencies()` and motif search through one consistent API, and results come back ready for pandas.

```python
from skbio import DNA
from skbio.diversity import alpha

print(DNA("ACGTACGTAC").kmer_frequencies(3))
print("gc:", round(DNA("ACCGGGTTTTA").gc_content(), 4))

# how varied is one sample, given counts of four species
print("shannon diversity:", round(alpha.shannon([4, 3, 2, 1]), 4))
```

Out:
```
{'ACG': 2, 'CGT': 2, 'GTA': 2, 'TAC': 2}
gc: 0.4545
shannon diversity: 1.2799
```

**Install it separately, and know what you are asking for.** It is deliberately
not in `requirements.txt`:

```bash
pip install scikit-bio
```

That pulls scipy, pandas, h5py, statsmodels and biom-format. On Windows,
biom-format had to **compile from source** during testing, which needs a C
compiler installed. Biopython arrives as a wheel in seconds; this does not. On
Linux or macOS it is usually smooth, and conda avoids the problem entirely.

Link: [scikit-bio](https://scikit.bio/) · [GitHub](https://github.com/scikit-bio/scikit-bio)

---

## Where to go next

- **Practice problems:** [Rosalind](https://rosalind.info/problems/locations/) — bioinformatics as coding exercises, graded automatically. Start at "Bioinformatics Stronghold".
- **The reference:** [Biopython Tutorial](https://biopython.org/docs/latest/Tutorial/index.html) — long, but the chapter list doubles as a map of the field.
- **Bulk downloads:** [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/) — easier than Entrez when you want whole genomes.
