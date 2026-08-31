# LEARN — Bioinformatics for people who already know Python

You know Python. You do not know biology. This file is the shortest path
from there to being useful.

Twenty points. Each one tells you **what it is**, **why you will hit it**,
and gives a snippet you can paste and run.

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

## Where to go next

- **Practice problems:** [Rosalind](https://rosalind.info/problems/locations/) — bioinformatics as coding exercises, graded automatically. Start at "Bioinformatics Stronghold".
- **The reference:** [Biopython Tutorial](https://biopython.org/docs/latest/Tutorial/index.html) — long, but the chapter list doubles as a map of the field.
- **Bulk downloads:** [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/) — easier than Entrez when you want whole genomes.
