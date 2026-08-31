# GLOSSARY — the ten words that block a Python developer

Papers and tool docs assume you know these. Nobody defines them, because
everyone around you learned them in a biology class you did not take.

Ten terms. Each one: a plain definition, why it matters to you, and a short
snippet that makes the word concrete instead of abstract.

Run the snippets from the repo root. Verified on Python 3.14.3 / Biopython 1.88.

## Contents

1. [Gene](#1-gene)
2. [Exon and intron](#2-exon-and-intron)
3. [Genome, transcriptome, proteome](#3-genome-transcriptome-proteome)
4. [SNP and variant](#4-snp-and-variant)
5. [Allele](#5-allele)
6. [Read](#6-read)
7. [Contig and scaffold](#7-contig-and-scaffold)
8. [Coverage and depth](#8-coverage-and-depth)
9. [Reference genome](#9-reference-genome)
10. [Annotation](#10-annotation)

---

## 1. Gene

A stretch of DNA that carries the instructions for one product, usually a
protein. Think of the genome as a very long file and a gene as one function
definition inside it.

Why you care: nearly every question you will be asked is really "what is
happening at this gene", so "find the region, then do something to it" is the
shape of most of your code.

```python
from Bio import SeqIO

genome = next(SeqIO.parse("data/sample.fasta", "fasta"))
gene = genome.seq[0:39]          # a gene is just a slice with a meaning attached
print("gene length:", len(gene))
print("product:", gene.translate())
```

```
gene length: 39
product: MAIVMGR*KGAR*
```

## 2. Exon and intron

In humans and other complex organisms, a gene is interrupted. **Exons** are the
parts that survive into the final message. **Introns** are cut out and thrown
away. The cutting is called splicing.

Why you care: the DNA sequence of a gene is *not* the sequence that makes the
protein. If you translate raw genomic DNA and get garbage, unremoved introns
are the usual reason.

```python
from Bio.Seq import Seq

gene = Seq("ATGGCC" + "GTAAGTCCTTAG" + "ATTGTAATGGGCCGCTGA")
#           exon 1     intron (removed)  exon 2

print("with intron :", gene.translate())
mrna = gene[:6] + gene[18:]       # splice the intron out
print("spliced     :", mrna.translate())
```

```
with intron : MAVSP*IVMGR*
spliced     : MAIVMGR*
```

Same DNA. Only the spliced version gives the real protein.

## 3. Genome, transcriptome, proteome

Three layers, three questions:

| Word | What it is | The question it answers |
|---|---|---|
| Genome | all the DNA | what *could* this cell do |
| Transcriptome | all the RNA present right now | what is it *actually* doing |
| Proteome | all the proteins present right now | what is doing the *work* |

Why you care: it tells you which file type you will be handed. Genome work
means FASTA and VCF. Transcriptome work means FASTQ and count tables.

```python
from Bio.Seq import Seq

dna = Seq("ATGGCCATTGTAATGGGCCGCTGA")
print("genome layer       :", dna)
print("transcriptome layer:", dna.transcribe())
print("proteome layer     :", dna.translate())
```

```
genome layer       : ATGGCCATTGTAATGGGCCGCTGA
transcriptome layer: AUGGCCAUUGUAAUGGGCCGCUGA
proteome layer     : MAIVMGR*
```

## 4. SNP and variant

A **variant** is any place where a sample's DNA differs from the reference. A
**SNP** (say "snip") is the simplest kind: one single base changed.

Why you care: "find the differences" is the core job of most clinical and
population work, and it reduces to comparing two strings position by position.

```python
reference = "ACGTACGTAC"
sample    = "ACGTATGTAC"

variants = [(i + 1, r, s) for i, (r, s) in enumerate(zip(reference, sample)) if r != s]
for pos, ref, alt in variants:
    print(f"position {pos}: {ref} -> {alt}")
```

```
position 6: C -> T
```

Real variant calling is harder than this, because reads carry errors and you
must decide whether a difference is real or noise. But the idea is this.

## 5. Allele

You carry two copies of most genes, one from each parent. The different versions
that can sit at one position are its **alleles**. If both copies match you are
*homozygous*; if they differ you are *heterozygous*.

Why you care: a variant is not just "present or absent" — it is present on one
copy or two. That third state is why VCF files store a genotype like `0/1` and
why a naive boolean model of variants is wrong.

```python
genotypes = {"sample1": "0/0", "sample2": "0/1", "sample3": "1/1"}

for name, gt in genotypes.items():
    a, b = gt.split("/")            # 0 = reference allele, 1 = alternate allele
    kind = "homozygous" if a == b else "heterozygous"
    print(f"{name}: {gt}  {kind}, alternate copies = {int(a) + int(b)}")
```

```
sample1: 0/0  homozygous, alternate copies = 0
sample2: 0/1  heterozygous, alternate copies = 1
sample3: 1/1  homozygous, alternate copies = 2
```

## 6. Read

Sequencing machines cannot read a whole chromosome. They read short fragments,
typically 50-300 bases each, and produce millions of them. Each fragment is a
**read**.

Why you care: this is the unit of raw data. A FASTQ file is a pile of reads, and
almost every early pipeline step is a loop over reads.

```python
from Bio import SeqIO

reads = list(SeqIO.parse("data/sample.fastq", "fastq"))
print("reads:", len(reads))
print("read length:", len(reads[0].seq))
print("total bases sequenced:", sum(len(r.seq) for r in reads))
```

```
reads: 2
read length: 36
total bases sequenced: 72
```

Two reads is a toy. A real run gives you tens of millions, which is why
streaming matters.

## 7. Contig and scaffold

If you have no reference genome, you must rebuild it from the reads. Overlapping
reads get merged into a **contig** — a continuous stretch. Contigs that you can
order and orient, but not fully join, get linked into a **scaffold**, with runs
of `N` standing in for the unknown gaps.

Why you care: `N` in an assembly is not a sequencing error. It is an honest "we
know something is here, we do not know what".

```python
contig_a = "ACGTACGTAC"
contig_b = "TTTTGGGGCC"
scaffold = contig_a + "N" * 20 + contig_b     # 20 bases of estimated gap

print("scaffold length:", len(scaffold))
print("known bases:", len(scaffold) - scaffold.count("N"))
print("gap fraction:", f"{scaffold.count('N') / len(scaffold):.0%}")
```

```
scaffold length: 40
known bases: 20
gap fraction: 50%
```

## 8. Coverage and depth

**Depth** at a position is how many reads cover it. **Coverage** usually means
the average depth across the whole target, written like `30x`.

Why you care: depth is the confidence number. A difference seen in 1 read out of
40 is probably an error; seen in 20 out of 40, it is probably a real
heterozygous variant. Low depth is the reason most doubtful calls are doubtful.

```python
from collections import Counter

reads = [(100, 136), (110, 146), (120, 156), (500, 536)]   # start, end per read
depth = Counter(pos for start, end in reads for pos in range(start, end))

print("bases covered:", len(depth))
print("max depth:", max(depth.values()))
print("depth at 125:", depth[125], "| at 505:", depth[505], "| at 400:", depth[400])
```

```
bases covered: 92
max depth: 3
depth at 125: 3 | at 505: 1 | at 400: 0
```

Position 400 returns 0 rather than raising, because `Counter` treats missing keys
as zero — which happens to be exactly right for depth.

## 9. Reference genome

One agreed sequence per species that everyone measures against. The current human
one is **GRCh38**; you will still meet the older **hg19**. Coordinates only mean
something relative to a named reference build.

Why you care: mixing coordinates from two builds silently gives wrong answers.
The positions look perfectly valid — they just point somewhere else. Always
record which build your numbers came from.

```python
variant = {"chrom": "chr1", "pos": 155_235_252, "ref": "C", "alt": "T",
           "build": "GRCh38"}                    # never store a position without this

print(f"{variant['chrom']}:{variant['pos']} {variant['ref']}>{variant['alt']}"
      f" ({variant['build']})")
```

```
chr1:155235252 C>T (GRCh38)
```

## 10. Annotation

The reference genome is just letters. **Annotation** is the separate layer that
says what those letters mean: this range is a gene, that range is an exon.
It arrives as GFF, GTF or BED, and it gets updated independently of the genome.

Why you care: your analysis is only as current as your annotation file. Two
people can use the same genome build, different annotation releases, and get
different gene counts from the same reads.

```python
import pandas as pd

annotation = pd.read_csv("data/regions.bed", sep="\t", header=None,
                         names=["chrom", "start", "end", "feature"])
annotation["length"] = annotation["end"] - annotation["start"]
print(annotation)
```

```
  chrom  start  end   feature  length
0  chr1    100  250  promoter     150
1  chr1    500  900     exon1     400
2  chr2     50  120     exon2      70
```

---

Next: the roadmap in [LEARN.md](LEARN.md) · proteins in [PROTEINS.md](PROTEINS.md)
