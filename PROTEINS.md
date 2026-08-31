# PROTEINS — sequence work meets 3D structure

Six short points, same shape as [LEARN.md](LEARN.md): what it is, why you care,
a snippet you can run.

DNA work is text processing. Protein work is text processing **plus geometry**,
because a protein's job depends on the shape it folds into. This file covers
both halves.

Run everything from the repo root. The structure examples use `data/sample.pdb`,
which is [1CRN](https://www.rcsb.org/structure/1CRN) — crambin, a small plant
protein of 46 amino acids, solved in 1981 and still the standard tiny test case.

Verified on Python 3.14.3 / Biopython 1.88.

## Contents

1. [Amino acid properties](#1-amino-acid-properties)
2. [Hydrophobicity](#2-hydrophobicity)
3. [Isoelectric point](#3-isoelectric-point)
4. [Reading a PDB file](#4-reading-a-pdb-file)
5. [Secondary structure](#5-secondary-structure)
6. [Domains and motifs](#6-domains-and-motifs)

---

## 1. Amino acid properties

Proteins use a 20-letter alphabet instead of DNA's 4. Each letter is one amino
acid, and each has real chemistry: size, charge, water-loving or water-fearing.
Files write them two ways — one letter (`A`) in sequences, three letters
(`Ala`) in structure files.

Why you care: you will need to convert between the two forms constantly, and
"molecular weight" is the first number anyone asks you for.

```python
from Bio.SeqUtils import seq1, seq3
from Bio.SeqUtils.ProtParam import ProteinAnalysis

print(seq1("AlaGlyLys"), "|", seq3("AGK"))

protein = ProteinAnalysis("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN")
print("weight:", round(protein.molecular_weight(), 1), "daltons")
print("most common:", sorted(protein.count_amino_acids().items(),
                             key=lambda kv: -kv[1])[:3])
```

```
AGK | AlaGlyLys
weight: 4736.4 daltons
most common: [('C', 6), ('T', 6), ('A', 5)]
```

Six cysteines (`C`) in 46 residues is unusually high. Cysteines form bridges
that lock a protein together, which is exactly why crambin is so stable.

Link: [Bio.SeqUtils.ProtParam](https://biopython.org/docs/latest/api/Bio.SeqUtils.ProtParam.html)

## 2. Hydrophobicity

Some amino acids avoid water, others attract it. When a protein folds, the
water-fearing ones hide in the middle and the water-loving ones face out. The
**Kyte-Doolittle** scale puts a number on this: positive means hydrophobic.

Why you care: a run of hydrophobic residues usually means the protein crosses a
membrane. This is how transmembrane regions get predicted, and it is just a
sliding window over a list of numbers — ordinary Python.

```python
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd     # Kyte-Doolittle scale

protein = ProteinAnalysis("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN")
print("GRAVY (whole protein):", round(protein.gravy(), 3))

profile = protein.protein_scale(kd, window=9, edge=1.0)   # sliding window of 9
print("windows:", len(profile))
print("most hydrophobic window starts at residue", profile.index(max(profile)) + 1,
      "score", round(max(profile), 2))
```

```
GRAVY (whole protein): 0.37
windows: 38
most hydrophobic window starts at residue 32 score 1.96
```

GRAVY is the average over the whole protein. The window version is the useful
one — it shows you *where* the hydrophobic patch is.

Link: [Hydrophilicity plot](https://en.wikipedia.org/wiki/Hydrophilicity_plot)

## 3. Isoelectric point

The pI is the pH at which a protein carries no net charge. Below its pI the
protein is positive, above it the protein is negative.

Why you care: it is the number that decides how you physically separate proteins
in the lab, and it is a quick sanity check that you translated the right frame.

```python
from Bio.SeqUtils.ProtParam import ProteinAnalysis

protein = ProteinAnalysis("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN")
print("pI:", round(protein.isoelectric_point(), 2))
print("charge at pH 7:", round(protein.charge_at_pH(7.0), 2))
```

```
pI: 5.73
charge at pH 7: -0.66
```

pI 5.73 means this protein is slightly negative at the pH of blood. That follows
from the sequence alone — no experiment needed.

Link: [Isoelectric point](https://en.wikipedia.org/wiki/Isoelectric_point)

## 4. Reading a PDB file

A **PDB** file lists the xyz coordinate of every atom in a solved structure.
Biopython models it as a strict hierarchy: Structure → Model → Chain → Residue →
Atom. Learn that chain and the whole module opens up.

Why you care: once you have coordinates, structural questions become geometry
questions, and geometry is just arithmetic.

```python
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

structure = PDBParser(QUIET=True).get_structure("crambin", "data/sample.pdb")

residues = [r for r in structure.get_residues() if r.id[0] == " "]   # drop water
print("chains:", len(list(structure.get_chains())),
      "residues:", len(residues),
      "atoms:", len(list(structure.get_atoms())))
print("sequence:", seq1("".join(r.get_resname() for r in residues)))
print("first CA at xyz:", [round(float(x), 2) for x in residues[0]["CA"].coord])
print("CA1 to CA2 distance:", round(float(residues[0]["CA"] - residues[1]["CA"]), 2), "angstroms")
```

```
chains: 1 residues: 46 atoms: 327
sequence: TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN
first CA at xyz: [16.97, 12.78, 4.34]
CA1 to CA2 distance: 3.79 angstroms
```

Two things worth memorising. `r.id[0] == " "` keeps only real amino acids and
throws away water and ions, which are also stored as residues. And subtracting
two atoms gives you the distance between them in angstroms — Biopython
overloads the minus operator for exactly this.

Link: [Biopython structural bioinformatics](https://biopython.org/docs/latest/Tutorial/chapter_pdb.html)

## 5. Secondary structure

Before a protein folds into its final shape, stretches of it form two repeating
local patterns: **alpha helices** (coils) and **beta sheets** (flat strands).
Everything else is loop. PDB files record which residues belong to which.

Why you care: helix and sheet positions are how you compare two structures, and
reading them straight from the file needs no extra tool.

```python
for line in open("data/sample.pdb"):
    if line.startswith("HELIX"):
        start, end = line[21:25], line[33:37]
    elif line.startswith("SHEET"):
        start, end = line[22:26], line[33:37]   # SHEET's start column differs by one
    else:
        continue
    print(f"{line[:5]} residues {start.strip()}-{end.strip()}")
```

```
HELIX residues 7-19
HELIX residues 23-30
SHEET residues 1-4
SHEET residues 32-35
```

Those columns are fixed-width, not delimited — the PDB format predates CSV
thinking, and HELIX and SHEET do not even use the same columns for the starting
residue. Read the spec, count the characters, and test against a real file.

If a file has no HELIX records, the depositors simply did not annotate them.
Computing secondary structure yourself needs the external **DSSP** program,
which Biopython wraps as `Bio.PDB.DSSP`.

You can also estimate the fractions from sequence alone, with no structure file:

```python
from Bio.SeqUtils.ProtParam import ProteinAnalysis

protein = ProteinAnalysis("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN")
helix, turn, sheet = protein.secondary_structure_fraction()
print(f"helix {helix:.0%}  turn {turn:.0%}  sheet {sheet:.0%}")
```

```
helix 15%  turn 33%  sheet 37%
```

Treat that as a rough guess, not a measurement. It counts which amino acids tend
to appear in each state; it does not know the actual fold.

Link: [PDB file format](https://www.wwpdb.org/documentation/file-format)

## 6. Domains and motifs

A **domain** is a chunk of protein that folds and works on its own, and the same
domain shows up across many different proteins. A **motif** is smaller — a short
pattern of residues that marks a function, written as a regular expression.

Why you care: this is the point where biology hands you a problem that is
literally `re`. PROSITE patterns translate almost directly into Python regex.

```python
import re

protein = "MAMAPRTEINSTRING"

# N-glycosylation site, PROSITE pattern N-{P}-[ST]:
# an N, then anything except P, then S or T
sites = [(m.start(1) + 1, m.group(1)) for m in re.finditer(r"(?=(N[^P][ST]))", protein)]
print("sites:", sites)
```

```
sites: [(10, 'NST')]
```

The lookahead `(?=...)` again, so overlapping hits are not missed. `{P}` in
PROSITE means "not P", which is `[^P]` in Python.

For real domain assignment you do not write the patterns yourself — you scan
against **Pfam** or **InterPro**, which hold tens of thousands of curated
models. The Python part of that job is submitting sequences and parsing results.

Link: [PROSITE patterns](https://prosite.expasy.org/prosuser.html) · [InterPro](https://www.ebi.ac.uk/interpro/)

---

Back to [LEARN.md](LEARN.md) · vocabulary in [GLOSSARY.md](GLOSSARY.md)
