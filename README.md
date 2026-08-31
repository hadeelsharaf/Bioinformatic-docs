# Python for Bioinformatics 

> **New to bioinformatics but comfortable with Python?** Start with **[LEARN.md](LEARN.md)** — 20 short points, each with a snippet you can run. Sample files are in `data/`.

Bioinformatics is the field that develops methods and software tools for understanding biological data. Units 6 and 7 in [this course](https://www.khanacademy.org/science/ap-biology) will help with understanding the *basics* of biology for this field.  

Next-generation sequencing (NGS) is one of the fundamental technological developments. Whole-genome sequencing (WGS), restriction site-associated DNA sequencing (RAD-Seq), ribonucleic acid sequencing (RNA-Seq), chromatin immunoprecipitation sequencing (ChIP-Seq), and several other technologies are routinely used to investigate important biological problems. These are called high-throughput (HT) sequencing technologies. [See this](https://htseq.readthedocs.io/en/master/index.html) for a python package to help with the HT sequencing. 

DNA in text files is represented as a string with sequence of [specific characters](https://www.bioinformatics.org/sms/iupac.html); so, knowing about the following topics will be helpful:

- File processing (txt and csv).
- String and [Regex](https://docs.python.org/3/library/re.html) functions.
- [BioPython](https://biopython.org/wiki/Documentation).

*Examples of Functions Bioinformatics:* 

All examples below use the same sequence:

```python
seq = 'ACCGGGTTTTA'
```

- Counting bases in a DNA sequence (Tetranucleotide Frequency): 

```python
seq.count('T')
```

```
4
```

```python
from collections import Counter

Counter(seq)
```

```
Counter({'T': 4, 'G': 3, 'A': 2, 'C': 2})
```

- Reverse Complement of DNA:

```python
def reverse_dna(dna):
    trans = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a'
    }
    complement = [trans.get(base, base) for base in dna]
    print(''.join(reversed(complement)))

reverse_dna(seq)
```

```
TAAAACCCGGT
```

   Note: the function above prints instead of returning, so its result cannot be reused. Returning `''.join(reversed(complement))` is more useful. BioPython does the same thing in one call: `Seq(seq).reverse_complement()`.

- Computing GC Content:
 A higher GC content level indicates a relatively higher melting temperature in molecular biology, and DNA sequences that encode proteins tend to be found in GC-rich regions.

```python
import re

GC_count = len(re.findall('[GC]', seq))
GC_percent = (GC_count * 100) / len(seq)
print(GC_percent)
```

```
45.45454545454545
```

- Transcribing DNA into mRNA: regions of DNA must be transcribed into a form of RNA called messenger RNA (mRNA).

```python
from Bio import Seq

Seq.transcribe(seq)
```

```
'ACCGGGUUUUA'
```

- Translating mRNA into Protein: mRNA makes protein. 

```python
rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'

from Bio import Seq

Seq.translate(rna)
```

```
'MAMAPRTEINSTRING*'
```

   Note on points 4,5: these functions can be done using string replacement and regex but using BioPython is the recommended approach.
- Finding Open Reading Frames ORF:  finding a region in DNA or RNA. 
   using regex: This region starts with M and ends with (\*).
   
  ![Regex diagram for the ORF pattern: M, then any run of non-star characters, then a star](./Aspose.Words.dbe14a37-9efa-4589-ae6f-146889bcecbe.006.png) 
  
   the following section is applied after a series of transcribing and translating steps

```python
import re

re.findall('(?=(M[^*]*)[*])', 'MAMAPR*MP*M')
```

```
['MAMAPR', 'MAPR', 'MP']
```

   The lookahead `(?=...)` lets the matches overlap, so a protein that starts at a later `M` is found too. This only scans one reading frame — see [LEARN.md point 12](LEARN.md#12-open-reading-frames-in-all-six-frames) for all six.

*Sequence file extensions:*

| Extension | Format |
| --- | --- |
| `.fasta` | FASTA |
| `.fa` | FASTA |
| `.fna` | FASTA (nucleotides) |
| `.faa` | FASTA (amino acids) |
| `.fq` | FASTQ |
| `.fastq` | FASTQ |

- To read or write to a file: 

```python
from Bio import SeqIO

# to read records in a fasta file
records = list(SeqIO.parse('data/sample.fasta', 'fasta'))
for record in records:
    print(record.id, len(record.seq))

# to write records to a fasta file
SeqIO.write(records, 'example.fasta', 'fasta')
```

```
seq1 39
seq2 36
seq3 23
3
```

   Do not open the file yourself in `'w'` mode first. `SeqIO.write` takes the
   filename directly, and opening in `'w'` mode truncates the file — if you then
   forget to write, the file is left empty.

For compressed a fastq files:

```python
# NEEDS YOUR OWN FILE: SRR003265.filt.fastq.gz is not in this repository.
import gzip
from Bio import SeqIO

# parsing records in compressed fastq
records = SeqIO.parse(gzip.open('./SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')

record = next(records)
```

   `'rt'` means read as text. Without the `t` you get bytes and the parser fails.
   There is a runnable version of this on the sample file in `data/` — see
   [LEARN.md point 7](LEARN.md#7-everything-is-gzipped-stream-never-read).
