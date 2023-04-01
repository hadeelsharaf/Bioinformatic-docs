# Python for Bioinformatics 

Bioinformatics is the field that develops methods and software tools for understanding biological data. Units 6 and 7 in [this course](https://www.khanacademy.org/science/ap-biology) will help with understanding the basics of biology for this field.  

Next-generation sequencing (NGS) is one of the fundamental technological developments. Whole-genome sequencing (WGS), restriction site-associated DNA sequencing (RAD-Seq), ribonucleic acid sequencing (RNA-Seq), chromatin immunoprecipitation sequencing (ChIP-Seq), and several other technologies are routinely used to investigate important biological problems. These are called high-throughput (HT) sequencing technologies. [See this](https://htseq.readthedocs.io/en/master/index.html) for a python package to help with the HT sequencing. 

DNA in text files is represented as a string with sequence of [specific characters](https://www.bioinformatics.org/sms/iupac.html); so, knowing about the following topics will be helpful:

- File processing (txt and csv).
- String and [Regex](https://docs.python.org/3/library/re.html) functions.
- [BioPython](https://biopython.org/wiki/Documentation).

*Examples of Functions Bioinformatics:* 

1) Counting bases in a DNA sequence (Tetranucleotide Frequency): 
   (Aspose.Words.dbe14a37-9efa-4589-ae6f-146889bcecbe.001.png)
1) Reverse Complement of DNA:
   (Aspose.Words.dbe14a37-9efa-4589-ae6f-146889bcecbe.002.png)
1) Computing GC Content: A higher GC content level indicates a relatively higher melting temperature in molecular biology, and DNA sequences that encode proteins tend to be found in GC-rich regions.
   (Aspose.Words.dbe14a37-9efa-4589-ae6f-146889bcecbe.003.png)
1) Transcribing DNA into mRNA: regions of DNA must be transcribed into a form of RNA called messenger RNA (mRNA).
   (Aspose.Words.dbe14a37-9efa-4589-ae6f-146889bcecbe.004.png)
1) Translating mRNA into Protein: mRNA makes protein. 
   (Aspose.Words.dbe14a37-9efa-4589-ae6f-146889bcecbe.005.png)
   Note on points 4,5: these functions can be done using string replacement and regex but using BioPython is the recommended approach.
1) Finding Open Reading Frames ORF:  finding a region in DNA or RNA. 
   using regex: This region starts with M and ends with (\*).
   (Aspose.Words.dbe14a37-9efa-4589-ae6f-146889bcecbe.006.png)
   the following section is applied after a series of transcribing and translating steps
   (Aspose.Words.dbe14a37-9efa-4589-ae6f-146889bcecbe.007.png)




*Sequence file extensions:*

(Aspose.Words.dbe14a37-9efa-4589-ae6f-146889bcecbe.008.png)

To read or write to a file: 

(Aspose.Words.dbe14a37-9efa-4589-ae6f-146889bcecbe.009.png)

For compressed a fastq files:
(Aspose.Words.dbe14a37-9efa-4589-ae6f-146889bcecbe.010.png)
