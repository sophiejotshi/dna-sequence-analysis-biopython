# DNA Sequence Analysis Using Biopython

## About the Project
This project is a basic bioinformatics pipeline written in Python to analyze DNA sequences in FASTA format. It was created as a mini project to understand how biological sequence data can be processed programmatically using Biopython.

The SARS-CoV-2 reference genome is used as the input sequence to demonstrate real biological data handling.


## What This Project Does
- Reads DNA sequences from a FASTA file
- Calculates DNA length
- Computes GC content and AT content
- Calculates nucleotide composition (%A, %T, %G, %C)
- Translates DNA into a protein sequence
- Counts stop codons in the reading frame
- Saves all results into a CSV file for further analysis


## Dataset Used
- SARS-CoV-2 reference genome  
- Accession ID: NC_045512.2 
- Source: NCBI Nucleotide Database


## Tools Used
- Python  
- Biopython  
- Pandas   


## How the Code Works
1. The FASTA file is read using Biopython’s 'SeqIO'.
2. Nucleotide counts and percentages are calculated.
3. The sequence is trimmed so its length is a multiple of three.
4. DNA is translated into a protein sequence until a stop codon is reached.
5. Stop codons are counted across the sequence.
6. All calculated features are stored in a CSV file.


## Biological Interpretation
- The genome length (~29.9 kb) matches the known size of SARS-CoV-2.
- The GC content reflects the AT-rich nature of RNA viruses.
- Protein translation stops early due to stop codons in the primary reading frame.
- The presence of many stop codons highlights the importance of proper ORF detection when analyzing viral genomes.


## Output
- dna_features_sars_cov_2.csv 
This file contains all calculated DNA features and can be used for visualization or further analysis in Excel, Power BI, or Python.



